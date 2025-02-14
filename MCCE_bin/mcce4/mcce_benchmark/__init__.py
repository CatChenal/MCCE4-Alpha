#!/usr/bin/env python

"""
Module: __ini__.py

Init module of mcce_benchmark.
"""
from datetime import datetime
from enum import Enum
import getpass
from importlib import resources
from inspect import signature
import logging
import os
from pathlib import Path
from shutil import which


# set virtual cores for NumExpr
# 64 == default, max virtual cores in NumExpr:
os.environ['NUMEXPR_MAX_THREADS'] = "64"
os.environ['NUMEXPR_NUM_THREADS'] = "32"


APP_NAME = "mcce4.mcce_benchmark"

# fail fast:
USER_MCCE = which("mcce")
if USER_MCCE is None:
    raise EnvironmentError(f"{APP_NAME}, init :: mcce executable not found.")
USER_MCCE = Path(USER_MCCE).parent


ENTRY_POINTS = {
    "setup": "bench_setup",
    "launch": "bench_batch",  # used by crontab :: launch 1 batch
    "analyze": "bench_analyze",
    "compare": "bench_compare",
}

# bench_setup sub-commands, also used throughout:
SUB0 = "pdbids"
SUB1 = "pkdb_pdbs"
SUB2 = "user_pdbs"
SUB3 = "launch"  # :: crontab job scheduler step of setup.
# full path of the launch EP for scheduling:
LAUNCHJOB = which(ENTRY_POINTS[SUB3])
DELTA = "\u0394"


class FILES(Enum):
    """Output file names; _PKL => dict to .pickle"""
    SETUP_OPTIONS = "bench_setup_options.txt"
    PRERUN_RPT = "prerun_report.md"
    # Next 2: bench setup cli args, saved at the same time;
    # if benchmark.log is deleted, the text file provides
    # the setup options in readeable form, while the pickle
    # file is used by analysis and comparison modules.
    CLI_ARGS_PKL = "cli_args.pickle"
    CLI_ARGS_TXT = "cli_args.txt"
    # Analysis files
    ALL_SUMCRG = "all_sum_crg.out"
    ALL_SUMCRG_TXT = "all_sum_crg.out.txt"
    ALL_SUMCRG_DIFF = "all_sumcrg_diff.txt"
    ALL_PKAS = "all_pK.out"
    ALL_PKAS_TXT = "all_pK.out.txt"
    ALL_PKAS_OOB = "all_pkas_oob.txt"  # out of bounds pKas
    JOB_PKAS_PKL = "job_pkas.pickle"    
    # not saved: EXPL_PKAS_PKL
    MATCHED_PKAS_TXT = "matched_pkas.txt"
    MATCHED_PKAS_STATS_PKL = "matched_pkas_stats.pickle"
    MATCHED_PKAS_STATS = "matched_pkas_stats.txt"  # saved pka_stats_dict["report"]
    MULTICONF_RES = "multiconformer_residues.txt"
    RESIDUES_STATS = "residues_stats.txt"
    RESIDUES_STATS_PKL = "residues_stats.pickle"
    CONF_COUNTS = "conf_counts.txt"
    CONFS_PER_RES = "confs_per_res.txt"
    CONFS_THRUPUT = "confs_throughput.txt"
    RES_COUNTS = "res_counts.txt"
    RES_OUTLIER = "outlier_residues.txt"
    RESID_OUTLIER = "outlier_resids.txt"
    RUN_TIMES = "run_times.txt"
    FIG_CONFS_TP = "confs_throughput.png"
    FIG_FIT_ALLPKS = "pkas_fit.png"
    FIG_FIT_PER_RES = "res_analysis.png"
    VERSIONS = "versions.txt"


RUNS_DIR = "runs"
ANALYZE_DIR = "analysis"
MCCE_EPS = 4  # default dielectric constant (epsilon) in MCCE
N_BATCH = 10  # number of jobs to maintain in the process queue
N_PDBS = 118  # UPDATE if change in packaged data!

MCCE_OUTPUTS = [
    "acc.atm",
    "acc.res",
    "entropy.out",
    "extra.tpl",
    "fort.38",
    "head1.lst",
    "head2.lst",
    "head3.lst",
    "mc_out",
    "name.txt",
    "new.tpl",
    "null",
    "pK.out",
    "respair.lst",
    "rot_stat",
    "run.log",
    "run.prm",
    "run.prm.record",
    "step0_out.pdb",
    "step1_out.pdb",
    "step2_out.pdb",
    "sum_crg.out",
    "vdw0.lst",
]


class Bench_Resources:
    """Immutable class to store package data paths and main constants."""
    __slots__ = (
        "_BENCH_DATA",
        "_BENCH_NGPB_TOOLS",
        "_BENCH_DB",
        "_BENCH_WT",
        "_BENCH_PROTS",
        "_BENCH_PDBS",
        "_DEFAULT_JOB",
        "_DEFAULT_JOB_SH",
        "_Q_BOOK",
        "_BENCH_Q_BOOK",
        "_BENCH_PH_REFS",
        "_BENCH_PARSE_PHE4",
    )

    def __init__(self, res_files=resources.files(f"{APP_NAME}.data")):
        self._BENCH_DATA = res_files
        self._BENCH_NGPB_TOOLS = self._BENCH_DATA.joinpath("ngpb_deps_paths.json")
        self._BENCH_DB = self._BENCH_DATA.joinpath("pkadbv1")
        self._BENCH_WT = self._BENCH_DB.joinpath("WT_pkas.csv")
        self._BENCH_PROTS = self._BENCH_DB.joinpath("proteins.tsv")
        self._BENCH_PDBS = self._BENCH_DB.joinpath(RUNS_DIR)
        self._DEFAULT_JOB = "default_run"
        self._DEFAULT_JOB_SH = self._BENCH_PDBS.joinpath(f"{self._DEFAULT_JOB}.sh")
        self._Q_BOOK = "book.txt"
        self._BENCH_Q_BOOK = self._BENCH_PDBS.joinpath(self._Q_BOOK)
        self._BENCH_PH_REFS = self._BENCH_DB.joinpath("refsets")
        self._BENCH_PARSE_PHE4 = self._BENCH_PH_REFS.joinpath("parse.e4")

    @property
    def BENCH_DATA(self):
        return self._BENCH_DATA

    @property
    def BENCH_NGPB_TOOLS(self):
        return self._BENCH_NGPB_TOOLS

    @property
    def BENCH_DB(self):
        return self._BENCH_DB

    @property
    def BENCH_PH_REFS(self):
        return self._BENCH_PH_REFS

    @property
    def BENCH_PARSE_PHE4(self):
        return self._BENCH_PARSE_PHE4

    @property
    def BENCH_WT(self):
        return self._BENCH_WT

    @property
    def BENCH_PROTS(self):
        return self._BENCH_PROTS

    @property
    def BENCH_PDBS(self):
        return self._BENCH_PDBS

    @property
    def Q_BOOK(self):
        return self._Q_BOOK

    @property
    def BENCH_Q_BOOK(self):
        return self._BENCH_Q_BOOK

    @property
    def DEFAULT_JOB(self):
        return self._DEFAULT_JOB

    @property
    def DEFAULT_JOB_SH(self):
        return self._DEFAULT_JOB_SH

    def __str__(self):
        return f"""
        BENCH_DATA = {str(self.BENCH_DATA)}
        BENCH_NGPB_TOOLS = {str(self.BENCH_NGPB_TOOLS)}
        BENCH_PH_REFS = {str(self.BENCH_PH_REFS)}
        BENCH_PARSE_PHE4 = {str(self.BENCH_PARSE_PHE4)}
        BENCH_DB = {str(self.BENCH_DB)}
        BENCH_WT = {str(self.BENCH_WT)}
        BENCH_PROTS = {str(self.BENCH_PROTS)}
        BENCH_PDBS = {str(self.BENCH_PDBS)}
        DEFAULT_JOB = {str(self.DEFAULT_JOB)}
        DEFAULT_JOB_SH = {str(self.DEFAULT_JOB_SH)}
        BENCH_Q_BOOK = {str(self.BENCH_Q_BOOK)}
        Q_BOOK = {str(self.Q_BOOK)}
        """


BENCH = Bench_Resources()


class Opts:
    def __init__(
        self, app_name: str = APP_NAME, cli_name: str = None, d_args: dict = None
    ):
        """Class Opts makes command line options available to all objects.
        Arguments:
        ---------
        app_name : str
            API name
        cli_name : str
            name of the command line tool or subcommand.
        d_args : dict
            cli args (Namespace) given as vars(args).
        """
        self.app_name = app_name
        self.cli_name = cli_name
        self.all = d_args

    def __str__(self):
        out = f"__{self.app_name}.{self.cli_name}__ - User options:\n"
        if not self.all:
            return out + "  (not set)\n"

        for k in self.all:
            if k == "func":
                sig = signature(self.all["func"])
                out = out + f"  {k}: {self.all['func'].__name__}{sig}\n"
            else:
                out = out + f"  {k}: {self.all[k]}\n"

        return out


cli_opts = Opts()


# Config for root logger:
logger = logging.getLogger()
logger.setLevel(logging.INFO)
# console/screen handler
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch_formatter = logging.Formatter("%(levelname)s - %(funcName)s: %(message)s")
ch.setFormatter(ch_formatter)
# file handlers
flog_format = "[%(asctime)s - %(levelname)s]: %(name)s, %(funcName)s:\n\t%(message)s"
fh = logging.FileHandler("benchmark.log", encoding="utf-8")
fh.setLevel(logging.INFO)
fh_formatter = logging.Formatter(flog_format)
fh.setFormatter(fh_formatter)
# add to logger
logger.addHandler(ch)
logger.addHandler(fh)


USER = getpass.getuser()


LOG_HDR = f"""
START\n{'-'*70}\n{datetime.now().strftime(format="%Y-%m-%d %H:%M:%S")} - {USER = }
APP DEFAULTS:
Globals:
{USER_MCCE = } :: mcce executable in use
{MCCE_EPS = } :: default protein epsilon
{N_BATCH = } :: default batch size for job submission
{N_PDBS = } : number of pdbs in the dataset
Saved setup options:
  CLI_ARGS_PKL: {FILES.CLI_ARGS_PKL.value}
  CLI_ARGS_TXT: {FILES.CLI_ARGS_TXT.value}
Default analysis output file names (fixed):
  VERSIONS: {FILES.VERSIONS.value}
  RUN_TIMES: {FILES.RUN_TIMES.value}
  ALL_PKAS: {FILES.ALL_PKAS.value}
  ALL_SUMCRG: {FILES.ALL_SUMCRG.value}
  ALL_SUMCRG_DIFF: {FILES.ALL_SUMCRG_DIFF.value}
  ALL_PKAS_OOB: {FILES.ALL_PKAS_OOB.value}
  CONF_COUNTS: {FILES.CONF_COUNTS.value}
  RES_COUNTS: {FILES.RES_COUNTS.value}
  MULTICONF_RES: {FILES.MULTICONF_RES.value}
  CONFS_PER_RES: {FILES.CONFS_PER_RES.value}
  CONFS_THRUPUT: {FILES.CONFS_THRUPUT.value}
  FIG_CONFS_TP: {FILES.FIG_CONFS_TP.value}
  JOB_PKAS_PKL: {FILES.JOB_PKAS_PKL.value}
  MATCHED_PKAS_TXT: {FILES.MATCHED_PKAS_TXT.value}
  MATCHED_PKAS_STATS_PKL: {FILES.MATCHED_PKAS_STATS_PKL.value}
  MATCHED_PKAS_STATS: {FILES.MATCHED_PKAS_STATS.value}
  RESIDUES_STATS: {FILES.RESIDUES_STATS.value}
  RESIDUES_STATS_PKL: {FILES.RESIDUES_STATS_PKL.value}
  RES_OUTLIER = {FILES.RES_OUTLIER.value}
  RESID_OUTLIER = {FILES.RESID_OUTLIER.value}
  FIG_FIT_ALLPKS = {FILES.FIG_FIT_ALLPKS .value}
  FIG_FIT_PER_RES = {FILES.FIG_FIT_PER_RES.value}
\n{'-'*70}
"""
