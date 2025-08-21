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
from pathlib import Path
from shutil import which
import sys


# TODO: To discuss whether this should rather be a setting in
# a scheduling script:
# # set virtual cores for NumExpr
# # 64 == default, max virtual cores in NumExpr:
# os.environ['NUMEXPR_MAX_THREADS'] = "64"
# os.environ['NUMEXPR_NUM_THREADS'] = "32"


APP_NAME = "mcce4.mcce_benchmark"

# fail fast:
USER_MCCE = which("mcce")
if USER_MCCE is None:
    #raise EnvironmentError(f"{APP_NAME}, init :: mcce executable not found.")
    sys.exit(f"EnvironmentError - {APP_NAME}, init :: mcce executable not found.")

USER_MCCE = Path(USER_MCCE).parent


ENTRY_POINTS = {
    "setup": "bench_setup",
    "launch": "bench_batch",  # used by crontab :: launch 1 batch of size 10
    "analyze": "bench_analyze",
    "compare": "bench_compare",
}

# bench_setup sub-commands; sub-commands SUB1, SUB2, SUB3 are the datasets folder names:
SUB0 = "pdbids"
SUB1 = "pkdb_v1"
SUB2 = "pkdb_vR"
SUB3 = "pkdb_snase"
SUB4 = "user_pdbs"
SUB5 = "launch" # :: crontab job scheduler step of setup.


# full path of the launch EP for scheduling:
LAUNCHJOB = which(ENTRY_POINTS[SUB5])
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
    MATCHED_PKAS_TXT = "matched_pkas.txt"
    MATCHED_PKAS_STATS_PKL = "matched_pkas_stats.pickle"
    MATCHED_PKAS_STATS = "matched_pkas_stats.txt"  # saved pka_stats_dict["report"]
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


# dataset-specific pdb count is now in BenchResources
N_PDBS = 156   # : size of the largest dataset, pkdb_vR 

RUNS_DIR = "runs"
ANALYZE_DIR = "analysis"
MCCE_EPS = 4  # default dielectric constant (epsilon) in MCCE
N_BATCH = 10  # number of jobs to maintain in the process queue
DEFAULT_JOB = "default_run"
DEFAULT_JOB_SH = f"{DEFAULT_JOB}.sh"
Q_BOOK = "book.txt"
BENCH_DATA = resources.files(f"{APP_NAME}.data")
datasets_dict = dict((d.name, d) for d in BENCH_DATA.glob("pkdb_*"))


class BenchResources:
    """Class to store package data paths and main constants for a selected dataset."""

    def __init__(self, dataset: str, data_dir=BENCH_DATA):
        # check given dataset corresponds to an existing folder:
        if not data_dir.joinpath(dataset).is_dir():
            sys.exit(f"Dataset: {dataset!r} is not setup.")
        
        self.BENCH_DATA = data_dir
        self.BENCH_DB = self.BENCH_DATA.joinpath(dataset)
        self.BENCH_WT = self.BENCH_DB.joinpath("pkas.csv")
        self.BENCH_PROTS = self.BENCH_DB.joinpath("proteins.tsv")
        self.N_PDBS = self.prots_count()
        self.BENCH_PDBS = self.BENCH_DB.joinpath(RUNS_DIR)
        #self.DEFAULT_JOB_SH = self.BENCH_PDBS.joinpath(f"{DEFAULT_JOB}.sh")
        self.DEFAULT_JOB_SH = self.BENCH_DB.joinpath(DEFAULT_JOB_SH)        
        self.BENCH_Q_BOOK = self.BENCH_PDBS.joinpath(Q_BOOK)
        if dataset == "pkdb_v1":
            self.BENCH_PH_REFS = self.BENCH_DB.joinpath("refsets")
            self.BENCH_PARSE_PHE4 = self.BENCH_PH_REFS.joinpath("parse.e4")

    def prots_count(self) -> int:
        with open(self.BENCH_PROTS) as fp:
            lines = fp.readlines()
        return len([ln for ln in lines if not ln.startswith("#")])

    def __str__(self):
        return f"""
        BENCH_DATA     = {str(self.BENCH_DATA)}
        BENCH_DB       = {str(self.BENCH_DB)}
        BENCH_WT       = {str(self.BENCH_WT)}
        BENCH_PROTS    = {str(self.BENCH_PROTS)}
        BENCH_PDBS     = {str(self.BENCH_PDBS)}
        DEFAULT_JOB_SH = {str(self.DEFAULT_JOB_SH)}
        BENCH_Q_BOOK   = {str(self.BENCH_Q_BOOK)}
        N_PDBS = {self.N_PDBS}
        """

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


USER = getpass.getuser()


LOG_HDR = f"""
START\n{'-'*70}\n{datetime.now().strftime(format="%Y-%m-%d %H:%M:%S")} - {USER = }
APP DEFAULTS:
Globals:
{USER_MCCE = } :: mcce executable in use
{MCCE_EPS = } :: default protein epsilon
{N_BATCH = } :: default batch size for job submission
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
