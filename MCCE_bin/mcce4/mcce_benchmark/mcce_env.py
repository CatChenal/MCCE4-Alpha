#!/usr/bin/env python

"""
Module: mcce_env.py

Modified version of ENV class from Stable-MCCE/bin/pdbio.py:
  - Simplified: Only 2 methods: load_runprm,  __str__;
  - Needs "rundir_path" Path as class parameter
  - Attributes:
    self.runprm: dict
    self.rundir: Path
"""

import logging
import os
from pathlib import Path
import shutil
import sys
from typing import Union
from mcce4.mcce_benchmark import BENCH, RUNS_DIR, SUB1
from mcce4.mcce_benchmark import io_utils as iou


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


DEFAULT_RUNPRM = "run.prm.record"


class ENV:
    def __init__(self, rundir_path: str, runprm_file: str = DEFAULT_RUNPRM) -> dict:
        self.rundir = Path(rundir_path)
        self.runprm = {}
        self.titr_bounds = None
        self.sumcrg_hdr = ""
        self.non_default_prm_file = None if runprm_file == DEFAULT_RUNPRM else runprm_file
        # populate self.runprm dict:
        self.load_runprm(runprm_file)
        if self.runprm and self.non_default_prm_file is None:
            self.get_titr_bounds()

    def load_runprm(self, runprm_file: str = DEFAULT_RUNPRM):
        """NOTE: If ENV is used for comparing two runs, only DEFAULT_RUNPRM
                 (run.prm.record) is a valid file.
        """
        fp = Path(self.rundir.joinpath(runprm_file))
        if not fp.exists():
            raise FileNotFoundError(f"Not found: {runprm_file} in {self.rundir}")

        with open(fp) as fin:
            lines = fin.readlines()

        for line in lines:
            entry_str = line.strip().split("#")[0]
            fields = entry_str.split()
            if len(fields) > 1:
                key_str = fields[-1]
                if key_str[0] == "(" and key_str[-1] == ")":
                    key = key_str.strip("()").strip()
                    # inconsistant output in run.prm.record:
                    if key == "EPSILON_PROT":
                        value = str(round(float(fields[0]), 1))
                    else:
                        value = fields[0]
                    self.runprm[key] = value

        return

    def get_titr_bounds(self):
        """
        Populate self.titr_bounds

        ph       "ph" for pH titration, "eh" for eh titration       (TITR_TYPE)
        0.0      Initial pH                                         (TITR_PH0)
        1.0      pH interval                                        (TITR_PHD)
        0.0      Initial Eh                                         (TITR_EH0)
        30.0     Eh interval (in mV)                                (TITR_EHD)
        15       Number of titration points                         (TITR_STEPS)
        """
        prec = 0
        if self.runprm["TITR_TYPE"] == "ph":
            b1 = float(self.runprm["TITR_PH0"])
            step = float(self.runprm["TITR_PHD"])
        else:
            b1 = float(self.runprm["TITR_EH0"])
            step = float(self.runprm["TITR_EHD"])

        b2 = [b1 + i * step for i in range(int(self.runprm["TITR_STEPS"]))]
        bounds = round(b1, prec), round(b2[-1], prec)
        self.titr_bounds = bounds

        return bounds

    def __str__(self):
        out = f"rundir: {self.rundir}\nrunprm dict:\n"
        for k in self.runprm:
            out = out + f"{k} : {self.runprm[k]}\n"
        return out


def valid_envs(env1: ENV, env2: ENV) -> tuple:
    """
    Return a 2-tuple:
        (bool, # True: valid,
         error/info message: includes differing keys).
    Step 6 params: excluded from diffing: run1 and run2 are still comparable if only one has run step 6
    """

    s6_keys = {
        "GET_HBOND_MATRIX",
        "GET_HBOND_NETWORK",
        "HBOND_ANG_CUTOFF",
        "HBOND_LOWER_LIMIT",
        "HBOND_UPPER_LIMIT",
        "MS_OUT",  # set in step4 if step6 has run
    }
    # For warning:
    path_keys = {"DELPHI_EXE", "MCCE_HOME"}
    all_keys = set(env1.runprm).union(set(env2.runprm))
    # not always there or nightmare to compare:
    all_keys.remove("EXTRA")
    all_keys.remove("RENAME_RULES")

    # populate diff dict:
    delta = {}
    msg = "OK"

    for k in all_keys:
        if k not in s6_keys:
            v1 = env1.runprm.get(k, None)
            v2 = env2.runprm.get(k, None)
            if v1 is None or v2 is None or v1 != v2:
                delta[k] = [("env1", v1), ("env2", v2)]

    n_delta = len(delta)
    if n_delta == 0:
        return True, msg

    if n_delta > 1:
        if set(delta.keys()) == path_keys:
            msg = "These path keys differ between the two run sets.\n"
            for k in delta:
                msg += f"{k}\t:\t{delta[k]}\n"
            return True, msg

        msg = "WARNING: A parcimonious comparison implies only one differing parameter between the two sets.\n"
        msg = msg + f"WARNING: Found {len(delta)} differing parameters:\n"
        for k in delta:
            msg += f"{k}\t:\t{delta[k]}\n"

    else:  # len == 1: check value of TITR_TYPE
        if delta.get("TITR_TYPE", None) is not None:
            msg = "ERROR: A valid comparison requires the two run sets to have the same TITR_TYPE.\n"
            msg = msg + f"ERROR: TITR_TYPE found: {delta}\n."
            return False, msg

    return True, msg


def get_ref_set(refset_name: str, subcmd: str = SUB1) -> Path:
    if subcmd != SUB1:
        msg = "'parse.e4' is the only reference available & applies to pH titrations setup with `bench_setup pkdb_pdbs"
        logger.error(msg)
        raise ValueError(msg)

    fp = Path(BENCH.BENCH_PH_REFS.joinpath(refset_name))

    return fp


def get_mcce_env_dir(
    bench_dir: str, subcmd: str = SUB1, is_refset: bool = False
) -> Union[Path, None]:
    """Return a path where to get run.prm.record.
    If is_refset is True then subcmd = SUB1 and bench_dir is the reference
    dataset name, e.g. 'parse.e4'."""

    if is_refset and subcmd != SUB1:
        raise ValueError(
            f"The reference dataset {bench_dir} is only available via {SUB1}."
        )

    bdir = Path(bench_dir)
    if is_refset:
        # then bench_dir is the name of a reference dataset
        bdir = get_ref_set(bench_dir, subcmd=subcmd)

    runs_dir = bdir.joinpath(RUNS_DIR)
    for fp in runs_dir.iterdir():
        # use first dir found
        if fp.is_dir():
            return fp

    return


def get_run_env(bench_dir: str, subcmd: str = SUB1, is_refset: bool = False) -> ENV:
    run_dir = get_mcce_env_dir(bench_dir, subcmd=subcmd, is_refset=is_refset)
    env = ENV(run_dir)
    return env


def validate_envs(
    bench_dir1: str, bench_dir2: str, subcmd: str = SUB1, dir2_is_refset: bool = False
) -> tuple:
    """
    Wrapper for fetching a run dir, instantiating the envs,
    and validating them.
    Return the 2-tuple from valid_envs: bool, msg.
    """
    env1 = get_run_env(bench_dir1, subcmd=subcmd)
    env2 = get_run_env(bench_dir2, subcmd=subcmd, is_refset=dir2_is_refset)

    return valid_envs(env1, env2)


def envkey_is_set(k: str) -> bool:
    return os.environ.get(k) is not None


def set_envkey(k: str, v: str):
    os.environ[k] = v
    return


def reset_envkey(k: str):
    if envkey_is_set(k):
        del os.environ[k]
    return


def export_path(
    path_id: str, path_value: str, run_cmd: bool = True
) -> Union[object, str]:
    """Run the command `export <path_id>=<path_value> to the user environment
    if run_cmd is True, else return the command sting.
    """
    cmd = f"export {path_id.upper()}={path_value}"
    if run_cmd:
        return iou.subprocess_run(cmd, check=True)
    return cmd


def set_LDLIB():
    """Export LD_LIBRARY_PATH if not found in the environment.
    This variable holds references to dynamically loaded libraries used by
    PBE solver 'ngpb'.
    """
    if "LD_LIBRARY_PATH" in os.environ:
        return

    ngpb_tools_d = iou.from_json(BENCH.BENCH_NGPB_TOOLS)
    lib_name = list(ngpb_tools_d.keys())[0]
    lib_paths = ngpb_tools_d[lib_name]

    result = export_path(lib_name, lib_paths)
    if isinstance(result, iou.subprocess.CalledProcessError):
        logger.error(result)

    set_envkey(lib_name, lib_paths)

    return


def get_user_env() -> str:
    """Return the env name from sys.prefix.
    Useful to check which is the activated env, if any.
    """
    user_prefix = sys.prefix
    if user_prefix == "/usr":
        # no conda env activated
        env = ""
    else:
        env = Path(user_prefix).name
        if "envs" not in user_prefix:
            if not env.endswith("conda3"):
                # /home/user/*conda3 if base
                raise EnvironmentError(
                    "You appear not to be using conda, which is required for scheduling."
                )
            else:
                env = "base"
    return env


def get_conda_paths() -> tuple:
    """Return the path to conda that is not in an env,
    which is what cron 'sees', presumably.
    Needed to build cmd '. <conda_path>/activate <env>' in crontab.
    TODO: Deprecate?
    """
    conda_path = Path(shutil.which("conda")).parent
    conda = str(conda_path)
    if conda_path.name == "bin":
        # base, done
        return (conda,)
    # else, assume "condabin", reset:
    return str(conda_path.parent.joinpath("bin")), conda


def check_openeye_env_activated():
    try:
        # no error if conda environment where file is found is activated:
        openeye_path = Path(shutil.which("oecheminfo.py"))
    except TypeError:
        logger.error(
            "OpenEye Toolkit not installed or conda environment not activated."
        )
        openeye_path = None

    return openeye_path


def get_conda_exe(raise_err: bool = False) -> str:
    """Return the path to the conda executable, if installed.
    Does not require an activated environment.
    """
    conda_exe = os.environ.get("CONDA_EXE")
    if raise_err:
        if conda_exe is None:
            raise EnvironmentError("Conda is required but not found.")
    return conda_exe


def get_conda_sh(conda_exe: str) -> str:
    return str(Path(conda_exe).parent.parent.joinpath("etc/profile.d/conda.sh"))


def get_openeye_license() -> str:
    oe_license = os.environ.get("OE_LICENSE")
    if oe_license is None:
        raise EnvironmentError("OpenEye license not found.")
    return oe_license


def get_openeye_env_bin(conda_exe: str) -> Path:
    """Find the conda env where oecheminfo.py is found.
    Does not require an activated environment.
    """
    conda_envs = Path(conda_exe).parent.parent.joinpath("envs")
    fp = None
    for fp in conda_envs.glob("*/bin/oecheminfo.py"):
        if fp.exists():
            # use the first one found
            break
    if fp is None:
        raise EnvironmentError("No conda environment found with OpenEye software.")

    return fp.parent


def get_zap_env() -> tuple:
    """Obtain the necessary information to setup the scheduler.
    Return a 5-tuple:
     - The path to the conda executable;
     - The path to the conda env bin;
     - The path to the conda activation script;
     - The path to the OpenEye software licence;
     - The name of the conda environment to be activated.
    """
    conda_exe = get_conda_exe(raise_err=True)
    conda_sh = get_conda_sh(conda_exe)
    conda_env_bin = get_openeye_env_bin(conda_exe)
    oe_env = conda_env_bin.parent.name
    conda_env_bin = str(conda_env_bin)

    try:
        oe_licence = get_openeye_license()
    except EnvironmentError:
        logger.critical(
            (
                "OpenEye software is installed, but the OE_LICENSE key was not found.\n"
                "Add the key to your environment; see: \n"
                "https://docs.eyesopen.com/toolkits/python/quickstart-python/license.html"
            ),
            exc_info=1,
        )

    return conda_exe, conda_env_bin, conda_sh, oe_licence, oe_env
