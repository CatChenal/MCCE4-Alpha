#!/usr/bin/env python

"""
Module: cli_utils.py

Contains functions used by all cli modules mostly
for args type/validation.
  - cli.py
  - pkanalysis.py
  - comparison.py
  - batch_submit.py
"""

import logging
from pathlib import Path
import shlex
from typing import Union
from mcce4.mcce_benchmark import N_PDBS


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


bench_default_setup_args = {
    "# Options for bench app:": ">>>>",
    "-job_name": "'default_run'",
    "-sentinel_file": "pK.out  # change if skipping steps",
    "#--redo_prerun": "",
    "#--launch": "",
    "#\n# next 2 options are mutually exclusive: choose one": ">>>>",
    "-n_pdbs": N_PDBS,
    "#-pdbids_file": "",
    "#\n# Options common to all mcce steps:": ">>>>",
    "-d": 4.0,
    "-u": "",
    "#\n# Options for step1:": ">>>>",
    "#--dry": "",
    "#--noter": "",
    "#--s1_norun": "",
    "#\n# Options for step2:": "   >>>",
    "-conf_making_level": 1,
    "-ftpl": "",
    "#--s2_norun": "",
    "#\n# Options for step3:": ">>>>",
    "-c": "1 99999",
    "-s": "delphi",
    "-t": "/tmp",
    "-p": 1,
    "-l": "",  # load step3 options from file
    "-ftpl": "",
    "-salt": 0.15,
    "#--refresh": "",
    "#--vdw": "",
    "#--fly": "",
    "#--debug": "",
    "#--s3_norun": "",
    "#\n# Options for step4:": ">>>>",
    "-titr_type": "ph",
    "-i": 0.0,
    "-interval": 1.0,
    "-n": 15,
    "#--ms": "",
    "#--s4_norun": "",
}


def write_default_setup_options_file(fp: Path):
    with open(fp, "w") as opts:
        opts.write(
            "# Note1: 'default_run' name must be changed if changing any default setting.\n"
        )
        opts.write(
            "# Note2: flag options start with '--' (default=False):: uncomment to turn on.\n#\n"
        )
        for k, v in bench_default_setup_args.items():
            if k.startswith("#--") or v == ">>>>":
                opts.write(f"{k}\n")
            else:
                opts.write(f"{k}={v}\n")
    return


def sh_split(arg_line):
    if arg_line.startswith("-c"):
        # shlex cannot handle options with multiple values from
        # a file: too many splits
        arg = "-c" + str([str(i) for i in arg_line.split("=")[1]])
        return arg
    for arg in shlex.split(arg_line, comments="#<"):
        yield arg


def args_int_or_float(x: str) -> Union[int, float]:
    if x.isnumeric():
        return int(x)
    elif "." in x:
        return float(x)
    else:
        raise TypeError(f"{x} is neither a valid integer nor a float.")


def arg_valid_dirpath(p: str) -> Union[None, Path]:
    """Return resolved path from the command line."""
    if not len(p):
        return None
    return Path(p).resolve()


def arg_valid_npdbs(n_pdbs: str) -> int:
    """Return validated number or 1."""
    try:
        n = abs(int(n_pdbs))
        if 0 < n <= N_PDBS:
            return n
        else:
            logger.warning(f"{n_pdbs= } not in (1, {N_PDBS}): reset to 1")
            return 1
    except ValueError:
        logger.warning(f"{n_pdbs= } not in (1, {N_PDBS}): reset to 1")
        return 1


def icase(s: str) -> str:
    """Return the correctly cased value for known strings."""
    if s.lower() == "pk.out":
        return "pK.out"
    elif s.lower() == "step2_out.pdb":
        return "step2_out.pdb"
    else:
        return s


def arg_valid_dir_or_file(p: str) -> Union[None, Path]:
    """Check if resolved path points to a dir or file."""
    if not len(p):
        return None
    pr = Path(p).resolve()
    if pr.is_dir():
        if len(list(pr.glob("*.pdb"))) == 0:
            logger.error(f"No pdbs in folder: {pr}.")
            raise ValueError(f"No pdbs in folder: {pr}.")
        return pr
    if pr.is_file() and pr.exists():
        return pr
    else:
        return None
