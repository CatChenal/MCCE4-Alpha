#!/usr/bin/env python

"""
Module: custom_sh.py

Functions for building a custom script when cli args are not all defaults.

Notes about step3.py in scritps:
* For custom scripts:
  The 'pound' variable (='' or '#') is a workaround for the disappearance of the `--norun`
  option, in the new step3.py. See:
  https://github.com/GunnerLab/MCCE4/issues/211

CHANGELOG:

* 11-24-24:
In the benchmark app, the `-u` option is used to obtain the names of the customized extra.tpl
or name.txt files in the bench_dir so that the app can automatically link them into each of the
proteins folders in bench_dir/runs.
If other key besides EXTRA or RENAME_RULES is given, a custom run.prm is created, then linked.

* 11-10-24:
 - Removed the obsolete -e option
"""

from argparse import Namespace
from enum import Enum
import logging
from pathlib import Path
from pprint import pformat
from mcce4.mcce_benchmark import RUNS_DIR
from mcce4.mcce_benchmark.io_utils import make_executable


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


# To test script is run inside a PDB folder:
RUN_SH_TEST_ECHO = """#!/bin/bash

echo "Using RUN_SH_TEST_ECHO as script: $PWD"
"""

# To test submit_script without running anything:
RUN_SH_NORUN = """#!/bin/bash

step1.py prot.pdb --norun
"""


# pseudo default, bypassing pbe solver
RUN_SH_PSEUDO = """#!/bin/bash

step1.py prot.pdb --dry
step2.py
step3.py --refresh
step4.py --xts

sleep 10
"""

SH_TEMPLATE = """#!/bin/bash

step1.py prot.pdb {wet}{noter}{d}{s1_norun}{u}{load_runprm}
step2.py {conf_making_level}{ftpl}{d}{s2_norun}{u}{load_runprm}
{pound}step3.py {d}{c}{s}{t}{p}{ftpl}{salt}{refresh}{vdw_relax}{fly}{skip_pb}{debug}{l}{u}{load_runprm}
step4.py --xts {titr_type}{i}{interval}{n}{ms}{s4_norun}{u}{load_runprm}

sleep 10
"""


class ScriptChoices(str, Enum):
    TEST_ECHO = RUN_SH_TEST_ECHO
    NORUN = RUN_SH_NORUN
    CUSTOM = SH_TEMPLATE


opt_flags = [
    "s1_norun",
    "s2_norun",
    "s3_norun",
    "s4_norun",
    "wet",
    "noter",
    "refresh",
    "fly",
    "debug",
    "skip_pb",
    "ms",
]


# maping of bench cli option names to mcce steps options:
cli_to_mcce_opt = {
    "conf_making_level": "l",
    "interval": "d",
    "s1_norun": "norun",
    "s2_norun": "norun",
    # "s3_norun": "norun",  # removed in new step3.py
    "s4_norun": "norun",
    "titr_type": "t",
    "wet": "dry",
}

# cli defaults per step; convenience dict to see where opts go:
defaults_per_step = {
    "s1": {
        "d": 4,
        "wet": False,
        "noter": False,
        "s1_norun": False,
        "u": "",
        "load_runprm": "",
    },
    "s2": {
        "d": 4,
        "conf_making_level": 1,
        "ftpl": "",
        "s2_norun": False,
        "u": "",
        "load_runprm": "",
    },
    "s3": {
        "d": 4,
        "c": [1, 99999],
        "s": "delphi",
        "t": "/tmp",
        "p": 1,
        "l": "",
        "ftpl": "",
        "salt": 0.15,
        "refresh": False,
        "vdw_relax": 0,
        "fly": False,
        "debug": False,
        "skip_pb": False,
        "s3_norun": False,
        "u": "",
        "load_runprm": "",
    },
    "s4": {
        "d": 4,
        "titr_type": "ph",
        "i": 0.0,
        "interval": 1.0,
        "n": 15,
        "ms": False,
        "s4_norun": False,
        "u": "",
        "load_runprm": "",
    },
}
# combined:
all_default_opts = {}
for S in defaults_per_step:
    all_default_opts.update(((k, v) for k, v in defaults_per_step[S].items()))


def cli_args_to_dict(sh_args: Namespace) -> dict:
    """Only return mcce steps args."""

    excluded_keys = [
        "subparser_name",
        "bench_dir",
        "n_pdbs",
        "pdbids_file",
        "redo_prerun",
        "pdbs_list",
        "sentinel_file",
        "job_name",
        "launch",
        "func",
        "help",
    ]
    d_args = {k: v for k, v in vars(sh_args).items() if k not in excluded_keys}
    return d_args


def all_opts_are_defaults(args: Namespace, debug: bool = False) -> bool:
    """Return True if sh_args are default in all the steps,
    else return False.
    Purpose: determine whether to write a custom script or
             use the default one.
    """
    # d_sh_args holds mcce steps options only
    d_sh_args = cli_args_to_dict(args)
    if debug:
        logger.info(pformat(d_sh_args, indent=2))

    is_default = True
    for opt in d_sh_args:
        is_default = is_default and d_sh_args[opt] == all_default_opts.get(opt)
        if debug:
            logger.info(
                f"\n{opt= }: {is_default= }\n\t{d_sh_args[opt]= } == {all_default_opts.get(opt)= }"
            )
        if not is_default:  # done
            break

    return is_default


def populate_custom_template(job_args: Namespace) -> str:
    """Return the custom template string filled with appropriate values."""

    if not vars(job_args):
        msg = "job_args cannot be None when using the CUSTOM template."
        logger.error(msg)
        raise ValueError(msg)

    d_args = cli_args_to_dict(job_args)
    # add missing keys needed by template from default opts:
    d_args.update(((k, v) for k, v in all_default_opts.items() if k not in d_args))

    d_all = {}
    # 'pound' only applies to step3:
    d_all["pound"] = ""

    # note: trailing spaces needed
    # flag options: don't output if False, except for 'wet'
    for fo in opt_flags:
        v = d_args.pop(fo)
        if fo == "wet":
            d_all["wet"] = "" if v else "--dry "
        elif fo.endswith("norun"):
            if fo == "s3_norun" and v:
                d_all["pound"] = "#"
            else:
                d_all[fo] = "--norun " if v else ""
        else:
            d_all[fo] = f"--{fo} " if v else ""

    # all remaining options with values:
    for k in d_args:
        v = d_args.get(k, "")
        if str(v) == str(all_default_opts[k]):
            d_all[k] = ""
        else:
            d_all[k] = f"-{cli_to_mcce_opt.get(k, k)} {v} "

    return ScriptChoices.CUSTOM.value.format(**d_all)


def write_run_script_from_template(
    bench_dir: str,
    job_name: str,
    script_template: ScriptChoices = ScriptChoices.CUSTOM,
    job_args: Namespace = None,
) -> None:
    """
    Write a custom shell script in <bench_dir>/runs/ to submit steps 1-4 when
    script_template is CUSTOM, or perform tests otherwise. job_args can be None for
    templates other than CUSTOM.
    Delete a pre-exisitng script with the same name.

    Args:
    script_template (ScriptChoices enum): one of TEST_ECHO, NORUN, CUSTOM (default)
    """
    bench_dir = Path(bench_dir)
    runs_dir = bench_dir.joinpath(RUNS_DIR)
    if not runs_dir.exists():
        msg = f"{bench_dir} does not have a {RUNS_DIR} subfolder!"
        logger.critical(msg)
        raise FileNotFoundError(msg)

    if not job_name:  # empty str"
        logger.error("'job_name' is required to have a value.")
        raise ValueError("'job_name' is required to have a value.")

    if job_name == "default_run":
        msg = "'job_name' cannot be 'default_run' for a custom script."
        logger.error(msg)
        raise ValueError(msg)

    if script_template is ScriptChoices.CUSTOM:
        if job_args is None:
            msg = "job_args cannot be None when using the CUSTOM template."
            logger.error(msg)
            raise ValueError(msg)

        sh_text = populate_custom_template(job_args)
    else:
        sh_text = script_template.value

    sh_path = runs_dir.joinpath(f"{job_name}.sh")
    logger.info(f"Creating custom script: {str(sh_path)}")
    if sh_path.exists():
        sh_path.unlink()
    with open(sh_path, "w") as fh:
        fh.write(sh_text)

    # make script executable:
    make_executable(sh_path)

    return
