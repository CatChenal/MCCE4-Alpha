#!/usr/bin/env python

"""
Module: run.py
Functions to launch mcce step1.py
"""

from argparse import Namespace
import logging
from pathlib import Path
import subprocess
from typing import Union

from mcce4.protinfo import RUN1_LOG


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


CUSTOM_S1_SH = """#!/bin/bash

step1.py prot.pdb {wet}{noter}{d}{e}{u}
sleep 1
"""

s1_defaults = {
    "wet": False,
    "noter": False,
    "d": 4,
    "e": "mcce",
    "u": "",
}


def cli_args_to_dict(sh_args: Namespace) -> dict:
    """Only return step1 args."""
    excluded_keys = ["pdb", "fetch"]
    d_args = {k: v for k, v in vars(sh_args).items() if k not in excluded_keys}

    return d_args


def populate_sh_template(job_args: Namespace) -> str:
    """Return the custom template string filled with values."""
    d_args = cli_args_to_dict(job_args)
    # add missing keys needed by template from s1_defaults:
    d_args.update(((k, v) for k, v in s1_defaults.items() if k not in d_args))

    d_all = {}
    # note: trailing spaces needed:
    # special cases:
    v = d_args.pop("wet")
    d_all["wet"] = "" if v else "--dry "
    v = d_args.pop("noter")
    d_all["noter"] = "--noter " if v else ""

    # all remaining options:
    for k in d_args:
        v = d_args.get(k, "")
        if str(v) == str(s1_defaults[k]):
            d_all[k] = ""
        else:
            d_all[k] = f"-{k} {v} "

    return CUSTOM_S1_SH.format(**d_all)


def write_script(dest_dirpath: Path, sh_txt: str):
    """Write an executable bash script to run mcce step1."""
    sh_path = dest_dirpath.joinpath("s1.sh")
    sh_path.write_text(sh_txt)
    # make executable:
    sh_path.chmod(0o755)

    return


def run_step1(pdb_dir: Path) -> Union[None, str]:
    """Run step1 in pdb_dir."""
    result = None
    try:
        proc = subprocess.Popen(
                f"{pdb_dir}/s1.sh",
                cwd=str(pdb_dir),
                close_fds=True,
                stdout=open(f"{pdb_dir}/{RUN1_LOG}", "w"),
                 )
        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            raise subprocess.CalledProcessError(proc.returncode, proc.args, stderr)
        
    except subprocess.CalledProcessError as e:
        result = e.stderr.decode()
        logger.error(f"  {result}")

    finally:
        return result


def already_softlinked(linked_fp: Path, parent_fp: Path) -> bool:
    """Return True if linked_fp name is same as parent_fp name.
    Return False if linked_fp is not a symlink or src and dest differ.
    """
    if not linked_fp.is_symlink():
        return False
    return linked_fp.readlink().name == parent_fp.name


def pdb_is_big(pdb_fp: Path) -> float:
    """Determine whether the pdb size is ge 1 MB.
    Used to adjust sleep time after running step1 so that run1.log
    is completely written.
    """
    size = round(pdb_fp.stat().st_size / (1024 * 1024), 1)
    is_big = size >= 1.
    logger.info(f"File Size is {size} (MB); {is_big = }")

    return size


def do_step1(pdb_fp: Path, args: Namespace) -> Union[None, str]:
    """Main function."""
    run_dir = pdb_fp.parent.resolve()
    result = None

    # setup prot.pdb as soft link:
    prot = run_dir.joinpath("prot.pdb")
    if not prot.exists():
        prot.symlink_to(pdb_fp.name)
    elif not already_softlinked(prot, pdb_fp):
        prot.unlink()
        prot.symlink_to(pdb_fp.name)

    sh_str = populate_sh_template(args)
    write_script(run_dir, sh_str)
    # launch step1:
    result = run_step1(run_dir)

    return result
