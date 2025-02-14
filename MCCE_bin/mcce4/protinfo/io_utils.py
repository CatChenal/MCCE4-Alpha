#!/usr/bin/env python

"""
Module: io_utils
"""

from functools import wraps, partial
import logging
from pathlib import Path
import subprocess
import shutil
from time import sleep
from typing import Tuple, Union
import requests


# reset libs logging to higher level so that no unnecessary message is logged
# when this module is used.
logging.getLogger("requests").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.WARNING)


logger = logging.getLogger(__name__)
#logger.setLevel(logging.INFO)


def retry(func=None, *, times=6):
    """Retry decorator."""
    if func is None:
        return partial(retry, times=times)

    @wraps(func)
    def wrapper(*args, **kwargs):
        attempt = 0
        while attempt < times:
            try:
                return func(*args, **kwargs)
            except Exception as ex:
                attempt += 1
                logger.debug(f"Exception {func}: {ex} (attempt: {attempt})")
        logger.error(f"No success after {attempt} times")
        return None

    return wrapper


class ENV:
    def __init__(self, rundir_path: str) -> dict:
        self.rundir = Path(rundir_path)
        self.runprm = {}
        # populate self.runprm dict:
        self.load_runprm()

    @retry(times=3)
    def load_runprm(self):
        fp = Path(self.rundir.joinpath("run.prm.record"))
        if not fp.exists():
            sleep(1)
            raise FileNotFoundError(f"Not found: run.prm.record in {self.rundir}")

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
                        value = round(float(fields[0]), 1)
                    else:
                        value = fields[0]
                self.runprm[key] = value

        return


def get_path_keys(pdbpath: Path) -> Union[dict, None]:
    """
    Return path keys from run.prm.record as a dict with keys:
    ["topology files","renaming file"].
    """
    rundir = pdbpath.parent
    env = ENV(rundir)
    return {
        "topologies": env.runprm["MCCE_HOME"] + "/param",
        "renaming file": env.runprm["RENAME_RULES"],
    }


def subprocess_run(
    cmd: str,
    capture_output: bool = True,
    check: bool = False,
    text: bool = True,
    shell: bool = True,
) -> Union[subprocess.CompletedProcess, subprocess.CalledProcessError]:
    """Wraps subprocess.run. Return CompletedProcess or err obj."""

    try:
        data = subprocess.run(
            cmd, capture_output=capture_output, check=check, text=text, shell=shell
        )
    except subprocess.CalledProcessError as e:
        data = e

    return data


def make_executable(sh_path: str) -> None:
    """Alternative to os.chmod(sh_path, stat.S_IXUSR): permission denied."""

    sh_path = Path(sh_path)
    cmd = f"chmod +x {str(sh_path)}"

    try:
        subprocess_run(cmd, capture_output=False, check=True)
    except subprocess.CalledProcessError:
        logger.exception("Error in subprocess cmd 'chmod +x'")
        raise

    return


def rcsb_download(pdb_fname: str) -> requests.Response:
    url_rscb = "https://files.rcsb.org/download/" + pdb_fname
    return requests.get(url_rscb, allow_redirects=True)


def rcsb_download_header(pdb_fname: str) -> requests.Response:
    url_rscb = "https://files.rcsb.org/header/" + pdb_fname
    return requests.get(url_rscb, allow_redirects=True,
                        headers = {"accept-encoding": "identity"})


def get_rcsb_pdb(pdbid: str) -> Union[Path, Tuple[None, str]]:
    """Given a pdb id, download the pdb file containing
    the biological assembly from rcsb.org.
    The file is downloaded with a pdb extension.
    """
    pdbid = pdbid.lower()
    pdb = pdbid + ".pdb"  # final pdb filename
    pdb1 = pdbid + ".pdb1"
    # Removed cif file download as this format cannot be converted 
    # since biopython was removed from dependencies.
    # The bioassembly and header are downloaded because the bioassembly file
    # will not have the remarks section, which mcce4.pdbio.py parses.

    content = None
    # list of bool to identify which pdb was saved:
    which_ba = [False, False]  # 0:  bio assembly, 1: pdb standard

    # try bio assembly:
    r1 = rcsb_download(pdb1)
    if r1.status_code < 400:
        which_ba[0] = True
        with open(pdb1, "wb") as fo:
            fo.write(r1.content)
    else:
        logger.warning(f"Could not download the bio assembly: {r1.reason}")
    
    if not which_ba[0]:
        # try standard format
        r2 = rcsb_download(pdb)
        if r2.status_code < 400:
            which_ba[1] = True
            with open(pdb, "wb") as fo:
                fo.write(r2.content)
        else:
            logger.warning(f"Could not download the pdb file: {r2.reason}")

        if not which_ba[1]:  # both False;
            return None, "Error: Could neither download the bio assembly or pdb file."
    else:
        # get the header
        r2 = rcsb_download_header(pdb)
        if r2.status_code < 400:
            with open(pdb, "w") as fo:
                fo.write(r2.text)
        else:
            logger.warning(f"Could not download the pdb file header, bioassembly with partial header used: {r2.reason}")
            # use bioassembly as pdb file:
            shutil.move(str(pdb1), str(pdb))
            return Path(pdb).resolve()

    cmd = ("sed -n '/MODEL/,$ p' " + pdb1 + " >> " + pdb + "; rm " + pdb1)

    o = subprocess_run(cmd, capture_output = True, check = True)
    if isinstance(o, subprocess.CalledProcessError):
        msg = "Could not prep the pdb file with cmd:\n" + cmd
        logger.critical(msg, exc_info=1)

    print("Download completed.")

    return Path(pdb).resolve()
