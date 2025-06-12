#!/usr/bin/env python

"""
Module: io_utils
"""

from functools import wraps, partial
import logging
from pathlib import Path
from time import sleep
from typing import Dict, Union


logger = logging.getLogger(__name__)


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
    def __init__(self, rundir_path: str) -> Dict:
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


def get_path_keys(pdbpath: Path) -> Union[Dict, None]:
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
