#!/usr/bin/env python

import logging
import shutil


# use a dedicated run.log file for step1 run by protinfo:
RUN1_LOG = "run1.log"
RPT = "protinfo.md"  # report name ends with RPT


USER_MCCE = shutil.which("mcce")
NO_MCCE_MSG = """The mcce executable was not found.
The protinfo report will not include any information or diagnostics
from MCCE step1.py."""

LOG_FILE = "protinfo.log"
fh = logging.FileHandler(LOG_FILE, delay=True, mode="w")
fh.setFormatter(logging.Formatter(fmt="[%(asctime)s - %(levelname)s]: %(name)s, %(funcName)s:\n\t%(message)s",
                                  datefmt="%Y-%m-%d %H:%M:%S"))
logging.basicConfig(handlers=[fh], level=logging.DEBUG)

logger = logging.getLogger()
