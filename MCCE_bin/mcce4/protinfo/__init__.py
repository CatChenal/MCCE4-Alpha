#!/usr/bin/env python

import logging
from pathlib import Path
import shutil


logging.basicConfig(level=logging.INFO, 
                    format="[ %(levelname)s ] %(name)s - %(funcName)s:\n  %(message)s")


USER_MCCE = shutil.which("mcce")
MCCE4 = Path(__file__).parent.parent.parent.parent


CLI_NAME = "protinfo"
# use a dedicated run.log file for step1 run by protinfo:
RUN1_LOG = "run1.log"
RPT = "protinfo.md"  # report name ends with RPT


class ERR:
    """Output error messages from fstrings."""

    def __init__(self, cli_name: str = CLI_NAME):
        self.cmd = cli_name

    def FETCH_EXISTING_FILE(self, pdb: Path) -> str:
        return f"""
The input pdb ({pdb.name}) already exists.
To use it, remove the --fetch flag; run:
    {self.cmd} {pdb.name}
To overwrite it with the downloaded biological assembly; run:
    {self.cmd} {pdb.stem} --fetch
"""

    def MISSING_FETCH_FLAG(self, pdbid: str) -> str:
        return f"""
The input pdb ({pdbid}) seems to be a pdbid. To download its
biological assembly, run this command:
    {self.cmd} {pdbid} --fetch
"""

    def CALL_NOT_IN_FILE_DIR(self) -> str:
        return f"Call {self.cmd} from where the pdb resides.\n"
