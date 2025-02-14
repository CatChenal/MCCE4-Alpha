#!/usr/bin/env python

import logging
from pathlib import Path
import shutil
import sys
from mcce4.constants import CANONICAL


logger = logging.getLogger()
logger.setLevel(logging.INFO)


APP_NAME = "tcms"

# console/screen handler
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch_formatter = logging.Formatter("%(levelname)s - %(message)s")
ch.setFormatter(ch_formatter)
# file handlers
flog_format = "[%(asctime)s - %(levelname)s]: %(name)s, %(funcName)s:\n\t%(message)s"
fh = logging.FileHandler(APP_NAME + ".log", encoding="utf-8")
fh.setLevel(logging.INFO)
fh_formatter = logging.Formatter(flog_format)
fh.setFormatter(fh_formatter)
# add to logger
logger.addHandler(ch)
logger.addHandler(fh)


MCCE_EXEC = shutil.which("mcce")
# FIX: Depends on the 'distribution'; error applies to Gunner Lab server installation.
# validate:
if "conda3" in MCCE_EXEC or "MCCE4" not in MCCE_EXEC:
    logger.error("Need to use MCCE4 without any conda environment.")
    sys.exit(1)


OUT_DIR = "tcms"
N_TOP = 5  # default value for 'n_top' option
OUT_TYPES = ["pdb", "gmx", "dcd"]
PARAM_DIR = Path(MCCE_EXEC).parent.parent.joinpath("param")
logger.debug(f"Default param folder: {PARAM_DIR}")


def get_mcce_filepaths(mcce_dir: Path, ph: float) -> tuple:
    """Returns a 3-tuple for head3.lst, step2_out.pdb and
    the msout file for the given ph if all found, else exits.
    """
    ok = True
    out = []
    for fname in ["head3.lst", "step2_out.pdb", "ms_out"]:
        fp = mcce_dir.joinpath(fname)
        ok = ok and fp.exists()
        if not ok:
            logger.error(f"Missing: {fp!r}")
            sys.exit(1)
        if fname == "ms_out":
            msout_fp = fp.joinpath(f"pH{ph:.0f}eH0ms.txt")
            if msout_fp.exists():
                out.append(msout_fp)
            else:
                msout_fp = fp.joinpath(f"pH{ph:.2f}eH0.00ms.txt")
                if msout_fp.exists():
                    out.append(msout_fp)
                else:
                    logger.error(f"Not found: msout file for pH: {ph}")
                    sys.exit(1)
        else:
            out.append(fp)

    return tuple(out)
