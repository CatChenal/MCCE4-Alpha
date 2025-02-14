#!/usr/bin/env python

"""
Module: step6.py

This program set up and run MCCE step 6, calculate Hydrogen Bond network.

Usage examples:

1. Run step 6 with default
    step6.py

2. Write run.prm for step 6, do not actually run step 6. (Dry run)
    step4.py --norun

3. Run step 6 with other customized parameters
    step6.py -u EXTRA=./extra.tpl
"""

import argparse
import os
import subprocess
from mccesteps import export_runprm
from mccesteps import record_runprm


def write_runprm(args):
    runprm = {}

    path = str(os.path.dirname(os.path.abspath(__file__)))
    base_path = os.path.dirname(path)
    # print(base_path)
    runprm["MCCE_HOME"] = base_path
    runprm["DO_ANALYSIS"] = "t"
    runprm["GET_HBOND_MATRIX"] = "t"
    runprm["EXTRA"] = "%s/extra.tpl" % base_path
    runprm["HBOND_LOWER_LIMIT"] = "1.2"
    runprm["HBOND_UPPER_LIMIT"] = "3.2"
    runprm["HBOND_ANG_CUTOFF"] = "90.0"
    runprm["GET_HBOND_NETWORK"] = "t"

    if args.u:
        fields = args.u.split(",")
        for field in fields:
            try:
                key, value = field.split("=")
                runprm[key] = value
            except ValueError:
                print(f"Each component format is 'KEY=VALUE'. Unrecognized: {field}.")

    export_runprm(runprm)
    record_runprm(runprm, "#STEP6")

    return


if __name__ == "__main__":

    # Get the command arguments
    helpmsg = (
        "Run step 6, hydrogen bond analysis, requires microstate output from step 4."
    )
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument(
        "--norun",
        default=False,
        action="store_true",
        help="Create run.prm but do not run step 6; default: %(default)s.",
    )
    parser.add_argument(
        "-e",
        metavar="/path/to/mcce",
        type=str,
        default="mcce",
        help="mcce executable location; default: %(default)s.",
    )
    parser.add_argument(
        "-u",
        metavar="Key=Value",
        default="",
        help="User customized variables; default: %(default)s.",
    )

    args = parser.parse_args()
    # print(args)

    write_runprm(args)
    if not args.norun:
        process = subprocess.Popen([args.e], stdout=subprocess.PIPE)
        for line in process.stdout:
            print(line.decode(), end="")
