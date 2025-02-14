#!/usr/bin/env python

"""
Module: step4.py

This program set up and run MCCE step 4.

Input:
 * energies/
 * head3.lst

Output:
 * fort.38
 * pK.out

Usage examples:

1. Run step 4 with default (pH titration from 0.0 to 14.0, without entropy correction)
    step4.py

2. Write run.prm for step 4, do not actually run step 4. (Dry run)
    step4.py --norun

3. Run step 4 at defined points
    step4.py -i 4.0 -d 1 -n 5

4. Run step 4 with entropy correction and 'ms_out' directory retention
    step4.py --xts --ms

5. Run step 4 using specific mcce executable
    step4.py -e /path/to/mcce

6. Run step 4 using Eh titration
    step4.py -t eh

7. Run step 4 with other customized parameters
    step4.py -u EXTRA=./extra.tpl
"""

import argparse
import os
import subprocess
from mccesteps import export_runprm
from mccesteps import record_runprm
from mccesteps import detect_runprm
from mccesteps import restore_runprm


def write_runprm(args):
    runprm = {}

    path = str(os.path.dirname(os.path.abspath(__file__)))
    base_path = os.path.dirname(path)
    # print(base_path)
    runprm["MCCE_HOME"] = base_path
    runprm["DO_MONTE"] = "t"
    runprm["EXTRA"] = "%s/extra.tpl" % base_path
    runprm["TITR_TYPE"] = args.t.lower()

    if runprm["TITR_TYPE"] == "ph":
        runprm["TITR_PH0"] = args.i
        runprm["TITR_EH0"] = "0.0"
    elif runprm["TITR_TYPE"] == "eh":
        runprm["TITR_PH0"] = "7.0"
        runprm["TITR_EH0"] = args.i

    runprm["TITR_PHD"] = args.d
    runprm["TITR_EHD"] = args.d
    runprm["TITR_STEPS"] = args.n
    runprm["BIG_PAIRWISE"] = "5.0"
    runprm["MONTE_SEED"] = "-1"
    runprm["MONTE_T"] = "298.15"
    runprm["MONTE_FLIPS"] = "3"
    runprm["MONTE_NSTART"] = "100"
    runprm["MONTE_NEQ"] = "300"
    runprm["MONTE_REDUCE"] = "0.001"
    runprm["MONTE_RUNS"] = "6"
    runprm["MONTE_NITER"] = "2000"
    runprm["MONTE_TRACE"] = "50000"
    runprm["NSTATE_MAX"] = "1000000"
    if args.xts:
        runprm["MONTE_TSX"] = "t"
    else:
        runprm["MONTE_TSX"] = "f"
    runprm["MFE_POINT"] = "f"
    runprm["MFE_CUTOFF"] = "-1.0"
    if args.ms:
        runprm["MS_OUT"] = "t"
    else:
        runprm["MS_OUT"] = "f"

    if args.u:
        fields = args.u.split(",")
        for field in fields:
            try:
                key, value = field.split("=")
                runprm[key] = value
            except ValueError:
                print(f"Each component format is 'KEY=VALUE'. Unrecognized: {field}.")

    if args.load_runprm:
        lines = open(args.load_runprm)
        for line in lines:
            entry_str = line.strip().split("#")[0]
            fields = entry_str.split()
            if len(fields) > 1:
                key_str = fields[-1]
                if key_str[0] == "(" and key_str[-1] == ")":
                    key = key_str.strip("()").strip()
                    value = fields[0]
                    runprm[key] = value

    # unconditionally force to run step 4 only in step4.py
    runprm["DO_PREMCCE"] = "f"
    runprm["DO_ROTAMERS"] = "f"
    runprm["DO_ENERGY"] = "f"
    runprm["DO_MONTE"] = "t"

    export_runprm(runprm)
    record_runprm(runprm, "#STEP4")

    return


if __name__ == "__main__":

    # Get the command arguments
    helpmsg = "Run mcce step 4, Monte Carlo sampling to simulate a titration."
    parser = argparse.ArgumentParser(description=helpmsg)

    parser.add_argument(
        "--norun",
        default=False,
        action="store_true",
        help="Create run.prm but do not run step 4; default: %(default)s.",
    )
    parser.add_argument(
        "-i",
        metavar="initial ph/eh",
        type=float,
        default="0.0",
        help="Initial pH/Eh of titration; default: %(default)s.",
    )
    parser.add_argument(
        "-d",
        metavar="interval",
        type=float,
        default="1.0",
        help="titration interval in pJ or mV; default: %(default)s.",
    )
    parser.add_argument(
        "-n",
        metavar="steps",
        type=int,
        default="15",
        help="number of steps of titration; default: %(default)s.",
    )
    parser.add_argument(
        "--xts",
        default=False,
        action="store_true",
        help="Enable entropy correction, default is false; default: %(default)s.",
    )
    parser.add_argument(
        "--ms",
        default=False,
        action="store_true",
        help="Enable microstate output; default: %(default)s.",
    )
    parser.add_argument(
        "-e",
        metavar="/path/to/mcce",
        type=str,
        default="mcce",
        help="mcce executable location; default: %(default)s.",
    )
    parser.add_argument(
        "-t",
        metavar="ph or eh",
        type=str,
        default="ph",
        help="titration type, pH or Eh; default: %(default)s.",
    )
    parser.add_argument(
        "-u",
        metavar="Key=Value",
        type=str,
        default="",
        help="User customized variables; default: %(default)s.",
    )
    parser.add_argument(
        "-load_runprm",
        metavar="prm_file",
        default="",
        help="Load additional run.prm file, overwrite default values.",
    )

    args = parser.parse_args()

    detected = detect_runprm()
    write_runprm(args)
    if not args.norun:
        process = subprocess.Popen([args.e], stdout=subprocess.PIPE)
        for line in process.stdout:
            print(line.decode(), end="")

    if detected:
        restore_runprm()
