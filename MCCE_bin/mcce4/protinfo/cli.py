#!/usr/bin/env python

__doc__ = """
Command line interface for the MCCE4 ProtInfo/protinfo tool, which gathers:
 * Info about the input protein headers (if any) from mcce4.pdbio.
 * Info from MCCE4 step1 run1.log & debug.log when step1 can be run.
Note: The file 'run1.log' is created by protinfo to gather the logging
      information emitted by step1.py.

This is the 'main' module for the cli, which calls the function
that outputs a single protein report: `get_single_pdb_report(args)`.

Options:
 1. pdb (required): a pdb file name or pdbid (assumed valid).
 2. --fetch (False if not used): If 'pdb' is a pdbid and flag is used,
    the biological assembly is downloaded from rcsb.org.
   Step1 options:
    --dry (False if not used): Keep water molecules.
    --noter (False if not used): Do not label terminal residues (for making ftpl).
    -d (4.0): Protein dielectric constant for delphi.
    -u (''): User selected, comma-separated KEY=var pairs from run.prm; e.g.:
           -u HOME_MCCE=/path/to/mcce_home,EXTRA=./extra.tpl.
    -e (mcce): mcce executable location.
    -h, --help  Show this help message and exit.
"""
from argparse import ArgumentParser, RawDescriptionHelpFormatter, Namespace
from functools import partial
import logging
from pathlib import Path
from pprint import pformat
import sys
import time
from typing import Tuple, Union

from mcce4.downloads import get_rcsb_pdb
from mcce4.io_utils import show_elapsed_time
from mcce4.protinfo import CLI_NAME, ERR, RPT, USER_MCCE
from mcce4.protinfo import parsers


logger = logging.getLogger(CLI_NAME)


NO_MCCE_MSG = """The mcce executable was not found.
The protinfo report will not include any information or diagnostics
from MCCE step1.py."""
if USER_MCCE is None:
    logger.warning(NO_MCCE_MSG)


USAGE = """
 >protinfo 1fat --fetch
 >protinfo 1fat.pdb --noter [+ other step1.py options, e.g. --wet, etc.]
"""


err = ERR()


def check_pdb_arg(input_pdb: str) -> Union[Path, str, Tuple[None, str]]:
    """Validate input_pdb str, which can be either a pdb id or a pdb file.
    The tuple output type is used to pass the error message.
    """
    pdb = Path(input_pdb).resolve()
    if not pdb.exists():
        if not pdb.suffix:
            # if no extension, assume pdbid:
            s = len(pdb.stem)
            if s < 4:
                return None, f"Invalid pdbid length: {s}; expected at least 4."

            return input_pdb

        return None, f"File not found: {pdb}"

    if not pdb.parent == Path.cwd():
        return None, err.CALL_NOT_IN_FILE_DIR()

    if pdb.stem == "prot":
        return None, "Input pdb cannot be named 'prot.pdb'."

    if pdb.suffix != ".pdb":
        return None, f"Not a valid extension: {pdb.suffix}"

    return pdb


def validate_pdb_inputs(args: Namespace) -> Union[Path, str]:
    """Validate args.pdb and args.fetch"""
    pdb = check_pdb_arg(args.pdb)
    # pdb is either a Path, a string or a tuple=(None, error_message):
    if isinstance(pdb, tuple):  # error:
        sys.exit(pdb[1])

    if isinstance(pdb, str):
        if not args.fetch:
            sys.exit(err.MISSING_FETCH_FLAG(pdb))
    else:
        if args.fetch:
            sys.exit(err.FETCH_EXISTING_FILE(pdb))

    return pdb


def get_pdb_rpt_climode(
    args: Union[Namespace, dict], do_checks: bool = True, do_fetch: bool = True
) -> Path:
    """Get info and save report for a single pdb.
    This function is called by the protinfo cli, re: 'climode'.
    Expected keys in args: pdb [Path, str], fetch [bool].
    Workflow:
      1. validate_pdb_inputs
      2. downloads.get_rcsb_pdb(pdb): Download the bio-assembly if pdb is id
      3. collect_info in two dicts:
         first: from pdbio.Structure, second: from run1.log protinfo.parser
      4. collect_info_lines from the two dicts
      5. save_report
    Returns:
      The validated pdb file path.
    Note:
      Output files are siloed into a prerun output folder, so that
      protinfo can be safely re-run in a mcce run folder where step 1 & 2 have
      already been run.
    """
    if isinstance(args, dict):
        args = Namespace(**args)

    if do_checks:
        pdb = validate_pdb_inputs(args)
    else:
        pdb = Path(args.pdb)

    if isinstance(pdb, Path):
        run_dir = pdb.parent.resolve()
    else:
        run_dir = Path.cwd()

    out_dir = run_dir.joinpath("prerun")
    out_dir.mkdir(exist_ok=True)

    if do_fetch:
        if not isinstance(pdb, Path):
            pdb = get_rcsb_pdb(pdb)
        else:
            if not pdb.exists():
                pdb = get_rcsb_pdb(pdb.stem)

        if not isinstance(pdb, Path):
            sys.exit("Could not download from rcsb.org.")

    inpdb_fp = out_dir.joinpath(pdb.name)
    if not inpdb_fp.exists():
        inpdb_fp.symlink_to(f"../{pdb.name}")

    logger.info(f"Processing {inpdb_fp.stem.upper()}.")
    prot_d, step1_d = parsers.collect_info(inpdb_fp, args)
    if args.save_dicts:
        prot_d_fp = out_dir.joinpath("prot_d.txt")
        prot_d_fp.write_text("# prot dict \n" + pformat(prot_d) + "\n")
        step1_d_fp = out_dir.joinpath("step1_d.txt")
        step1_d_fp.write_text("# step1 dict \n" + pformat(step1_d) + "\n")

    #logger.info(f"Command line options used:\n{pformat(vars(args))}")

    parsers.write_report(inpdb_fp, prot_d, step1_d)

    return inpdb_fp


# Define alternate 'get_pdb_rpt_climode' function with preset flags
# to False. For use when the pdb already exists, as in the bench app.
get_pdb_rpt = partial(get_pdb_rpt_climode, do_checks=False, do_fetch=False)


def arg_valid_pdb_len(p: str) -> Union[None, str]:
    """Return None if pdb is empty str."""
    if not len(p):
        return None
    if not p.endswith(".pdb") and len(p) < 4:
        # pdbid given for fetching: at least 4 chars
        return None
    return p


def pi_parser():
    p = ArgumentParser(
        prog=f"{CLI_NAME}",
        description=__doc__,
        formatter_class=RawDescriptionHelpFormatter,
        usage=USAGE,
        epilog="""Report issues here:
        https://github.com/GunnerLab/MCCE4/issues""",
    )
    p.add_argument(
        "pdb",
        type=arg_valid_pdb_len,
        help="""A pdb file name (in the current directory) or
        a pdbid (assumed valid).""",
    )
    p.add_argument(
        "--fetch",
        default=False,
        action="store_true",
        help="Download the biological assembly of given pdbid from rcsb.org.",
    )
    p.add_argument(
        "--save_dicts",
        default=False,
        action="store_true",
        help="""Save the structure and step1 dicts as text files. 
        Enables reusing the parsed data independently of the report; default: %(default)s.""",
    )

    s1 = p.add_argument_group("s1", "step1 options")
    # step1.py prot.pdb {dry}{noter}{d}{u}{e}
    s1.add_argument(
        "-d",
        metavar="epsilon",
        type=int,
        default=4,
        help="protein dielectric constant for delphi; %(default)s.",
    )
    s1.add_argument(
        "-e",
        metavar="/path/to/mcce",
        type=str,
        default="mcce",
        help="mcce executable location; default: %(default)s.",
    )
    s1.add_argument(
        "-u",
        metavar="Key=Value",
        type=str,
        default="",
        help="""
        User selected, comma-separated KEY=var pairs from run.prm; default: %(default)s.
        Example: -u HOME_MCCE=/path/to/mcce_home,H2O_SASCUTOFF=0.05,EXTRA=./extra.tpl
        Note: No space after a comma!""",
    )
    s1.add_argument(
        "--wet",
        default=False,
        action="store_true",
        help="Keep water molecules; %(default)s.",
    )
    s1.add_argument(
        "--noter",
        default=False,
        action="store_true",
        help="Do not label terminal residues (for making ftpl); %(default)s.",
    )

    return p


def prot_info_cli(argv=None):
    """Cli 'main' function: produces a Markdown report for a single pdb."""
    # print(f"Start of protinfo. Details logged to {LOG_FILE!s}")

    cli_parser = pi_parser()
    args = cli_parser.parse_args(argv)

    start_t = time.time()
    prerun_pdb = get_pdb_rpt_climode(args)
    show_elapsed_time(start_t=start_t, writer=logger.info,
                      info="Setup, data collection & report creation")

    rpt_fp = Path(f"{prerun_pdb.stem}_{RPT}")
    if rpt_fp.exists():
        print(rpt_fp.read_text())

    return


if __name__ == "__main__":
    prot_info_cli(sys.argv)
