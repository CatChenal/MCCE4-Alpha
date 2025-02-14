#!/usr/bin/env python

"""
Module: cli.py

Command line interface for tcrgms: Tautomeric topN Charge Microsates.
This tool outputs
  - A file listing the topN tautomeric charge microstates, along with
    their related properties: energy (E), net charge (sum_crg), count,
    and occupancy (occ).
  - The <topN> pdb files of each charge state.

Usage:
If called inside a mcce output folder (not a requirement) at pH7 (default):
  tcrgms
Otherwise:
  tcrgms path/to/mcce_dir -ph 4

If called at the command line, the --overwrite flag is set via user prompt.
"""
from argparse import ArgumentParser, Namespace, RawDescriptionHelpFormatter
import logging
import pandas as pd
from pathlib import Path
import shutil
import sys
import time
from typing import Tuple, Union
from mcce4.constants import IONIZABLE_RES
from mcce4.io_utils import mf
from mcce4 import ms_analysis as msa
from mcce4.tcrgms import process_topn as proctop
from mcce4.tcrgms import (
    get_mcce_filepaths,
    PARAM_DIR,
    APP_NAME,
    N_TOP,
    OUT_DIR,
    OUT_TYPES,
)


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


USAGE = """
tcms :: Tautomeric topN Charge Microsates.
This tool outputs:
  - Two files listing the topN tautomeric charge microstates, along with
    their related properties: energy (E), net charge (sum_crg), count,
    and occupancy (occ); one is an extended verion of the other.
  - The <topN> files of each charge state with final format depending on the `output_type` option:
      'pdb' (default): standard pdb format;
      'gmx': Gromacs pdb format;
      'dcd': pseudo trajectory file made up of the topN charge microstates (at most).

Usage:
If called inside a mcce output folder (not a requirement) at pH7 (default) & n_top=5:
  tcms input.pdb

Otherwise:
  tcms path/to/mcce_dir/input.pdb -ph 4 n_top 10
"""


ERR_PROT_USED = """
File 'prot.pdb' is not soft-linked to the actual input pdb.
Either give the path of the actual input pdb at the command line;
Or soft-link it to the actual pdb, e.g.: ln -s <input.pdb> prot.pdb,
then re-run the same command.
"""


def arg_valid_pdb(p: str) -> Path:
    """Return resolved path from the command line."""
    fp = Path(p).resolve()
    if not fp.exists():
        logger.error(mf("Path not found: {!r}", fp))
        sys.exit(1)
    if fp.suffix != ".pdb":
        logger.error("File extension not '.pdb'.")
        sys.exit(1)
    if fp.name == "prot.pdb":
        if fp.is_symlink():
            inpdb_fp = fp.parent.joinpath(fp.readlink())
            logger.info(mf("Input pdb: {}", inpdb_fp))
            return inpdb_fp
        else:
            logger.error(ERR_PROT_USED)
            sys.exit(1)
    return fp


def tc_parser():
    p = ArgumentParser(
        prog=f"{APP_NAME}",
        usage=USAGE,
        formatter_class=RawDescriptionHelpFormatter,
        epilog="""Report issues & feature requests here:
        https://github.com/GunnerLab/MCCE4/issues""",
    )
    p.add_argument(
        "inputpdb_filepath",
        type=arg_valid_pdb,
        help="The pdb filepath in a mcce output dir.",
    )
    p.add_argument(
        "-ph",
        default=7.0,
        type=float,
        help="pH point at which the charge microstates are retrieved; Default: %(default)s.",
    )
    p.add_argument(
        "-n_top",
        default=N_TOP,
        type=int,
        help="Number of most favorable charge microstates to return; Default: %(default)s.",
    )
    p.add_argument(
        "-output_type",
        default="pdb",
        type=str,
        choices=OUT_TYPES,
        help="""
        Format of the charge microstates output file(s):
          'pdb': standard pdb format;
          'gmx': Gromacs, no waters;
          'dcd': files bundled as a pseudo trajectory;
        Default: %(default)s.""",
    )
    p.add_argument(
        "--dry",
        default=False,
        action="store_true",
        help="Output files without waters; Default: %(default)s.",
    )
    p.add_argument(
        "--overwrite",
        default=False,
        action="store_true",
        help="Overwrite existing output files; Default: %(default)s.",
    )

    return p


def get_output_dirname(ph: float, n_top: int) -> str:
    return OUT_DIR + f"_ph{ph:.2f}_top{n_top}"


def display_ionizable_res(df_fp: Path):
    df = pd.read_csv(df_fp, sep="\t")
    msk = df["info"] == "totals"
    mski = df["residues"].str[:3].isin(IONIZABLE_RES)
    new = pd.concat([df[mski], df[msk]]).replace(pd.NA, " ")
    print(new.to_string())


def run_tcrgms(
    args: Union[Namespace, dict], tool_prompt=False
) -> Tuple[list, Tuple[int, int]]:
    """Main function: get topN charge ms with outputs: (at most) n pdbs & 2 result files."""
    logger.info("Starting the processing of Tautomeric Charge MicroStates...\n")

    if isinstance(args, dict):
        args = Namespace(**args)

    mcce_dir = args.inputpdb_filepath.parent
    inpdb = args.inputpdb_filepath.name
    outname = get_output_dirname(args.ph, args.n_top)
    output_dir = mcce_dir.joinpath(outname)
    if output_dir.exists():
        if tool_prompt:
            ans = input(
                (
                    "The %s output folder already exists in %s.\n"
                    "Would you like to proceed and rewrite this folder? y/n "
                    % (outname, str(mcce_dir))
                )
            )
            if ans.strip().lower() == "y":
                args.overwrite = True
            else:
                args.overwrite = False

        if not args.overwrite:
            logger.info(
                "The %s output folder already exists in %s.\n & overwrite option is False. Exiting."
                % (outname, str(mcce_dir))
            )
            sys.exit()

        shutil.rmtree(output_dir)

    output_dir.mkdir()

    args.output_type = args.output_type.lower()
    if args.output_type == "dcd" or args.output_type == "gmx":
        args.output_type = "pdb"
        logger.warning(
            "This output type is not yet implemented. Using standard pdb default."
        )

    # get the dry status from the original run: without this info in the condition,
    # waters would be output even in a dry run (via conformers)
    dry_run = proctop.get_run_dry_opt(mcce_dir)
    dry = dry_run or args.dry or args.output_type == "gmx"

    msg = f"""
    Input options:
      Input pdb: {args.inputpdb_filepath};
      pH point: {args.ph:.2f};
      Top N crgms to process: {args.n_top = };
      Output type: {args.output_type};
      Remove waters? {args.dry};
      Overwrite existing files? {args.overwrite};
    Output folder: {output_dir};
    """
    logger.info(msg)

    h3_fp, step2_fp, msout_fp = get_mcce_filepaths(mcce_dir, args.ph)

    # MSOut, conformer microstates class instantiation
    start_t = time.time()
    Confs = msa.get_confs_collection(h3_fp)
    conformers = Confs.conformers

    mso = msa.MSout(msout_fp)
    elapsed = time.time() - start_t

    N_free = len(mso.free_residues)
    N_fixed = len(mso.fixed_iconfs)

    msg = f"""
    Conformers & microstates info:
    len(conformers) = {Confs.N};
    msa.MSout.free_residues: {N_free:,};
    msa.MSout.fixed_iconfs: {N_fixed:,};
    Total conformer microstates: {mso.N_ms:,};
    Unique conformer microstates: {mso.N_uniq:,}.
    """
    logger.info(msg)
    logger.info(
        f"Loading conformers & msa.MSout took: {elapsed:,.2f} s ({elapsed/60:,.2f} min).\n"
    )

    logger.info("Getting fixed res charges...")
    # fixedres_crg_dict, keys=resid :: {'VALA0002': (4, 0), 'PHEA0003': (5, 0), ...
    # 1st output, net_charge not used (computed in df):
    _, fixedres_crg_dict = msa.fixed_res_crg(
        conformers, mso.fixed_iconfs, no_trailing_underscore=True
    )

    logger.info("Getting topN charge microstates...")
    start_t = time.time()
    CRGMS = msa.Charge_Microstates(mso.microstates, conformers)
    elapsed = time.time() - start_t
    logger.info(
        f"Loading msa.Charge_Microstates took: {elapsed:,.2f} s ({elapsed/60:,.2f} min)."
    )
    logger.info(
        f"Unique charge microstates: {len(CRGMS.confms_by_crg_stateid.keys()):,}\n"
    )

    logger.info("Getting topN dict...")
    topN_lst = CRGMS.get_topN_cms(n_top=args.n_top)
    if not topN_lst:
        logger.info("NO DATA: All the cms have occupancies below the 1% threshold.")
        sys.exit()

    # Will be used along with the 'fixed residues dict' to create a df.
    raw_topN_crg_d = proctop.get_topN_dict(topN_lst, conformers, dry=dry)
    # {'confs': {'NTR01A0001': (0, 0), 'LYS+1A0001': (3, 1), ...
    if not dry:
        topN_crg_d = proctop.process_waters(raw_topN_crg_d)
        del raw_topN_crg_d
    else:
        topN_crg_d = raw_topN_crg_d

    logger.info("Updating dict with tautomers if any...")
    tauto_conflist = proctop.multi_tauto_conflist(PARAM_DIR)
    topN_crg_d = proctop.update_topNcrg_d_tauto(
        topN_crg_d, tauto_conflist, conformers, args.n_top
    )

    logger.info("Outputing non-canonically charged residues, if any...")
    nc_found, nc_fp = proctop.get_non_canonical(topN_crg_d, output_dir)

    logger.info("Creating working dataframe...")
    top_df = proctop.get_crg_vec_df(fixedres_crg_dict, topN_crg_d, mso.N_ms)

    logger.info("Writing the mcce pdbs...")
    # format common keys for final pdb REMARK 250:
    remark_args = {"INPDB": inpdb, "T": f"{mso.T:.2f}", "PH": f"{mso.pH:.2f}"}

    proctop.write_tcrgms_pdbs(
        output_dir,
        inpdb,
        top_df,
        conformers,
        remark_args,
        n_top=args.n_top,
        name_pref=proctop.S2PREF,
    )

    logger.info("Converting the mcce pdbs...")
    proctop.mccepdbs_to_pdbs(output_dir, step2_fp, rm_prefix=proctop.S2PREF)

    logger.info(
        f"File creation over: {args.n_top} tautomeric charge microstates saved to pdbs in {str(output_dir)!r}"
    )

    logger.info("Writing the 'master' and 'user' files from df...")
    # save df to 2 tsv files: 1 'master': all cols; 2: 'user': no 'iconf_' cols:
    proctop.topNdf_to_tsv(output_dir, top_df, n_top=args.n_top)
    tsv_fp = output_dir.joinpath(f"top{args.n_top}_{APP_NAME}.tsv")
    logger.info(
        (
            f"File of at most {args.n_top} tautomeric charge microstates vectors saved as {str(tsv_fp)!r}\n"
            f"\tMaster file is: {outname}/top{args.n_top}_master_{APP_NAME}.tsv."
        )
    )
    if tsv_fp.exists():
        logger.info(f"\nIonizable residues in {str(tsv_fp)!r}:\n")
        display_ionizable_res(tsv_fp)

    if nc_found:
        logger.info("Non canonical residues found: " + str(nc_fp))
        logger.info(nc_fp.read_text())
    return


def tcrgms_cli(argv=None, tool_prompt=False):
    """Cli 'main' function: Obtain the top N charge microstates from the mcce run dir given
    by the input pdb path and outputs, their pdb files."""

    cli_parser = tc_parser()
    args = cli_parser.parse_args(argv)
    run_tcrgms(args, tool_prompt=tool_prompt)

    return


if __name__ == "__main__":
    tcrgms_cli(sys.argv, tool_prompt=True)
