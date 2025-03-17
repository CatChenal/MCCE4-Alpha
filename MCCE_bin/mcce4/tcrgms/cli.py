#!/usr/bin/env python

"""
Module: cli.py

Command line interface for tcrgms: Tautomeric topN Charge Microsates.

This tool outputs:
  - A file listing the topN tautomeric charge microstates, along with
    their related properties (totals): energy (E), net charge (sum_crg), count,
    and occupancy (occ).
  - The <topN> step2_out.pdb, .pqr and .pdb files of each charge state.
  - The output folder name follow this format: tcms_ph<X as float>_top<N>
Note: The number of charge ms returned is at most N.

The function that processes the entire pipeline is `run_tcrgms`.

Usage:
If called inside a mcce output folder (not a requirement) at pH7 (default):
  > tcrgms
Otherwise:
  > tcrgms path/to/mcce_dir -ph 4 -residue_kinds ASP GLU HEM

Note: If called at the command line, the --overwrite flag is set via user prompt.
"""
from argparse import ArgumentParser, Namespace, RawDescriptionHelpFormatter
import logging
import pandas as pd
from pathlib import Path
import shutil
import sys
import time
from typing import Tuple, Union

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
    IONIZABLE_RES,
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
      'gmx': Gromacs pdb format.

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
    """Return path from the command line."""
    fp = Path(p)
    if not fp.exists():
        logger.error(mf("Path not found: {!r}", fp))
        sys.exit(1)
    if fp.suffix != ".pdb":
        logger.error("File extension not '.pdb'.")
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
        help="The pdb filepath in a mcce run dir.",
    )
    p.add_argument(
        "-ph",
        default="7",
        type=str,
        help="pH point (e.g.: 7, 7.5), at which the charge microstates are retrieved; Default: %(default)s.",
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
        Format of the charge microstates output file(s) beside step2_out and pqr:
          'pdb': standard pdb format;
          'gmx': Gromacs, no waters.
        Default: %(default)s.""",
    )
    p.add_argument(
        "-residue_kinds",
        nargs="?",
        default=IONIZABLE_RES,
        help="Filter mcce residues (including cofactors) with these kinds, e.g. ASP GLU HEM; Default: %(default)s.",
    )
    p.add_argument(
        "--wet",
        default=False,
        action="store_true",
        help="Output files with waters; Default: %(default)s.",
    )
    p.add_argument(
        "--overwrite",
        default=False,
        action="store_true",
        help="Overwrite existing output files; Default: %(default)s.",
    )
    p.add_argument(
        "--all_res",
        default=False,
        action="store_true",
        help="Overwrite default filtering by ionizable residues (include all res); Default: %(default)s.",
    )

    return p


def get_output_dirname(ph: Union[int, float], n_top: int) -> str:
    if type(ph) is int:
        return OUT_DIR + f"_ph{ph}_top{n_top}"
    
    return OUT_DIR + f"_ph{ph:.2f}_top{n_top}"


def extend_residue_kinds(res_kinds: list) -> list:
    """Return the IONIZABLES_RES list augmented with new kinds in 'res_kinds'
    in the same order as in IONIZABLES_RES, i.e.:
    acid, base, polar, N term, C term, followed by user-provided cofactors or groups.
    """
    if not res_kinds:
        return []

    userlst = [res.upper() for res in res_kinds]
    ioniz_set = set(IONIZABLE_RES)
    sym_diff = ioniz_set.symmetric_difference(userlst)
    new_res = sym_diff.difference(ioniz_set)
    if new_res:
        return IONIZABLE_RES + sorted(new_res)
    else:
        return IONIZABLE_RES


def display_residue_kinds(df_fp: Path, residue_kinds: list):
    """Display the content of df_fp file filtered by residue_kinds list."""
    df = pd.read_csv(df_fp, sep="\t")
    msk = df["info"] == "totals"
    mski = df["residues"].str[:3].isin(residue_kinds)
    new = pd.concat([df[mski], df[msk]]).replace(pd.NA, " ")
    print(new.to_string() + "\n")


def run_tcrgms(args: Union[Namespace, dict], tool_prompt=False) -> Tuple[list, Tuple[int, int]]:
    """Main function (processing pipeline runner): get topN charge ms with outputs:
    (at most) n step2.pdbs, standard.pdbs, s2.pqr & 2 result files and a summary.txt.
    """
    logger.info("Starting the processing of Tautomeric Charge MicroStates...\n")

    if isinstance(args, dict):
        args = Namespace(**args)

    mcce_dir = args.inputpdb_filepath.parent
    inpdb = args.inputpdb_filepath.name

    # to get the output name with same ph format as input
    if "." in args.ph:
        args.ph = float(args.ph)
    else:
        args.ph = int(args.ph)
    outname = get_output_dirname(args.ph, args.n_top)
    output_dir = mcce_dir.joinpath(outname)
    if output_dir.exists():
        if tool_prompt:
            ans = input(
                (
                    "The %s output folder already exists in %s.\n"
                    "Would you like to proceed and rewrite this folder? y/n " % (outname, str(mcce_dir))
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
    if args.output_type != "pdb":
        args.output_type = "pdb"
        logger.warning("This output type is not yet implemented. Using standard pdb default.")

    # does not matter how the initial run was done, if dry here, HOH will be removed
    dry = not args.wet or args.output_type == "gmx"

    msg = f"""
    Input options:
      Input pdb: {args.inputpdb_filepath};
      pH point: {args.ph:.2f};
      Top N crgms to process: {args.n_top = };
      Output type: {args.output_type};
      Keep waters? {args.wet};
      Overwrite existing files? {args.overwrite};
    Output folder: {output_dir}
    """
    logger.info(msg)

    h3_fp, step2_fp, msout_fp = get_mcce_filepaths(mcce_dir, args.ph)

    # Load conformer microstates, instantiate MSOut class
    start_t = time.time()
    conformers = msa.read_conformers(h3_fp)

    mso = msa.MSout(msout_fp)
    elapsed = time.time() - start_t
    logger.info(f"Loading conformers & msa.CMSout took: {elapsed:,.2f} s ({elapsed/60:,.2f} min).\n")
    logger.info(mso)

    logger.info("Getting fixed res charges...")
    # fixedres_crg_dict, keys=resid :: {'VALA0002': (4, 0), 'PHEA0003': (5, 0), ...
    # The 1st output (net_charge) is not used at this point:
    _, fixedres_crg_dict = msa.fixed_res_crg(conformers, mso.fixed_iconfs, no_trailing_underscore=True)

    logger.info("Getting topN charge microstates...")
    start_t = time.time()

    residue_kinds = args.residue_kinds
    if residue_kinds != IONIZABLE_RES:
        # extend to ionizable res + new kinds, or just sort:
        residue_kinds = extend_residue_kinds(residue_kinds)

    if args.all_res:
        CRGMS = msa.Charge_Microstates(mso.microstates, conformers)
    else:
        CRGMS = msa.Charge_Microstates(mso.microstates, conformers, residue_kinds)

    elapsed = time.time() - start_t
    logger.info(f"Loading msa.Charge_Microstates took: {elapsed:,.2f} s ({elapsed/60:,.2f} min).")
    logger.info(f"Unique charge microstates: {len(CRGMS.confms_by_crg_stateid.keys()):,}\n")

    logger.info("Getting topN dict...")
    topN_lst = CRGMS.get_topN_cms(n_top=args.n_top)
    if not topN_lst:
        logger.info("NO DATA: All the cms have occupancies below the 1% threshold.")
        sys.exit()

    # output dict will be used along with the 'fixed residues dict' to create a df.
    raw_topN_crg_d = proctop.get_topN_dict(topN_lst, conformers, dry=dry)
    # {'confs': {'NTR01A0001': (0, 0), 'LYS+1A0001': (3, 1), ...
    if not dry:
        # remove DM
        topN_crg_d = proctop.process_waters(raw_topN_crg_d)
        del raw_topN_crg_d
    else:
        topN_crg_d = raw_topN_crg_d

    logger.info("Updating dict with tautomers if any...")
    tauto_conflist = proctop.multi_tauto_conflist(PARAM_DIR)
    topN_crg_d = proctop.update_topNcrg_d_tauto(topN_crg_d, tauto_conflist, conformers, args.n_top)

    logger.info("Creating working dataframe...")
    top_df = proctop.get_crg_vec_df(fixedres_crg_dict, topN_crg_d)

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

    logger.info(f"File creation over: {args.n_top} tautomeric charge microstates saved to pdbs in {str(output_dir)!r}")

    logger.info("Writing the 'master' and 'user' files from df...")
    # save df to 2 tsv files: 1 'master': all cols; 2: 'user': no 'iconf_' cols:
    if args.all_res:
        proctop.topNdf_to_tsv(output_dir, top_df, n_top=args.n_top)
    else:
        proctop.topNdf_to_tsv(output_dir, top_df, n_top=args.n_top, res_kinds=residue_kinds)

    tsv_fp = output_dir.joinpath(f"top{args.n_top}_{APP_NAME}.tsv")
    logger.info(
        (
            f"File of at most {args.n_top} tautomeric charge microstates vectors saved as {str(tsv_fp)!r}\n"
            f"\tMaster file is: {outname}/top{args.n_top}_master_{APP_NAME}.tsv."
        )
    )

    logger.info("Writing the summary file...")
    # get summary from the user file:
    proctop.write_summary_file(tsv_fp, args.n_top)

    if args.all_res:  # display a subset:
        logger.info(f"\nIonizable residues in {str(tsv_fp.name)!r}:\n")
        display_residue_kinds(tsv_fp, residue_kinds)

    msg = "Summary file contents:\n\n" + tsv_fp.parent.joinpath("summary.txt").read_text()
    logger.info(msg)

    return


def cli_main(argv=None, tool_prompt=False):
    """Cli 'main' function: Obtain the top N charge microstates from the mcce run dir
    given by the input pdb path and output their step2_out, pqr and pdb files, along
    a summary file, and a .tsv file of cms state vectors & totals.
    """
    cli_parser = tc_parser()
    args = cli_parser.parse_args(argv)
    run_tcrgms(args, tool_prompt=tool_prompt)

    return


if __name__ == "__main__":
    cli_main(sys.argv, tool_prompt=True)
