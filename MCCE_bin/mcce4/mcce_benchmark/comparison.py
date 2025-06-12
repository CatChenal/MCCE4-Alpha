#!/usr/bin/env python

"""
Cli end point for comparison of two sets of runs. `compare`

Cli parser with options:
  -dir1: path to run set 1
  -dir2: path to run set 2
  -o: path of output folder
  --user_pdbs: Absence means 'pkdb_pdbs'
               Flag enabling the creation of the analysis output files
               in each of the sets if analysis folder not found.
               Thus, this switch enables by-passing the 'intra set'
               analysis files generation using `bench_analyze <sub-command>`.
  --dir2_is_refset: Flag presence indicate dir2 holds the NAME of a reference dataset,
               currently 'parse.e4' for pH titrations.
  The flags are mutually exclusive.

USAGE:
 >bench_compare -dir1 <d1> -dir2 <d2> -diff_point <ph/eh pt> -o <output dir> [+ 2 flags]

Flag options:
1. With 'user_pdbs' flag: the 2 sets were created with 'bench_setup user_pdbs':
   >bench_compare -dir1 <path to set 1> -dir2 <path to set 2> -o <comp> --user_pdbs

2. With 'pkdb_pdbs' flag: the 2 sets were created with 'bench_setup pkdb_pdbs':
   The analysis for these runs outputs pKa stats against the experimental values in pKDBv1.
   >bench_compare -dir1 <path to set 1> -dir2 <path to set 2> -o <comp> --pkdb_pdbs

3. With 'dir2_is_refset' flag: indicates that dir2 is a refset name;
   If used, --pkdb_pdbs must also be present.
   >bench_compare -dir1 <d1> dir2 parse.e4 -o <comp> --dir2_is_refset --pkdb_pdbs
"""

from argparse import ArgumentParser, RawDescriptionHelpFormatter, Namespace
import logging
import pandas as pd
from pathlib import Path
from typing import Union
import sys
from mcce4.io_utils import files_diff
from mcce4.mcce_benchmark import BENCH, ENTRY_POINTS, SUB1, SUB2, FILES, cli_opts
from mcce4.mcce_benchmark import ANALYZE_DIR, RUNS_DIR
from mcce4.mcce_benchmark import cli_utils as clu, io_utils as iou
from mcce4.mcce_benchmark.cleanup import clear_folder
from mcce4.mcce_benchmark import mcce_env as mcenv, pkanalysis, plots


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def get_allsumcrg_diff(
    sc1: Path,
    sc2: Path,
    titr_type: str = "ph",
    titr_point: Union[int, float] = 7,
    save_path: str = None,
    return_df: bool = True,
) -> Union[None, pd.DataFrame]:
    """Return the charge difference (sc1 - sc2) at point 'titr_point' from two
    collated sum_crg.out files.
    Args:
      - sc1 (Path): Path to the first all_smcrg.out file.
      - sc2 (Path): Path to the second (control) all_smcrg.out file.
      - titr_type (str, "ph"): Titration type.
      - titr_point (int or float, 7): Titration point for filtered output.
      - save_path (str, None): Save the unfiltered df with given path if not None.
      - return_df (bool, True): Default is True because in comparison.py
        the plots are created right after the diffing operation.
    Note:
      Due to current header formating for sum_crg.out, titr_point should be int
      for pH titrations (even though mcce allows factional pH points), if not e.g.,
      7.5 is recast to 8.
    """
    if save_path is not None:
        diff_fp = Path(save_path).parent
    else:
        diff_fp = Path.cwd()
    sc_diff = files_diff(sc1, sc2, out_dir=diff_fp, return_df=True, collated=True)
    if sc_diff is None:
        logger.error("Could not diff all_sum_crg.out files.")
        return

    titr = sc_diff.columns[0][-2:]
    titr_type = titr_type.lower()
    if titr.lower() != titr_type:
        logger.error("Mismatched titr types in sum_crg files viz runprm.")
        raise TypeError("Mismatched titr types in sum_crg files viz runprm.")

    prec = 0
    if titr.lower() != "ph":
        prec = 1
    titr_point_recast = f"{titr_point:.{prec}f}"
    if titr_point != titr_point_recast:
        logger.warning(
            f"titr_point {titr_point} recast to {titr_point_recast} as per titr type precision."
        )

    # split the pdb:resid col:
    cols = sc_diff.columns.to_list()
    col1, col2 = cols[0].split(":")
    sc_diff[[col1, col2]] = sc_diff[cols[0]].str.split(pat=":", expand=True)
    sc_diff.drop(columns=cols[0], inplace=True)
    sc_diff = sc_diff[[col1, col2] + cols[1:]]
    # save entire df:
    if save_path is not None:
        Path(save_path).write_text(sc_diff.to_string(index=False) + "\n")

    missing_point = False
    try:
        sc_diff = sc_diff[[col1, col2, str(titr_point_recast)]]
    except KeyError:
        logger.error(
            (
                f"Titration point: {titr_point_recast} not found.\n"
                f"Diffed all_sum_crg columns are: {sc_diff.columns.to_list()}."
            )
        )
        missing_point = True

    if return_df and not missing_point:
        return sc_diff
    return


def compare_runs(args: Union[dict, Namespace]):

    if isinstance(args, dict):
        args = Namespace(**args)

    if args.dir2_is_refset and args.subcmd != SUB1:
        logger.error(f"The reference dataset {args.dir2} is only available via {SUB1}.")
        raise TypeError(
            f"The reference dataset {args.dir2} is only available via {SUB1}."
        )

    # check: analysis was run?
    analyze1 = args.dir1.joinpath(ANALYZE_DIR)
    if not analyze1.exists():
        logger.error(f"Analysis folder missing for {args.dir1}.")
        raise ValueError(f"Analysis folder missing for {args.dir1}.")

    if args.dir2_is_refset:
        ref_path = mcenv.get_ref_set(args.dir2)
        analyze2 = ref_path.joinpath(ANALYZE_DIR)
    else:
        args.dir2 = Path(args.dir2)
        analyze2 = args.dir2.joinpath(ANALYZE_DIR)
    if not analyze2.exists():
        logger.error(f"Analysis folder missing for {args.dir2}.")
        raise ValueError(f"Analysis folder missing for {args.dir2}.")

    ok, msg = mcenv.validate_envs(
        args.dir1, args.dir2, subcmd=args.subcmd, dir2_is_refset=args.dir2_is_refset
    )
    if not ok:
        logger.error(f"Runs failed validation:\n{msg}")
        raise ValueError(f"Runs failed validation:\n{msg}")
    if ok and msg != "OK":
        logger.warning(f"Runs validation warning:\n{msg}")

    # conf making levels:
    lev1, lev2 = 1, 1
    lev1 = iou.get_setup_args_vals(args.dir1, "conf_making_level")
    if lev1:
        lev1 = lev1[0]
    lev2 = iou.get_setup_args_vals(args.dir2, "conf_making_level")
    if lev2:
        lev2 = lev2[0]

    # get levels for plot:
    if lev1 != lev2:
        logger.warning(f"Differing conformer making levels: {lev1 = }, {lev2 = }")
        level = lev1, lev2
    else:
        level = lev1

    out_dir = Path(args.o)
    if not out_dir.exists():
        out_dir.mkdir()
    else:
        clear_folder(out_dir)
        logger.info(f"Cleared comparison output folder: {out_dir}")

    # 1. get collated sum_crg.out diff:
    logger.info("Calculating sum_crg diff file.")
    sc1 = analyze1.joinpath(FILES.ALL_SUMCRG.value)
    sc2 = analyze2.joinpath(FILES.ALL_SUMCRG.value)

    titr_type = iou.get_setup_args_vals(args.dir1, "titr_type")
    if titr_type:
        if titr_type[0] is not None:
            titr_type = titr_type[0]
        else:
            titr_type = "ph"
            logger.warning("titr_type not found in setup_args; reset to 'ph'.")
    else:
        titr_type = "ph"
        logger.warning("titr_type not found in setup_args; reset to 'ph'.")

    txt_fp = out_dir.joinpath(FILES.ALL_SUMCRG_DIFF.value)
    # default: collated sum_crg diff at ph 7
    all_sumcrg_df = get_allsumcrg_diff(
        sc1,
        sc2,
        titr_type=titr_type,
        titr_point=args.diff_point,
        save_path=txt_fp,
        return_df=True,
    )
    if all_sumcrg_df is None:
        logger.warning(
            (
                "Could not diff the collated sumcrg files, "
                "or the titration point was not found."
            )
        )
    else:
        sets_names = (analyze1.parent.name, analyze2.parent.name)
        for pdbid in all_sumcrg_df.Diff.unique():
            sumcrgfig_fp = out_dir.joinpath(f"sumcrg_diff_{pdbid}.png")
            plots.plot_sumcrg_diff(
                all_sumcrg_df, pdbid, sets_names, out_fp=sumcrgfig_fp
            )

    # 2. get pkas from serialized dicts from all_pks1, all_pks2 & match pkas:
    logger.info("Matching the pkas and saving pk values & diff to file.")
    d1 = iou.from_pickle(analyze1.joinpath(FILES.JOB_PKAS_PKL.value))
    d2 = iou.from_pickle(analyze2.joinpath(FILES.JOB_PKAS_PKL.value))
    matched_pkas = pkanalysis.match_pks(d1, d2)

    # 'pretty' file:
    matched_fp = out_dir.joinpath(FILES.MATCHED_PKAS_TXT.value)
    pkanalysis.matched_pks_txt(
        matched_fp,
        matched_pkas,
        sets_names=(analyze1.parent.name, analyze2.parent.name),
    )

    # 3. matched pkas stats
    logger.info("Calculating the pkas & residues stats.")
    matched_df = iou.matched_pks_to_df(matched_fp)
    d_stats, res_stats_df = pkanalysis.matched_pks_stats(
        matched_df, titr_type=titr_type, level=level
    )
    logger.info(d_stats["report"])
    # dict to pickle:
    iou.to_pickle(d_stats, out_dir.joinpath(FILES.MATCHED_PKAS_STATS_PKL.value))

    # save res_stats_df:
    res_stats_fp = out_dir.joinpath(FILES.RESIDUES_STATS.value)
    res_stats_fp.write_text(
        res_stats_df.to_string(index=False, float_format="{:.2f}".format) + "\n"
    )
    # save res_stats_df to dict, to be reused in res plot:
    res_stats_d = res_stats_df.set_index("res").to_dict(orient="index")
    iou.to_pickle(res_stats_d, out_dir.joinpath(FILES.RESIDUES_STATS_PKL.value))

    if d_stats["fit"][0] is None:
        logger.info("Data could not be fitted: no plot generated.")
    else:
        logger.info("Plotting pkas fit -> pic.")
        save_to = out_dir.joinpath(FILES.FIG_FIT_ALLPKS.value)
        plots.plot_pkas_fit(matched_fp, d_stats, out_fp=save_to, comparison=True)

    # 4. get figure for matched residues analysis:
    logger.info("Plotting matched residues analysis -> pic.")
    plots.plot_res_analysis(
        matched_fp,
        res_stats_d,
        level=level,
        out_fp=out_dir.joinpath(FILES.FIG_FIT_PER_RES.value),
    )

    # 5. Outlier residue analysis
    logger.info("Doing outlier residue analysis.")
    # all res:
    _ = pkanalysis.res_outlier_count(matched_fp, grp_by="res")
    # by res type:
    _ = pkanalysis.res_outlier_count(matched_fp, grp_by="resid")

    return


# ........................................................................
CLI_NAME = ENTRY_POINTS["compare"]  # as per pyproject.toml entry point

DESC = """
Description:
Compare two sets of runs, ~ A/B testing; (convention: B is 'reference', i.e.: A - B).
IMPORTANT: 'bench_analyze' must be run on each set before comparing them.

Options:
  -dir1: path to run set 1
  -dir2: path to run set 2
  -o:    path to output dir
  --user_pdbs: Absence means 'pkdb_pdbs'
               Flag enabling the creation of the analysis output files
               in each of the sets if analysis folder not found.
               Thus, this switch enables by-passing the 'intra set'
               analysis files generation using `bench_analyze <sub-command>`.

  --dir2_is_refset: Flag presence indicates that dir2 value is a refset name;
               If used, --user_pdbs must NOT be present.

  (mce) >bench_compare -dir1 <d1> -dir2 parse.e4 --dir2_is_refset -o output/dir/path
"""


def compare_parser():
    """Cli arguments parser for use in benchmark comparison."""

    p = ArgumentParser(
        prog=f"{CLI_NAME} ",
        description=DESC,
        formatter_class=RawDescriptionHelpFormatter,
        epilog="""
        Post an issue for all errors and feature requests at:
        https://github.com/CatChenal/MCCE4/issues
        """,
    )

    p.add_argument(
        "-dir1",
        required=True,
        type=clu.arg_valid_dirpath,
        help="""Path to run set 1.""",
    )
    p.add_argument(
        "-dir2",
        required=True,
        type=str,
        help="""Path to run set 2 OR name of a reference set.""",
    )
    p.add_argument(
        "-diff_point",
        default=7.0,
        type=float,
        help="""
        Titration point at which the sum_crg diff is calculated;
        default for pH titrations: %(default)s.""",
    )
    p.add_argument(
        "-o",
        required=True,
        type=clu.arg_valid_dirpath,
        help="""Path to comparison results folder.""",
    )
    # cannot have --user_pdbs & --dir2_is_refset together:
    mutex = p.add_mutually_exclusive_group()
    mutex.add_argument(
        "--user_pdbs",
        default=False,
        action="store_true",
        help="""
        Flag enabling the creation of the analysis output files on
        runs setup with user's pdbs in each of the sets if their analysis
        folder is not found.
        Thus, this switch enables the by-passing of the
        bench_analyze <sub-command> step.
        """,
    )
    mutex.add_argument(
        "--dir2_is_refset",
        default=False,
        action="store_true",
        help="""
        Flag presence indicate dir2 holds the NAME of a reference dataset, currently 'parse.e4'.
        """,
    )

    return p


def compare_cli(argv=None):
    """Command line interface for MCCE benchmarking comparison entry point.
    """
    cli_parser = compare_parser()
    args = cli_parser.parse_args(argv)

    cli_opts.cli_name = CLI_NAME
    cli_opts.all = vars(args)
    logger.info(cli_opts)

    # OK to compare?
    for d in [args.dir1, args.dir2]:
        if isinstance(d, Path):
            bench = iou.Pathok(d)
            pct = iou.pct_completed(bench.joinpath(RUNS_DIR, BENCH.Q_BOOK))
            if pct < 1.0:
                logger.info(
                    f"Runs not 100% processed in {d}, try again later; {pct = :.1%}"
                )
                return

    # Add subcmd to args to hold inferred comp type:
    subcmd = SUB1
    if args.user_pdbs:
        subcmd = SUB2
    setattr(args, "subcmd", subcmd)
    env = mcenv.get_run_env(Path(args.dir1), subcmd=subcmd)
    # Add titr to args set:
    setattr(args, "titr", env.runprm["TITR_TYPE"])

    compare_runs(args)

    return


if __name__ == "__main__":
    compare_cli(sys.argv[1:])
