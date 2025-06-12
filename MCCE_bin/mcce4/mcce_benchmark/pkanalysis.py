#!/usr/bin/env python

"""
Cli end point for comparison of two sets of runs. `compare`.

Cli parser with 2 sub-commands with same options:
  -dir1 <path/to/set/of/runs1>; must exists
  -dir2 <path/to/set/of/runs2>; must exists

1. pkdb_pdbs (SUB1)
2. user_pdbs (SUB2)
"""

from argparse import ArgumentParser, RawDescriptionHelpFormatter, Namespace
from functools import partial
import logging
from pathlib import Path
from typing import Tuple, Union
import sys
import numpy as np
import pandas as pd
from mcce4.io_utils import mccepdbline_positions, mcfile2df, mf
from mcce4.mcce_benchmark import (
    BENCH,
    ENTRY_POINTS,
    SUB1,
    SUB2,
    FILES,
    cli_opts,
    ANALYZE_DIR,
    RUNS_DIR,
)
from mcce4.mcce_benchmark import io_utils as iou, mcce_env as mcenv, plots
from mcce4.mcce_benchmark.cleanup import clear_folder


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


round2 = partial(round, ndigits=2)


DELTA = "\u0394"


def get_mcce_version(runs_dir: str) -> None:
    """MCCE version(s) from runs/*/run.log files -> VERSIONS file."""
    runs_fp = iou.Pathok(runs_dir, raise_err=False)
    if not runs_fp:
        return None
    cmd = (
        f"grep -m1 'Version' {str(runs_fp)}/*/run.log | awk -F: '/Version/ "
        + "{print $2 $3}' | sort -u"
    )
    out = iou.subprocess_run(cmd)
    if out is iou.subprocess.CalledProcessError:
        logger.error("Error fetching Version.")
        return
    msg = "MCCE Version(s) found in run.log files:\n"
    for v in [o.strip() for o in out.stdout.splitlines()]:
        msg += f"\t{v}\n"
    runs_fp.parent.joinpath(ANALYZE_DIR, FILES.VERSIONS.value).write_text(msg)

    return


def collate_files(bench_dir: Union[Path, str], fname: str) -> None:
    """Assumed: all same titr_type."""
    if fname not in ["pK.out", "sum_crg.out"]:
        logger.error(
            "Collation of files applicable to pK.out or sum_crg.out in <bench_dir>/runs subfolders."
        )
        sys.exit(1)

    bench = Path(bench_dir)
    # output dir
    analyze = bench.joinpath(ANALYZE_DIR)
    if not analyze.exists():
        analyze.mkdir()

    globbed = list(bench.joinpath(RUNS_DIR).glob(f"*/{fname}"))
    if not globbed:
        logger.error(f"No {fname} files found.")
        sys.exit(1)

    # output files
    pk_file = fname == "pK.out"
    if pk_file:
        fout_fp = bench.joinpath(ANALYZE_DIR, FILES.ALL_PKAS.value)
        txt_fp = bench.joinpath(ANALYZE_DIR, FILES.ALL_PKAS_TXT.value)
    else:
        fout_fp = bench.joinpath(ANALYZE_DIR, FILES.ALL_SUMCRG.value)
        txt_fp = bench.joinpath(ANALYZE_DIR, FILES.ALL_SUMCRG_TXT.value)

    hdr0, titr = "", ""
    with open(fout_fp, "w") as fout:
        for fp in globbed:
            p = fp.parent.name
            with open(fp) as fin:
                for j, line in enumerate(fin.readlines()):
                    if j == 0 and not hdr0:
                        hdr0 = line
                        titr = hdr0[:4]
                        _ = fout.write(f"   PDB:{hdr0.strip()}\n")
                    else:
                        if line.startswith(titr):
                            continue
                        if line.endswith("----------\n"):
                            continue
                        if pk_file and "too sharp" in line:
                            line = line.replace("too sharp", "toosharp")
                        _ = fout.write(f"{p}:" + line)

    # if all_pK.out: add note column & process oob values:
    allout_df = mcfile2df(fout_fp)
    # split the pdb:resid col:
    cols = allout_df.columns.to_list()
    col1, col2 = cols[0].split(":")
    allout_df[[col1, col2]] = allout_df[cols[0]].str.split(pat=":", expand=True)
    allout_df.drop(columns=cols[0], inplace=True)
    allout_df = allout_df[[col1, col2] + cols[1:]]
    # write the 'pretty'/split version:
    txt_fp.write_text(allout_df.to_string(index=False) + "\n")

    logger.info(f"Created {str(fout_fp)!r} and {str(txt_fp)!r}.")
    return


def get_all_pks_df(fpath: str, reduced_ok=True) -> Union[pd.DataFrame, None]:
    """Load <bench_dir>/analysis/all_pK.out into a pandas.DataFrame.

    Args:
      - fpath (str): bench_dir or a file path.
      - reduced_ok (bool, True): extract and save NA values & return a reduced df.
    Return:
      A pandas.DataFrame or None.
    """

    p = Path(fpath)
    if p.is_dir():
        allfp = p.joinpath(ANALYZE_DIR, FILES.ALL_PKAS_TXT.value)
    else:
        if p.is_file() and p.name == FILES.ALL_PKAS_TXT.value:
            allfp = p
        else:
            logger.error(
                f"Function loads {FILES.ALL_PKAS_TXT.value} into a df, not {p.name}."
            )
            sys.exit(1)

    # Load all_pK.out:
    allout_df = mcfile2df(allfp)
    if reduced_ok:
        save_fp = allfp.parent.joinpath(FILES.ALL_PKAS_OOB.value)
        if not save_fp.exists():  # not created
            oob_df = allout_df.loc[allout_df["pKa/Em"].isna()]
            if oob_df.shape[0]:
                oob_df.to_csv(Path(save_fp), sep="\t")
        allout_df = allout_df.loc[~allout_df["pKa/Em"].isna()].drop(columns="note")

    return allout_df


def save_all_run_times(pdbs_dir: str, overwrite: bool = True) -> None:
    """Return mcce step times from run.log saved to a tab-separated file."""

    pdbs = iou.Pathok(pdbs_dir)
    # output dir
    analyze = pdbs.parent.joinpath(ANALYZE_DIR)
    if not analyze.exists():
        analyze.mkdir()

    RUNS = '"/' + RUNS_DIR + '/"'  # len(RUNS) :: 6
    cmd = (
        "awk '/Total time/ {idx=index(FILENAME," + RUNS + ")+6; len=length(FILENAME); "
    )
    cmd = cmd + "print substr(FILENAME, idx, len-idx-7), $4, $(NF-1)}' "
    cmd = cmd + str(pdbs) + "/*/run.log | sed 's/MC:/step4/'"
    out = iou.subprocess_run(cmd)
    if out is iou.subprocess.CalledProcessError:
        logger.error("Error fetching run times.")
        return

    fp = analyze.joinpath(FILES.RUN_TIMES.value)
    if fp.exists():
        if overwrite:
            fp.unlink()
        else:
            return

    time_list = ["  ".join(line.split()) + "\n" for line in out.stdout.splitlines()]
    with open(fp, "w") as o:
        o.writelines(["PDB   step  seconds\n"])
        o.writelines(time_list)

    return


def get_step2_count(step2_out_path: str, kind: str) -> int:
    """Return the count of items given by `kind`, one of ['res', 'confs']
    from a step2_out.pdb file.
    """

    kind = kind.lower()
    if kind not in ["res", "confs"]:
        logger.error(f"kind must be one of ['res','confs']; Given: {kind}.")
        raise ValueError(f"kind must be one of ['res','confs']; Given: {kind}.")

    if kind == "res":
        cmd = "awk '!/CTR|NTR|HETATM/ {print substr($0,22,5)}' "
    else:
        cmd = "awk '/ATOM|HETATM/!/_000/ {print substr($0,22,5)}' "
    # complete cmd with common part:
    cmd = cmd + step2_out_path + "| uniq | wc -l"
    data = iou.subprocess_run(cmd)
    if isinstance(data, iou.subprocess.SubprocessError):
        logger.exception(f"Error fetching {kind} count.")
        raise data

    data = data.stdout.strip()
    if not data:
        logger.info(
            f"No count from step2_out.pdb in {Path(step2_out_path).parent.name}"
        )
        return 0

    return int(data)


def save_all_counts(pdbs_dir: str, kind: str, overwrite: bool = True) -> None:
    """Save the count of items given by `kind` from step2_out.pdb in all subfolders
    of pdbs_dir to a tab-separated file; format: DIR \t n.
    """
    pdbs = iou.Pathok(pdbs_dir)

    kind = kind.lower()
    if kind not in ["res", "confs"]:
        logger.error(f"'kind' must be one of ['res','confs']; Given: {kind}.")
        raise ValueError(f"'kind' must be one of ['res','confs']; Given: {kind}.")

    # output dir
    analyze = pdbs.parent.joinpath(ANALYZE_DIR)
    if not analyze.exists():
        analyze.mkdir()

    fname = FILES.CONF_COUNTS.value
    if kind == "res":
        fname = FILES.RES_COUNTS.value

    fp = analyze.joinpath(fname)
    if fp.exists():
        if overwrite:
            fp.unlink()
        else:
            return

    count_list = []
    for dp in pdbs.glob("./*/step2_out.pdb"):
        N = get_step2_count(str(dp), kind=kind)
        count_list.append(f"{dp.parent.name}  {N}\n")

    with open(fp, "w") as o:
        o.writelines([f"PDB  {kind}\n"])
        o.writelines(count_list)

    return


def save_confs_per_res(pdbs_dir: str, overwrite: bool = True) -> None:
    """Save conf per res to txt."""
    pdbs = iou.Pathok(pdbs_dir)
    # output dir
    analyze = pdbs.parent.joinpath(ANALYZE_DIR)
    if not analyze.exists():
        analyze.mkdir()

    # confs file:
    save_all_counts(pdbs, kind="confs", overwrite=overwrite)
    txt_count = analyze.joinpath(FILES.CONF_COUNTS.value)
    # res file:
    save_all_counts(pdbs, kind="res", overwrite=overwrite)
    txt_res = analyze.joinpath(FILES.RES_COUNTS.value)

    df_res = iou.txt2df(txt_res, header=0).set_index("PDB")
    df_res.sort_index(inplace=True)
    df_confs = iou.txt2df(txt_count, header=0).set_index("PDB")
    df_confs.sort_index(inplace=True)

    df = df_confs.merge(df_res, on="PDB")
    df["confs_per_res"] = round2(df.confs / df.res)

    # final output:
    txt_fin = analyze.joinpath(FILES.CONFS_PER_RES.value)
    if txt_fin.exists() and overwrite:
        txt_fin.unlink()
    txt_fin.write_text(df.to_string() + "\n")

    return


def save_confs_throughput(pdbs_dir: str, overwrite: bool = True) -> pd.DataFrame:
    """
    Obtain and save the average time & conformer throughput per step in a txt
    file, FILES.CONFS_THRUPUT.
    """
    pdbs = iou.Pathok(pdbs_dir)
    # output dir
    analyze = pdbs.parent.joinpath(ANALYZE_DIR)
    if not analyze.exists():
        analyze.mkdir()

    # times file:
    txt_time = analyze.joinpath(FILES.RUN_TIMES.value)
    if txt_time.exists() and overwrite:
        txt_time.unlink()
    save_all_run_times(pdbs, overwrite=overwrite)
    df_time = iou.txt2df(txt_time, header=0).set_index("PDB")
    df_time.sort_index(inplace=True)

    # confs file:
    txt_count = analyze.joinpath(FILES.CONF_COUNTS.value)
    if txt_count.exists() and overwrite:
        txt_count.unlink()
    save_all_counts(pdbs, kind="confs", overwrite=overwrite)
    df_confs = iou.txt2df(txt_count, header=0).set_index("PDB")
    df_confs.sort_index(inplace=True)

    df = df_confs.merge(df_time, on="PDB")
    logger.info(f"merged: {df.columns.to_list()}")
    df["confs_per_sec"] = round2(df.confs / df.seconds)
    df["confs_per_min"] = round2(df.confs_per_sec * 60)

    # final output:
    txt_fin = analyze.joinpath(FILES.CONFS_THRUPUT.value)
    if txt_fin.exists() and overwrite:
        txt_fin.unlink()

    gp = df[df.seconds > 0].groupby(by="step").aggregate("mean")  # , as_index=True
    txt_fin.write_text(gp.to_string() + "\n")

    return


def job_pks_to_dict(book_fpath: str) -> dict:
    """
    Uses the 'book' file to retrieve COMPLETED jobs, together with all_pK.out
    instead of iterating over subfolders (which may not be there.)
    Canonical dir struc: book_fpath points to <bench_dir>/runs/book.txt
    Uses <bench_dir>/analysis/all_pK.out.
    """

    book_fp = iou.Pathok(book_fpath)
    completed_dirs = iou.get_book_dirs_for_status(book_fp)  # default 'c'
    all_out_fp = book_fp.parent.parent.joinpath(ANALYZE_DIR, FILES.ALL_PKAS_TXT.value)
    # reduced df:
    all_df = get_all_pks_df(all_out_fp)
    c_resid, c_pk = all_df.columns[1:3]
    calc_pks = {}
    for pdbid in completed_dirs:
        pk_df = all_df.loc[all_df.PDB == pdbid]  # filter for this dir
        for r, row in pk_df.iterrows():
            resid = row[c_resid]
            if resid.startswith("NTG"):
                resid = "NTR" + resid[3:]
            key = (pdbid, resid)
            calc_pks[key] = row[c_pk]

    return calc_pks


def expl_pks_masterfile_to_df(
    infile: Path = BENCH.BENCH_WT, drop_na: bool = False
) -> pd.DataFrame:
    """Load an experimental pKas master file into a pandas.DataFrame.
    Currently, there is only one file for wild-type proteins:
    /data/pkadbv1/WT_pkas.csv'.
    """

    def to_float(value: str):
        value = value.strip()
        try:
            new = float(value)
        except (ValueError, TypeError):
            new = pd.NA
            logger.debug(f"Error converting {value} -> NA")
        return new

    def to_upper(s: str) -> str:
        return s.upper()

    expl_df = pd.read_csv(
        infile,
        usecols=["PDB ID", "Res Name", "Chain", "Res ID", "Expt. pKa"],
        comment="#",
        converters={
            "PDB ID": to_upper,
            "Res Name": to_upper,
            "Expt. pKa": to_float,
        },
    )
    if drop_na:
        expl_df.dropna(how="any", inplace=True)
    return expl_df.sort_values(by=["PDB ID", "Res ID"], ignore_index=True)


def expl_pks_to_dict(expl_type: str = "WT") -> dict:
    """Uses packaged resource proteins file specified by expl_type.
    Args:
      - expl_type (str, "WT"): Which protein type experimental data to use.
    Note:
      Currently, only "WT" is available as expl_type.
    Return:
      - A dict with 2-tuple keys, e.g.: ('135L', 'GLU-A0007_').
    """
    res_to_mcce = {
        "ARG": "ARG+",
        "ASP": "ASP-",
        "CYS": "CYS-",
        "GLU": "GLU-",
        "HIS": "HIS+",
        "LYS": "LYS+",
        "N-TERM": "NTR+",
        "C-TERM": "CTR-",
        "TYR": "TYR-",
    }
    if expl_type.upper() != "WT":
        logger.error("Only WT experimental pkas file is available; in use.")
        expl_type = "WT"
    pkas_df = expl_pks_masterfile_to_df(drop_na=True)

    expl_pks = {}
    for _, row in pkas_df.iterrows():
        key = (
            row["PDB ID"],
            f'{res_to_mcce[row["Res Name"]]}{row["Chain"]}{int(row["Res ID"]):04d}_',
        )
        expl_pks[key] = row["Expt. pKa"]

    return expl_pks


def match_pks(calc_pks: dict, expl_pks: dict) -> list[Tuple]:
    """Return a list of 4-tuples:
      id=<pdb>/<res>, calculated pka, experimental pka, difference (calculated - experimental).
    Convention: second dict is the reference/control.
    """
    calculated_ids = list(set([key[0] for key in calc_pks]))
    pkas = []
    for key in expl_pks:
        if key[0] not in calculated_ids:
            continue

        if key in calc_pks:
            calc_pka = calc_pks[key]
        elif key[1][3] == "-":
            if key[1][:3] == "TYR":
                calc_pka = 14.0
            else:
                calc_pka = 0.0
        elif key[1][3] == "+":
            calc_pka = 14.0
        else:
            logger.error(f"Parsing error of job pKas for {key}")

        pkas.append(
            (
                "{}/{}".format(*key),
                round2(calc_pka),
                round2(expl_pks[key]),
                round2(calc_pka - expl_pks[key]),
            )
        )

    return pkas


def matched_pks_txt(fpath: str, matched_pkas: list, sets_names: tuple) -> None:
    """Write a list of 4-tuples (as in a matched pkas list) to a txt file.
    Args:
    fpath (str): Output file path
    matched_pkas (list): List of 4-tuples: key, pka_set1, pka_set2, delta (1 - 2).
    sets_names (tuple): A 2-tuple to set names of columns 2 & 3 in the output file.
    """
    fp = Path(fpath)
    cols = f"key,{sets_names[0].lower()},{sets_names[1].lower()},delta".split(",")
    mdf = pd.DataFrame(matched_pkas, columns=cols)
    cols = mdf.columns.to_list()
    col1, col2 = ["pdb", "resid"]
    mdf[[col1, col2]] = mdf[cols[0]].str.split(pat="/", expand=True)
    mdf.drop(columns=cols[0], inplace=True)
    mdf = mdf[[col1, col2] + cols[1:]]
    fp = fp.parent.joinpath(FILES.MATCHED_PKAS_TXT.value)
    iou.df2txt(mdf, fp, replace=True)

    return


def matched_pks_stats(
    matched_df: pd.DataFrame,
    titr_type: str = "ph",
    prec: int = 2,
    level: Union[int, Tuple[int, int]] = None,
) -> Tuple[dict, pd.DataFrame]:
    """Perform a fit of the data and gather statistics/info.
    Also save the report to FILES.MATCHED_PKAS_STATS in analysis dir.
    Returns:
      A 2-tuple (dict, pandas.DataFrame).
        Dictionary:
                {"N": N matched,
                "titr_vals": [pKa | Em], id of titration values,
                "level": Tuple[(level number, name)]
                # LLS fit data if converged, else (None, error message)
                "fit": (b, m, rmse, R^2),
                "data_stats": (rmsd, corr),
                "bounds": bounds used for within proportion calc,
                "within": absolute delta pKs within bounds: [(b, pct within b),..],
                "res_table": rmsd info by res,
                "report": pkas_stats}.
        Pandas.DataFrame: holds the residues stats.
    Note:
      The report key includes the residues stats.
      The matched_df cols should be : PDB,resid,<set1>,<set2>,delta
    """
    def get_rmsd(delta):
        """Used for aggregating with pandas."""
        return np.sqrt(np.mean(delta**2))

    d_out = dict()
    N = matched_df.shape[0]
    d_out["N"] = N

    matched_df["res"] = matched_df.resid.apply(lambda x: x[:3])
    matched_df.drop(["pdb", "resid"], axis=1, inplace=True)
    # cols: <set1>, <set2>, delta
    # col1:A :: calc; col2:B :: ref (calc or expl)
    A_name, B_name = matched_df.columns[:2].to_list()
    A = matched_df.iloc[:, 0]
    B = matched_df.iloc[:, 1]

    # add delta absolute value:
    matched_df["delta_abs"] = matched_df.delta.abs()
    matched_df.drop([A_name, B_name], axis=1, inplace=True)

    # get res stats:
    res_grp = (
        matched_df.groupby("res", as_index=False)
        .agg(
            n=pd.NamedAgg(column="res", aggfunc="size"),
            delta_mean=pd.NamedAgg(column="delta", aggfunc="mean"),
            delta_abs_mean=pd.NamedAgg(column="delta_abs", aggfunc="mean"),
            delta_std=pd.NamedAgg(column="delta", aggfunc="std"),
            delta_rmsd=pd.NamedAgg(column="delta", aggfunc=get_rmsd),
        )
        .round(prec)
    )

    # add overall stats row:
    delta = matched_df["delta"]
    res_grp.loc[res_grp.shape[0]] = {
        "res": "All",
        "n": N,
        "delta_mean": round(delta.mean(), prec),
        "delta_abs_mean": round(matched_df["delta_abs"].mean(), prec),
        "delta_std": round(delta.std(), prec),
        "delta_rmsd": round(get_rmsd(delta), prec),
    }

    # convert to text for later assignment to 'report' key:
    fmt = "{:>." + str(prec) + "f}"
    grp_txt = res_grp.to_string(index=False, float_format=fmt.format)

    if titr_type == "ph":
        titr_d = 1
        which_val = "pKa"
    else:
        titr_d = 30  # mV, default interval for eH titr
        which_val = "Em"
    d_out["titr_vals"] = which_val

    if B_name == "expl":
        B_name = "pKDB"

    levtupl = iou.get_level_tup(level)
    d_out["level"] = levtupl

    # initialize report text:
    txt = f"Residues stats for {A.name!r} "

    multi_valued = isinstance(levtupl[0], tuple)
    if not multi_valued:
        if levtupl[0]:
            txt = txt + f"v. {B_name!r} @ {levtupl[1]!r} level:\n"
        else:
            txt = txt + f"v. {B_name!r}:\n"
    else:
        l1 = levtupl[0][0]
        l2 = levtupl[1][0]
        if not l1 and not l2:
            txt = f"v. {B_name!r}:\n"
        else:
            txt = f"(L: {levtupl[0][1]}) v. {B_name!r} (L: {levtupl[1][1]}):\n"

    txt = txt + f"Number of {which_val}s matched: {N:,}\n"

    rmsd = get_rmsd(A - B)
    corr = plots.correl(A, B)
    d_out["data_stats"] = rmsd, corr

    # get fit results:
    converged, pfit, rmse, r_sqrd, err_msg = plots.fit_series(A, B)

    if converged:
        b, m = pfit.convert()
        txt = (
            txt
            + f"Fit line: y = {m:.{prec}f}.x + {b:.{prec}f}, R^2: {r_sqrd:.{prec}f}\n"
        )

        d_out["fit"] = (round2(b), round2(m), round2(rmse), round2(r_sqrd))
    else:
        txt = txt + f"No fit line ({err_msg})\n"
        # The 'fit' key will be a 2-tuple:
        d_out["fit"] = (False, f"({err_msg})")

    txt = txt + f"RMSD: {rmsd:.{prec}f}, Corr = {corr:.{prec}f}\n"

    comp_bounds = [3 * titr_d, 2 * titr_d, 1.5 * titr_d, titr_d, 0.5 * titr_d]
    d_out["bounds"] = comp_bounds
    txt = txt + f"Absolute {DELTA}pK proportion within these bounds:\n"
    within_tups = []
    for b in comp_bounds:
        v = matched_df.delta_abs[matched_df.delta_abs.le(b)].count() / N
        within_tups.append((round(b, 1), round(v,2)))
        txt = txt + f"  abs({DELTA}pK) within {b:>3} units: {v:>7.{prec}%}\n"
    d_out["within"] = within_tups

    # add stats grouped by res in own key:
    d_out["res_table"] = grp_txt
    d_out["report"] = txt + "\n" + grp_txt + "\n"

    return d_out, res_grp


def res_outlier_count(
    matched_fp: str,
    grp_by: str = "res",
    replace: bool = False,
    bounds: tuple = (0, 14),
    return_df: bool = False,
) -> Union[None, pd.DataFrame]:
    """Return counts per residue type for diff(col1-col2) > 3,
    and pKa values beyond titration bounds in a df.
    Save df to FILES.RES_OUTLIER or FILES.RESID_OUTLIER in the
    parent folder of matched_fp.
    Args:
    matched_fp (str): file path of matched pkas file;
    grp_by (str, "res"): one of ["res", "resid"];
    replace (bool, False): To overwrite existing file;
    bounds (tuple, (0,14)): Default titration bounds.
    """
    if grp_by not in ["res", "resid"]:
        logger.error(f"grp_by not in ['res','resid']: {grp_by}")
        raise ValueError(f"grp_by not in ['res','resid']: {grp_by}")

    matched_df = iou.matched_pks_to_df(matched_fp)  # :: pdb, resid, set1, set2, delta
    N = matched_df.shape[0]
    cols_pks = matched_df.columns[2:4].to_list()
    idx_name = f"{grp_by.upper()} | {cols_pks[0]} v. {cols_pks[1]}"

    matched_df["res"] = matched_df.resid.str[:3]
    matched_df.sort_values(by="res", inplace=True)
    matched_df = matched_df[["res", "resid", cols_pks[0], cols_pks[1]]]

    if grp_by == "res":
        matched_df[["res", "resi"]] = matched_df.resid.str.split("[-|+]", expand=True)
        matched_df.drop(columns=["resi", "resid"], inplace=True)

    set1 = matched_df.iloc[:, -2]
    set2 = matched_df.iloc[:, -1]

    matched_df["delta"] = abs(set1 - set2)
    matched_df["Delta over 3"] = matched_df.delta > 3.0
    matched_df["Out of bounds"] = (abs(set1 - bounds[0]) < 0.01) | (
        abs(set1 - bounds[1]) < 0.01
    )
    matched_df.drop(columns=[set1.name, set2.name, "delta"], inplace=True)
    msk = matched_df["Out of bounds"]  # == True
    gp_oob = matched_df[msk].groupby(grp_by).count()
    gp_oob.drop(columns="Delta over 3", inplace=True)
    msk = matched_df["Delta over 3"]  # == True
    gp_del3 = matched_df[msk].groupby(grp_by).count()
    gp_del3.drop(columns="Out of bounds", inplace=True)

    out_df = (
        gp_oob.merge(gp_del3, how="left", on=grp_by).replace({np.nan: 0}).astype(int)
    )
    oob_name = f"Out of bounds {bounds}"
    out_df.rename(columns={"Out of bounds": oob_name}, inplace=True)
    out_df.index.name = idx_name
    pcts = [f"{s/N:.1%}" for s in out_df.sum()]
    out_df.loc[f"pct (N:{N})"] = pcts

    if grp_by == "res":
        outlier_fp = Path(matched_fp).parent.joinpath(FILES.RES_OUTLIER.value)
    else:
        out_df.drop(columns=["res_x", "res_y"], inplace=True)
        outlier_fp = Path(matched_fp).parent.joinpath(FILES.RESID_OUTLIER.value)
    logging.info(f"\n{out_df.to_string()}\n")

    if outlier_fp.exists() and replace:
        outlier_fp.unlink()
    outlier_fp.write_text(out_df.to_string() + "\n")
    if return_df:
        return out_df
    return


def analyze_runs(bench_dir: Path, subcmd: str):
    """Create all analysis output files."""
    bench = iou.Pathok(bench_dir)
    pdbs = bench.joinpath(RUNS_DIR)
    book_fp = pdbs.joinpath(BENCH.Q_BOOK)

    analyze = bench.joinpath(ANALYZE_DIR)
    if not analyze.exists():
        analyze.mkdir()
    else:
        clear_folder(analyze)

    get_mcce_version(pdbs)

    logger.info("Collating pK.out and sum_crg.out files.")
    collate_files(bench, fname="pK.out")
    collate_files(bench, fname="sum_crg.out")

    logger.info("Calculating conformers and residues counts into file.")
    save_all_counts(pdbs, kind="confs", overwrite=True)
    save_all_counts(pdbs, kind="res", overwrite=True)
    save_all_run_times(pdbs, overwrite=True)
    save_confs_per_res(pdbs)

    logger.info("Calculating conformers thoughput into file.")
    save_confs_throughput(pdbs)

    logger.info("Plotting conformers throughput per step -> pic.")
    # get conf-making level from serialized setup args:
    lev = iou.get_setup_args_vals(bench, "conf_making_level")
    if lev:
        if lev[0] is not None:
            lev = lev[0]
        else:
            lev = 0  # non-existing level
    else:
        lev = 0

    # get number of completed runs:
    n_complete = len(iou.get_book_dirs_for_status(book_fp))
    thruput_df = iou.txt2df(analyze.joinpath(FILES.CONFS_THRUPUT.value), header=0)
    plots.plot_conf_thrup(
        thruput_df,
        n_complete,
        analyze.parent.name,
        level=lev,
        out_fp=analyze.joinpath(FILES.FIG_CONFS_TP.value),
    )

    logger.info("Getting calculated pKas to dict.")
    # effective calculated pkas for all completed runs:
    calc_pkas = job_pks_to_dict(book_fp)
    iou.to_pickle(calc_pkas, analyze.joinpath(FILES.JOB_PKAS_PKL.value))

    env = mcenv.get_run_env(bench, subcmd=subcmd)
    titr_type = env.runprm["TITR_TYPE"]
    titr_bounds = env.get_titr_bounds()

    if subcmd == SUB1:
        logger.info("Getting experimental pKas to dict.")
        expl_pkas = expl_pks_to_dict()
        # calc_pkas: done

        logger.info("Matching the pkas and saving list to file.")
        matched_pkas = match_pks(calc_pkas, expl_pkas)
        # 'pretty' file:
        matched_fp = analyze.joinpath(FILES.MATCHED_PKAS_TXT.value)
        matched_pks_txt(matched_fp, matched_pkas, sets_names=(bench.name, "expl"))

        res_outlier_count(matched_fp, bounds=titr_bounds)

        logger.info("Calculating the matched pkas stats & saving dict.")
        matched_df = iou.matched_pks_to_df(matched_fp)
        d_stats, res_stats_df = matched_pks_stats(
            matched_df, titr_type=titr_type, level=lev
        )
        # save pk stats dict:
        iou.to_pickle(d_stats, analyze.joinpath(FILES.MATCHED_PKAS_STATS_PKL.value))
        # save res_stats_df to dict, to be reused in res plot:
        res_stats_d = res_stats_df.set_index("res").to_dict(orient="index")
        iou.to_pickle(res_stats_d, analyze.joinpath(FILES.RESIDUES_STATS_PKL.value))

        # save res_stats_df to text:
        res_stats_fp = analyze.joinpath(FILES.RESIDUES_STATS.value)
        res_stats_fp.write_text(
            res_stats_df.to_string(index=False, float_format="{:.2f}".format) + "\n"
        )

        # log & save mini report:
        rpt = d_stats["report"]
        logger.info(rpt)
        rpt_fp = analyze.joinpath(FILES.MATCHED_PKAS_STATS.value)
        rpt_fp.write_text(rpt)

        if not d_stats["fit"][0]:
            logger.info("Plotting pkas data without fit -> pic.")
        else:
            logger.info("Plotting pkas fit -> pic.")

        plots.plot_pkas_fit(
            matched_fp, d_stats, out_fp=analyze.joinpath(FILES.FIG_FIT_ALLPKS.value)
        )

        logger.info("Plotting residues analysis -> pic.")
        plots.plot_res_analysis(
            matched_fp,
            res_stats_d,
            lev,
            out_fp=analyze.joinpath(FILES.FIG_FIT_PER_RES.value),
        )

    return


# ................................................................................
CLI_NAME = ENTRY_POINTS["analyze"]  # as per pyproject.toml entry point
DESC = """
Description:
Create analysis output files for a set of runs in <bench_dir>/analysis.

The main command is {} along with one of 2 sub-commands:
- Sub-command 1: {}: analyze pKas against pKaDBv1;
- Sub-command 2: {}: analyze pKas for user's pdbs;

Output files: listed in benchmark.info.
"""


def pkdb_pdbs_analysis(args: Union[dict, Namespace]) -> None:
    """Processing tied to sub-command 1: pkdb_pdbs."""
    if isinstance(args, dict):
        args = Namespace(**args)
    analyze_runs(args.bench_dir, SUB1)
    return


def user_pdbs_analysis(args: Union[dict, Namespace]) -> None:
    """Processing tied to sub-command 2: user_pdbs."""
    if isinstance(args, dict):
        args = Namespace(**args)
    analyze_runs(args.bench_dir, SUB2)
    return


def icase(s: str) -> str:
    """Return the correctly cased value for known strings."""
    if s.lower() == "pk.out":
        return "pK.out"
    elif s.lower() == "step2_out.pdb":
        return "step2_out.pdb"
    else:
        return s


def arg_valid_dirpath(p: str):
    """Return resolved path from the command line."""
    if not len(p):
        return None
    return Path(p).resolve()


def analyze_parser():
    """Cli arguments parser with sub-commands for use in benchmark analysis."""
    p = ArgumentParser(
        description=DESC.format(CLI_NAME, SUB1, SUB2),
        formatter_class=RawDescriptionHelpFormatter,
        epilog="""
        Post an issue for all errors and feature requests at:
        https://github.com/CatChenal/MCCE4/issues
        """,
    )

    subparsers = p.add_subparsers(
        required=True,
        title=f"{CLI_NAME} sub-commands",
        dest="subparser_name",
        description="Sub-commands of MCCE benchmarking analysis cli.",
        help=f"""
        The 2 choices for the benchmarking process:
            1) Analyze dataset of completed runs viz pKaDBv1: {SUB1}
            2) Analyze dataset of completed runs of user-provided pdbs : {SUB2}
        """
    )

    sub1 = subparsers.add_parser(
        SUB1,
        formatter_class=RawDescriptionHelpFormatter,
        help="""
        Sub-command for analyzing a benchmarking set against the pKaDBv1
        using the same dataset and structure: <bench_dir>/runs folder."
        """,
    )
    sub1.add_argument(
        "-bench_dir",
        required=True,
        type=arg_valid_dirpath,
        help="""The user's directory where the /runs subfolder is setup.
        """,
    )
    sub1.set_defaults(func=pkdb_pdbs_analysis)

    sub2 = subparsers.add_parser(
        SUB2,
        formatter_class=RawDescriptionHelpFormatter,
        help="""
        Sub-command for analyzing a benchmarking set of user's pdbs.
        """,
    )
    sub2.add_argument(
        "-bench_dir",
        required=True,
        type=arg_valid_dirpath,
        help="""The user's directory where the /runs subfolder is setup.
        """,
    )
    sub2.set_defaults(func=user_pdbs_analysis)
    return p


def analyze_cli(argv=None):
    """
    Command line interface for MCCE benchmarking analysis entry point.
    """
    cli_parser = analyze_parser()
    args = cli_parser.parse_args(argv)

    cli_opts.cli_name = CLI_NAME
    cli_opts.all = vars(args)
    logger.info(cli_opts)

    # OK to analyze?
    bench = iou.Pathok(args.bench_dir)
    book_fp = bench.joinpath(RUNS_DIR, BENCH.Q_BOOK)
    logger.info(f"book_fp: {book_fp}")
    pct = iou.pct_completed(book_fp)
    if pct is None:
        logger.critical((f"Something seems wrong with the book file:\n{str(book_fp)}\nTry:\n"
                         " grep '\b[ce]\b' runs/book.txt | wc -l; cat runs/book.txt | wc -l \n"
                         "in the bench_dir."
                         )
        )
        sys.exit(1)

    if pct < 1.0:
        logger.info(f"Runs not 100% processed, try again later; {pct = :.1%}")
        return

    args.func(args)

    return


if __name__ == "__main__":
    analyze_cli(sys.argv[1:])
