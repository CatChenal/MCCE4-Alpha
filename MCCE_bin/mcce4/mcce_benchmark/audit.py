#!/usr/bin/env python

"""Module: audit.py
Contains functions to query and manage data.


Main functions:
--------------
Note: proteins.tsv should be considered the ground truth.

def list_complete_runs(benchmarks_dir:str, like_runs:bool=False) -> list:
def cp_completed_runs(src_dir:str, dest_dir:str) -> None:
def proteins_df(prot_tsv_file:Path=BENCH.BENCH_PROTS, return_excluded:bool=None) -> pd.DataFrame:
def get_usable_prots(prot_tsv_file:Path=BENCH.BENCH_PROTS) -> list:
def valid_pdb(pdb_dir:str, return_name:bool = False) -> Union[bool, Path, None]:
def list_all_valid_pdbs(pdbs_dir:Path = BENCH.BENCH_PDBS) -> tuple:
def list_all_valid_pdbs_dirs(pdbs_dir:Path = BENCH.BENCH_PDBS) -> tuple:
def multi_model_pdbs(pdbs_dir:Path = BENCH.BENCH_PDBS) -> Union[np.ndarray, None]:
def reset_multi_models(pdbs_dir:Path = BENCH.BENCH_PDBS, debug:bool = False) -> list:
def update_proteins_multi(proteins_file:Path = BENCH.BENCH_PROTS):
def rewrite_book_file(book_file:Path) -> None:
def pdb_list_from_book(book_file:Path = Path(BENCH.Q_BOOK)) -> list:
def pdb_list_from_runs_folder(pdbs_dir:Path = BENCH.BENCH_PDBS) -> list:
def prots_symdiff_runs(prot_tsv_file:Path=BENCH.BENCH_PROTS,
def update_data(prot_tsv_file:Path=BENCH.BENCH_PROTS,
def same_pdbs_book_vs_runs() -> bool:
def pdb_list_from_experimental_pkas(pkas_file:Path=BENCH.BENCH_WT) -> list:
def proteins_to_tsv(prot_file:str) -> list:
"""

import logging
import numpy as np
import pandas as pd
from pathlib import Path
import subprocess
import shutil
from typing import Union
from mcce4.mcce_benchmark import BENCH, RUNS_DIR, SUB1, SUB2
from mcce4.mcce_benchmark.pkanalysis import analyze_runs, expl_pks_masterfile_to_df


logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)


MULTI_ACTIVE_MSG = """Multi-model folder {!r} contains multiple 'active' pdbs with
{!r} being the second one.
The only 'active' pdb (to be selected as 'prot.pdb'), must be the one
listed in the 'Use' colummn in the 'proteins.tsv' file.
The function audit.reset_multi_models() must be re-run to fix the problem.
"""


def redo_refset_analysis():
    """Rerun pkanalysis in the reference dataset."""
    analyze_runs(BENCH.BENCH_PARSE_PHE4, SUB1)
    print("redo_refset_analysis over.")
    return


def list_complete_runs(benchmarks_dir: str, sentinel_file: str = "pK.out") -> list:
    """Return a list of folders that contain sentinel_file (pK.out default).
    Search in benchmarks_dir/RUNS_DIR subfolders.
    """
    search_dir = Path(benchmarks_dir).joinpath(RUNS_DIR)
    complete = list(fp.parent for fp in search_dir.glob(f"./*/{sentinel_file}"))

    return complete


def cp_completed_runs(src_dir: str, dest_dir: str) -> None:
    complete = list_complete_runs(src_dir)
    for fp in complete:
        dest = Path(dest_dir).joinpath(RUNS_DIR, fp.name)
        if dest.exists():
            shutil.copytree(
                fp, dest, dirs_exist_ok=True, ignore=shutil.ignore_patterns("prot.pdb*")
            )
            logger.info("Copied:", dest)
    return


def proteins_df(
    prot_tsv_file: Path = BENCH.BENCH_PROTS, return_excluded: bool = None
) -> pd.DataFrame:
    """
    Load data/pkadbv1/proteins.tsv into a pandas.DataFrame.
    Args:
    return_excluded (bool, None): If None: return entire set;
    If True: return df of commented out entries; If False, return the
    non-excluded entries
    and see the reasons), else return df of "good to go"  proteins.
    """
    df = pd.read_csv(prot_tsv_file, sep="\t")
    df.sort_values(by="PDB", inplace=True)
    if return_excluded is None:
        return df
    if return_excluded:
        return df[df.Use.isna()]
    else:
        return df[~df.Use.isna()]


def get_usable_prots(prot_tsv_file: Path = BENCH.BENCH_PROTS) -> list:
    """
    Return a list of uncommented pdb ids from proteins.tsv.
    """
    df = proteins_df(return_excluded=False)
    return df.PDB.to_list()


def valid_pdb(pdb_dir: str, return_name: bool = False) -> Union[bool, Path, None]:
    """Return whether 'pdb_dir' contains a valid pdb (bool, default), or its
    Path if 'return_name'=True if valid, else None.
    Used by list_all_valid_pdbs.
    """
    pdb_dir = Path(pdb_dir)
    # single model pdb
    pdb = pdb_dir.joinpath(f"{pdb_dir.name.lower()}.pdb")
    ok = pdb.exists()
    if ok:  # found
        if return_name:
            return pdb
        else:
            return ok

    # multi-model protein: main pdb was renamed with .full extension
    files = [
        fp for fp in pdb_dir.glob("*.[pdb full]*") if not fp.name.startswith("model")
    ]
    found_full = False
    found_active = False
    active = []
    for f in files:
        if f.suffix == ".full":
            found_full = True
            continue
        if f.name.startswith(f"{pdb_dir.name.lower()}_"):
            found_active = True
            if return_name:
                pdb = f
            # to check if several are 'valid'
            active.append(True)
            if len(active) > 1:
                logger.error(MULTI_ACTIVE_MSG.format(pdb_dir.name, f.name))
                found_active = False
                break

    # FIX: 'and' => will fail if all multi prots are removed:
    valid = found_full and found_active
    if return_name:
        if not valid:
            return None
        return pdb
    else:
        return valid


def list_all_valid_pdbs(pdbs_dir: Path = BENCH.BENCH_PDBS) -> tuple:
    """Return a list ["PDB/pdb[_*].pdb", ..] of valid pdb.
    Return a 2-tuple of lists: (valid_folders, invalid_folders), with
    each list = ["PDB/pdb[_*].pdb", ..].
    For managing packaged data.
    """
    pdbs_dir = Path(pdbs_dir)
    if not pdbs_dir.is_dir():
        raise FileNotFoundError(f"Not a dir or does not exists: {pdbs_dir}.")

    valid = []
    invalid = []
    for fp in pdbs_dir.glob("*"):
        if fp.is_dir() and not fp.name.startswith("."):
            p = valid_pdb(fp, return_name=True)
            if p is None:
                invalid.append(fp.name)
            else:
                valid.append(f"{fp.name}/{p.name}")
    valid.sort()
    invalid.sort()
    logger.info(f"{len(valid) = }; {len(invalid) = }")
    if len(invalid):
        logger.warning(f"Invalid pdbs: {invalid}")

    return valid, invalid


def list_all_valid_pdbs_dirs(pdbs_dir: Path = BENCH.BENCH_PDBS) -> tuple:
    """Check that all subfolders of RUNS_DIR contain a pdb file with
    the same name.
    Return a 2-tuple: (valid_folders, invalid_folders).
    For managing packaged data.
    """
    pdbs_dir = Path(pdbs_dir)
    if not pdbs_dir.is_dir():
        raise FileNotFoundError(f"Not a dir or does not exists: {pdbs_dir}.")

    valid = []
    invalid = []
    for fp in pdbs_dir.glob("*"):
        if fp.is_dir() and not fp.name.startswith("."):
            v = valid_pdb(fp)
            if v:
                valid.append(fp.name)
            else:
                invalid.append(fp.name)
    valid.sort()
    invalid.sort()
    logger.info(f"Valid folders: {len(valid)}; Invalid folders: {len(invalid)}")

    return valid, invalid


def multi_model_pdbs(pdbs_dir: Path = BENCH.BENCH_PDBS) -> Union[np.ndarray, None]:
    """
    Query RUNS_DIR for pdb with multiple models.
    Return dir/pdb name in a numpy array or None.
    """
    multi_models = None
    query_path = pdbs_dir.joinpath("*/")
    n_parts = len(query_path.parts) + 2
    f1 = n_parts - 2
    f2 = n_parts - 1
    prt_fields = "{print $" + str(f1) + '","$' + str(f2) + "}"
    cmd = f"grep '^MODEL' {query_path}/*.pdb|sed -e '/model/d; s/:MODEL/\//g'|gawk -F'/|:' '{prt_fields}'|uniq"  # noqa W605
    try:
        data = subprocess.check_output(
            cmd, stderr=subprocess.STDOUT, check=True, text=True, shell=True
        )
        data = data.splitlines()
        if data:
            multi_models = np.array([line.split(",") for line in data])
    except subprocess.CalledProcessError as e:
        logger.exception(f"Error in subprocess cmd.\nException: {e}")
        raise

    return multi_models


def reset_multi_models(pdbs_dir: Path = BENCH.BENCH_PDBS, debug: bool = False) -> list:
    """Use multi model entries in 'data/pkadbv1/proteins.tsv' to select the
    model<x>.pdbs corresponding to the proteins 'Use' column, if found.

    Rename pdbid.pdb -> pdbid.pdb.full
    n = matched model number from Use.split(".")[1]
    new_n = Use.replace(".", "")
    Rename model{n:02}.pdb -> pdbid_{new_n}.pdb

    Should be re-run every time the 'Use' column in 'data/pkadbv1/proteins.tsv' is
    changed for one or more multi-model proteins.
    For managing packaged data.
    """
    prots_df = proteins_df()
    multi = prots_df[prots_df.Model == "multi"]

    missing_data = []

    for i, ro in multi.iterrows():
        pdb = ro.PDB.lower()
        chain, n = ro.Use.split(".")
        n = int(n.strip())

        mdir = pdbs_dir.joinpath(ro.PDB)
        pdb_path = mdir.joinpath(f"{pdb}.pdb")
        path_full = mdir.joinpath(f"{pdb}.pdb.full")
        use_prot = mdir.joinpath(f"{pdb}_{chain}{n}.pdb")
        modl_prot = mdir.joinpath(f"model{n:02}.pdb")

        if path_full.exists():
            # possibly already processed, check if protein in use matches:
            if use_prot.exists():
                # matched, done:
                continue
            # else, check the model pdb
            if not modl_prot.exists():
                print("Could not find {modl_prot} to rename as {use_prot}.")
                missing_data.append(modl_prot)
                continue

            # rename previously used prots:
            pdb_glob = f"{mdir.name.lower()}_*.pdb"
            used_prots = list(mdir.glob(pdb_glob))
            if used_prots:
                for p in used_prots:
                    if debug:
                        print(f"shutil.move({p}, {p}.x)")
                    else:
                        _ = shutil.move(p, f"{p}.x")

            if debug:
                print(f"shutil.copy({modl_prot}, {use_prot})")
            else:
                _ = shutil.copy(modl_prot, use_prot)

        else:
            if not pdb_path.exists():
                print("Expected pdb not found: {pdb_path}.")
                missing_data.append(pdb_path)
                continue
            if debug:
                print(f"shutil.move({pdb_path}, {path_full})")
            else:
                _ = shutil.move(pdb_path, path_full)

            if not modl_prot.exists():
                logger.error("Could not find {modl_prot} to rename as {use_prot}.")
                missing_data.append(modl_prot)
                continue

    return missing_data


def update_proteins_multi(proteins_file: Path = BENCH.BENCH_PROTS):
    """Update 'data/pkadbv1/proteins.tsv' Model column from
    list of multi-model proteins.
    For managing packaged data.
    """
    multi_models = multi_model_pdbs()
    if multi_models is None:
        print("No multi model pdbs.")
        return

    prots_df = proteins_df(proteins_file)
    for pdb in multi_models[:, 0]:
        prots_df.loc[prots_df.PDB == pdb, "Model"] = "multi"

    prots_df.to_csv(BENCH.BENCH_PROTS, index=False, sep="\t")

    return


def rewrite_book_file(book_file: Path, subcmd: str = SUB1) -> None:
    """Re-write RUNS_DIR/book file with validated entries if subcmd is
    'pkdb_pdbs', else with the directories names found in book_file parent folder.
    """
    if subcmd == SUB1:
        valid, invalid = list_all_valid_pdbs_dirs(book_file.parent)
    else:
        valid = [d.name for d in list(book_file.parent.glob("./*/"))
                 if (d.is_dir() and d.name.isupper())]
    with open(book_file, "w") as book:
        book.writelines([f"{v}\ti\n" for v in valid])

    return


def pdb_list_from_book(book_file: Path = Path(BENCH.Q_BOOK)) -> list:
    pdbs = []
    with open(book_file) as book:
        for line in book:
            rawtxt = line.strip().split("#")[0]
            pdbs.append(rawtxt.split()[0])

    return pdbs


def pdb_list_from_runs_folder(pdbs_dir: Path = BENCH.BENCH_PDBS) -> list:
    """pdbs_dir: folder with one PDBID folder for each pdbid.pdb file"""
    pdbs_dirs = [
        fp.name
        for fp in pdbs_dir.glob("*")
        if fp.is_dir() and not fp.name.startswith(".")
    ]
    pdbs_dirs.sort()

    return pdbs_dirs


def prots_symdiff_runs(
    prot_tsv_file: Path = BENCH.BENCH_PROTS, pdbs_dir: Path = BENCH.BENCH_PDBS
) -> tuple:
    """Get the symmetric difference btw the list of usable proteins
    and the runs/subfolders list.
    Return a tuple of lists: extra_dirs, missing_dirs.
    Package data management.
    """
    # list from "ground truth" file, proteins.tsv:
    curated_ok = get_usable_prots(prot_tsv_file)

    # list from folder setup: would differ if a change occured
    # without related change in proteins.tsv
    dir_list = pdb_list_from_runs_folder(pdbs_dir)

    s1 = set(curated_ok)
    s2 = set(dir_list)

    extra = s1.symmetric_difference(s2)
    extra_dirs = []
    missing_dirs = []

    for x in extra:
        if x not in s1:
            # print(f"Extra dir: {x}")
            extra_dirs.append(x)
        else:
            # print(f"Missing dir for: {x}")
            missing_dirs.append(x)

    return extra_dirs, missing_dirs


def update_data(
    prot_tsv_file: Path = BENCH.BENCH_PROTS, pdbs_dir: Path = BENCH.BENCH_PDBS
) -> None:
    """Delete extra subfolders from runs/ when corresponding pdb is not
    in proteins.tsv.
    """
    # Note: if missing dirs: no pdb => will need to be downloaded again and
    # processed (mannually) as per jmao.

    extra_dirs, missing_dirs = prots_symdiff_runs(prot_tsv_file, pdbs_dir)
    if not extra_dirs:
        logger.info("No extra dirs.")
    else:
        for x in extra_dirs:
            dx = pdbs_dir.joinpath(x)
            shutil.rmtree(dx)
        logger.info("Removed extra dirs.")

    book = pdbs_dir.joinpath(BENCH.Q_BOOK)
    rewrite_book_file(book)
    logger.info("Wrote fresh book file.")

    return


def same_pdbs_book_vs_runs() -> bool:
    """
    Compares the list of pdbs in the Q_BOOK with the list
    obtained from the runs/ folder.
    For managing packaged data.
    """
    book_pbs = pdb_list_from_book()
    pdbs = pdb_list_from_runs_folder()
    same = len(book_pbs) == len(pdbs)
    if not same:
        logger.warning(
            f"The lists differ in lengths:\n\t{len(book_pbs) = }; {len(pdbs) = }"
        )
        return same

    df = pd.DataFrame(zip(book_pbs, pdbs), columns=["book", "runs_dir"])
    comp = df[df.book != df.runs_dir]
    same = len(comp) == 0
    if not same:
        logger.warning(f"The lists differ in data:\n{comp}")

    return same


def to_float(value):
    """Conversion function to be used in pd.read_csv.
    Return NA on conversion failure.
    """
    try:
        new = float(value)
    except TypeError:
        new = np.nan

    return new


def pdb_list_from_experimental_pkas(pkas_file: Path = BENCH.BENCH_WT) -> list:
    """Parses valid pKa values from an experimental pKa file and return
    their pdb ids in a list.
    """
    if pkas_file.name != BENCH.BENCH_WT.name:
        raise TypeError(
            "Only wild type proteins listed in 'WT_pkas.csv' are currently considered."
        )

    pks_df = expl_pks_masterfile_to_df(drop_na=True)
    pdbs = pks_df["PDB ID"].unique().tolist()

    return sorted(pdbs)


##############################################################
# transformation of legacy files


def proteins_to_tsv(prot_file: str) -> list:
    """Transform initial text file to a tab separated file (.tsv)
    with additional column 'Model' to indicate whether pdb is single- or
    multi-modelled. No duplicates.
    """

    prots = Path(prot_file)
    tsv = prots.parent.joinpath("proteins.tsv")
    uniq = []
    with open(prots) as p, open(tsv, "w") as out:
        for i, line in enumerate(p):
            if i == 0:
                out.write("\t".join(line.split()) + "\tModel\n")
                continue

            pdb, fields = line.split("   ", maxsplit=1)
            if pdb not in uniq:
                uniq.append(pdb)
                if line.startswith("#"):
                    out.write(f"{pdb}\t{fields.strip()}\t\t\n")
                else:
                    other = fields.split("   ", maxsplit=1)
                    out.write(
                        f"{pdb}\t{other[0].strip()}\t{other[1].strip()}\tsingle\n"
                    )
    return
