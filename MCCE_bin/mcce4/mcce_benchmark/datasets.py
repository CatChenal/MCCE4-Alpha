#!/usr/bin/env python
"""Module: datasets.py

Functions:
----------

"""
from argparse import ArgumentParser, RawDescriptionHelpFormatter, Namespace
from functools import partial
import os
from pathlib import Path
import re
import sys

import pandas as pd

from mcce4 import downloads
from mcce4.mcce_benchmark import BenchResources, BENCH_DATA, RUNS_DIR, datasets_dict
from mcce4.mcce_benchmark import io_utils as iou, audit
from mcce4.protinfo.cli import get_pdb_rpt_climode


read_tsv = partial(pd.read_csv, sep="\t") #, converters={"PDB": str})
dataset_names = list(datasets_dict.keys())


def get_pkdb(args: Namespace) -> Path:
    pkdb = BENCH_DATA.joinpath(args.dataset_name)
    if not pkdb.exists():
        sys.exit(f"Error, dir not found: {pkdb!r}; Available: {dataset_names}.")
    return pkdb


def get_runs_dir(args: Namespace, create: bool=True) -> Path:
    pkdb = BENCH_DATA.joinpath(args.dataset_name)
    runs_dir = pkdb.joinpath(RUNS_DIR)
    if not runs_dir.exists():
        if not create:
            msg = (f"Error, directory not found: {runs_dir}\n"
                    "In order to further prep a dataset, its runs folder must"
                    "first be setup with its (bioassembly) pdbs."
                    )
            sys.exit(msg)

        runs_dir.mkdir()
    
    return runs_dir


def get_proteins_file(args: Namespace) -> Path:
    prots = BENCH_DATA.joinpath(args.dataset_name, args.proteins_file)
    if not prots.exists():
        sys.exit(f"Error, dir not found: {prots!r}.")
    return prots


#TODO: setup_dataset function
# SUB0 = "setup_dataset"
# -name -pka_file
# add folder in mcce4/mcce_benchmark/data called args.name
# copy args.pkas_file as pkas.csv into pkdb folder (assumed to have useable column names)
# add runs dir
# create the proteins file
# get the bioassemblies
# get the prerun reports
# finalize proteins file with data in protinfo dicts
# update pkas file with excluded as per prerun


def pkdb_r_prep(csv_in, csv_out):
    """Process the (assumed) intact, dowloaded csv file of PKAD-R to
    obtain a clean version."
    """
    #TODO: Adapt cleanup fn from notebook
    return


def create_proteins_file_from_pkas_file(dataset):
    """
    Create unique pdb list -> proteins file
    Also save proteins file with Chain to verify if it matched the one
    used in the bioassembly.
    """
    prots_chn = "proteins_chain.tsv"
    bench = BenchResources(dataset=dataset)
    pk_df = pd.read_csv(bench.BENCH_WT)
    
    show_cols = ["PDB ID","Chain","Source", "Res#", "Expt. method"]
    pdbids = pk_df[show_cols].sort_values(by="PDB ID")
    pdbids = pdbids.drop_duplicates()
    pdbids = pdbids.rename(columns={"PDB ID": "PDB"})
    pdbids["PDB"] = pdbids["PDB"].astype("O")
    pdbids.to_csv(bench.BENCH_DB.joinpath(prots_chn), index=False, sep="\t")

    p_uniq = pdbids.drop("Chain", axis=1)
    p_uniq = p_uniq.drop_duplicates()
    p_uniq["Protein_type"] = "None"
    p_uniq["Molecule"] = "None"
    p_uniq["Biounits"] = "None"
    p_uniq["Res#"] = 0
    p_uniq["Reason_excluded"] = "None"
    p_uniq["Method"] = "None"
    p_uniq["Bioassembly_id"] = 1

    # remove entries that are both included and excluded;
    # case where a row (res) was excluded in pkas.csv:
    ids = p_uniq.PDB.tolist()
    for p in ids:
        xp = f"#{p}"
        ok = [p_uniq.PDB == p]
        xmsk = p_uniq.PDB == xp
        x = p_uniq[xmsk]
        if len(ok) and len(x):  # drop erroneous # (came from a single row)
            p_uniq = p_uniq.drop(x.index)

    p_uniq = p_uniq[["PDB","Protein_type", "Biounits","Res#","Reason_excluded",
                     "Method", "Bioassembly_id", "Source"]]
    p_uniq["Res#"] = p_uniq["Res#"].astype(int)
    p_uniq["PDB"] = p_uniq["PDB"].astype("O")
    p_uniq = p_uniq.sort_values(by="PDB")
    p_uniq = p_uniq.drop_duplicates()
    p_uniq.to_csv(bench.BENCH_PROTS, index=False, sep="\t")
    print(f"Unique PDB in {dataset}: {len(p_uniq)}")

    return


def download_assemblies(args):
    pkdb = get_pkdb(args)
    prot_tsv_path = get_proteins_file(args)
    assemblies_dir = get_runs_dir(args)

    df = read_tsv(prot_tsv_path)
    msk = df.PDB.str.startswith("#")
    df = df[~msk]
    
    os.chdir(str(assemblies_dir))
   
    failed_downloads = []
    for i, ro in df.iterrows():
        pdbid = ro.PDB
        bio = int(ro.Bioassembly_id)
        pdb_dir = Path(ro.PDB)
        pdb_dir.mkdir(exist_ok=True)

        os.chdir(str(pdb_dir))
        if not Path(f"{pdbid.lower()}.pdb").exists():
            print(i, "Missing, to download:", pdbid)
            try:
                out = downloads.get_rcsb_pdb(pdbid, bioassembly_id=bio)
            except Exception:
                failed_downloads.append(pdbid)
                continue
        os.chdir("../")

    print("Tot folders:", i, "; FAILED DOWNLOADS:\n", failed_downloads)
    
    return


def run_protinfo(args):
    pkdb = get_pkdb(args)
    prot_tsv_path = get_proteins_file(args)
    runs_dir = get_runs_dir(args, create=False)

    df = read_tsv(prot_tsv_path)
    msk = df.PDB.str.startswith("#")
    df = df[~msk]

    pi_args = {
    "pdb": "xxx.pdb",
    "d": 4,
    "e": "mcce",
    "noter": False,
    "u": "",
    "wet": False,
    "fetch": False,
    "save_dicts": True,
    }

    os.chdir(str(runs_dir))

    failed_preruns = []
    for i, pdbid in enumerate(df["PDB"]):
        pdb_dir = Path(pdbid)
        pdb_dir.mkdir(exist_ok=True)
        os.chdir(str(pdb_dir))
        
        pi_args["pdb"] = f"{pdbid.lower()}.pdb"
        try:
            get_pdb_rpt_climode(pi_args, do_fetch=False)
        except (FileNotFoundError, ZeroDivisionError):
            failed_preruns.append(pdbid)
    
        os.chdir("../")

    print("FAILED RUNS\n", failed_preruns)
    
    return


def txt2dict(txt_fp: Path) -> dict:
    """Open a txt file assumed to contain a printed dict
    and evaluate the text as such.
    """
    if not txt_fp.exists():
        print(f"No saved dict to txt file found in {txt_fp.parent!s}.")
        return None
    
    d = eval(txt_fp.read_text())
    if not isinstance(d, dict):
        print(f"Could not recover dict from txt file {txt_fp!s}")
        return None

    return d


def get_res_tot(seqres: str) -> str:
    """Get the residues total from the Seqres Species value.
    """
    # Any totals for DNS, RNA?
    oxy, doxy = None, None
    try:
        doxy = seqres.index("; dNTPs")
    except ValueError:
        pass
    try:
        oxy = seqres.index("; NTPs")
    except ValueError:
        pass

    if (doxy, oxy) == (None, None):
        return seqres.split("Total: ")[1]
    if oxy is None:
        nuc = doxy
    elif doxy is None:
        nuc = oxy
    elif doxy < oxy:
        nuc = doxy
    else:
        nuc = oxy

    return seqres[:nuc].split("Total: ")[1]


def prerun_data_dict(pkdb: Path, pdb: str) -> dict:
    """
    Process specific keys from protinfo saved dicts into a dict
    that will be used to update a specific row in proteins.tsv,
    which is identified by pdb.
    """
    # Load saved dict from its txt file
    pre_dir =  pkdb.joinpath("runs", pdb, "prerun")
    prot_d_fp = pre_dir.joinpath("prot_d.txt")
    prot_d = txt2dict(prot_d_fp)
    if not prot_d:
        # structure could not be parsed, no s1 either
        return None
    
    A = "\u212B"
    cols = ["Protein_type", "Function", "Biounits", "Res#", "Reason_excluded",
            "Method", "Resolution"]
    # preset out dict with None to preserve order:
    out = {k: None for k in cols}
    out["Protein_type"] = ""
    out["Reason_excluded"] = ""  # reset as string

    no_s1 = False
    tot_res = 0

    pdb_struc = prot_d["PDB.Structure"]
    if not pdb_struc.get("HEADERLESS"):
        # reasons for exclusion: unusable, multimodels (no step1):
        no_s1 = (pdb_struc.get("UNUSABLE") is not None
                 or pdb_struc.get("MultiModels") is not None
                 )
        if pdb_struc.get("UNUSABLE"):
            out["Reason_excluded"] += "Unusable: caveat or obsolete; "
        if pdb_struc.get("MultiModels"):
            out["Reason_excluded"] += "Multiple models; "

        out["Protein_type"] = pdb_struc.get("Molecule")
        out["Function"] = pdb_struc.get("Function")

        #print(f"{pdb = }: {pdb_struc.get('Biounits') = }")
        out["Biounits"] = pdb_struc.get("Biounits")
        if pdb_struc.get("Seqres Species"):
            tot_res = int(get_res_tot(pdb_struc["Seqres Species"]))
        out["Res#"] = tot_res
        out["Method"] = pdb_struc.get("Method")
        resol = pdb_struc.get("Resolution")
        if resol:
            if "ANGSTROMS." in resol:
                out["Resolution"] = resol.replace("ANGSTROMS.", A)
            else:
                out["Resolution"] = "N.A."  # not applicable
    
    if no_s1:
        return out

    s1_d_fp = pre_dir.joinpath("step1_d.txt")
    s1_d = txt2dict(s1_d_fp)
    if not s1_d is None:
        tpls = []
        pdb_s1 = s1_d["MCCE.Step1"]
        # process no tpl ("Labeling" key) -> Reason_excluded column
        if pdb_s1.get("Labeling") is not None:
            if pdb_s1["Labeling"]:
                # 1st item is a header, skip:
                tpls = [tpl.split(":", maxsplit=1)[0] for tpl in pdb_s1['Labeling'][1].strip().split("; ")]
                if tpls:
                    out["Reason_excluded"] += "no tpl: " + ", ".join(tpls) + ";"

    return out


def update_proteins_per_prerun(args):
    """
    Get info from the prerun reports to update proteins.tsv with.
    Note:
      Function MUST be run after running `protinfo` on a dataset.
    """
    pkdb = get_pkdb(args)
    runs_dir = get_runs_dir(args, create=False)
    prot_tsv_path = get_proteins_file(args)
 
    df = read_tsv(prot_tsv_path)
    msk = df.PDB.str.startswith("#")
    updated = False

    for i, ro in df[~msk].iterrows():
        pdb = ro.PDB
        upd_d = prerun_data_dict(pkdb, pdb)
        if upd_d:
            updated |= True
            for k in upd_d:
                df.loc[i, k] = upd_d[k]

                if k == "Reason_excluded" and len(upd_d[k]):
                    if not args.mark_excluded:
                        # always exclude prots
                        # multi models, unusable, with nucleic acids, unknown ligands:
                        if re.findall(r"_D?[AGTC]|Unusable|UNL|Multiple", upd_d[k]):
                            if not pdb.startswith("#"):
                                df.loc[i, "PDB"] = f"#{pdb}"
                        continue
                    if not pdb.startswith("#"):
                        df.loc[i, "PDB"] = f"#{pdb}"

    if updated:
        df.to_csv(prot_tsv_path, index=False, sep="\t")

    print("Proteins update over.")

    return


def update_pkfile_excluded_prots(args):
    """
    Update the pkas file with excluded pdbs from proteins.tsv.
    To be run after updating proteins.tsv with protinfo dicts.
    """
    prot_file = BenchResources(args.dataset_name).BENCH_PROTS
    xprots = audit.proteins_df(prot_file, return_excluded=True)
    if xprots.empty:
        print(f"There are no excluded proteins in {prot_file!s}.")
        return
    excluded = [prot[1:] for prot in xprots.PDB.tolist()]
   
    pka_file = BenchResources(args.dataset_name).BENCH_WT
    pkas = pd.read_csv(pka_file)
    # get the pdbids to exclude via the comment char, #:
    msk = pkas["PDB ID"].isin(excluded)
    for i in pkas[msk].index:
        pdb = pkas.loc[i,"PDB ID"]
        if not pdb.startswith("#"):
            pkas.loc[i,"PDB ID"] = f"#{pdb}"
    
    pkas.to_csv(pka_file, index=False)
    print("Updated pkas file with excluded proteins.")

    return


CLI_NAME = "bench_datasetprep"
SUB1 = "get_biopdbs"
SUB2 = "get_protinfo"
SUB3 = "postprerun_proteins"
SUB4 = "postprerun_pkas"

def cli_parser():
    p = ArgumentParser(
        formatter_class=RawDescriptionHelpFormatter,
        epilog="""
        Post an issue for all errors and feature requests at:
        https://github.com/GunnerLab/MCCE4/issues
        """,
    )
    # common parser
    cp = ArgumentParser(add_help=False)
    cp.add_argument(
        "dataset_name",
        type=str,
        help=f"""The name the experimental dataset for which to obtain the pdbs;
        One of: {dataset_names}.
        """,
    )
    cp.add_argument(
        "-proteins_file",
        default="proteins.tsv",
        type=str,
        help=f"The name of the dataset proteins file; default: proteins.tsv.",
    )
    subparsers = p.add_subparsers(
        required=True,
        title=f"{CLI_NAME} sub-commands",
        dest="subparser_name",
        description="Sub-commands for preping a dataset.",
        help=f"""
          1) Get the (bioassembly) pdbs: {SUB1}
          2) Run protinfo on the pdbs: {SUB2}
          3) Update proteins file with protinfo data: {SUB3}
          4) Update the pkas file with excluded pdbs: {SUB4}
        """,
    )

    sub1 = subparsers.add_parser(
        SUB1,
        formatter_class=RawDescriptionHelpFormatter,
        parents=[cp],
    )
    sub1.set_defaults(func=download_assemblies)

    sub2 = subparsers.add_parser(
        SUB2,
        formatter_class=RawDescriptionHelpFormatter,
        parents=[cp],
    )
    sub2.set_defaults(func=run_protinfo)

    sub3 = subparsers.add_parser(
        SUB3,
        formatter_class=RawDescriptionHelpFormatter,
        parents=[cp],
    )
    sub3.add_argument(
        "--mark_excluded",
        default=False,
        action="store_true",
        help="""Exclude a pdb based on the 'Reason_excluded' column in proteins.tsv.
        This is False by default as the snase dataset would have only 1 protein left."""
    )
    sub3.set_defaults(func=update_proteins_per_prerun)

    sub4 = subparsers.add_parser(
        SUB4,
        formatter_class=RawDescriptionHelpFormatter,
        parents=[cp],
    )
    sub4.set_defaults(func=update_pkfile_excluded_prots)

    return p


def assemblies_cli(argv=None):
    """
    Command line interface for preping a dataset.
    """
    p = cli_parser()
    args = p.parse_args(argv)
    if args.dataset_name not in datasets_dict:
        sys.exit(f"Unknown dataset name: {args.dataset_name}")
    if args.subparser_name == SUB2:
        # add option"
        args.save_dicts = True
    args.func(args)

    return

if __name__ == "__main__":
     assemblies_cli(sys.argv[:1])
