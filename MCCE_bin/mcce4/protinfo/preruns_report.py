#!/usr/bin/env python3

"""Module: preruns_report.py

To create a reduced status report on a collection of prerun (protinfo)
folders with dictionaries, which are saved when the --save_dicts option
is given to the cli or in args to the cli.get_pdb_rpt function.
"""

from collections import defaultdict
from pathlib import Path

import pandas as pd


PRERUN_RPT = "preruns_report.tsv"


def get_lst_size(ro) -> int:
    """Function for df.apply
    """
    if ro["MissingTpl"] == 0:
        return 0
    if isinstance(ro["MissingTpl"], list):
        n = len(ro["MissingTpl"])
    else:
        n = 1
        
    return n


def extract_struc_info(struc_d: dict) -> dict:
    """Output dict keys:
    ["Name", "Function", "Cofactors", "TotalRes, "Models", "Unusable"]
    """
    pdb_struc = {}
    # Nmae: '1B5E Dcmp Hydroxymethylase From T4'
    name = struc_d["Name"].split(maxsplit=1)[1].strip()  
    pdb_struc["Name"] = name
    pdb_struc["Function"] = struc_d["PDB.Structure"].get("Function")
    
    cofac = struc_d["PDB.Structure"].get("Cofactors")
    if cofac is None:
        pdb_struc["Cofactors"] = None
    else:
        colst = []
        for k in cofac.keys():
            val = cofac[k]
            colst.append(f"{k}: {val[1]} ({val[0]})")
        if len(colst) == 1:
            pdb_struc["Cofactors"] = colst[0]
        else:
            pdb_struc["Cofactors"] = colst

    res = struc_d["PDB.Structure"].get("Model 1 Residues")
    if res is None:
        pdb_struc["TotalRes"] = 0
    else:
        # to sum 'Total: ' values
        # ### Model 1 Residues:
        #  - A:
        # 'RESIDUES': 'A: 12, C: 8, D: 7, E: ... S: 10, T: 7, V: 6, W: 6, Y: 3; Total: 1>
        #     'Ratio: 29.5%'
        sum_res = 0
        for chn in res.keys():
            res_info = res[chn].get("RESIDUES")
            if res_info is None:
                continue
            sum_res += int(res_info.split("; ")[1].split(":")[-1])
        pdb_struc["TotalRes"] = sum_res

    mdls = struc_d["PDB.Structure"].get("Models")
    if mdls is None:
        pdb_struc["Models"] = None
    else:
        pdb_struc["Models"] = int(mdls)

    unusable = struc_d["PDB.Structure"].get("UNUSABLE")
    if unusable is None:
        pdb_struc["Unusable"] = False
    else:
        pdb_struc["Unusable"] = unusable
    
    return pdb_struc


def extract_s1_info(s1_d: dict) -> dict:
    """Output dict keys:
    ["MissingTpl", "TplInfo", "Step1"] (last one for Status)
    """
    pdb_s1 = {}
    no_data2 = {"MissingTpl": None, "TplInfo": None, "Step1": None}

    step1_data = s1_d.get("MCCE.Step1")
    if step1_data is None:
        pdb_s1.update(no_data2)
        return pdb_s1

    if not isinstance(step1_data, dict):
        pdb_s1.update({"MissingTpl": None, "TplInfo": None, "Step1": step1_data})
        return pdb_s1
    
    labeling_val = s1_d["MCCE.Step1"].get("Labeling")
    if labeling_val is None:
        pdb_s1.update({"MissingTpl": 0, "TplInfo": None, "Step1": "OK"})
    else:
        created, info = [], []
        for tpls in labeling_val[1:]:
            for tpl in tpls.split("; ")[:-1]:
                if not tpl:
                    continue
                tpl_created, tpl_info = tpl.split(": ")
                #print(f"{tpl = }:: {tpl_created}, {tpl_info}")
                created.append(tpl_created)
                info.append(tpl_info)
        if not created:
            pdb_s1.update({"MissingTpl": 0, "TplInfo": None, "Step1": "OK"})
        else:
            if len(created) == 1:
                pdb_s1.update({"MissingTpl": created[0], "TplInfo": info[0], "Step1": "OK"})
            else:
                pdb_s1.update({"MissingTpl": created, "TplInfo": info, "Step1": "OK"})

    status = s1_d["MCCE.Step1"].get("Status")
    if status is not None:
        pdb_s1.update({"Step1": status})
        
    return pdb_s1


def get_runs_reduced_info_dict(runs_dir: Path) -> dict:
    """
    - runs_dir: Path to proteins folders dir.
    """
    all_d = defaultdict(dict)
    no_data1 = {"Name":None, "Function":None, "Cofactors":None,
                "TotalRes":None, "Models":None, "Unusable":None}
    no_data2 = {"MissingTpl": None, "TplInfo": None, "Step1": None}

    for d in runs_dir.iterdir():
        # if not d.stem.isupper():
        #     continue
        if not d.joinpath("prerun").is_dir():
            continue

        pdb = d.stem
        struc = d.joinpath("prerun", "prot_d.txt")
        if not struc.exists():
            all_d[pdb].update({"Name":None, "Function":None,
                               "Cofactors":None, "TotalRes":None,
                               "Models":None, "Unusable":"Not prot_d.txt"})
        else:
            struc_d = eval(struc.read_text())
            if isinstance(struc_d, dict):
                d_struc = extract_struc_info(struc_d)
                all_d[pdb].update(d_struc)
            else:
                if struc_d is not None:
                    print(f"Eval error: struc_d not dict & not None: {struc_d}")
                all_d[pdb].update(no_data1)

        # mcce step1 data:
        s1 = d.joinpath("prerun", "step1_d.txt")
        if not s1.exists():
            all_d[pdb].update({"MissingTpl": None, "TplInfo": None,
                               "Step1": "No step1_d.txt"})
        else:
            s1_d = eval(s1.read_text())
            if isinstance(s1_d, dict):
                d_s1 = extract_s1_info(s1_d)
                all_d[pdb].update(d_s1)
            else:
                if s1_d is not None:
                    print(f"Eval error: s1_d not dict or None: {s1_d}")
                all_d[pdb].update(no_data2)

    return all_d


def get_preruns_report(runs_dir: Path) -> pd.DataFrame:
    """
    - runs_dir: Path to proteins folders dir.
    """
    all_d = get_runs_reduced_info_dict(runs_dir)
    df = pd.DataFrame.from_dict(all_d, orient="index")
    df.index.name = "PDBID"
    df = df.reset_index()
    rpt_fp = runs_dir.joinpath(PRERUN_RPT)
    df.to_csv(rpt_fp, index=False, sep="\t")
    print(f"Pre-runs report: {rpt_fp!s}")

    return
