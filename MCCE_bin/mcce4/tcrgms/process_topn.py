#!/usr/bin/env python

"""
Module: process_topn.py

Hold processing functions to obtain a tautomeric, topN charge
microstate, tcrgms, and produce these outputs:
  - A csv file containing the topN vectors of resid and their charges;
  - PDB files (n=topN), using the conformers listed in each vector.
"""
from copy import deepcopy
from datetime import datetime
import logging
from pathlib import Path
from pprint import pformat
import subprocess
from typing import Any, TextIO, Tuple, Union
import numpy as np
import pandas as pd
from mcce4.tcrgms import APP_NAME, CANONICAL
from mcce4.tcrgms import s2_to_pdb as mc2pdb


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


# additional col per crg vector to hold iconf:
N_ADDED_COLS_PER_TOP_VEC = 1
S2PREF = "s2_"


conf2tauto = {
    "ASP01": "OD1",
    "ASP02": "OD2",
    "CTR01": "HO",
    "CTR02": "HXT",
    "GLU01": "OE1",
    "GLU02": "OE2",
    "HIS01": "NE2",  # "HIS01":"HSE"
    "HIS02": "ND1",  # "HIS02":"HSD"
}


def subprocess_run(
    cmd: str, capture_output=True, check: bool = False, text=True, shell=True
) -> Union[subprocess.CompletedProcess, subprocess.CalledProcessError]:
    """Wraps subprocess.run. Return CompletedProcess or err obj."""

    try:
        data = subprocess.run(
            cmd, capture_output=capture_output, check=check, text=text, shell=shell
        )
    except subprocess.CalledProcessError as e:
        data = e

    return data


def get_run_dry_opt(run_dir=Path) -> bool:
    """Get the dry status of the original run from run.prm.record."""
    dry = True
    runprm = run_dir.joinpath("run.prm.record")
    if not runprm.exists():
        logger.warning("File run.prm.record not found: dry run assumed.")
        return dry
    data = subprocess_run("grep H2O_SASCUTOFF " + str(runprm))
    if isinstance(data, subprocess.CompletedProcess):
        dry = float(data.stdout.split()[0]) == -0.01

    return dry


def get_non_canonical(
    topN_crg_d: dict, out_dir: Path
) -> Tuple[bool, Union[Path, None]]:
    """Save the non-canonical residues list to file if any are found
    and return (True, fpath) else return (False, None).
    """
    non_canonical = []

    for k in topN_crg_d:
        for ck in topN_crg_d[k]["confs"]:
            if ck.startswith("HOH"):
                continue

            _, crg = topN_crg_d[k]["confs"][ck]
            defres = CANONICAL.get(ck[:3])
            if defres is None:
                if ck[:5] == "HIS01":
                    non_canonical.append([k, ck[:3] + ck[5:], crg])
                continue
            if crg != defres:
                non_canonical.append([k, ck[:3] + ck[5:], crg])

    if non_canonical:
        nc_fp = out_dir.joinpath("non_canonical_res.txt")
        nc = pd.DataFrame(non_canonical, columns=["cms", "res", "crg|tauto"]).to_string(
            index=False
        )
        nc_fp.write_text(nc + "\n")
        return True, nc_fp
    else:
        return False, None


def multi_tauto_conflist(param_dir: Union[str, Path]) -> dict:
    """Return a dict with 3-char resid as key and the list from
    a ftpl CONFLIST lines as value.
    """
    tauto_conflist = {}
    if isinstance(param_dir, Path):
        param_dir = str(param_dir)

    # collect all the conflist lines if there are more than 2 confs:
    tp_kind = "ftpl"
    if "0" in param_dir:
        tp_kind = "tpl"

    cmd = (
        """awk -v n=3 '/CONFLIST/ {$1=""; $3=""; if (NF-n >= 3) {print $0}}'"""
        f" {param_dir}/*.{tp_kind}"
    )
    out = subprocess_run(cmd)
    if out is subprocess.CalledProcessError:
        # FIX: Maybe abort instead?
        logger.error(
            (
                "Error fetching tautomeric conformers from ftpl files."
                "Output of function is an empty dict."
            )
        )
        return tauto_conflist

    for line in out.stdout.splitlines():
        k, v = line.split(":")
        tauto_conflist[k.strip()] = v.strip()

    return tauto_conflist


def get_topN_dict(
    topN_lst: list,
    conformers: list,
    no_trailing_underscore: bool = True,
    dry: bool = False,
) -> dict:
    """
    Initial implementation.
    Format of confermers ids w/o trailing underscore if no_trailing_underscore is True (default),
    which matches dict output of ms_analysis.fixed_res_crg with no_trailing_underscore set True
    by cli.
    """
    topN_crg_d = {}
    upto = -4 if no_trailing_underscore else -3

    for i, ro in enumerate(topN_lst, start=1):
        dd = {}
        for ic in ro[0][0]:  # state
            conf = conformers[ic]
            if dry and conf.confid.startswith("HOH"):
                continue
            dd[conf.confid[:upto]] = conf.iconf, int(conf.crg)
        topN_crg_d[i] = {
            "confs": dd,
            "E": round(ro[0][1], 0),
            "sum_crg": -99,
            "count": ro[1],
            "occ": ro[2],
        }

    return topN_crg_d


def process_waters(raw_topN_crg_d: dict) -> dict:
    """
    Remove dummy waters.
    Add NA entries for missing HOHs from the water set the to obtain the same number of waters
    in all cms.
    Return the updated copy of the input dict.
    """
    topN_crg_d = deepcopy(raw_topN_crg_d)
    hoh_set = set()

    # 1st pass: populate the waters set & remove dummies:
    for k in raw_topN_crg_d:
        for ck in raw_topN_crg_d[k]["confs"]:
            if not ck.startswith("HOH"):
                continue
            if ck.startswith("HOHDM"):
                _ = topN_crg_d[k]["confs"].pop(ck)
                continue
            hoh_set.add(ck)
    # 2nd pass: add an entry for waters not in the cms:
    for k in topN_crg_d:
        hoh_confs = set([ck for ck in topN_crg_d[k]["confs"] if ck.startswith("HOH")])
        for h in hoh_set.difference(hoh_confs):
            topN_crg_d[k]["confs"].update({h: (np.nan, np.nan)})

    return topN_crg_d


def crg_d_with_tauto(crg_conf_d: dict, tauto_conflist: dict, conformers: list) -> dict:
    """Return the tautomeric form of a conformer in the input dict with
       all confids transformed into resid.
    Notes:
      crg_conf_d source:  topN_crg_d[n]["confs"]: {'ASP-1A0119':(247, -1),..}
    Returns:
      The updated input dict where the crg value is either the original crg or
      the tautomer identifier, if any.
    """
    for k in crg_conf_d:
        res = k[:3]
        conflist = tauto_conflist.get(res)
        if conflist is None:
            continue
        # get confid:
        v = crg_conf_d[k]  # (iconf, crg)
        confid = conformers[v[0]].confid[:5]
        tauto = conf2tauto.get(confid)
        if tauto is None:
            continue
        # update crg in input list:
        crg_conf_d[k] = v[0], tauto

    return crg_conf_d


def update_topNcrg_d_tauto(
    topNcrg_d: dict, tauto_conflist: dict, conformers: list, n_top: int
) -> dict:
    """Update each of the 'confs' subkey dicts with tautomer id."""
    for k in topNcrg_d:
        topNcrg_d[k]["confs"] = crg_d_with_tauto(
            topNcrg_d[k]["confs"], tauto_conflist, conformers
        )

    return topNcrg_d


def get_df_shape_cols(n_fixed: int, topN_crg_d: dict) -> tuple:
    """Return a shape tuple given the inputs and defined columns."""
    n_top = len(topN_crg_d)
    cols = (
        ["residues"]
        # N_ADDED_COLS_PER_TOP_VEC = 1, source:
        + " ".join(f"iconf_{i} {i}" for i in range(1, n_top + 1)).split()
        + ["info"]
    )
    shape = (
        n_fixed + len(topN_crg_d[1]["confs"]) + len(topN_crg_d[1].keys()) - 1,
        n_top + n_top * N_ADDED_COLS_PER_TOP_VEC + 2,
    )

    return shape, cols


def get_preset_df(
    shape: tuple, colnames: list = None, fill_val: Any = "x", fill_type: type = str
) -> pd.DataFrame:
    """Create an empty pandas.DataFrame of given shape with cells
    filled with fill_val.
    """
    a = np.empty(shape, dtype=fill_type)
    data = np.full_like(a, fill_value=fill_val, dtype=fill_type)
    if colnames is None:
        return pd.DataFrame(data)

    return pd.DataFrame(data, columns=colnames)


def topNdf_to_tsv(output_dir: Path, top_df: pd.DataFrame, n_top: int):
    """Save topN_df to 2 tab-separated-values (.tsv) files in output_dir:
    1. The input df is saved as-is to 'top<n>_master_tcgrms.tsv';
    2. The df w/o 'iconf_<n>' columns is saved to 'top<n>_tcgrms.tsv'.
       This file is also referred to as 'the user file'.
    """
    tsv_fp = output_dir.joinpath(f"top{n_top}_master_{APP_NAME}.tsv")
    top_df.to_csv(tsv_fp, sep="\t", index=False)
    tsv_fp = output_dir.joinpath(f"top{n_top}_{APP_NAME}.tsv")
    keep_cols = [c for c in top_df.columns if not c.startswith("iconf_")]
    top_df[keep_cols].to_csv(tsv_fp, sep="\t", index=False)

    return


def get_chain(df: pd.DataFrame, n_exclude_last_rows: int = 4) -> list:
    out = df.residues.str[3]
    out[-n_exclude_last_rows:] = ""

    return out


def get_crg_vec_df(fixedres_crg_d: dict, topN_crg_d: dict, n_ms: int) -> pd.DataFrame:
    """
    Creates a charge vector DataFrame.
    Returns a "master df" with cols (e.g. if n_top=3):
      [residues, iconf_1, 1, iconf_2, 2, iconf_3, 3, info]
      Columns "iconf_n" are not output in the topN user file;
      they are in the 'master df' file, topN_master.tsv.

    REM: Structure of topN_crg_d at start:
        topN_crg_d[1] = {"confs": dd,    # dict
                         "E": E,
                         "sum_crg": -99,  # not yet computed
                         "count": count,
                         "occ": occ
                         }
    """
    n_fixed = len(fixedres_crg_d)
    shape, cols = get_df_shape_cols(n_fixed, topN_crg_d)
    df = get_preset_df(shape, colnames=cols)

    for r in range(1, len(topN_crg_d) + 1):
        sum_crg_fx = 0
        sum_crg_fr = 0
        # add fixed res
        for i, k in enumerate(fixedres_crg_d):
            if r == 1:  # do only once
                df.loc[i, "residues"] = k
                df.loc[i, "info"] = "fixed"
            iconf, crg = fixedres_crg_d[k]
            df.loc[i, f"iconf_{r}"] = iconf
            df.loc[i, str(r)] = crg
            sum_crg_fx += int(crg)
        i += 1
        # add free res+
        for j, kk in enumerate(topN_crg_d[r]["confs"]):
            tup = topN_crg_d[r]["confs"][kk]
            # match fixed res format:
            df.loc[i + j, "residues"] = kk[:3] + kk[5:]
            df.loc[i + j, f"iconf_{r}"] = tup[0]
            df.loc[i + j, str(r)] = tup[1]
            df.loc[i + j, "info"] = "free"
            try:
                sum_crg_fr += int(tup[1])
            except (TypeError, ValueError):
                # assumed: neutral tautomers!
                pass
        q = i + j + 1
        # add totals
        freekeys = list(topN_crg_d[r].keys())
        freekeys.remove("confs")
        for o, ok in enumerate(freekeys):
            df.loc[q + o, "info"] = "totals"
            df.loc[q + o, "residues"] = ok
            if ok == "sum_crg":
                df.loc[q + o, str(r)] = sum_crg_fx + sum_crg_fr
            else:
                if ok == "occ":
                    df.loc[q + o, str(r)] = f"{topN_crg_d[r][ok]:.2%}"
                else:
                    df.loc[q + o, str(r)] = f"{topN_crg_d[r][ok]:,.0f}"
    # sort:
    msk = df["info"] == "totals"
    df = pd.concat([df[~msk].sort_values(by="iconf_1"), df[msk]], ignore_index=True)
    df["chain"] = get_chain(df)

    return df


def get_pdb_remark(remark_data: dict):
    """
    Return a REMARK 250 header to prepend the final pdb with.
    Args:
      remark_data (dict): Keys: INPDB, T, PH, E, SUM_CRG, COUNT, OCC;
      the values are assumed to be strings reflecting appropriate formats.

    > REMARK 250 is mandatory if other than X-ray, NMR, neutron, or electron study.
    [Ref]: https://www.wwpdb.org/documentation/file-format-content/format33/remarks1.html#REMARK%20250
    """

    R250 = """REMARK 250
REMARK 250 EXPERIMENTAL DETAILS
REMARK 250   EXPERIMENT TYPE               : MCCE simulation; MCCE4 {TOOL} tool.
REMARK 250   DATE OF DATA COLLECTION       : {DATE}
REMARK 250   REMARK: DATE OF DATA COLLECTION is the creation date of this pdb.
REMARK 250 EXPERIMENTAL CONDITIONS
REMARK 250   SIMULATION INPUT PDB          : {INPDB}
REMARK 250   TEMPERATURE                   : {T} (K)
REMARK 250   PH                            : {PH}
REMARK 250 CHARGE MICROSTATE INFORMATION
REMARK 250   ENERGY                        : {E} (kcal/mol)
REMARK 250   NET CHARGE                    : {SUM_CRG}
REMARK 250   COUNT                         : {COUNT}
REMARK 250   OCCUPANCY                     : {OCC}
REMARK 250 REMARK:
REMARK 250  This pdb was created from a tautomeric charge microstate vector
REMARK 250  extracted by {TOOL}, a MCCE4 tool.
REMARK 250
"""
    return R250.format(
        TOOL=APP_NAME, DATE=datetime.today().strftime("%d-%b-%y"), **remark_data
    )


def confs_to_pdb(
    step2_fh: TextIO, selected_confids: dict, output_pdb: str, remark_data: dict
) -> None:
    """Read step2_out coordinate line for each conformer in list `selected_confs`
    and creates a pdb file in step2 format.

    Args:
    step2_fh (TextIO): File handle of 'step2_out.pdb'.
    selected_confids (dict): A microstate's dict of conformer ids.
        Note: dict is for matching confIDs, values are all 1 (i.e. True & unused).
    output_pdb (str): Output pdb file_path.
    """
    remark = get_pdb_remark(remark_data)
    with open(output_pdb, "w") as out:
        out.write(remark)
        for line in step2_fh:
            if len(line) < 82:
                continue
            if line[80:82] == "BK":
                out.write(line)
                continue
            confID = line[17:20] + line[80:82] + line[21:26] + "_" + line[27:30]
            # ATOM     14  CB  LYS A0001_001   1.180   5.987  12.487   2.000       0.000      01O000M000
            # 012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
            #          1         2         3         4         5         6         7         8

            # standard pdb line:                  x       y       z  occ  B_iso
            # ATOM     14  CB  LYS A   8       5.987  11.822  24.520  0.57 13.03           C   <
            # 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
            # 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
            # 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
            # 55 - 60        Real(6.2)     occupancy    Occupancy.
            # 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
            # 77 - 78        LString(2)    element      Element symbol, right-justified.
            # 79 - 80        LString(2)    charge       Charge  on the atom.

            if selected_confids.get(confID) is not None:
                out.write(line)

    return


def write_tcrgms_pdbs(
    out_dir: Path,
    inpdb_fname: str,
    tcrgms_df: pd.DataFrame,
    conformers: list,
    remark_args: dict,
    n_top: int,
    name_pref: str = S2PREF,
) -> None:
    """
    Write the topN tcms pdbs from the master df iconf_N columns.
    """
    # split df viz "totals" rows:
    msk = tcrgms_df["info"] == "totals"
    df_tots = tcrgms_df[msk].set_index("residues")
    df = tcrgms_df[~msk]

    # get effective topN size from df:
    n_datacols = int((len(df.columns) - 3) / 2)

    pdb_name = Path(inpdb_fname).stem
    step2_fh = open(out_dir.parent.joinpath("step2_out.pdb"))

    for i in range(1, n_datacols + 1):
        si = str(i)
        R250_d = dict((k, v) for k, v in remark_args.items())
        R250_d["E"] = df_tots.loc["E", si]
        R250_d["SUM_CRG"] = df_tots.loc["sum_crg", si]
        R250_d["COUNT"] = df_tots.loc["count", si]
        R250_d["OCC"] = df_tots.loc["occ", si]

        pdb_out = out_dir.joinpath(f"{name_pref}tcms{i}_{pdb_name}.pdb")

        # TODO?: use conf.occ instead of 1, in case its needed later:
        # confids_dict = dict((conformers[cid].confid, conformers[cid].occ) for cid in df[f"iconf_{i}"])
        confids_dict = dict(
            (conformers[cid].confid, 1) for cid in df[f"iconf_{i}"] if cid is not np.nan
        )
        confs_to_pdb(step2_fh, confids_dict, pdb_out, R250_d)

        step2_fh.seek(0)

    step2_fh.close()

    return


def mccepdbs_to_pdbs(out_dir: Path, s2_path: Path, rm_prefix: str = None) -> None:
    """Convert MCCE formatted pdbs in `out_dir` to pdb format.

    Args:
    out_dir (Path): Folder path for mcce pdbs
    s2_path (Path): Path to step2_out.pdb
    rm_prefix (str): If not None, the output file name won't have this prefix.
    """
    s2_rename_dict, chn_dict = mc2pdb.get_mcce_pdb_res_dict(s2_path)
    logger.debug(pformat(s2_rename_dict, sort_dicts=False))

    remove_prefix = rm_prefix is not None and rm_prefix
    glob_str = "*.pdb"
    if remove_prefix:
        glob_str = f"{rm_prefix}*.pdb"

    for fp in out_dir.glob(glob_str):
        converter = mc2pdb.Mcce2PDBConverter(
            fp, s2_rename_dict, chn_dict, output_dir=out_dir
        )
        out_name = f"{fp.stem}.pdb"
        if remove_prefix:
            out_name = f"{fp.stem.removeprefix(rm_prefix)}.pdb"
        converter.mcce_to_pdb(out_name)

    return
