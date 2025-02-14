#!/usr/bin/env python

"""
Module: info.py

Holds functions to gather information about a protein obtained
from protinfo.bio_parser and protinfo.log_parser:
 - Info about the input protein from mcce4.pdbio
 - Info from MCCE step1 run1.log, debug.log when step1 can be run

Main steps:
  1. validate_pdb_inputs
  2. queries.get_rcsb_pdb(pdb):
     download bio-assembly if pdb is id
  3. collect_info in two dicts:
     first: from pdbio, second: from run1.log parser
  4. collect_info_lines from the two dicts
"""


from argparse import Namespace
from collections import defaultdict
import logging
from pathlib import Path
from pprint import pformat
from typing import Tuple, Union
from mcce4.protinfo import USER_MCCE, RPT
from mcce4.protinfo import parsers, run

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def hetero_count_per_chain(atomlines: list) -> Union[dict, None]:
    """Given a list of pdb coordinates lines, parse the hetero atoms to
    return a dictionary holding the count of free cofactors and waters,
    per chain (key).
    """
    # get unique hetero species per chain:
    hetero_set = defaultdict(set)
    for line in atomlines:
        if not line.startswith("HETATM"):
            continue
        # name = line[17:20]; chn = line[21]; seq = int(line[22:26])
        hetero_set[line[21]].add((line[17:20], int(line[22:26])))
    if not hetero_set:
        return None
    # get species count per chain:
    tots_per_chain = defaultdict(dict)
    for c in hetero_set:
        cntr = defaultdict(int)
        for val in hetero_set[c]:
            cntr[val[0]] += 1
        tots_per_chain[c]=  dict(cntr)
        
    return dict(tots_per_chain)


def get_cofactors_change(pdb_heteros: dict, step1_heteros: dict) -> list:
    """Return a list of changes in cofactors counts between the starting structure
    dictionary entry key prot_d["PDB.Structure"]["Model 1 Free Cofactors & Waters"] passed in pdb_heteros
    and step1_d["MCCE.Step1"]["Free Cofactors"] passed in step1_heteros.
    """
    # check what was removed:
    diff = []
    for chn in pdb_heteros:
        for k in pdb_heteros[chn]:
            if step1_heteros.get(chn) is None:
                continue
            if step1_heteros[chn].get(k) is None:
                diff.append(f"Removed all {pdb_heteros[chn][k]} {k} in {chn}.")
            else:
                # get non-zero count difference
                k_rem = step1_heteros[chn][k]
                d = pdb_heteros[chn][k] - k_rem
                if d:
                    diff.append(f"Removed {d} {k} in {chn}; {k_rem} remaining.")

    return diff


def update_s1_dict_cofactors_change(pdb: Path, prot_d: dict, step1_d: dict):
    """Update step1_d["MCCE.Step1"]["Free Cofactors"] list with cofactor changes if any.
    """
    if prot_d["PDB.Structure"].get("Model 1 Free Cofactors & Waters") is None:
        return

    pdb_heteros = prot_d["PDB.Structure"]["Model 1 Free Cofactors & Waters"]

    # get cofactors count from step1 pdb:
    s1_pdb = pdb.parent.joinpath("step1_out.pdb")
    s1_heteros = hetero_count_per_chain(s1_pdb.read_text().splitlines())
    if s1_heteros is None:
        return

    diff = get_cofactors_change(pdb_heteros, s1_heteros)
    if diff:
        step1_d_heteros = []
        idx = 0
        # initial list
        step1_d_heteros = step1_d["MCCE.Step1"]["Free Cofactors"].copy()
        if step1_d_heteros[idx].startswith("Total deleted cofactors"):
            idx = 1
        changes = " ".join(h1 for h1 in diff)
        if step1_d_heteros[idx] != changes:
            step1_d_heteros.insert(idx, changes)
            step1_d["MCCE.Step1"]["Free Cofactors"] = step1_d_heteros

    return


def collect_info(pdb: Path, args: Namespace) -> Tuple[dict, Union[dict, None]]:
    """Return at least one dict holding info from mcce4.pdbio.
    The second dict is None when step1 cannot be run, otherwise
    it contains info from the parsed run1.log file.
    Args:
      pdb (Path): File path of validated pdb.
      args (argparse.Namespace): Cli arguments (for creating step1 script).
    Returns:
      A 2-tuple of dicts when step1 can run, else (dict1, None).
    """
    step1_d = None

    # structural info:
    prot_d = parsers.info_input_prot(pdb)

    is_multi = prot_d["PDB.Structure"].get("MultiModels") is not None
    is_malformed = prot_d["PDB.Structure"].get("Malformed PDB") is not None
    DO_STEP1 = USER_MCCE is not None and not is_multi and not is_malformed

    if DO_STEP1:
        result = run.do_step1(pdb, args)
        if result is None:  # no error message
            step1_d = parsers.info_s1_log(pdb)
            # update cofactors section changes if any:
            update_s1_dict_cofactors_change(pdb, prot_d, step1_d)
        else:
            step1_d = {"MCCE.Step1": f"Error in run script 's1.sh': {result}."}

    return prot_d, step1_d


def write_report(pdb: Path, prot_d: dict, s1_d: Union[dict, None]):
    """Write the prerun report for pdb in its parent folder using the
    information passed in the parsers dicts, prot_d and step1_d.
    Args:
      pdb (Path): the pdb filepath; expected: path of pdb in prerun folder.
      prot_d (dict): The dictionary of sections from mcce4.pdbio.
      s1_d ([dict, None]): The dictionary of sections from step1 log parser.
    """
    name = prot_d.pop("Name")
    if s1_d is None:
        dict_lst = [prot_d]
    else:
        dict_lst = [prot_d, s1_d]

    # save report in the call folder (e.g. ../prerun/.)
    rpt_fp = pdb.parent.parent.joinpath(f"{pdb.stem}_{RPT}")
    with open(rpt_fp, "w") as rpt:
        rpt.write(f"---\n# {name}\n")

        for i, subd in enumerate(dict_lst):
            # h2: section hdrs, PDB.Structure or MCCE.Step1
            h2 = list(subd.keys())[0]
            rpt.write(f"## {h2}\n")

            for k in subd[h2]:
                if not subd[h2][k]:
                    continue

                if i == 0:
                    if isinstance(subd[h2][k], (str, int, float)):
                        rpt.write(f"### {k}: {subd[h2][k]}\n")
                        continue

                    rpt.write(f"### {k}:\n")
                    if isinstance(subd[h2][k], list):
                        for val in subd[h2][k]:
                            rpt.write(f"  - {val}\n")
                    elif isinstance(subd[h2][k], dict):
                        for kk, vv in subd[h2][k].items():
                            if isinstance(vv, (dict, tuple)):
                                out = pformat(vv, sort_dicts=True, compact=True, width=160)[1:-1]
                                rpt.write(f"  - {kk}:\n {out}\n")
                            else:
                                rpt.write(f"  - {kk}: {vv}\n")
                else:
                    rpt.write(f"### {k}:\n")
                    for val in subd[h2][k]:
                        if k == "Distance Clashes":
                            if isinstance(val, str) and val.startswith("Clashes"):
                                rpt.write(f"<details><summary>{val}</summary>\n\n")
                            elif isinstance(val, str) and val.endswith("end_clash"):
                                rpt.write("\n</details>\n")
                            else:
                                rpt.write(f"- {val}\n")

                        elif k == "Labeling":
                            if val.startswith("Generic") or val.startswith(
                                "Unloadable"
                            ):
                                rpt.write(
                                    f"<strong><font color='red'>{val}:</font></strong>  \n"
                                )
                            elif val.startswith("Likely"):
                                val = val.replace(
                                    "Likely cause",
                                    "<strong><font color='red'>Likely cause</font></strong>",
                                )
                                rpt.write(f"{val}  \n")
                            else:
                                rpt.write(f"{val}\n")
                        else:
                            if isinstance(val, tuple):
                                ter, lst = val
                                rpt.write(
                                    f" - <strong>{ter}</strong>: {', '.join(lst)}\n"
                                )
                            elif isinstance(val, list):
                                ter, lst = val
                                rpt.write(
                                    f" - <strong>{ter}</strong>: {', '.join(x for x in lst)}\n"
                                )
                            elif isinstance(val, dict):
                                for kk, vv in val.items():
                                    rpt.write(f"  - {kk}: {vv}\n")
                            else:
                                rpt.write(f"  - {val}\n")

                rpt.write("\n")

    return
