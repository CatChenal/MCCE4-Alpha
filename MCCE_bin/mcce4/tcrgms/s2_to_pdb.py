#!/usr/bin/env python

"""
Module: s2_to_pdb.py

Contains functions and convert class to convert a pdb
in mcce step2_out format to a regular pdb.
"""

from collections import defaultdict
import logging
from pathlib import Path
from typing import Tuple, Union
from mcce4.constants import ALL_RES, COMMON_HETATMS


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


TER = ["NTR", "CTR"]
MCCE_FIELDS = 11


def get_mcce_pdb_res_dict(step2_path: Path) -> Tuple[dict, dict]:
    """Return a 2-tuple of dicts for res with multiple entries,
    i.e. terminal residues; and to detect chain change.
    Sample output, dict1: {'A0001': ['NTR', 'LYS'],
                           'A0129': ['LEU', 'CTR']}.
    Note: Order matters!
    """
    ter_dict = defaultdict(list)
    chn_dict = defaultdict(int)
    with open(step2_path) as fh:
        for i, line in enumerate(fh.readlines(), start=1):
            _, _, _, res, conf, *_ = line.split()
            res = res.upper()
            res_id = conf.split("_")[0]
            if res in ALL_RES:
                chn_dict[res_id[0]] = i
            # order needed (see L75), can't use set():
            if res not in ter_dict[res_id]:
                ter_dict[res_id].append(res)

    return dict((k, v) for k, v in ter_dict.items() if len(v) > 1), dict(chn_dict)


class Mcce2PDBConverter:
    """A class for converting MCCE PDB files to PDB file."""

    def __init__(self,
                 mccepdb_path: Path,
                 s2rename_dict: dict,
                 chain_dict: dict,
                 output_dir: Path):
        """
        Attributes:
          - mccepdb_path (Path): The path to the MCCE PDB file.
          - s2rename_dict (dict): Dict for renaming step2_out.pdb TER residues;
                              1st output of: get_mcce_pdb_res_dict(step2_path)
          - chain_dict (dict): Dict keeps track of the change of chains if any;
                              2nd output of: get_mcce_pdb_res_dict(step2_path)
          - output_dir (Path): Dir path for the converted file.
        """
        self.mccepdb_path = mccepdb_path
        self.rename_dict = s2rename_dict
        # values of chain_dict is the step2 line number (starting at 1) 
        self.chain_dict = chain_dict
        self.output_dir = output_dir

    def get_pdb_line(self, mcce_line: str):
        """Return the standard pdb line."""
        try:
            # check the line from step2_out.pdb has correct format:
            rec, idx, atm, res, conf, x, y, z, rad, crg, hist = mcce_line.split()
        except ValueError:
            logger.error(
                f"This mcce line is invalid (expected {MCCE_FIELDS} fields):\n\t{mcce_line}"
            )
            return

        res_index, conf_number = conf.split("_")
        res = res.upper()
        if res == "NTG":
            res = "GLY"
        else:
            # {'A0001': ['NTR', 'LYS'], 'A0129': ['LEU', 'CTR']}
            if self.rename_dict.get(res_index) is not None:
                if res in TER:  # list has 2 items
                    ix = int(res == "NTR")
                    res = self.rename_dict[res_index][ix]
            else:
                if res.replace("_", "") in COMMON_HETATMS:
                    rec = "HETATM"

        res_index = f"{res_index[0]}{res_index[1:].lstrip('0'):>4}"
        align = "^" if len(atm) < 3 else ">"
        # element = atm[0]
        pline = (
            f"{rec:<6s}{idx:>5s} {atm:{align}4s} {res:<3s} {res_index:>4s}{' ':4s}"
            f"{float(x):>8.3f}{float(y):>8.3f}{float(z):>8.3f}  1.00{' ':17s}{atm[0]:<2s}\n"
        )

        return pline

    def mcce_to_pdb(self, out_name: str):
        """
        Converts an MCCE PDB file to PDB file.
        Args:
            out_name (str): Converted PDB filename.
        """
        # lines from pdb w/step2 format + remark:
        with open(self.mccepdb_path) as mcfp:
            mcce_lines = mcfp.readlines()

        n_remarks = 0
        prev_chn, line_chn = None, None
        chains = list(self.chain_dict.keys())
        n_chains = len(chains)
        multichains = n_chains > 1
        if multichains:
            c = 0
            prev_chn, line_chn = chains[c], self.chain_dict[chains[c]]

        out_fp = self.output_dir.joinpath(out_name)
        with open(out_fp, "w") as ofp:
            for i, mcce_line in enumerate(mcce_lines, start=1):
                if mcce_line.startswith("REMARK"):
                    n_remarks = i
                    # nothing to convert:
                    ofp.write(mcce_line)
                    continue

                new_line = self.get_pdb_line(mcce_line)
                ofp.write(new_line)
                if multichains:
                    if i == line_chn + n_remarks:
                        ofp.write("TER\n")
                        c += 1 
                        try:
                            prev_chn, line_chn = chains[c], self.chain_dict[chains[c]]
                        except IndexError:
                            break
            if not multichains:
                ofp.write("TER\n")
        return
