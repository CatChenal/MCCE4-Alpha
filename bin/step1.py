#!/usr/bin/env python

"""
Module: step1.py

This module is a command line interface to set up and run MCCE step 1.

Input:
 * PDB file

Output:
 * run.prm
 * run.prm.record
 * step1_out.pdb
 * list_rot.gold
 * head1.lst

Usage examples:
0. Prep
    cd rundir
   Optional: Soft-link the input protein, else 'prot.pdb' in examples
             is a placeholder for the pdb
    ln -s my_protein.pdb prot.pdb

1. Run step 1
    step1.py prot.pdb

2. Write run.prm for step 1, do not actually run step 1. (No run)
    step1.py prot.pdb --norun

3. Run step 1 with specified mcce
    step1.py prot.pdb -e /path/to/mcce

4. Run step 1 with other customized parameters
    step1.py prot.pdb -u HOME_MCCE=/path/to/mcce_home,H2O_SASCUTOFF=0.05
"""

import argparse
import math
import os
from pathlib import Path
import subprocess
import sys
import re
from collections import defaultdict, Counter
from mccesteps import export_runprm
from mccesteps import record_runprm
from mccesteps import detect_runprm
from mccesteps import restore_runprm
from conserve_het import label_het

# ligands and renaming rules.
Possible_ligands = {"HIS": "HIL", "MET": "MEL", "CYS": "CYL"}
Possible_receptors = ["HEM", "HEC", "HEB", "CLA"]

# distance between heavy atoms smaller than this is considered bonded
BOND_threshold = 2.7


class Atom:
    def __init__(self):
        self.icount = 0
        self.name = ""
        self.element = ""  # two-letter code for element name
        self.resname = ""
        self.chainid = ""
        self.seqnum = 0
        self.icode = ""
        self.xyz = ()
        return

    def loadline(self, line):
        self.icount = int(line[6:11])
        self.name = line[12:16]
        if self.name[:2] == " D":
            self.name = " H" + self.name[2:]
        self.resname = line[17:20]
        self.chainid = line[21]
        self.seqnum = int(line[22:26])
        self.icode = line[26]
        self.xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        self.resid = (self.resname, self.chainid, self.seqnum, self.icode)

        if len(self.name.strip()) < 4:
            self.element = self.name[:2]
        else:
            self.element = " H"

        return

    def printme(self):
        return "ATOM %6d %4s %3s %c%4d%c   %8.3f%8.3f%8.3f\n" % (
            self.icount,
            self.name,
            self.resname,
            self.chainid,
            self.seqnum,
            self.icode,
            self.xyz[0],
            self.xyz[1],
            self.xyz[2],
        )


class Residue:
    def __init__(self):
        self.resid = ()
        self.atoms = []


def ddvv(v1, v2):
    dx = v1[0] - v2[0]
    dy = v1[1] - v2[1]
    dz = v1[2] - v2[2]
    return dx * dx + dy * dy + dz * dz


def shortest_d(res1, res2):
    ddmin = 1.0e10
    for atom1 in res1.atoms:
        if atom1.element == " H":
            continue
        for atom2 in res2.atoms:
            if atom2.element == " H":
                continue
            dd = ddvv(atom1.xyz, atom2.xyz)
            if ddmin > dd:
                ddmin = dd
    return math.sqrt(ddmin)


def fix_format(fname):
    # make sure NT atoms appear before CTR atoms.
    # truncate ending special characters from Windows.
    NTR_atoms = [" CA ", " N  ", " HA ", " H  ", " H2 ", " H3 "]
    CTR_atoms = [" C  ", " O  ", " HO ", " OXT", " HXT"]

    pdblines = open(fname).readlines()

    new_pdblines = []
    cur_resid = ""
    ntr_atoms = []
    mid_atoms = []
    ctr_atoms = []
    for line in pdblines:
        line = line.rstrip() + "\n"
        if line[:6] == "ATOM  " or line[:6] == "HETATM":
            # Remove alternative location
            if line[16] == "A":
                line = line[:16] + " " + line[17:]
            elif line[16] != " ":
                continue

            resid = line[17:27]
            if resid != cur_resid:
                new_pdblines += ntr_atoms
                new_pdblines += mid_atoms
                new_pdblines += ctr_atoms
                ntr_atoms = []
                mid_atoms = []
                ctr_atoms = []
                cur_resid = resid
            if line[12:16] in NTR_atoms:
                ntr_atoms.append(line)
            elif line[12:16] in CTR_atoms:
                ctr_atoms.append(line)
            else:
                mid_atoms.append(line)
        else:
            new_pdblines.append(line)

    # last one
    new_pdblines += ntr_atoms
    new_pdblines += mid_atoms
    new_pdblines += ctr_atoms
    
    return new_pdblines


def group_residues(lines):
    residues = []
    atoms = []
    for line in lines:
        atom = Atom()
        atom.loadline(line)
        atoms.append(atom)

    resids = []
    for atom in atoms:
        if atom.resid in resids:
            ires = resids.index(atom.resid)
            residues[ires].atoms.append(atom)
        else:
            residue = Residue()
            residue.resid = atom.resid
            residue.atoms.append(atom)
            residues.append(residue)
            resids.append(atom.resid)

    return residues


def rename_ligand(res):
    #print("   Renaming amino acid %s to ligand %s" % (str(res.resid), Possible_ligands[res.resid[0]]))
    print("   Distance below bond threshold, renaming amino acid ", end='')
    for item in res.resid:
        print(str(item) + " ", end='')
    print("to ligand ", end='')
    print(Possible_ligands[res.resid[0]])
    for atom in res.atoms:
        if atom.resname in Possible_ligands:
            atom.resname = Possible_ligands[atom.resname]


def identify_ligands(lines):
    newlines = []

    lines = [x.strip() for x in lines if x[:6] == "ATOM  " or x[:6] == "HETATM"]
    residues = group_residues(lines)

    receptors = []
    for res in residues:
        resname = res.resid[0]
        if resname in Possible_receptors:
            receptors.append(res)

    ligands = []
    for res in residues:
        resname = res.resid[0]
        if resname in Possible_ligands:
            ligands.append(res)

    for receptor in receptors:
        for ligand in ligands:
            if shortest_d(receptor, ligand) < BOND_threshold:
                rename_ligand(ligand)

    for res in residues:
        for atom in res.atoms:
            newlines.append(atom.printme())

    return newlines


def write_runprm(args):
    runprm = {}

    path = str(os.path.dirname(os.path.abspath(__file__)))
    base_path = os.path.dirname(path)
    # print(base_path)
    runprm["INPDB"] = "step0_out.pdb"
    runprm["DO_PREMCCE"] = "t"
    runprm["MCCE_HOME"] = base_path
    runprm["MINIMIZE_SIZE"] = "t"
    if args.noter:
        runprm["TERMINALS"] = "f"
    else:
        runprm["TERMINALS"] = "t"
    runprm["CLASH_DISTANCE"] = "2.0"
    runprm["H2O_SASCUTOFF"] = "0.05"
    runprm["IGNORE_INPUT_H"] = "t"
    runprm["RENAME_RULES"] = "%s/name.txt" % base_path
    runprm["EPSILON_PROT"] = args.d

    if args.dry:
        runprm["H2O_SASCUTOFF"] = "-0.01"

    if args.u:
        fields = args.u.split(",")
        for field in fields:
            try:
                key, value = field.split("=")
                runprm[key] = value
            except ValueError:
                print(f"Each component format is 'KEY=VALUE'. Unrecognized: {field}.")

    if args.load_runprm:
        lines = open(args.load_runprm)
        for line in lines:
            entry_str = line.strip().split("#")[0]
            fields = entry_str.split()
            if len(fields) > 1:
                key_str = fields[-1]
                if key_str[0] == "(" and key_str[-1] == ")":
                    key = key_str.strip("()").strip()
                    value = fields[0]
                    runprm[key] = value

    # unconditionally force to run step 1 only in step1.py
    runprm["DO_PREMCCE"] = "t"
    runprm["DO_ROTAMERS"] = "f"
    runprm["DO_ENERGY"] = "f"
    runprm["DO_MONTE"] = "f"

    export_runprm(runprm)
    record_runprm(runprm, "#STEP1")

    return

def insert_termini_warning(text: str) -> str:
    """Insert a warning when multiple termini were labelled, which
    happens when a chain has breaks and termini labelling is turned on.
    The warning directs the user to turn off labelling by setting parameter
    'TERMINALS' to 'f'.
    """
    pattern = re.compile(r"(   Identify NTR and CTR\.\.\..*?Done)", re.DOTALL)
    match = pattern.search(text)
    try:
        block = match.group(1).splitlines()[1:-1]
    except AttributeError:
        return text
    # number of labelled ter should be 2 per 
    ch  = len(set([line[20] for line in block]))
    if len(block) > ch * 2:
        warning = ("      Multiple NTR/CTR residues indicate breaks in the chain(s). "
                   "Turn off labelling of termini by setting parameter 'TERMINALS' to 'f'.\n"
                   )
        insertion_point = match.end(1) - 7  # 7 == len("   Done")
        return text[:insertion_point] + warning + text[insertion_point:]
    
    return text


def filter_and_condense_missing_atoms(log_text):
    # This function takes all the "Missing Heavy Atom" messages
    # generated in the C code, removes the messages associated 
    # with NTR or CTR. If multiple atoms are removed from a 
    # conformer, they are reformatted to be on one line.

    # Collect cleaned entries
    cleaned_lines = []
    missing_atoms = defaultdict(list)

    # Pattern to match all "Missing heavy atom ..." lines
    pattern = re.compile(
        r'(Missing heavy atom\s+(\S+)\s+of conf\s+(\S+)\s+in "([^"]+)")'
    )

    # Store lines to remove
    lines_to_remove = set()

    for match in pattern.finditer(log_text):
        full_line, atom, conf, residue = match.groups()

        # Flag line for removal
        lines_to_remove.add(full_line)

        # Skip unwanted conformers
        if 'NTR' in conf or 'CTR' in conf or conf.endswith('BK'):
            continue

        key = (conf, residue)
        missing_atoms[key].append(atom)

    # Generate new condensed lines
    for (conf, residue), atoms in missing_atoms.items():
        atoms_str = ", ".join(atoms)
        cleaned_lines.append(
                f"   Missing heavy atoms for {conf} in \"{residue}\":  {atoms_str}."
        )

    # Replace old section in original text
    new_lines = []
    for line in log_text.splitlines():
        if not any(line.strip().startswith(to_remove) for to_remove in lines_to_remove):
            new_lines.append(line)
    new_log_text = "\n".join(new_lines)

    # Insert condensed lines before final "Missing heavy atoms detected" message
    insert_point = new_log_text.find("Missing heavy atoms detected.")
    if insert_point != -1:
        before = new_log_text[:insert_point].rstrip()
        after = new_log_text[insert_point:]
        final_text = before + "\n" + "\n".join(cleaned_lines) +"\n\n   " + after
    else:
        final_text = new_log_text + "\n" + "\n".join(cleaned_lines)

    return final_text

if __name__ == "__main__":

    # Get the command arguments
    helpmsg = "Run mcce step 1, premcce to format PDB file to MCCE PDB format."
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument(
        "--norun",
        default=False,
        action="store_true",
        help="Create run.prm but do not run step 1; default: %(default)s.",
    )
    parser.add_argument(
        "--noter",
        default=False,
        action="store_true",
        help="Do not label terminal residues (for making ftpl); default: %(default)s.",
    )
    parser.add_argument(
        "-e",
        metavar="/path/to/mcce",
        default="mcce",
        help="mcce executable location; default: %(default)s.",
    )
    parser.add_argument(
        "-u",
        metavar="Key=Value",
        default="",
        help="User customized variables; default: %(default)s.",
    )
    parser.add_argument(
        "-d",
        metavar="epsilon",
        default="4.0",
        help="protein dielectric constant for PBE solvers; default: %(default)s.",
    )
    parser.add_argument(
        "-load_runprm",
        metavar="prm_file",
        default="",
        help="Load additional run.prm file, overwrite default values."
    )
    parser.add_argument(
        "--dry",
        default=False,
        action="store_true",
        help="Delete all water molecules; default: %(default)s.",
    )

    parser.add_argument("prot", metavar="pdb", nargs=1)
    args = parser.parse_args()

    detected = detect_runprm()
    write_runprm(args)

    if not args.norun:
        new_pdblines = fix_format(args.prot[0])
        print("Preprocessing input pdb file, identifying ligands ...")
        new_pdblines = identify_ligands(new_pdblines)
        print("Done, writing to step0_out.pdb ...")
        with open("step0_out.pdb", "w") as fo:
            fo.writelines(new_pdblines)

        process = subprocess.Popen([args.e], close_fds=True, stdout=subprocess.PIPE)
        process_str = ""
        for line in process.stdout:
            # print(line.decode(), end="") # unformatted "Missing heavy atom" code
            process_str += line.decode()

        # inserts termini warning if mult NTR in a chain
        process_str = insert_termini_warning(process_str)
        # re-formats Missing heavy atom output
        processed_str = filter_and_condense_missing_atoms(process_str)
        print(processed_str)

        if not Path("step1_out.pdb").exists():
            sys.exit("[Step1.py Error]: Output not found: step1_out.pdb.")

        # convert to HETATM
        lines = open("step1_out.pdb").readlines()
        newlines = label_het(lines)
        open("step1_out.pdb", "w").writelines(newlines)

    if detected:
        restore_runprm()
