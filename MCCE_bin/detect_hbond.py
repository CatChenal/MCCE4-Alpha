#!/usr/bin/env python
"""
This script detects hydrogen bonds between atoms from two conformers and list them in file hah.txt.
This script is self-contained, with no external dependencies and can be run independently.

The input file is step2_out.pdb, the output file is hah.txt.

Usage: detect_hbond.py [step2_out.pdb]
"""

import argparse
import numpy as np
import os
import sys

# Set up argument parsing
parser = argparse.ArgumentParser(description="Detect hydrogen bonds between atoms from two conformers.")
parser.add_argument("-inpdb", type=str, default="step2_out.pdb", help="MCCE step2_out formatted input pdbfile")
args = parser.parse_args()

# Constants
# The parameters are roughly based on https://ctlee.github.io/BioChemCoRe-2018/h-bond/
# DNEAR = 2.0     # Min Distance cutoff for hydrogen bond between heavy atoms
# DFAR  = 4.0     # Max Distance cutoff for hydrogen bond between heavy atoms
# ANGCUT = 100    # Angle cutoff for hydrogen bond between heavy atoms, 180 is ideal, 100 is the smallest angle allowed

# To match old program, we use the following parameters
DNEAR    = 2.0  # Min Distance cutoff for hydrogen bond between heavy atoms
DFAR     = 4.0  # Max Distance cutoff for hydrogen bond between heavy atoms
ANGCUT   = 90   # Angle cutoff for hydrogen bond between heavy atoms, 180 is ideal, 100 is the smallest angle allowed
ANGBLOCK = 60   # Angle cutoff for blocking atoms, 180 is ideal, 90 is the smallest angle allowed for not blocking

# Minimum charge for hydrogen bond donor/acceptor, my best guess
MIN_NCRG = -0.2  # Minimum heavy atom charge for hydrogen bond donor/acceptor
MIN_HCRG =  0.2  # Minimum H atom charge for hydrogen bond donor/acceptor


# Setup input and output filenames
if args.inpdb == "step2_out.pdb":
   in_file = "step2_out.pdb"
   if not os.path.isfile(in_file):
      print(f"Error: Default input file '{in_file}' not found in the working directory.")
      sys.exit(1)
   out_file = "hah.txt"
   blocking_file = "blocking.txt"
else:
    in_file = args.inpdb
    if not os.path.isfile(in_file):
      print(f"Error: Default input file '{in_file}' not found in the working directory.")
      sys.exit(1)
    input_file_name = os.path.splitext(os.path.basename(args.inpdb))[0]  # Extract filename without extension
    out_file = f"{input_file_name}.txt"
    blocking_file = f"{input_file_name}_blocking_file.txt"

# Optional input parameter files
class Atom:
    def __init__(self):
        self.atom_name = ""
        self.res_name = ""
        self.res_seq = 0
        self.chain_id = ""
        self.iCode = ""
        self.xyz = (0.0, 0.0, 0.0)
        self.charge = 0.0
        self.radius = 0.0
        self.confNum = 0
        self.confType = ""
        self.conn12 = []

    def loadline(self, line):
        self.atom_name = line[12:16]
        self.res_name = line[17:20]
        self.res_seq = int(line[22:26])
        self.chain_id = line[21]
        self.iCode = line[26]
        self.xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        self.charge = float(line[62:74])
        self.confNum = int(line[27:30])
        self.confType = "%3s%2s" % (self.res_name, line[80:82])
        self.confID = "%5s%c%04d%c%03d" % (self.confType, self.chain_id, self.res_seq, self.iCode, self.confNum)


    def printline(self):
        print("%s \"%4s\" %8.3f %8.3f %8.3f %8.3f" % (self.confID, self.atom_name, self.x, self.y, self.z, self.charge))       


def dist(atom1, atom2):
    return ((atom1.xyz[0] - atom2.xyz[0])**2 + (atom1.xyz[1] - atom2.xyz[1])**2 + (atom1.xyz[2] - atom2.xyz[2])**2)**0.5

def is_H(atom_name):
    atom = atom_name.strip()
    return (len(atom) < 4 and atom_name[1] == "H") or (len(atom) == 4 and atom_name[0] == "H")

def avv(v1, v2):
    # angle between two vectors
    return np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))

def detect_hbond(pdb_file):
    """
    Detect hydrogen bonds between atoms from two conformers and list them in file hah.txt.

    Solution:
    1. Read the input pdb file
    2. Group atoms into conformers
    3. Make intra 12 connectivity of atoms based on atom distance
    4. Identify potential hydrogen bond donors and acceptors based on charge, including backbone
    5. Pick two atoms from poetntial hydrogen bond donors and acceptors
    6. For each atom pair, exclude those within the conformer
    7.     Check distance between donor and acceptor
           Check D-H-A angle
           Check if connected atoms on Acceptor is blocking the D--H--A path by using D--?-A angle
           If all passed, record a potential hydrogen bond
    8. Write the hydrogen bond list to hah.txt
    """

    # Read the input pdb file
    atoms = []
    lines = open(pdb_file).readlines()
    for line in lines:
        if line.startswith("ATOM  ") or line.startswith("HETATM"):
            atom = Atom()
            atom.loadline(line)
            atoms.append(atom)

    # Group atoms into conformers. conformers is a dictionary with confID as key and a list of atoms as value.
    conformers = {}
    for atom in atoms:
        if atom.confID not in conformers:
            conformers[atom.confID] = []
        conformers[atom.confID].append(atom)

    # make intra 12 connectivity of atoms based on atom distance
    for confID in conformers:
        atoms_inconf = conformers[confID]
        for i in range(len(atoms_inconf)-1):
            for j in range(i+1, len(atoms_inconf)):
                atom1 = atoms_inconf[i]
                atom2 = atoms_inconf[j]
                if dist(atom1, atom2) < 2.0:
                    atom1.conn12.append(atom2)
                    atom2.conn12.append(atom1)

    # Identify potential hydrogen bond donors and acceptors based on charge
    donors_acceptors = []
    for confID in conformers:
        atoms_inconf = conformers[confID]
        for atom in atoms_inconf:
            # These are conditions for an atom to be a donor or acceptor (heavy atom)
            # 1. atom is not H
            # 2. atom charge is less than -0.2
            if not is_H(atom.atom_name) and atom.charge < MIN_NCRG:
                donors_acceptors.append(atom)

    # Pick two atoms as potential hydrogen bond donors and acceptors
    hbond_records = []
    blocking_records = []
    for donor in donors_acceptors:
        for acceptor in donors_acceptors:
            if donor.confID[:3]+donor.confID[5:11] != acceptor.confID[:3]+acceptor.confID[5:11]:  # exclude atoms from the same residue
                distance = dist(donor, acceptor)
                if DNEAR < distance < DFAR:
                    for h in donor.conn12:
                        if is_H(h.atom_name) and h.charge > MIN_HCRG:
                            Vhd = (donor.xyz[0] - h.xyz[0], donor.xyz[1] - h.xyz[1], donor.xyz[2] - h.xyz[2])
                            Vha = (acceptor.xyz[0] - h.xyz[0], acceptor.xyz[1] - h.xyz[1], acceptor.xyz[2] - h.xyz[2])
                            angle = avv(Vhd, Vha) * 180 / np.pi
                            if angle > ANGCUT: # this means the D--H--A angle is good
                                # Check if connected atoms on Acceptor is blocking the D--H--A path by using H--A-? angle
                                blocking = False
                                for q in acceptor.conn12:  # since donor and acceptor are in different conformers, we don't need to check if q is donor
                                    Vqa = (acceptor.xyz[0] - q.xyz[0], acceptor.xyz[1] - q.xyz[1], acceptor.xyz[2] - q.xyz[2])
                                    blocking_angle = avv(Vqa, Vha) * 180 / np.pi 
                                    if blocking_angle < ANGBLOCK:  # this means q is blocking the D--H--A path
                                        blocking = True
                                        blocking_records.append((h, acceptor, q, blocking_angle))
                                        break
                                if not blocking:
                                    distance_h2a = dist(h, acceptor)
                                    hbond_records.append((donor, h, acceptor, distance_h2a, angle))

    # Write the hydrogen bond list to hah.txt
    lines = []
    for record in hbond_records:
        donor, h, acceptor, distance, angle = record
        line = "%s  %s  %s~%s--%s%5.2f %6d\n" % (donor.confID, acceptor.confID, donor.atom_name, h.atom_name, acceptor.atom_name, distance, angle)
        lines.append(line)

    # Write the blocking list to blocking.txt
    lines_blocking = []
    for record in blocking_records:
        h, acceptor, q, blocking_angle = record
        line = "%s  %s  %s--%s~%s %3d\n" % (h.confID, acceptor.confID, h.atom_name, acceptor.atom_name, q.atom_name, blocking_angle)
        lines_blocking.append(line)

    open(out_file, "w").writelines(lines)
    open(blocking_file, "w").writelines(lines_blocking)

    
    return

if __name__ == "__main__":
    print("Detecting hydrogen bonds...")
    print("Criteria:")
    print("  Donor-Acceptor (heavy atoms) Distance: %5.2f - %5.2f" % (DNEAR, DFAR))
    print("  D-H-A Angle (Angle >= this to qualify H bond): %5.2f" % (ANGCUT))
    print("  H--A-? Blocking Angle (Angle <= this to block H bond): %5.2f" % (ANGBLOCK))
    print("  Heavy Atom Charge for H bond Donor/Acceptor: < %5.2f" % (MIN_NCRG))
    print("  H Atom Charge for H bond Donor: > %5.2f" % (MIN_HCRG))
    print()
    print("Input file:    %s" % (in_file))
    print("Output file:   %s" % (out_file))
    print("Blocking file: %s" % (blocking_file))
    print()
    print("Questions about step 6 hah.txt outpot:")
    print("  Step 6 excludes hydrogen bonds with the backbone atoms")
    print("  Step 6 doesn't check if the potential hydrogen bonds are blocked by 3rd atom")
    print()
    detect_hbond(args.inpdb)



