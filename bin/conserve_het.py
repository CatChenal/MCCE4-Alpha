#!/usr/bin/env python
"""
Conserve HETATM in mcce pdb file
"""
import sys

Amino_acids = {"ALA",
               "ARG",
               "ASN",
               "ASP",
               "CTR",
               "CYD",
               "CYL",
               "CYS",
               "GLN",
               "GLY",
               "GLU",
               "HIL",
               "HIS",
               "ILE",
               "LEU",
               "LYS",
               "MET",
               "MEL",
               "NTG",
               "NTR",
               "PHE",
               "PRO",
               "SER",
               "THR",
               "TRP",
               "TYR",
               "VAL"}

def label_het(pdblines):
    newlines = []
    for line in pdblines:
        if line[:6] == "ATOM  " and line[17:20] not in Amino_acids:
            newline = "HETATM" + line[6:]
            newlines.append(newline)
        else:
            newlines.append(line)
    return newlines


if __name__ == "__main__":
    fname = sys.argv[1]
    lines = open(fname).readlines()
    newlines = label_het(lines)
    sys.stdout.writelines(newlines)

    