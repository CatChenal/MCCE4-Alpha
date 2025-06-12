#!/usr/bin/env python
"""
Created on Apr 30 05:00:00 2025

@author: Gehan Ranepura
"""

import argparse
from openeye import oechem

# Set up argument parsing
parser = argparse.ArgumentParser(description="Process an OEB file and print atom details.")
parser.add_argument("oeb_file", nargs="?", default="oeassigncharges.oeb.gz", help="Path to the OEB file (default: %(default)s)")
args = parser.parse_args()

# Open the OEB file
ifs = oechem.oemolistream(args.oeb_file)

# Iterate through the molecules in the file
for mol in ifs.GetOEGraphMols():
    # Retrieve the molecule's title (typically the residue or molecule name)
    molecule_name = mol.GetTitle() if mol.GetTitle() else "Unnamed Molecule"
    print(f"Molecule: {molecule_name}")

    # Iterate through the atoms in the molecule
    for atom in mol.GetAtoms():
        atom_name = atom.GetName() if atom.GetName() else "Unnamed"
        symbol = oechem.OEGetAtomicSymbol(atom.GetAtomicNum())
        charge = atom.GetPartialCharge()
        print(f"  Atom Name: {atom_name:>4} | Symbol: {symbol:>2} | Charge: {charge:>6.3f}")
    print("-" * 60)

