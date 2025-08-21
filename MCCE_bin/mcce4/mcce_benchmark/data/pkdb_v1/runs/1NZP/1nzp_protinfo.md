---
# 1NZP :: Solution Structure Of The Lyase Domain Of Human Dna Polymerase Lambda
## PDB.Structure
### Function: DNA BINDING PROTEIN/TRANSFERASE
### First Release: 19-FEB-03
### Method: SOLUTION NMR
### Resolution: NOT APPLICABLE.
### Molecule: DNA POLYMERASE LAMBDA
### Seqres Species: Residues: A:87; Total: 87
### Missing:
  - Residues:
 'A': [('MET', 241), ('Count', 1)]

### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 10, ARG: 3, ASN: 3, ASP: 2, CYS: 1, GLN: 4, GLU: 6, GLY: 5, HIS: 5, ILE: 7, LEU: 8, LYS: 10, MET: 1, PHE: 1, PRO: 3, SER: 7, THR: 3, TRP: 1, TYR: 3, VAL: 3', 'Total: 86', 'Ionizable: 30',
              'Ratio: 34.9%')

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "ALA A 242"
 - <strong>CTR</strong>: "HIS A 327"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.52: " CA  NTR A 242" to " CB  ALA A 242"
- d= 1.85: "HG21 ILE A 256" to "HD13 LEU A 260"
- d= 1.62: "HG22 ILE A 256" to "HG22 ILE A 313"
- d= 1.86: " O   LYS A 259" to "HG23 VAL A 262"

</details>

