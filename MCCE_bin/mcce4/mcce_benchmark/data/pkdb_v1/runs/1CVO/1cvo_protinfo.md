---
# 1CVO :: The Solution Structure Of Cardiotoxin V From Naja Naja Atra
## PDB.Structure
### Function: CYTOTOXIN
### First Release: 28-APR-93
### Method: SOLUTION NMR
### Resolution: NOT APPLICABLE.
### Molecule: CARDIOTOXIN V
### Seqres Species: Residues: A:62; Total: 62
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 3, ARG: 1, ASN: 5, ASP: 2, CYS: 8, GLN: 1, GLU: 1, GLY: 2, HIS: 1, ILE: 1, LEU: 7, LYS: 11, PHE: 4, PRO: 5, SER: 2, THR: 4, TYR: 2, VAL: 2', 'Total: 62', 'Ionizable: 26',
              'Ratio: 41.9%')

### Disulfides:
  - CYS A  3 -- 2.06 Å --> CYS A  22
  - CYS A  15 -- 2.02 Å --> CYS A  40
  - CYS A  44 -- 2.03 Å --> CYS A  55
  - CYS A  56 -- 2.07 Å --> CYS A  61

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "LEU A   1"
 - <strong>CTR</strong>: "ASN A  62"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05

### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A  62":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.52: " CA  NTR A   1" to " CB  LEU A   1"

</details>

