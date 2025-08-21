---
# 1KXI :: Structure Of Cytotoxin Homolog Precursor
## PDB.Structure
### Function: CYTOTOXIN
### First Release: 29-AUG-96
### Method: X-RAY DIFFRACTION
### Resolution: 2.19 ANGSTROMS.
### Molecule: CARDIOTOXIN V
### Seqres Species: Residues: A:62, B:62; Total: 124
### Total waters: 73
### Biounits: DIMERIC: A, B;
### Models: 1
### Chains: A, B
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 3, ARG: 1, ASN: 5, ASP: 2, CYS: 8, GLN: 1, GLU: 1, GLY: 2, HIS: 1, ILE: 1, LEU: 7, LYS: 11, PHE: 4, PRO: 5, SER: 2, THR: 4, TYR: 2, VAL: 2', 'Total: 62', 'Ionizable: 26',
              'Ratio: 41.9%')
  - B:
 'RESIDUES': ('ALA: 3, ARG: 1, ASN: 5, ASP: 2, CYS: 8, GLN: 1, GLU: 1, GLY: 2, HIS: 1, ILE: 1, LEU: 7, LYS: 11, PHE: 4, PRO: 5, SER: 2, THR: 4, TYR: 2, VAL: 2', 'Total: 62', 'Ionizable: 26',
              'Ratio: 41.9%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 43
  - B:
 'HOH': 30

### Disulfides:
  - CYS A  3 -- 2.02 Å --> CYS A  22
  - CYS A  15 -- 2.02 Å --> CYS A  40
  - CYS A  44 -- 2.02 Å --> CYS A  55
  - CYS A  56 -- 2.03 Å --> CYS A  61
  - CYS B  3 -- 2.02 Å --> CYS B  22
  - CYS B  15 -- 2.02 Å --> CYS B  40
  - CYS B  44 -- 2.02 Å --> CYS B  55
  - CYS B  56 -- 2.04 Å --> CYS B  61

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "LEU A   1", "LEU B   1"
 - <strong>CTR</strong>: "ASN A  62", "ASN B  62"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 73.

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   1" to " CB  LEU A   1"
- d= 1.54: " CA  NTR B   1" to " CB  LEU B   1"

</details>

