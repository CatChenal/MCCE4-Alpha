---
# 1TQO :: Cryogenic Crystal Structure Of Staphylococcal Nuclease Variant Truncated Delta+Phs I92E;
## PDB.Structure
### Function: HYDROLASE
### First Release: 17-JUN-04
### Method: X-RAY DIFFRACTION
### Resolution: 2.00 ANGSTROMS.
### Molecule: THERMONUCLEASE
### Seqres Species: Residues: A:138; Total: 138
### Total waters: 85
### Missing:
  - Residues:
 'A': [('ALA', 1), ('THR', 2), ('SER', 3), ('THR', 4), ('LYS', 5), ('LYS', 6), ('ASN', 144),
       ('Count', 7)]

### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 13, ARG: 5, ASN: 6, ASP: 7, GLN: 5, GLU: 13, GLY: 9, HIS: 2, ILE: 4, LEU: 12, LYS: 18, MET: 4, PHE: 4, PRO: 4, SER: 2, THR: 7, TRP: 1, TYR: 7, VAL: 8', 'Total: 131', 'Ionizable: 52',
              'Ratio: 39.7%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 85

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "LEU A   7"
 - <strong>CTR</strong>: "ASP A 143"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 85.

### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 143":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   7" to " CB  LEU A   7"

</details>

