---
# 2OVO :: The Crystal And Molecular Structure Of The Third Domain Of Silver Pheasant Ovomucoid (Omsvp3);
## PDB.Structure
### Function: PROTEINASE INHIBITOR (KAZAL)
### First Release: 11-JUN-85
### Method: X-RAY DIFFRACTION
### Resolution: 1.50 ANGSTROMS.
### Molecule: OVOMUCOID THIRD DOMAIN
### Seqres Species: Residues: A:56; Total: 56
### Total waters: 31
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 4, ARG: 1, ASN: 5, ASP: 2, CYS: 6, GLU: 3, GLY: 4, HIS: 1, LEU: 4, LYS: 4, MET: 1, PHE: 2, PRO: 3, SER: 5, THR: 4, TYR: 3, VAL: 4', 'Total: 56', 'Ionizable: 20',
              'Ratio: 35.7%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 31

### Disulfides:
  - CYS A  8 -- 2.00 Å --> CYS A  38
  - CYS A  16 -- 2.03 Å --> CYS A  35
  - CYS A  24 -- 1.99 Å --> CYS A  56

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "LEU A   1"
 - <strong>CTR</strong>: "CYS A  56"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 31.

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   1" to " CB  LEU A   1"
- d= 2.00: " SG  CYS A   8" to " SG  CYS A  38"
- d= 1.99: " SG  CYS A  24" to " SG  CYS A  56"

</details>

