---
# 135L :: X-Ray Structure Of Monoclinic Turkey Egg Lysozyme At 1.3 Angstroms Resolution;
## PDB.Structure
### Function: HYDROLASE(O-GLYCOSYL)
### First Release: 10-JUN-93
### Method: X-RAY DIFFRACTION
### Resolution: 1.30 ANGSTROMS.
### Molecule: TURKEY EGG WHITE LYSOZYME
### Seqres Species: Residues: A:129; Total: 129
### Total waters: 114
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 13, ARG: 10, ASN: 14, ASP: 6, CYS: 8, GLN: 1, GLU: 2, GLY: 13, HIS: 2, ILE: 6, LEU: 9, LYS: 7, MET: 2, PHE: 2, PRO: 2, SER: 10, THR: 7, TRP: 6, TYR: 4, VAL: 5', 'Total: 129', 'Ionizable: 39',
              'Ratio: 30.2%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 114

### Disulfides:
  - CYS A  6 -- 2.00 Å --> CYS A 127
  - CYS A  30 -- 2.02 Å --> CYS A 115
  - CYS A  64 -- 2.03 Å --> CYS A  80
  - CYS A  76 -- 1.99 Å --> CYS A  94

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "LYS A   1"
 - <strong>CTR</strong>: "LEU A 129"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 114.

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.49: " CA  NTR A   1" to " CB  LYS A   1"
- d= 1.99: " SG  CYS A  76" to " SG  CYS A  94"

</details>

