---
# 1LYS :: X-Ray Structure Of A Monoclinic Form Of Hen Egg-White Lysozyme Crystallized At 313K. Comparison Of Two Independent Molecules;
## PDB.Structure
### Function: HYDROLASE(O-GLYCOSYL)
### First Release: 03-DEC-93
### Method: X-RAY DIFFRACTION
### Resolution: 1.72 ANGSTROMS.
### Molecule: HEN EGG WHITE LYSOZYME
### Seqres Species: Residues: A:129, B:129; Total: 258
### Total waters: 215
### Biounits: MONOMERIC: A,B;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 12, ARG: 11, ASN: 14, ASP: 7, CYS: 8, GLN: 3, GLU: 2, GLY: 12, HIS: 1, ILE: 6, LEU: 8, LYS: 6, MET: 2, PHE: 3, PRO: 2, SER: 10, THR: 7, TRP: 6, TYR: 3, VAL: 6', 'Total: 129', 'Ionizable: 38',
              'Ratio: 29.5%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 98

### Disulfides:
  - CYS A  6 -- 1.99 Å --> CYS A 127
  - CYS A  30 -- 2.00 Å --> CYS A 115
  - CYS A  64 -- 2.01 Å --> CYS A  80
  - CYS A  76 -- 2.00 Å --> CYS A  94
  - CYS B  6 -- 1.97 Å --> CYS B 127
  - CYS B  30 -- 2.03 Å --> CYS B 115
  - CYS B  64 -- 2.01 Å --> CYS B  80
  - CYS B  76 -- 2.02 Å --> CYS B  94

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "LYS A   1"
 - <strong>CTR</strong>: "LEU A 129"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 98.

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.54: " CA  NTR A   1" to " CB  LYS A   1"
- d= 1.99: " SG  CYS A   6" to " SG  CYS A 127"
- d= 2.00: " SG  CYS A  76" to " SG  CYS A  94"

</details>

