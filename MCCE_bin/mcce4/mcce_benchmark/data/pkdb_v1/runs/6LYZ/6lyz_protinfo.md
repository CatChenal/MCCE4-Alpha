---
# 6LYZ :: Real-Space Refinement Of The Structure Of Hen Egg-White Lysozyme
## PDB.Structure
### Function: HYDROLASE
### First Release: 01-FEB-75
### Method: X-RAY DIFFRACTION
### Resolution: 2.00 ANGSTROMS.
### Molecule: HEN EGG WHITE LYSOZYME
### Seqres Species: Residues: A:129; Total: 129
### Total waters: 101
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 12, ARG: 11, ASN: 14, ASP: 7, CYS: 8, GLN: 3, GLU: 2, GLY: 12, HIS: 1, ILE: 6, LEU: 8, LYS: 6, MET: 2, PHE: 3, PRO: 2, SER: 10, THR: 7, TRP: 6, TYR: 3, VAL: 6', 'Total: 129', 'Ionizable: 38',
              'Ratio: 29.5%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 101

### Disulfides:
  - CYS A  6 -- 2.03 Å --> CYS A 127
  - CYS A  30 -- 2.04 Å --> CYS A 115
  - CYS A  64 -- 2.01 Å --> CYS A  80
  - CYS A  76 -- 2.11 Å --> CYS A  94

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "LYS A   1"
 - <strong>CTR</strong>: "LEU A 129"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 101.

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.52: " CA  NTR A   1" to " CB  LYS A   1"
- d= 1.70: " NH1 ARG A   5" to " O   TRP A 123"
- d= 1.77: " NH2 ARG A   5" to " O   ARG A 125"
- d= 1.74: " NH1 ARG A  45" to " NH1 ARG A  68"
- d= 1.47: " O   LYS A  97" to " CG  ASP A 101"
- d= 1.01: " O   LYS A  97" to " OD1 ASP A 101"
- d= 1.70: " O   LYS A  97" to " OD2 ASP A 101"

</details>

