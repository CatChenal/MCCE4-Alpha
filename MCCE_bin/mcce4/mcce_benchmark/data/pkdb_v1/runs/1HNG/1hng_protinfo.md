---
# 1HNG :: Crystal Structure At 2.8 Angstroms Resolution Of A Soluble Form Of The Cell Adhesion Molecule Cd2;
## PDB.Structure
### Function: T LYMPHOCYTE ADHESION GLYCOPROTEIN
### First Release: 10-AUG-94
### Method: X-RAY DIFFRACTION
### Resolution: 2.80 ANGSTROMS.
### Molecule: CD2
### Seqres Species: Residues: A:176, B:176; Total: 352
### Missing:
  - Residues:
 'A': [('ARG', 1), ('Count', 1)], 'B': [('ARG', 1), ('Count', 1)]

### Biounits: DIMERIC: A, B;
### Models: 1
### Chains: A, B
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 8, ARG: 10, ASN: 12, ASP: 9, CYS: 4, GLN: 5, GLU: 14, GLY: 11, HIS: 2, ILE: 8, LEU: 18, LYS: 13, MET: 6, PHE: 5, PRO: 5, SER: 11, THR: 13, TRP: 4, TYR: 5, VAL: 12', 'Total: 175', 'Ionizable: 57',
              'Ratio: 32.6%')
  - B:
 'RESIDUES': ('ALA: 8, ARG: 10, ASN: 12, ASP: 9, CYS: 4, GLN: 5, GLU: 14, GLY: 11, HIS: 2, ILE: 8, LEU: 18, LYS: 13, MET: 6, PHE: 5, PRO: 5, SER: 11, THR: 13, TRP: 4, TYR: 5, VAL: 12', 'Total: 175', 'Ionizable: 57',
              'Ratio: 32.6%')

### Disulfides:
  - CYS A 110 -- 2.01 Å --> CYS A 174
  - CYS A 117 -- 2.00 Å --> CYS A 157
  - CYS B 110 -- 1.99 Å --> CYS B 174
  - CYS B 117 -- 2.00 Å --> CYS B 157

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "ASP A   2", "ASP B   2"
 - <strong>CTR</strong>: "GLU A 176", "GLU B 176"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.56: " CA  NTR A   2" to " CB  ASP A   2"
- d= 2.00: " SG  CYS A 117" to " SG  CYS A 157"
- d= 1.55: " CA  NTR B   2" to " CB  ASP B   2"
- d= 1.99: " SG  CYS B 110" to " SG  CYS B 174"
- d= 2.00: " SG  CYS B 117" to " SG  CYS B 157"

</details>

