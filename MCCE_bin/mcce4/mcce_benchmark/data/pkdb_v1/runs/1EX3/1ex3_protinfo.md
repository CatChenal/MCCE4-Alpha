---
# 1EX3 :: Crystal Structure Of Bovine Chymotrypsinogen A (Tetragonal)
## PDB.Structure
### Function: HYDROLASE
### First Release: 28-APR-00
### Method: X-RAY DIFFRACTION
### Resolution: 3.00 ANGSTROMS.
### Molecule: CHYMOTRYPSINOGEN A
### Seqres Species: Residues: A:245; Total: 245
### Total waters: 42
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 22, ARG: 4, ASN: 14, ASP: 9, CYS: 10, GLN: 10, GLU: 5, GLY: 23, HIS: 2, ILE: 10, LEU: 19, LYS: 14, MET: 2, PHE: 6, PRO: 9, SER: 28, THR: 23, TRP: 8, TYR: 4, VAL: 23', 'Total: 245', 'Ionizable: 48',
              'Ratio: 19.6%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 42

### Disulfides:
  - CYS A  1 -- 1.95 Å --> CYS A 122
  - CYS A  42 -- 2.09 Å --> CYS A  58
  - CYS A 136 -- 2.10 Å --> CYS A 201
  - CYS A 168 -- 1.95 Å --> CYS A 182
  - CYS A 191 -- 1.99 Å --> CYS A 220

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "CYS A   1"
 - <strong>CTR</strong>: "ASN A 245"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 42.

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.55: " CA  NTR A   1" to " CB  CYS A   1"
- d= 1.95: " SG  CYS A   1" to " SG  CYS A 122"
- d= 1.95: " SG  CYS A 168" to " SG  CYS A 182"
- d= 1.99: " SG  CYS A 191" to " SG  CYS A 220"

</details>

