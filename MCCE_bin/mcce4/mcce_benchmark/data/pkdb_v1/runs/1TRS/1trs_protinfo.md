---
# 1TRS :: The High-Resolution Three-Dimensional Solution Structures Of The Oxidized And Reduced States Of Human Thioredoxin;
## PDB.Structure
### Function: ELECTRON TRANSPORT
### First Release: 10-MAY-94
### Method: SOLUTION NMR
### Resolution: NOT APPLICABLE.
### Molecule: THIOREDOXIN
### Seqres Species: Residues: A:105; Total: 105
### Total waters: 7
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 11, ASN: 3, ASP: 7, CYS: 2, GLN: 5, GLU: 10, GLY: 5, HIS: 1, ILE: 4, LEU: 6, LYS: 12, MET: 2, PHE: 9, PRO: 3, SER: 7, THR: 5, TRP: 1, TYR: 1, VAL: 11', 'Total: 105', 'Ionizable: 33',
              'Ratio: 31.4%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 7

### Disulfides:
  - CYS A  32 -- 2.02 Å --> CYS A  35

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "MET A   1"
 - <strong>CTR</strong>: "VAL A 105"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 7.

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   1" to " CB  MET A   1"

</details>

