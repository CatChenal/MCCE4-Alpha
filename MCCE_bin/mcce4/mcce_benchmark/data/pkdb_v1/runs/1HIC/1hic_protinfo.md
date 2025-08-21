---
# 1HIC :: The Nmr Solution Structure Of Hirudin(1-51) And Comparison With Corresponding Three-Dimensional Structures Determined Using The; Complete 65-Residue Hirudin Polypeptide Chain;
## PDB.Structure
### Function: HIRUDIN
### First Release: 30-APR-92
### Method: SOLUTION NMR
### Resolution: NOT APPLICABLE.
### Molecule: HIRUDIN VARIANT
### Seqres Species: Residues: A:51; Total: 51
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ASN: 4, ASP: 2, CYS: 6, GLN: 4, GLU: 4, GLY: 8, HIS: 1, ILE: 1, LEU: 3, LYS: 3, PRO: 2, SER: 4, THR: 4, TYR: 1, VAL: 4', 'Total: 51', 'Ionizable: 17',
              'Ratio: 33.3%')

### Disulfides:
  - CYS A  6 -- 2.07 Å --> CYS A  14
  - CYS A  16 -- 2.16 Å --> CYS A  28
  - CYS A  22 -- 2.09 Å --> CYS A  39

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "VAL A   1"
 - <strong>CTR</strong>: "HIS A  51"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05

### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A  51":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.54: " CA  NTR A   1" to " CB  VAL A   1"
- d= 1.82: " O   GLU A  17" to "HD21 ASN A  20"

</details>

