---
# 1ANS :: Three-Dimensional Structure In Solution Of Neurotoxin Iii From The Sea Anemone Anemonia Sulcata;
## PDB.Structure
### Function: TOXIN
### First Release: 09-JUN-94
### Method: SOLUTION NMR
### Resolution: NOT APPLICABLE.
### Molecule: NEUROTOXIN III
### Seqres Species: Residues: A:27; Total: 27
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ARG: 1, ASN: 1, CYS: 6, GLN: 1, GLU: 1, GLY: 5, LYS: 1, PRO: 4, SER: 2, TRP: 2, TYR: 2, VAL: 1', 'Total: 27', 'Ionizable: 11',
              'Ratio: 40.7%')

### Disulfides:
  - CYS A  3 -- 2.03 Å --> CYS A  17
  - CYS A  4 -- 2.02 Å --> CYS A  11
  - CYS A  6 -- 2.02 Å --> CYS A  22

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "ARG A   1"
 - <strong>CTR</strong>: "VAL A  27"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   1" to " CB  ARG A   1"

</details>

