---
# 1EGF :: Solution Structure Of Murine Epidermal Growth Factor Determined By Nmr Spectroscopy And Refined By Energy Minimization With Restraints;
## PDB.Structure
### Function: GROWTH FACTOR
### First Release: 01-OCT-91
### Method: SOLUTION NMR
### Resolution: NOT APPLICABLE.
### Molecule: EPIDERMAL GROWTH FACTOR
### Seqres Species: Residues: A:53; Total: 53
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ARG: 4, ASN: 3, ASP: 4, CYS: 6, GLN: 1, GLU: 2, GLY: 6, HIS: 1, ILE: 2, LEU: 4, MET: 1, PRO: 2, SER: 6, THR: 2, TRP: 2, TYR: 5, VAL: 2', 'Total: 53', 'Ionizable: 22',
              'Ratio: 41.5%')

### Disulfides:
  - CYS A  6 -- 1.92 Å --> CYS A  20
  - CYS A  14 -- 1.53 Å --> CYS A  31
  - CYS A  33 -- 1.86 Å --> CYS A  42

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "ASN A   1"
 - <strong>CTR</strong>: "ARG A  53"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   1" to " CB  ASN A   1"
- d= 1.92: " SG  CYS A   6" to " SG  CYS A  20"
- d= 1.53: " SG  CYD A  14" to " SG  CYS A  31"
- d= 1.85: " SG  CYS A  33" to " SG  CYS A  42"

</details>

