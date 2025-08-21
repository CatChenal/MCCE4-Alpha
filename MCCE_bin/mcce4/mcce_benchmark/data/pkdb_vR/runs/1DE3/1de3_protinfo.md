---
# 1DE3 :: Solution Structure Of The Cytotoxic Ribonuclease Alpha-Sarcin
## PDB.Structure
### Function: HYDROLASE
### First Release: 12-NOV-99
### Method: SOLUTION NMR
### Resolution: NOT APPLICABLE.
### Molecule: RIBONUCLEASE ALPHA-SARCIN
### Seqres Species: Residues: A:150; Total: 150
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 5, ARG: 4, ASN: 11, ASP: 11, CYS: 4, GLN: 3, GLU: 6, GLY: 13, HIS: 8, ILE: 4, LEU: 9, LYS: 20, PHE: 6, PRO: 13, SER: 9, THR: 11, TRP: 2, TYR: 8, VAL: 3', 'Total: 150', 'Ionizable: 61',
              'Ratio: 40.7%')

### Disulfides:
  - CYS A  6 -- 2.04 Å --> CYS A 148
  - CYS A  76 -- 2.04 Å --> CYS A 132

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "ALA A   1"
 - <strong>CTR</strong>: "HIS A 150"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.54: " CA  NTR A   1" to " CB  ALA A   1"
- d= 1.37: "HG23 THR A   3" to "HD23 LEU A  24"
- d= 2.00: "HG21 ILE A  69" to "HG21 ILE A 135"
- d= 1.81: "HD22 LEU A  94" to "HG22 ILE A 123"
- d= 1.96: "HD12 LEU A  95" to " O   TYR A 124"
- d= 1.89: " OE2 GLU A  96" to "HD11 ILE A 123"
- d= 1.88: "HD12 ILE A 135" to "HD22 LEU A 145"

</details>

