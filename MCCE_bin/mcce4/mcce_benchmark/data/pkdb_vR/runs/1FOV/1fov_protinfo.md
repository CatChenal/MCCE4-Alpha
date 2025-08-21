---
# 1FOV :: Glutaredoxin 3 From Escherichia Coli In The Fully Oxidized Form
## PDB.Structure
### Function: ELECTRON TRANSPORT
### First Release: 29-AUG-00
### Method: SOLUTION NMR
### Resolution: NOT APPLICABLE.
### Molecule: GLUTAREDOXIN 3
### Seqres Species: Residues: A:82; Total: 82
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 8, ARG: 5, ASN: 2, ASP: 6, CYS: 2, GLN: 3, GLU: 5, GLY: 7, HIS: 2, ILE: 6, LEU: 8, LYS: 6, MET: 1, PHE: 2, PRO: 4, SER: 4, THR: 4, TYR: 4, VAL: 3', 'Total: 82', 'Ionizable: 30',
              'Ratio: 36.6%')

### Disulfides:
  - CYS A  11 -- 2.02 Å --> CYS A  14

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "ALA A   1"
 - <strong>CTR</strong>: "LYS A  82"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   1" to " CB  ALA A   1"
- d= 1.80: "HD21 ASN A   2" to " OD1 ASP A  58"
- d= 1.90: " O   PRO A  32" to "HD21 ASN A  36"
- d= 1.87: " O   THR A  51" to "HE21 GLN A  54"

</details>

