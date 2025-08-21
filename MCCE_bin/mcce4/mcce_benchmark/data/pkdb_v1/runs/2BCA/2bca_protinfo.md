---
# 2BCA :: High-Resolution Solution Structure Of Calcium-Loaded Calbindin D9K
## PDB.Structure
### Function: CALCIUM-BINDING PROTEIN
### First Release: 18-AUG-93
### Method: SOLUTION NMR
### Resolution: NOT APPLICABLE.
### Molecule: CALBINDIN D9K
### Seqres Species: Residues: A:76; Total: 76
### Missing:
  - Residues:
 'A': [('MET', 0), ('Count', 1)]

### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 2, ASN: 2, ASP: 4, GLN: 4, GLU: 13, GLY: 6, ILE: 2, LEU: 12, LYS: 10, PHE: 5, PRO: 3, SER: 6, THR: 2, TYR: 1, VAL: 3', 'Total: 75', 'Ionizable: 28',
              'Ratio: 37.3%')

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "LYS A   1"
 - <strong>CTR</strong>: "GLN A  75"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05

### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A  75":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.54: " CA  NTR A   1" to " CB  LYS A   1"
- d= 1.86: "HD21 ASN A  21" to " OE1 GLN A  22"
- d= 1.80: "HE21 GLN A  22" to " OE1 GLU A  60"
- d= 1.93: " O   LYS A  29" to "HE21 GLN A  33"
- d= 1.82: "HD21 ASN A  56" to " OD2 ASP A  58"
- d= 1.80: " O   PHE A  63" to "HE21 GLN A  67"
- d= 1.82: " O   LYS A  71" to "HE21 GLN A  75"

</details>

