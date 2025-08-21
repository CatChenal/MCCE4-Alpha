---
# 1XOB :: Thioredoxin (Reduced Dithio Form), Nmr, 20 Structures
## PDB.Structure
### Function: ELECTRON TRANSPORT
### First Release: 28-NOV-95
### Method: SOLUTION NMR
### Resolution: NOT APPLICABLE.
### Molecule: THIOREDOXIN, REDUCED DITHIO FORM
### Seqres Species: Residues: A:108; Total: 108
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 12, ARG: 1, ASN: 4, ASP: 11, CYS: 2, GLN: 3, GLU: 5, GLY: 9, HIS: 1, ILE: 9, LEU: 13, LYS: 10, MET: 1, PHE: 4, PRO: 5, SER: 3, THR: 6, TRP: 2, TYR: 2, VAL: 5', 'Total: 108', 'Ionizable: 32',
              'Ratio: 29.6%')

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "SER A   1"
 - <strong>CTR</strong>: "ALA A 108"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05

### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 108":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   1" to " CB  SER A   1"
- d= 1.86: " OD1 ASP A  20" to "HD22 ASN A  83"
- d= 1.73: " OE1 GLU A  30" to "HD22 ASN A  59"
- d= 1.82: " O   ASP A  47" to "HE21 GLN A  50"
- d= 1.84: " OD1 ASN A  59" to "HE21 GLN A  62"
- d= 1.90: " O   ASN A  59" to "HD21 ASN A  63"
- d= 1.89: "HE21 GLN A  98" to " OE2 GLU A 101"
- d= 1.81: " O   PHE A 102" to "HD21 ASN A 106"

</details>

