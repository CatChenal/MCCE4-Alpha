---
# 2PPT :: Crystal Structure Of Thioredoxin-2
## PDB.Structure
### Function: OXIDOREDUCTASE
### First Release: 30-APR-07
### Method: X-RAY DIFFRACTION
### Resolution: 1.92 ANGSTROMS.
### Molecule: THIOREDOXIN-2
### Seqres Species: Residues: A:155, B:155; Total: 310
### Cofactors:
  -  ZN:
 'ZINC ION', 3

### Total cofactors: 3
### Total waters: 127
### Missing:
  - Residues:
 'A': [('GLY', -2), ('SER', -1), ('HIS', 0), ('ALA', 150), ('ARG', 151), ('ALA', 152), ('Count', 6)],
 'B': [('GLY', -2), ('SER', -1), ('HIS', 0), ('MET', 1), ('MET', 2), ('GLY', 3), ('ALA', 4), ('LYS', 5), ('MET', 6), ('ALA', 150), ('ARG', 151), ('ALA', 152), ('Count', 12)]

### Biounits: DIMERIC: A, B;
### Models: 1
### Chains: A, B
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 28, ARG: 12, ASN: 1, ASP: 6, CYS: 6, GLN: 7, GLU: 4, GLY: 17, HIS: 3, ILE: 8, LEU: 15, LYS: 7, MET: 4, PHE: 5, PRO: 10, SER: 3, THR: 4, TRP: 2, VAL: 7', 'Total: 149', 'Ionizable: 38',
              'Ratio: 25.5%')
  - B:
 'RESIDUES': ('ALA: 27, ARG: 12, ASN: 1, ASP: 6, CYS: 6, GLN: 7, GLU: 4, GLY: 16, HIS: 3, ILE: 8, LEU: 15, LYS: 6, MET: 1, PHE: 5, PRO: 10, SER: 3, THR: 4, TRP: 2, VAL: 7', 'Total: 143', 'Ionizable: 37',
              'Ratio: 25.9%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 76, '_ZN': 1
  - B:
 'HOH': 51, '_ZN': 2

### Links:
  - SG CYS A 14 -- 2.32 Å --> ZN  ZN A 300
  - SG CYS A 17 -- 2.46 Å --> ZN  ZN A 300
  - SG CYS A 34 -- 2.36 Å --> ZN  ZN A 300
  - SG CYS A 37 -- 2.30 Å --> ZN  ZN A 300
  - SG CYS B 14 -- 2.30 Å --> ZN  ZN B 301
  - SG CYS B 17 -- 2.15 Å --> ZN  ZN B 301
  - SG CYS B 34 -- 2.16 Å --> ZN  ZN B 301
  - SG CYS B 37 -- 2.05 Å --> ZN  ZN B 301
  - OE2 GLU B 58 -- 1.85 Å --> ZN  ZN B 302
  - NE2 HIS B 111 -- 2.12 Å --> ZN  ZN B 302

## MCCE.Step1
### Renamed:
  - "ZN    ZN A 300" to "ZN   _ZN A 300"
  - "ZN    ZN B 301" to "ZN   _ZN B 301"
  - "ZN    ZN B 302" to "ZN   _ZN B 302"

### Termini:
 - <strong>NTR</strong>: "MET A   1", "ALA B   7"
 - <strong>CTR</strong>: "GLY A 149", "GLY B 149"

### Free Cofactors:
  - Removed all 76 HOH in A. Removed all 51 HOH in B.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 127.

### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 149":   OXT
  -    Missing heavy atoms for CTR01 in "CTR B 149":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   1" to " CB  MET A   1"
- d= 1.53: " CA  NTR B   7" to " CB  ALA B   7"
- d= 1.85: " OE2 GLU B  58" to "ZN   _ZN B 302"

</details>

