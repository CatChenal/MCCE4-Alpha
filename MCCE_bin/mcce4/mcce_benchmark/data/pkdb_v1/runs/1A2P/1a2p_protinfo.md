---
# 1A2P :: Barnase Wildtype Structure At 1.5 Angstroms Resolution
## PDB.Structure
### Function: RIBONUCLEASE
### First Release: 07-JAN-98
### Method: X-RAY DIFFRACTION
### Resolution: 1.50 ANGSTROMS.
### Molecule: BARNASE
### Seqres Species: Residues: A:110, B:110, C:110; Total: 330
### Cofactors:
  -  ZN:
 'ZINC ION', 3

### Total cofactors: 3
### Total waters: 415
### Missing:
  - Residues:
 'A': [('ALA', 1), ('GLN', 2), ('Count', 2)],
 'B': [('ALA', 1), ('GLN', 2), ('Count', 2)],
 'C': [('ALA', 1), ('GLN', 2), ('Count', 2)]

### Biounits: MONOMERIC: A,B,C;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 7, ARG: 6, ASN: 6, ASP: 9, GLN: 3, GLU: 3, GLY: 10, HIS: 2, ILE: 8, LEU: 7, LYS: 8, PHE: 4, PRO: 3, SER: 9, THR: 9, TRP: 3, TYR: 7, VAL: 4', 'Total: 108', 'Ionizable: 35',
              'Ratio: 32.4%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 134, '_ZN': 1

### Links:
  - ND1 HIS A 18 -- 2.36 Å --> ZN  ZN A 112
  - OE2 GLU A 60 -- 2.21 Å --> ZN  ZN A 112
  - ZN  ZN A 112 -- 2.25 Å --> O  HOH A 227
  - ND1 HIS B 18 -- 2.21 Å --> ZN  ZN B 112
  - OE1 GLU B 60 -- 2.25 Å --> ZN  ZN B 112
  - OE2 GLU B 60 -- 2.32 Å --> ZN  ZN B 112
  - ZN  ZN B 112 -- 2.42 Å --> O  HOH B 260
  - ND1 HIS C 18 -- 2.06 Å --> ZN  ZN C 112
  - OE1 GLU C 60 -- 2.32 Å --> ZN  ZN C 112
  - OE2 GLU C 60 -- 2.12 Å --> ZN  ZN C 112
  - NZ LYS C 62 -- 2.09 Å --> ZN  ZN C 112
  - ZN  ZN C 112 -- 1.99 Å --> O  HOH C 139

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE ZN A 112', 'HIS A  18', 'GLU A  60', 'LYS A  62', 'HOH A 227']
  - AC2: ['BINDING SITE FOR RESIDUE ZN B 112', 'HIS B  18', 'GLU B  60', 'LYS B  62', 'HOH B 260']
  - AC3: ['BINDING SITE FOR RESIDUE ZN C 112', 'HIS C  18', 'GLU C  60', 'LYS C  62', 'HOH C 139']

## MCCE.Step1
### Renamed:
  - "ZN    ZN A 112" to "ZN   _ZN A 112"

### Termini:
 - <strong>NTR</strong>: "VAL A   3"
 - <strong>CTR</strong>: "ARG A 110"

### Free Cofactors:
  - Removed all 134 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 134.

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   3" to " CB  VAL A   3"

</details>

