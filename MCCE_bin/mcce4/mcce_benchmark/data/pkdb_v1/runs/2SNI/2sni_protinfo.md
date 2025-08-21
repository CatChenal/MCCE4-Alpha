---
# 2SNI :: Structural Comparison Of Two Serine Proteinase-Protein Inhibitor Complexes. Eglin-C-Subtilisin Carlsberg And Ci-2-Subtilisin Novo;
## PDB.Structure
### Function: COMPLEX (PROTEINASE/INHIBITOR)
### First Release: 05-SEP-88
### Method: X-RAY DIFFRACTION
### Resolution: 2.10 ANGSTROMS.
### Molecule: SUBTILISIN NOVO, CHYMOTRYPSIN INHIBITOR 2
### Seqres Species: Residues: E:275, I:83; Total: 358
### Cofactors:
  -  CA:
 'CALCIUM ION', 2

### Total cofactors: 2
### Total waters: 168
### Missing:
  - Residues:
 'I': [('SER', 1), ('SER', 2), ('VAL', 3), ('GLU', 4), ('LYS', 5), ('LYS', 6), ('PRO', 7), ('GLU', 8), ('GLY', 9), ('VAL', 10), ('ASN', 11), ('THR', 12), ('GLY', 13), ('ALA', 14), ('GLY', 15), ('ASP', 16), ('ARG', 17), ('HIS', 18), ('ASN', 19),
       ('Count', 19)]

### Biounits: DIMERIC: E, I;
### Models: 1
### Chains: E, I
### Model 1 Residues:
  - E:
 'RESIDUES': ('ALA: 37, ARG: 2, ASN: 18, ASP: 10, GLN: 11, GLU: 4, GLY: 33, HIS: 6, ILE: 13, LEU: 15, LYS: 11, MET: 5, PHE: 3, PRO: 14, SER: 37, THR: 13, TRP: 3, TYR: 10, VAL: 30', 'Total: 275', 'Ionizable: 43',
              'Ratio: 15.6%')
  - I:
 'RESIDUES': ('ALA: 3, ARG: 4, ASN: 1, ASP: 4, GLN: 2, GLU: 7, GLY: 3, ILE: 6, LEU: 6, LYS: 6, MET: 1, PHE: 1, PRO: 4, SER: 1, THR: 3, TRP: 1, TYR: 1, VAL: 10', 'Total: 64', 'Ionizable: 22',
              'Ratio: 34.4%')

### Model 1 Free Cofactors & Waters:
  - E:
 'HOH': 123, '_CA': 2
  - I:
 'HOH': 45

### Links:
  - NE2 GLN E  2 -- 2.38 Å --> CA  CA E 276
  - OD1 ASP E 41 -- 2.43 Å --> CA  CA E 276
  - OD2 ASP E 41 -- 2.56 Å --> CA  CA E 276
  - O  LEU E 75 -- 2.34 Å --> CA  CA E 276
  - ND2 ASN E 77 -- 2.43 Å --> CA  CA E 276
  - O  ILE E 79 -- 2.35 Å --> CA  CA E 276
  - O  VAL E 81 -- 2.31 Å --> CA  CA E 276
  - O  GLY E 169 -- 2.85 Å --> CA  CA E 277
  - O  TYR E 171 -- 3.01 Å --> CA  CA E 277
  - O  VAL E 174 -- 2.83 Å --> CA  CA E 277
  - O  GLU E 195 -- 3.03 Å --> CA  CA E 277
  - OD2 ASP E 197 -- 2.81 Å --> CA  CA E 277
  - CA  CA E 277 -- 2.92 Å --> O  HOH E 391

### Sites:
  - ACT: ['catalytic site', 'ASP E  32', 'HIS E  64', 'SER E 221']
  - IO1: ['ion binding site', 'GLN E   2', 'ASP E  41', 'LEU E  75', 'ASN E  77', 'ILE E  79', 'VAL E  81']
  - IO2: ['ion binding site', 'GLY E 169', 'TYR E 171', 'VAL E 174', 'GLU E 195', 'ASP E 197', 'HOH E 391']
  - RSB: ['inhibitor reactive site', 'MET I  59', 'GLU I  60']
  - AC1: ['BINDING SITE FOR RESIDUE CA E 276', 'GLN E   2', 'ASP E  41', 'LEU E  75', 'ASN E  77', 'ILE E  79', 'VAL E  81']
  - AC2: ['BINDING SITE FOR RESIDUE CA E 277', 'GLY E 169', 'TYR E 171', 'VAL E 174', 'ASP E 197', 'HOH E 391']

## MCCE.Step1
### Renamed:
  - "CA    CA E 276" to "CA   _CA E 276"
  - "CA    CA E 277" to "CA   _CA E 277"

### Termini:
 - <strong>NTR</strong>: "ALA E   1", "LEU I  20"
 - <strong>CTR</strong>: "GLN E 275", "GLY I  83"

### Free Cofactors:
  - Removed all 123 HOH in E.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 168.
  - Species and properties with assigned default values in debug.log:

  - _CA+2: ['VDW_RAD', 'VDW_EPS']


### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.52: " CA  NTR E   1" to " CB  ALA E   1"
- d= 1.53: " CA  NTR I  20" to " CB  LEU I  20"

</details>

