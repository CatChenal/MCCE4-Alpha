---
# 4EQP :: Crystal Structure Of Staphylococcal Nuclease Variant Delta+Phs I72D At Cryogenic Temperature;
## PDB.Structure
### Function: HYDROLASE
### First Release: 19-APR-12
### Method: X-RAY DIFFRACTION
### Resolution: 1.35 ANGSTROMS.
### Molecule: THERMONUCLEASE
### Seqres Species: Residues: A:143; Total: 143
### Cofactors:
  -  CA:
 'CALCIUM ION', 1
  - THP:
 "THYMIDINE-3',5'-DIPHOSPHATE", 1

### Total cofactors: 2
### Total waters: 121
### Missing:
  - Residues:
 'A': [('ALA', 1), ('THR', 2), ('SER', 3), ('THR', 4), ('LYS', 5), ('LYS', 6), ('GLU', 142), ('ASP', 143), ('ASN', 144), ('ALA', 145), ('ASP', 146), ('SER', 147), ('GLY', 148), ('GLN', 149),
       ('Count', 14)]

### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 13, ARG: 5, ASN: 6, ASP: 7, GLN: 5, GLU: 11, GLY: 9, HIS: 2, ILE: 4, LEU: 12, LYS: 18, MET: 4, PHE: 4, PRO: 4, SER: 2, THR: 7, TRP: 1, TYR: 7, VAL: 8', 'Total: 129', 'Ionizable: 50',
              'Ratio: 38.8%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 121, 'THP': 1, '_CA': 1

### Links:
  - OD2 ASP A 21 -- 2.82 Å --> CA  CA A 201
  - OD1 ASP A 40 -- 2.74 Å --> CA  CA A 201
  - O  THR A 41 -- 2.88 Å --> CA  CA A 201
  - OE2 GLU A 43 -- 3.02 Å --> CA  CA A 201
  - CA  CA A 201 -- 3.12 Å --> O5P THP A 202
  - CA  CA A 201 -- 2.81 Å --> O  HOH A 334
  - CA  CA A 201 -- 2.88 Å --> O  HOH A 341

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE CA A 201', 'ASP A  21', 'ASP A  40', 'THR A  41', 'GLU A  43', 'THP A 202', 'HOH A 334', 'HOH A 341']
  - AC2: ['BINDING SITE FOR RESIDUE THP A 202', 'ARG A  35', 'ASP A  40', 'LYS A  84', 'TYR A  85', 'ARG A  87', 'LEU A  89', 'TYR A 113', 'TYR A 115', 'LYS A 127', 'CA A 201', 'HOH A 308', 'HOH A 312', 'HOH A 315', 'HOH A 323', 'HOH A 332', 'HOH A 334', 'HOH A 335', 'HOH A 348', 'HOH A 349', 'HOH A 359', 'HOH A 365', 'HOH A 371', 'HOH A 379']

## MCCE.Step1
### Renamed:
  - "CA    CA A 201" to "CA   _CA A 201"

### Termini:
 - <strong>NTR</strong>: "LEU A   7"
 - <strong>CTR</strong>: "SER A 141"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
THP: https://pubchem.ncbi.nlm.nih.gov/#query=THP&tab=substance; 

### Free Cofactors:
  - Removed all 121 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 121.
  - Species and properties with assigned default values in debug.log:

  - _CA+2: ['VDW_RAD', 'VDW_EPS']

  - THPBK: ['VDW_RAD', 'VDW_EPS']


### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 141":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.52: " CA  NTR A   7" to " CB  LEU A   7"

</details>

