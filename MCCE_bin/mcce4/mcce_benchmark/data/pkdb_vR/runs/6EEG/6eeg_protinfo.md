---
# 6EEG :: Crystal Structure Of Staphylcoccal Nuclease Variant Delta+Phs T62E At Cryogenic Temperature;
## PDB.Structure
### Function: HYDROLASE
### First Release: 14-AUG-18
### Method: X-RAY DIFFRACTION
### Resolution: 1.95 ANGSTROMS.
### Molecule: THERMONUCLEASE
### Seqres Species: Residues: A:143; Total: 143
### Cofactors:
  - THP:
 "THYMIDINE-3',5'-DIPHOSPHATE", 1
  -  CA:
 'CALCIUM ION', 1

### Total cofactors: 2
### Total waters: 82
### Missing:
  - Residues:
 'A': [('ALA', 1), ('THR', 2), ('SER', 3), ('THR', 4), ('LYS', 5), ('LYS', 6), ('GLU', 142), ('ASP', 143), ('ASN', 144), ('ALA', 145), ('ASP', 146), ('SER', 147), ('GLY', 148), ('GLN', 149),
       ('Count', 14)]

### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 13, ARG: 5, ASN: 6, ASP: 6, GLN: 5, GLU: 12, GLY: 9, HIS: 2, ILE: 5, LEU: 12, LYS: 18, MET: 4, PHE: 4, PRO: 4, SER: 2, THR: 6, TRP: 1, TYR: 7, VAL: 8', 'Total: 129', 'Ionizable: 50',
              'Ratio: 38.8%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 82, 'THP': 1, '_CA': 1

### Links:
  - OD2 ASP A 21 -- 2.77 Å --> CA  CA A 202
  - OD1 ASP A 40 -- 2.72 Å --> CA  CA A 202
  - O  THR A 41 -- 2.88 Å --> CA  CA A 202
  - OE2 GLU A 43 -- 3.02 Å --> CA  CA A 202
  - O5P THP A 201 -- 3.11 Å --> CA  CA A 202
  - CA  CA A 202 -- 2.93 Å --> O  HOH A 304
  - CA  CA A 202 -- 2.75 Å --> O  HOH A 319

### Sites:
  - AC1: ['binding site for residue THP A 201', 'ARG A  35', 'ASP A  40', 'LYS A  84', 'TYR A  85', 'ARG A  87', 'LEU A  89', 'TYR A 113', 'TYR A 115', 'LYS A 127', 'CA A 202', 'HOH A 304', 'HOH A 309', 'HOH A 314', 'HOH A 316', 'HOH A 323', 'HOH A 330', 'HOH A 349', 'HOH A 351', 'HOH A 362', 'HOH A 366', 'HOH A 374']
  - AC2: ['binding site for residue CA A 202', 'ASP A  21', 'ASP A  40', 'THR A  41', 'GLU A  43', 'THP A 201', 'HOH A 304', 'HOH A 319']

## MCCE.Step1
### Renamed:
  - "CA    CA A 202" to "CA   _CA A 202"

### Termini:
 - <strong>NTR</strong>: "LEU A   7"
 - <strong>CTR</strong>: "SER A 141"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
THP: https://pubchem.ncbi.nlm.nih.gov/#query=THP&tab=substance; 

### Free Cofactors:
  - Removed all 82 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 82.
  - Species and properties with assigned default values in debug.log:

  - THPBK: ['VDW_RAD', 'VDW_EPS']

  - _CA+2: ['VDW_RAD', 'VDW_EPS']


### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 141":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.55: " CA  NTR A   7" to " CB  LEU A   7"

</details>

