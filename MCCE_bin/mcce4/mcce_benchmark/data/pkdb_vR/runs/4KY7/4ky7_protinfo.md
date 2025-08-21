---
# 4KY7 :: Crystal Structure Of Staphylococcal Nuclease Variant Delta+Phs V74D At Cryogenic Temperature;
## PDB.Structure
### Function: HYDROLASE
### First Release: 28-MAY-13
### Method: X-RAY DIFFRACTION
### Resolution: 1.55 ANGSTROMS.
### Molecule: THERMONUCLEASE
### Seqres Species: Residues: A:143; Total: 143
### Cofactors:
  - THP:
 "THYMIDINE-3',5'-DIPHOSPHATE", 1
  -  CA:
 'CALCIUM ION', 1

### Total cofactors: 2
### Total waters: 132
### Missing:
  - Residues:
 'A': [('ALA', 1), ('THR', 2), ('SER', 3), ('THR', 4), ('ASP', 143), ('ASN', 144), ('ALA', 145), ('ASP', 146), ('SER', 147), ('GLY', 148), ('GLN', 149),
       ('Count', 11)]

### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 13, ARG: 5, ASN: 6, ASP: 7, GLN: 5, GLU: 12, GLY: 9, HIS: 2, ILE: 5, LEU: 12, LYS: 20, MET: 4, PHE: 4, PRO: 4, SER: 2, THR: 7, TRP: 1, TYR: 7, VAL: 7', 'Total: 132', 'Ionizable: 53',
              'Ratio: 40.2%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 132, 'THP': 1, '_CA': 1

### Links:
  - OD2 ASP A 21 -- 2.79 Å --> CA  CA A 202
  - OD1 ASP A 40 -- 2.62 Å --> CA  CA A 202
  - O  THR A 41 -- 2.80 Å --> CA  CA A 202
  - O4P THP A 201 -- 3.12 Å --> CA  CA A 202
  - CA  CA A 202 -- 2.89 Å --> O  HOH A 341
  - CA  CA A 202 -- 2.72 Å --> O  HOH A 353
  - CA  CA A 202 -- 3.11 Å --> O  HOH A 423

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE THP A 201', 'ARG A  35', 'ASP A  40', 'LYS A  84', 'TYR A  85', 'ARG A  87', 'LEU A  89', 'TYR A 113', 'TYR A 115', 'LYS A 127', 'CA A 202', 'HOH A 328', 'HOH A 338', 'HOH A 348', 'HOH A 351', 'HOH A 353', 'HOH A 357', 'HOH A 368', 'HOH A 381', 'HOH A 401', 'HOH A 410', 'HOH A 423']
  - AC2: ['BINDING SITE FOR RESIDUE CA A 202', 'ASP A  21', 'ASP A  40', 'THR A  41', 'THP A 201', 'HOH A 341', 'HOH A 353']

## MCCE.Step1
### Renamed:
  - "CA    CA A 202" to "CA   _CA A 202"

### Termini:
 - <strong>NTR</strong>: "LYS A   5"
 - <strong>CTR</strong>: "GLU A 142"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
THP: https://pubchem.ncbi.nlm.nih.gov/#query=THP&tab=substance; 

### Free Cofactors:
  - Removed all 132 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 132.
  - Species and properties with assigned default values in debug.log:

  - THPBK: ['VDW_RAD', 'VDW_EPS']

  - _CA+2: ['VDW_RAD', 'VDW_EPS']


### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 142":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   5" to " CB  LYS A   5"

</details>

