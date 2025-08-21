---
# 4HMI :: Crystal Structure Of Staphylococcal Nuclease Variant Delta+Phs V99K At Cryogenic Temperature;
## PDB.Structure
### Function: HYDROLASE
### First Release: 18-OCT-12
### Method: X-RAY DIFFRACTION
### Resolution: 1.75 ANGSTROMS.
### Molecule: THERMONUCLEASE
### Seqres Species: Residues: A:143; Total: 143
### Cofactors:
  -  CA:
 'CALCIUM ION', 1
  - THP:
 "THYMIDINE-3',5'-DIPHOSPHATE", 1

### Total cofactors: 2
### Total waters: 90
### Missing:
  - Residues:
 'A': [('ALA', 1), ('THR', 2), ('SER', 3), ('THR', 4), ('LYS', 5), ('GLU', 142), ('ASP', 143), ('ASN', 144), ('ALA', 145), ('ASP', 146), ('SER', 147), ('GLY', 148), ('GLN', 149),
       ('Count', 13)]

### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 13, ARG: 5, ASN: 6, ASP: 6, GLN: 5, GLU: 11, GLY: 9, HIS: 2, ILE: 5, LEU: 12, LYS: 20, MET: 4, PHE: 4, PRO: 4, SER: 2, THR: 7, TRP: 1, TYR: 7, VAL: 7', 'Total: 130', 'Ionizable: 51',
              'Ratio: 39.2%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 90, 'THP': 1, '_CA': 1

### Links:
  - OD2 ASP A 21 -- 2.91 Å --> CA  CA A 201
  - OD1 ASP A 40 -- 2.74 Å --> CA  CA A 201
  - O  THR A 41 -- 2.92 Å --> CA  CA A 201
  - OE2 GLU A 43 -- 2.89 Å --> CA  CA A 201
  - CA  CA A 201 -- 3.11 Å --> O5P THP A 202
  - CA  CA A 201 -- 2.81 Å --> O  HOH A 322
  - CA  CA A 201 -- 2.84 Å --> O  HOH A 331

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE CA A 201', 'ASP A  21', 'ASP A  40', 'THR A  41', 'GLU A  43', 'THP A 202', 'HOH A 322', 'HOH A 331']
  - AC2: ['BINDING SITE FOR RESIDUE THP A 202', 'ARG A  35', 'ASP A  40', 'LYS A  84', 'TYR A  85', 'ARG A  87', 'LEU A  89', 'TYR A 113', 'TYR A 115', 'LYS A 127', 'CA A 201', 'HOH A 305', 'HOH A 306', 'HOH A 309', 'HOH A 311', 'HOH A 313', 'HOH A 322', 'HOH A 329', 'HOH A 342', 'HOH A 358', 'HOH A 366', 'HOH A 369', 'HOH A 372']

## MCCE.Step1
### Renamed:
  - "CA    CA A 201" to "CA   _CA A 201"

### Termini:
 - <strong>NTR</strong>: "LYS A   6"
 - <strong>CTR</strong>: "SER A 141"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
THP: https://pubchem.ncbi.nlm.nih.gov/#query=THP&tab=substance; 

### Free Cofactors:
  - Removed all 90 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 90.
  - Species and properties with assigned default values in debug.log:

  - _CA+2: ['VDW_RAD', 'VDW_EPS']

  - THPBK: ['VDW_RAD', 'VDW_EPS']


### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 141":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.54: " CA  NTR A   6" to " CB  LYS A   6"

</details>

