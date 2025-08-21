---
# 2RBM :: Crystal Structure Of Staphylococcal Nuclease Variant Delta+Phs I72K At Cryogenic Temperature;
## PDB.Structure
### Function: HYDROLASE
### First Release: 19-SEP-07
### Method: X-RAY DIFFRACTION
### Resolution: 1.90 ANGSTROMS.
### Molecule: THERMONUCLEASE
### Seqres Species: Residues: A:143; Total: 143
### Cofactors:
  -  CA:
 'CALCIUM ION', 1
  - PO4:
 'PHOSPHATE ION', 1
  - THP:
 "THYMIDINE-3',5'-DIPHOSPHATE", 1

### Total cofactors: 3
### Total waters: 128
### Missing:
  - Residues:
 'A': [('ALA', 1), ('THR', 2), ('SER', 3), ('THR', 4), ('LYS', 5), ('GLU', 142), ('ASP', 143), ('ASN', 144), ('ALA', 145), ('ASP', 146), ('SER', 147), ('GLY', 148), ('GLN', 149),
       ('Count', 13)]

### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 13, ARG: 5, ASN: 6, ASP: 6, GLN: 5, GLU: 11, GLY: 9, HIS: 2, ILE: 4, LEU: 12, LYS: 20, MET: 4, PHE: 4, PRO: 4, SER: 2, THR: 7, TRP: 1, TYR: 7, VAL: 8', 'Total: 130', 'Ionizable: 51',
              'Ratio: 39.2%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 128, 'PO4': 1, 'THP': 1, '_CA': 1

### Links:
  - OD2 ASP A 21 -- 2.81 Å --> CA  CA A 150
  - OD1 ASP A 40 -- 2.72 Å --> CA  CA A 150
  - O  THR A 41 -- 2.81 Å --> CA  CA A 150
  - CA  CA A 150 -- 2.87 Å --> O  HOH A 161
  - CA  CA A 150 -- 2.85 Å --> O  HOH A 164

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE CA A 150', 'ASP A  21', 'ASP A  40', 'THR A  41', 'GLU A  43', 'THP A 152', 'HOH A 161', 'HOH A 164', 'HOH A 281']
  - AC2: ['BINDING SITE FOR RESIDUE PO4 A 151', 'GLU A  52', 'MET A  65', 'HOH A 267']
  - AC3: ['BINDING SITE FOR RESIDUE THP A 152', 'ARG A  35', 'ASP A  40', 'ASP A  83', 'LYS A  84', 'TYR A  85', 'ARG A  87', 'LEU A  89', 'TYR A 113', 'TYR A 115', 'LYS A 127', 'CA A 150', 'HOH A 155', 'HOH A 159', 'HOH A 160', 'HOH A 161', 'HOH A 176', 'HOH A 214', 'HOH A 225', 'HOH A 228', 'HOH A 281']

## MCCE.Step1
### Renamed:
  - "CA    CA A 150" to "CA   _CA A 150"

### Termini:
 - <strong>NTR</strong>: "LYS A   6"
 - <strong>CTR</strong>: "SER A 141"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
PO4: https://pubchem.ncbi.nlm.nih.gov/#query=PO4&tab=substance; THP: https://pubchem.ncbi.nlm.nih.gov/#query=THP&tab=substance; 

### Free Cofactors:
  - Removed all 128 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 128.
  - Species and properties with assigned default values in debug.log:

  - _CA+2: ['VDW_RAD', 'VDW_EPS']

  - PO4BK: ['VDW_RAD', 'VDW_EPS']

  - THPBK: ['VDW_RAD', 'VDW_EPS']


### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 141":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.55: " CA  NTR A   6" to " CB  LYS A   6"

</details>

