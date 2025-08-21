---
# 1TT2 :: Cryogenic Crystal Structure Of Staphylococcal Nuclease Variant Truncated Delta+Phs I92K;
## PDB.Structure
### Function: HYDROLASE
### First Release: 21-JUN-04
### Method: X-RAY DIFFRACTION
### Resolution: 1.85 ANGSTROMS.
### Molecule: THERMONUCLEASE
### Seqres Species: Residues: A:138; Total: 138
### Cofactors:
  -  CA:
 'CALCIUM ION', 1
  - DMS:
 'DIMETHYL SULFOXIDE', 2
  - THP:
 "THYMIDINE-3',5'-DIPHOSPHATE", 1
  - GOL:
 'GLYCERIN; PROPANE-1,2,3-TRIOL', 1

### Total cofactors: 5
### Total waters: 116
### Missing:
  - Residues:
 'A': [('ALA', 1), ('THR', 2), ('SER', 3), ('THR', 4), ('LYS', 5), ('GLU', 142), ('ASP', 143), ('ASN', 144),
       ('Count', 8)]

### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 13, ARG: 5, ASN: 6, ASP: 6, GLN: 5, GLU: 11, GLY: 9, HIS: 2, ILE: 4, LEU: 12, LYS: 20, MET: 4, PHE: 4, PRO: 4, SER: 2, THR: 7, TRP: 1, TYR: 7, VAL: 8', 'Total: 130', 'Ionizable: 51',
              'Ratio: 39.2%')

### Model 1 Free Cofactors & Waters:
  - A:
 'DMS': 2, 'GOL': 1, 'HOH': 116, 'THP': 1, '_CA': 1

### Links:
  - OD1 ASP A 21 -- 2.48 Å --> CA  CA A 502
  - OD1 ASP A 40 -- 2.45 Å --> CA  CA A 502
  - O  THR A 41 -- 2.56 Å --> CA  CA A 502
  - O  HOH A 216 -- 2.69 Å --> CA  CA A 502
  - O  HOH A 217 -- 2.61 Å --> CA  CA A 502
  - O  HOH A 218 -- 2.35 Å --> CA  CA A 502
  - O4P THP A 501 -- 3.09 Å --> CA  CA A 502

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE CA A 502', 'ASP A  21', 'ASP A  40', 'THR A  41', 'HOH A 216', 'HOH A 217', 'HOH A 218', 'THP A 501']
  - AC2: ['BINDING SITE FOR RESIDUE DMS A 402', 'GLY A  79', 'GLN A  80', 'GLN A 106', 'TYR A 115', 'LYS A 116']
  - AC3: ['BINDING SITE FOR RESIDUE DMS A 403', 'PHE A  76']
  - AC4: ['BINDING SITE FOR RESIDUE THP A 501', 'ARG A  35', 'ASP A  40', 'LYS A  71', 'LYS A  84', 'TYR A  85', 'ARG A  87', 'LEU A  89', 'TYR A 113', 'TYR A 115', 'HOH A 201', 'HOH A 204', 'HOH A 206', 'HOH A 207', 'HOH A 217', 'HOH A 219', 'HOH A 220', 'HOH A 221', 'HOH A 249', 'HOH A 250', 'CA A 502']
  - AC5: ['BINDING SITE FOR RESIDUE GOL A 401', 'ALA A  12', 'LEU A  14', 'LYS A  70', 'ILE A  72', 'TYR A  85', 'GLY A  86', 'HOH A 279']

## MCCE.Step1
### Renamed:
  - "CA    CA A 502" to "CA   _CA A 502"

### Termini:
 - <strong>NTR</strong>: "LYS A   6"
 - <strong>CTR</strong>: "SER A 141"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
DMS: https://pubchem.ncbi.nlm.nih.gov/#query=DMS&tab=substance; THP: https://pubchem.ncbi.nlm.nih.gov/#query=THP&tab=substance; GOL: https://pubchem.ncbi.nlm.nih.gov/#query=GOL&tab=substance; 

### Free Cofactors:
  - Removed all 116 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 116.
  - Species and properties with assigned default values in debug.log:

  - _CA+2: ['VDW_RAD', 'VDW_EPS']

  - DMSBK: ['VDW_RAD', 'VDW_EPS']

  - THPBK: ['VDW_RAD', 'VDW_EPS']

  - GOLBK: ['VDW_RAD', 'VDW_EPS']


### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 141":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   6" to " CB  LYS A   6"

</details>

