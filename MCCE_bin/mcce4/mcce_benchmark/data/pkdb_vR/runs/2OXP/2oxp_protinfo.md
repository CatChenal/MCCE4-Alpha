---
# 2OXP :: Crystal Structure Of Staphylococcal Nuclease Mutant V66D/P117G/H124L/S128A;
## PDB.Structure
### Function: HYDROLASE
### First Release: 20-FEB-07
### Method: X-RAY DIFFRACTION
### Resolution: 2.00 ANGSTROMS.
### Molecule: THERMONUCLEASE
### Seqres Species: Residues: A:149; Total: 149
### Cofactors:
  - THP:
 "THYMIDINE-3',5'-DIPHOSPHATE", 1

### Total cofactors: 1
### Total waters: 77
### Missing:
  - Residues:
 'A': [('ALA', 1), ('THR', 2), ('SER', 3), ('THR', 4), ('LYS', 5), ('LYS', 6), ('LYS', 45), ('HIS', 46), ('PRO', 47), ('LYS', 48), ('LYS', 49), ('GLY', 50), ('GLU', 142), ('ASP', 143), ('ASN', 144), ('ALA', 145), ('ASP', 146), ('SER', 147), ('GLY', 148), ('GLN', 149),
       ('Count', 20)]

### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 13, ARG: 5, ASN: 5, ASP: 7, GLN: 5, GLU: 11, GLY: 9, HIS: 2, ILE: 5, LEU: 12, LYS: 18, MET: 4, PHE: 3, PRO: 4, SER: 2, THR: 8, TRP: 1, TYR: 7, VAL: 8', 'Total: 129', 'Ionizable: 50',
              'Ratio: 38.8%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 77, 'THP': 1

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE THP A 300', 'ARG A  35', 'LYS A  70', 'LYS A  71', 'LYS A  84', 'TYR A  85', 'ARG A  87', 'LEU A  89', 'TYR A 113', 'TYR A 115', 'HOH A 220', 'HOH A 223', 'HOH A 248', 'HOH A 256', 'HOH A 260', 'HOH A 266']

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "LEU A   7", "VAL A  51"
 - <strong>CTR</strong>: "THR A  44", "SER A 141"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
THP: https://pubchem.ncbi.nlm.nih.gov/#query=THP&tab=substance; 

### Free Cofactors:
  - Removed all 77 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 77.
  - Species and properties with assigned default values in debug.log:

  - THPBK: ['VDW_RAD', 'VDW_EPS']


### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A  44":   OXT
  -    Missing heavy atoms for CTR01 in "CTR A 141":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   7" to " CB  LEU A   7"
- d= 1.54: " CA  NTR A  51" to " CB  VAL A  51"

</details>

