---
# 3E5S :: Crystal Structure Of Staphylococcal Nuclease Variant Delta+Phs L103K At Cryogenic Temperature;
## PDB.Structure
### Function: HYDROLASE
### First Release: 14-AUG-08
### Method: X-RAY DIFFRACTION
### Resolution: 2.00 ANGSTROMS.
### Molecule: THERMONUCLEASE
### Seqres Species: Residues: A:143, B:143; Total: 286
### Cofactors:
  - THP:
 "THYMIDINE-3',5'-DIPHOSPHATE", 2
  -  CA:
 'CALCIUM ION', 2

### Total cofactors: 4
### Total waters: 151
### Missing:
  - Residues:
 'A': [('ALA', 1), ('THR', 2), ('SER', 3), ('THR', 4), ('LYS', 5), ('LYS', 6), ('GLU', 142), ('ASP', 143), ('ASN', 144), ('ALA', 145), ('ASP', 146), ('SER', 147), ('GLY', 148), ('GLN', 149), ('Count', 14)],
 'B': [('ALA', 1), ('THR', 2), ('SER', 3), ('THR', 4), ('LYS', 5), ('LYS', 6), ('GLU', 142), ('ASP', 143), ('ASN', 144), ('ALA', 145), ('ASP', 146), ('SER', 147), ('GLY', 148), ('GLN', 149), ('Count', 14)]

### Biounits: MONOMERIC: A,B;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 13, ARG: 5, ASN: 6, ASP: 6, GLN: 5, GLU: 11, GLY: 9, HIS: 2, ILE: 5, LEU: 11, LYS: 19, MET: 4, PHE: 4, PRO: 4, SER: 2, THR: 7, TRP: 1, TYR: 7, VAL: 8', 'Total: 129', 'Ionizable: 50',
              'Ratio: 38.8%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 76, 'THP': 1, '_CA': 1

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE THP A 150', 'ARG A  35', 'LYS A  84', 'TYR A  85', 'ARG A  87', 'TYR A 113', 'TYR A 115', 'LYS A 127', 'CA A 151', 'HOH A 153', 'HOH A 157', 'HOH A 159', 'HOH A 160', 'HOH A 183', 'HOH A 189', 'HOH A 208']
  - AC2: ['BINDING SITE FOR RESIDUE CA A 151', 'ASP A  21', 'ARG A  35', 'ASP A  40', 'GLU A  43', 'THP A 150', 'HOH A 156', 'HOH A 208', 'HOH A 224']
  - AC3: ['BINDING SITE FOR RESIDUE THP B 150', 'ARG B  35', 'LYS B  84', 'TYR B  85', 'ARG B  87', 'LEU B  89', 'TYR B 113', 'TYR B 115', 'LYS B 127', 'CA B 151', 'HOH B 153', 'HOH B 159', 'HOH B 166', 'HOH B 172', 'HOH B 183', 'HOH B 196', 'HOH B 205', 'HOH B 226']
  - AC4: ['BINDING SITE FOR RESIDUE CA B 151', 'ASP B  21', 'ASP B  40', 'THP B 150', 'HOH B 155', 'HOH B 167']

## MCCE.Step1
### Renamed:
  - "CA    CA A 151" to "CA   _CA A 151"

### Termini:
 - <strong>NTR</strong>: "LEU A   7"
 - <strong>CTR</strong>: "SER A 141"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
THP: https://pubchem.ncbi.nlm.nih.gov/#query=THP&tab=substance; 

### Free Cofactors:
  - Removed all 76 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 76.
  - Species and properties with assigned default values in debug.log:

  - THPBK: ['VDW_RAD', 'VDW_EPS']

  - _CA+2: ['VDW_RAD', 'VDW_EPS']


### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 141":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   7" to " CB  LEU A   7"

</details>

