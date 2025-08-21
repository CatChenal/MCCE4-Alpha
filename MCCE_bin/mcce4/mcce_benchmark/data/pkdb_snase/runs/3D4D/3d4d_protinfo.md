---
# 3D4D :: Crystal Structure Of Staphylococcal Nuclease Variant Delta+Phs Y91E At Cryogenic Temperature;
## PDB.Structure
### Function: HYDROLASE
### First Release: 14-MAY-08
### Method: X-RAY DIFFRACTION
### Resolution: 2.10 ANGSTROMS.
### Molecule: THERMONUCLEASE
### Seqres Species: Residues: A:143, B:143; Total: 286
### Cofactors:
  -  CA:
 'CALCIUM ION', 2
  - THP:
 "THYMIDINE-3',5'-DIPHOSPHATE", 2

### Total cofactors: 4
### Total waters: 183
### Missing:
  - Residues:
 'A': [('ALA', 1), ('THR', 2), ('SER', 3), ('THR', 4), ('LYS', 5), ('LYS', 6), ('GLU', 142), ('ASP', 143), ('ASN', 144), ('ALA', 145), ('ASP', 146), ('SER', 147), ('GLY', 148), ('GLN', 149), ('Count', 14)],
 'B': [('ALA', 1), ('THR', 2), ('SER', 3), ('THR', 4), ('LYS', 5), ('LYS', 6), ('GLU', 142), ('ASP', 143), ('ASN', 144), ('ALA', 145), ('ASP', 146), ('SER', 147), ('GLY', 148), ('GLN', 149), ('Count', 14)]

### Biounits: MONOMERIC: A,B;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 13, ARG: 5, ASN: 6, ASP: 6, GLN: 5, GLU: 12, GLY: 9, HIS: 2, ILE: 5, LEU: 12, LYS: 18, MET: 4, PHE: 4, PRO: 4, SER: 2, THR: 7, TRP: 1, TYR: 6, VAL: 8', 'Total: 129', 'Ionizable: 49',
              'Ratio: 38.0%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 81, 'THP': 1, '_CA': 1

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE CA A 150', 'ASP A  21', 'ASP A  40', 'THR A  41', 'GLU A  43']
  - AC2: ['BINDING SITE FOR RESIDUE CA B 150', 'ASP B  21', 'ASP B  40', 'THR B  41', 'GLU B  43']
  - AC3: ['BINDING SITE FOR RESIDUE THP A 151', 'ARG A  35', 'ASP A  40', 'LYS A  84', 'TYR A  85', 'ARG A  87', 'LEU A  89', 'TYR A 113', 'TYR A 115', 'LYS A 127']
  - AC4: ['BINDING SITE FOR RESIDUE THP B 151', 'ARG B  35', 'ASP B  40', 'LYS B  84', 'TYR B  85', 'ARG B  87', 'LEU B  89', 'TYR B 113', 'TYR B 115', 'LYS B 127']

## MCCE.Step1
### Renamed:
  - "CA    CA A 150" to "CA   _CA A 150"

### Termini:
 - <strong>NTR</strong>: "LEU A   7"
 - <strong>CTR</strong>: "SER A 141"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
THP: https://pubchem.ncbi.nlm.nih.gov/#query=THP&tab=substance; 

### Free Cofactors:
  - Removed all 81 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 81.
  - Species and properties with assigned default values in debug.log:

  - _CA+2: ['VDW_RAD', 'VDW_EPS']

  - THPBK: ['VDW_RAD', 'VDW_EPS']


### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 141":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   7" to " CB  LEU A   7"

</details>

