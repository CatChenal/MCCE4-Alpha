---
# 3EVQ :: Crystal Structure Of Staphylococcal Nuclease Variant Delta+Phs L25E At Cryogenic Temperature;
## PDB.Structure
### Function: HYDROLASE
### First Release: 13-OCT-08
### Method: X-RAY DIFFRACTION
### Resolution: 2.00 ANGSTROMS.
### Molecule: THERMONUCLEASE
### Seqres Species: Residues: A:143; Total: 143
### Cofactors:
  - THP:
 "THYMIDINE-3',5'-DIPHOSPHATE", 1
  -  CA:
 'CALCIUM ION', 1

### Total cofactors: 2
### Total waters: 55
### Missing:
  - Residues:
 'A': [('ALA', 1), ('THR', 2), ('SER', 3), ('THR', 4), ('LYS', 5), ('LYS', 6), ('GLU', 142), ('ASP', 143), ('ASN', 144), ('ALA', 145), ('ASP', 146), ('SER', 147), ('GLY', 148), ('GLN', 149),
       ('Count', 14)]

### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 13, ARG: 5, ASN: 6, ASP: 6, GLN: 5, GLU: 12, GLY: 9, HIS: 2, ILE: 5, LEU: 11, LYS: 18, MET: 4, PHE: 4, PRO: 4, SER: 2, THR: 7, TRP: 1, TYR: 7, VAL: 8', 'Total: 129', 'Ionizable: 50',
              'Ratio: 38.8%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 55, 'THP': 1, '_CA': 1

### Links:
  - OD2 ASP A 21 -- 2.95 Å --> CA  CA A 151
  - OD1 ASP A 40 -- 2.79 Å --> CA  CA A 151
  - O  THR A 41 -- 2.79 Å --> CA  CA A 151
  - CA  CA A 151 -- 2.86 Å --> O  HOH A 160
  - CA  CA A 151 -- 2.66 Å --> O  HOH A 163

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE THP A 150', 'ARG A  35', 'ASP A  40', 'LYS A  84', 'TYR A  85', 'ARG A  87', 'LEU A  89', 'TYR A 113', 'TYR A 115', 'LYS A 127', 'CA A 151', 'HOH A 152', 'HOH A 154', 'HOH A 163', 'HOH A 172', 'HOH A 183', 'HOH A 184', 'HOH A 185', 'HOH A 198']
  - AC2: ['BINDING SITE FOR RESIDUE CA A 151', 'ASP A  21', 'ASP A  40', 'THR A  41', 'GLU A  43', 'THP A 150', 'HOH A 160', 'HOH A 163']

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
  - Removed all 55 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 55.
  - Species and properties with assigned default values in debug.log:

  - THPBK: ['VDW_RAD', 'VDW_EPS']

  - _CA+2: ['VDW_RAD', 'VDW_EPS']


### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 141":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.54: " CA  NTR A   7" to " CB  LEU A   7"

</details>

