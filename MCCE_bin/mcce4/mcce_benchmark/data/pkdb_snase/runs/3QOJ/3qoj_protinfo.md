---
# 3QOJ :: Cryogenic Structure Of Staphylococcal Nuclease Variant D+Phs/V23K
## PDB.Structure
### Function: HYDROLASE
### First Release: 10-FEB-11
### Method: X-RAY DIFFRACTION
### Resolution: 1.60 ANGSTROMS.
### Molecule: THERMONUCLEASE
### Seqres Species: Residues: A:143; Total: 143
### Cofactors:
  -  CA:
 'CALCIUM ION', 1
  - THP:
 "THYMIDINE-3',5'-DIPHOSPHATE", 1

### Total cofactors: 2
### Total waters: 111
### Missing:
  - Residues:
 'A': [('ALA', 1), ('THR', 2), ('SER', 3), ('THR', 4), ('LYS', 5), ('LYS', 6), ('GLU', 142), ('ASP', 143), ('ASN', 144), ('ALA', 145), ('ASP', 146), ('SER', 147), ('GLY', 148), ('GLN', 149),
       ('Count', 14)]

### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 13, ARG: 5, ASN: 6, ASP: 6, GLN: 5, GLU: 11, GLY: 9, HIS: 2, ILE: 5, LEU: 12, LYS: 19, MET: 4, PHE: 4, PRO: 4, SER: 2, THR: 7, TRP: 1, TYR: 7, VAL: 7', 'Total: 129', 'Ionizable: 50',
              'Ratio: 38.8%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 111, 'THP': 1, '_CA': 1

### Links:
  - OD2 ASP A 21 -- 2.86 Å --> CA  CA A 150
  - OD1 ASP A 40 -- 2.66 Å --> CA  CA A 150
  - O  THR A 41 -- 2.88 Å --> CA  CA A 150
  - OE1 GLU A 43 -- 3.01 Å --> CA  CA A 150
  - O  HOH A 45 -- 2.90 Å --> CA  CA A 150
  - CA  CA A 150 -- 3.17 Å --> O5P THP A 151
  - CA  CA A 150 -- 2.78 Å --> O  HOH A 181

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE CA A 150', 'ASP A  21', 'ASP A  40', 'THR A  41', 'GLU A  43', 'HOH A  45', 'THP A 151', 'HOH A 181']
  - AC2: ['BINDING SITE FOR RESIDUE THP A 151', 'ARG A  35', 'ASP A  40', 'HOH A  48', 'LYS A  84', 'TYR A  85', 'ARG A  87', 'TYR A 113', 'TYR A 115', 'LYS A 127', 'CA A 150', 'HOH A 152', 'HOH A 156', 'HOH A 173', 'HOH A 181', 'HOH A 192', 'HOH A 215', 'HOH A 220', 'HOH A 227', 'HOH A 229', 'HOH A 247']

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
  - Removed all 111 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 111.
  - Species and properties with assigned default values in debug.log:

  - _CA+2: ['VDW_RAD', 'VDW_EPS']

  - THPBK: ['VDW_RAD', 'VDW_EPS']


### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 141":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.54: " CA  NTR A   7" to " CB  LEU A   7"

</details>

