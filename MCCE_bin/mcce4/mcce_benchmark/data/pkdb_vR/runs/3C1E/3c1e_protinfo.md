---
# 3C1E :: Crystal Structure Of Staphylococcal Nuclease Variant Delta+Phs L125K At Cryogenic Temperature;
## PDB.Structure
### Function: HYDROLASE
### First Release: 23-JAN-08
### Method: X-RAY DIFFRACTION
### Resolution: 1.90 ANGSTROMS.
### Molecule: THERMONUCLEASE
### Seqres Species: Residues: A:143; Total: 143
### Cofactors:
  -  CA:
 'CALCIUM ION', 1
  - THP:
 "THYMIDINE-3',5'-DIPHOSPHATE", 1

### Total cofactors: 2
### Total waters: 103
### Missing:
  - Residues:
 'A': [('ALA', 1), ('THR', 2), ('SER', 3), ('THR', 4), ('LYS', 5), ('ASP', 143), ('ASN', 144), ('ALA', 145), ('ASP', 146), ('SER', 147), ('GLY', 148), ('GLN', 149),
       ('Count', 12)]

### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 13, ARG: 5, ASN: 6, ASP: 6, GLN: 5, GLU: 12, GLY: 9, HIS: 2, ILE: 5, LEU: 11, LYS: 20, MET: 4, PHE: 4, PRO: 4, SER: 2, THR: 7, TRP: 1, TYR: 7, VAL: 8', 'Total: 131', 'Ionizable: 52',
              'Ratio: 39.7%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 103, 'THP': 1, '_CA': 1

### Links:
  - OD2 ASP A 21 -- 2.73 Å --> CA  CA A 150
  - OD1 ASP A 40 -- 2.72 Å --> CA  CA A 150
  - O  THR A 41 -- 2.83 Å --> CA  CA A 150
  - CA  CA A 150 -- 2.97 Å --> O  HOH A 157
  - CA  CA A 150 -- 2.73 Å --> O  HOH A 236

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE CA A 150', 'ASP A  21', 'ASP A  40', 'THR A  41', 'GLU A  43', 'HOH A 157', 'HOH A 236']
  - AC2: ['BINDING SITE FOR RESIDUE THP A 151', 'ARG A  35', 'LEU A  37', 'ASP A  40', 'LYS A  84', 'TYR A  85', 'ARG A  87', 'LEU A  89', 'TYR A 113', 'TYR A 115', 'LYS A 127', 'HOH A 157', 'HOH A 162', 'HOH A 184', 'HOH A 198', 'HOH A 229', 'HOH A 234', 'HOH A 239', 'HOH A 240', 'HOH A 243', 'HOH A 254']

## MCCE.Step1
### Renamed:
  - "CA    CA A 150" to "CA   _CA A 150"

### Termini:
 - <strong>NTR</strong>: "LYS A   6"
 - <strong>CTR</strong>: "GLU A 142"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
THP: https://pubchem.ncbi.nlm.nih.gov/#query=THP&tab=substance; 

### Free Cofactors:
  - Removed all 103 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 103.
  - Species and properties with assigned default values in debug.log:

  - _CA+2: ['VDW_RAD', 'VDW_EPS']

  - THPBK: ['VDW_RAD', 'VDW_EPS']


### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 142":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.55: " CA  NTR A   6" to " CB  LYS A   6"

</details>

