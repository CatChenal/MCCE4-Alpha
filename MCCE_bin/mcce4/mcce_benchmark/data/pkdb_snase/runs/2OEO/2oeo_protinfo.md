---
# 2OEO :: Cryogenic Crystal Structure Of Staphylococcal Nuclease Variant Truncated Delta+Phs I92D;
## PDB.Structure
### Function: HYDROLASE
### First Release: 30-DEC-06
### Method: X-RAY DIFFRACTION
### Resolution: 2.00 ANGSTROMS.
### Molecule: STAPHYLOCOCCAL THERMONUCLEASE
### Seqres Species: Residues: A:138; Total: 138
### Cofactors:
  -  CA:
 'CALCIUM ION', 1
  - THP:
 "THYMIDINE-3',5'-DIPHOSPHATE", 1

### Total cofactors: 2
### Total waters: 59
### Missing:
  - Residues:
 'A': [('ALA', 1), ('THR', 2), ('SER', 3), ('THR', 4), ('LYS', 5), ('GLU', 142), ('ASP', 143), ('ASN', 144),
       ('Count', 8)]
  - ResAtoms:
 'A': ['LYS_6=(CG,CD,CE,NZ)']

### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 17, ARG: 5, ASN: 6, ASP: 7, GLN: 5, GLU: 11, GLY: 9, HIS: 2, ILE: 4, LEU: 12, LYS: 14, MET: 4, PHE: 4, PRO: 4, SER: 3, THR: 7, TRP: 1, TYR: 7, VAL: 8', 'Total: 130', 'Ionizable: 46',
              'Ratio: 35.4%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 59, 'THP': 1, '_CA': 1

### Links:
  - OD2 ASP A 21 -- 2.78 Å --> CA  CA A 301
  - OD1 ASP A 40 -- 2.64 Å --> CA  CA A 301
  - O  THR A 41 -- 2.67 Å --> CA  CA A 301
  - O  HOH A 235 -- 2.93 Å --> CA  CA A 301
  - O  HOH A 236 -- 2.88 Å --> CA  CA A 301

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE CA A 301', 'ASP A  21', 'ASP A  40', 'THR A  41', 'HOH A 235', 'HOH A 236']
  - AC2: ['BINDING SITE FOR RESIDUE THP A 300', 'ARG A  35', 'LEU A  36', 'LEU A  37', 'ASP A  40', 'ASP A  83', 'LYS A  84', 'TYR A  85', 'ARG A  87', 'LEU A  89', 'TYR A 113', 'TYR A 115', 'HOH A 208', 'HOH A 210', 'HOH A 212', 'HOH A 235', 'HOH A 239', 'HOH A 241']

## MCCE.Step1
### Renamed:
  - "CA    CA A 301" to "CA   _CA A 301"

### Termini:
 - <strong>NTR</strong>: "LYS A   6"
 - <strong>CTR</strong>: "SER A 141"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
THP: https://pubchem.ncbi.nlm.nih.gov/#query=THP&tab=substance; 

### Free Cofactors:
  - Removed all 59 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 59.
  - Species and properties with assigned default values in debug.log:

  - _CA+2: ['VDW_RAD', 'VDW_EPS']

  - THPBK: ['VDW_RAD', 'VDW_EPS']


### Missing Heavy Atoms:
  -    Missing heavy atoms for LYS01 in "LYS A   6":   CG ,  CD ,  CE ,  NZ 
  -    Missing heavy atoms for CTR01 in "CTR A 141":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.52: " CA  NTR A   6" to " CB  LYS A   6"

</details>

