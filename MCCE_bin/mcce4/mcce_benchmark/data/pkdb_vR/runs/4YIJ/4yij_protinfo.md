---
# 4YIJ :: Crystal Structure Of Staphylcoccal Nuclease Variant Delta+Phs A109E At Cryogenic Temperature;
## PDB.Structure
### Function: HYDROLASE
### First Release: 02-MAR-15
### Method: X-RAY DIFFRACTION
### Resolution: 1.64 ANGSTROMS.
### Molecule: NUCLEASE A
### Seqres Species: Residues: A:143, B:143; Total: 286
### Cofactors:
  -  CA:
 'CALCIUM ION', 2
  - THP:
 "THYMIDINE-3',5'-DIPHOSPHATE", 2

### Total cofactors: 4
### Total waters: 276
### Missing:
  - Residues:
 'A': [('ALA', 1), ('THR', 2), ('SER', 3), ('THR', 4), ('GLU', 142), ('ASP', 143), ('ASN', 144), ('ALA', 145), ('ASP', 146), ('SER', 147), ('GLY', 148), ('GLN', 149), ('Count', 12)],
 'B': [('ASP', 143), ('ASN', 144), ('ALA', 145), ('ASP', 146), ('SER', 147), ('GLY', 148), ('GLN', 149), ('Count', 7)]

### Biounits: MONOMERIC: A,B;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 12, ARG: 5, ASN: 6, ASP: 6, GLN: 5, GLU: 12, GLY: 9, HIS: 2, ILE: 5, LEU: 12, LYS: 20, MET: 4, PHE: 4, PRO: 4, SER: 2, THR: 7, TRP: 1, TYR: 7, VAL: 8', 'Total: 131', 'Ionizable: 52',
              'Ratio: 39.7%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 138, 'THP': 1, '_CA': 1

### Links:
  - OD1 ASP A 19 -- 2.87 Å --> CA  CA A 201
  - OD2 ASP A 19 -- 2.36 Å --> CA  CA A 201
  - OD1 ASP A 21 -- 2.35 Å --> CA  CA A 201
  - CA  CA A 201 -- 2.36 Å --> O  HOH A 350
  - CA  CA A 201 -- 2.40 Å --> O  HOH A 354
  - CA  CA A 201 -- 2.35 Å --> O  HOH A 365
  - CA  CA A 201 -- 2.40 Å --> O  HOH A 378
  - OD1 ASP B 19 -- 2.83 Å --> CA  CA B 201
  - OD2 ASP B 19 -- 2.36 Å --> CA  CA B 201
  - OD1 ASP B 21 -- 2.36 Å --> CA  CA B 201
  - CA  CA B 201 -- 2.31 Å --> O  HOH B 312
  - CA  CA B 201 -- 2.36 Å --> O  HOH B 344
  - CA  CA B 201 -- 2.41 Å --> O  HOH B 353
  - CA  CA B 201 -- 2.40 Å --> O  HOH B 360

### Sites:
  - AC1: ['binding site for residue CA A 201', 'ASP A  19', 'ASP A  21', 'HOH A 350', 'HOH A 354', 'HOH A 365', 'HOH A 378']
  - AC2: ['binding site for residue THP A 202', 'ARG A  35', 'ASP A  40', 'LYS A  84', 'TYR A  85', 'ARG A  87', 'LEU A  89', 'TYR A 113', 'HOH A 310', 'HOH A 315', 'HOH A 350', 'HOH A 359', 'HOH A 363', 'HOH A 364', 'HOH A 365', 'HOH A 366', 'HOH A 385', 'HOH A 405', 'HOH A 419', 'ARG B 105']
  - AC3: ['binding site for residue CA B 201', 'ASP B  19', 'ASP B  21', 'HOH B 312', 'HOH B 344', 'HOH B 353', 'HOH B 360']
  - AC4: ['binding site for residue THP B 202', 'ARG A 105', 'ALA B   1', 'ARG B  35', 'ASP B  40', 'LYS B  84', 'TYR B  85', 'ARG B  87', 'LEU B  89', 'TYR B 113', 'TYR B 115', 'HOH B 306', 'HOH B 307', 'HOH B 312', 'HOH B 315', 'HOH B 344', 'HOH B 365', 'HOH B 366', 'HOH B 377', 'HOH B 392']

## MCCE.Step1
### Renamed:
  - "CA    CA A 201" to "CA   _CA A 201"

### Termini:
 - <strong>NTR</strong>: "LYS A   5"
 - <strong>CTR</strong>: "SER A 141"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
THP: https://pubchem.ncbi.nlm.nih.gov/#query=THP&tab=substance; 

### Free Cofactors:
  - Removed all 138 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 138.
  - Species and properties with assigned default values in debug.log:

  - _CA+2: ['VDW_RAD', 'VDW_EPS']

  - THPBK: ['VDW_RAD', 'VDW_EPS']


### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 141":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.52: " CA  NTR A   5" to " CB  LYS A   5"

</details>

