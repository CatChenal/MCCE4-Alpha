---
# 1PRW :: Crystal Structure Of Bovine Brain Ca++ Calmodulin In A Compact Form
## PDB.Structure
### Function: METAL BINDING PROTEIN
### First Release: 20-JUN-03
### Method: X-RAY DIFFRACTION
### Resolution: 1.70 ANGSTROMS.
### Molecule: CALMODULIN
### Seqres Species: Residues: A:149; Total: 149
### Cofactors:
  - ACE:
 'ACETYL GROUP', 1
  - M3L:
 'N-TRIMETHYLLYSINE', 1
  -  CA:
 'CALCIUM ION', 4

### Total cofactors: 6
### Total waters: 172
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 11, ARG: 6, ASN: 6, ASP: 17, GLN: 6, GLU: 21, GLY: 11, HIS: 1, ILE: 8, LEU: 9, LYS: 7, MET: 9, PHE: 8, PRO: 2, SER: 4, THR: 12, TYR: 2, VAL: 7', 'Total: 147', 'Ionizable: 54',
              'Ratio: 36.7%')

### Model 1 Free Cofactors & Waters:
  - A:
 'ACE': 1, 'HOH': 172, 'M3L': 1, '_CA': 4

### Links:
  - C  ACE A  0 -- 1.33 Å --> N  ALA A  1
  - C  GLU A 114 -- 1.33 Å --> N  M3L A 115
  - C  M3L A 115 -- 1.33 Å --> N  LEU A 116
  - OD1 ASP A 20 -- 2.24 Å --> CA  CA A 377
  - OD1 ASP A 22 -- 2.43 Å --> CA  CA A 377
  - OD2 ASP A 22 -- 3.40 Å --> CA  CA A 377
  - OD1 ASP A 24 -- 2.30 Å --> CA  CA A 377
  - O  THR A 26 -- 2.42 Å --> CA  CA A 377
  - OE2 GLU A 31 -- 2.49 Å --> CA  CA A 377
  - OE1 GLU A 31 -- 2.50 Å --> CA  CA A 377
  - OD1 ASP A 56 -- 2.25 Å --> CA  CA A 378
  - OD1 ASP A 58 -- 2.38 Å --> CA  CA A 378
  - OD1 ASN A 60 -- 2.20 Å --> CA  CA A 378
  - O  THR A 62 -- 2.37 Å --> CA  CA A 378
  - OE1 GLU A 67 -- 2.43 Å --> CA  CA A 378
  - OE2 GLU A 67 -- 2.54 Å --> CA  CA A 378
  - OD1 ASP A 93 -- 2.31 Å --> CA  CA A 380
  - OD1 ASP A 95 -- 2.36 Å --> CA  CA A 380
  - OD1 ASN A 97 -- 2.41 Å --> CA  CA A 380
  - O  TYR A 99 -- 2.27 Å --> CA  CA A 380
  - OE2 GLU A 104 -- 2.47 Å --> CA  CA A 380
  - OE1 GLU A 104 -- 2.47 Å --> CA  CA A 380
  - OD1 ASP A 129 -- 2.30 Å --> CA  CA A 379
  - OD1 ASP A 131 -- 2.42 Å --> CA  CA A 379
  - OD1 ASP A 133 -- 2.40 Å --> CA  CA A 379
  - O  GLN A 135 -- 2.30 Å --> CA  CA A 379
  - OE1 GLU A 140 -- 2.39 Å --> CA  CA A 379
  - OE2 GLU A 140 -- 2.61 Å --> CA  CA A 379
  - O  HOH A 178 -- 2.43 Å --> CA  CA A 378
  - O  HOH A 200 -- 2.54 Å --> CA  CA A 377
  - O  HOH A 203 -- 2.23 Å --> CA  CA A 380
  - O  HOH A 291 -- 2.39 Å --> CA  CA A 379

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE CA A 377', 'ASP A  20', 'ASP A  22', 'ASP A  24', 'THR A  26', 'GLU A  31', 'HOH A 200']
  - AC2: ['BINDING SITE FOR RESIDUE CA A 378', 'ASP A  56', 'ASP A  58', 'ASN A  60', 'THR A  62', 'GLU A  67', 'HOH A 178']
  - AC3: ['BINDING SITE FOR RESIDUE CA A 379', 'ASP A 129', 'ASP A 131', 'ASP A 133', 'GLN A 135', 'GLU A 140', 'HOH A 291']
  - AC4: ['BINDING SITE FOR RESIDUE CA A 380', 'ASP A  93', 'ASP A  95', 'ASN A  97', 'TYR A  99', 'GLU A 104', 'HOH A 203']

## MCCE.Step1
### Renamed:
  - "CA    CA A 377" to "CA   _CA A 377"
  - "CA    CA A 378" to "CA   _CA A 378"
  - "CA    CA A 379" to "CA   _CA A 379"
  - "CA    CA A 380" to "CA   _CA A 380"

### Termini:
 - <strong>CTR</strong>: "LYS A 148"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
ACE: https://pubchem.ncbi.nlm.nih.gov/#query=ACE&tab=substance; M3L: https://pubchem.ncbi.nlm.nih.gov/#query=M3L&tab=substance; 

### Load Structure:
  -    Warning: load_pdb(): "Duplicate backbone atom ignored, " CH3 ACE A   0".
  -    Warning: load_pdb(): "Duplicate backbone atom ignored, " C   ACE A   0".
  -    Warning: load_pdb(): "Duplicate backbone atom ignored, " O   ACE A   0".

### Free Cofactors:
  - Removed all 172 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 172.
  - Species and properties with assigned default values in debug.log:

  - ACEBK: ['VDW_RAD', 'VDW_EPS']

  - M3LBK: ['VDW_RAD', 'VDW_EPS']

  - _CA+2: ['VDW_RAD', 'VDW_EPS']


### Distance Clashes:
<details><summary>Clashes found</summary>

- No clash found.

</details>

