---
# 3ICB :: The Refined Structure Of Vitamin D-Dependent Calcium-Binding Protein From Bovine Intestine. Molecular Details, Ion Binding, And; Implications For The Structure Of Other Calcium-Binding Proteins;
## PDB.Structure
### Function: CALCIUM-BINDING PROTEIN
### First Release: 09-SEP-86
### Method: X-RAY DIFFRACTION
### Resolution: 2.30 ANGSTROMS.
### Molecule: CALCIUM-BINDING PROTEIN
### Seqres Species: Residues: A:75; Total: 75
### Cofactors:
  -  CA:
 'CALCIUM ION', 2
  - SO4:
 'SULFATE ION', 1

### Total cofactors: 3
### Total waters: 36
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 2, ASN: 2, ASP: 4, GLN: 4, GLU: 13, GLY: 5, ILE: 2, LEU: 12, LYS: 10, PHE: 5, PRO: 4, SER: 6, THR: 2, TYR: 1, VAL: 3', 'Total: 75', 'Ionizable: 28',
              'Ratio: 37.3%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 36, 'SO4': 1, '_CA': 2

### Links:
  - O  ALA A 14 -- 2.16 Å --> CA  CA A 76
  - O  GLU A 17 -- 2.39 Å --> CA  CA A 76
  - O  ASP A 19 -- 2.38 Å --> CA  CA A 76
  - O  GLN A 22 -- 2.22 Å --> CA  CA A 76
  - OE1 GLU A 27 -- 2.51 Å --> CA  CA A 76
  - OE2 GLU A 27 -- 2.27 Å --> CA  CA A 76
  - OD1 ASP A 54 -- 2.38 Å --> CA  CA A 77
  - OD1 ASN A 56 -- 2.45 Å --> CA  CA A 77
  - OD1 ASP A 58 -- 2.29 Å --> CA  CA A 77
  - O  GLU A 60 -- 2.26 Å --> CA  CA A 77
  - OE1 GLU A 65 -- 2.29 Å --> CA  CA A 77
  - OE2 GLU A 65 -- 2.61 Å --> CA  CA A 77
  - CA  CA A 76 -- 2.54 Å --> O  HOH A 80
  - CA  CA A 77 -- 2.25 Å --> O  HOH A 79

### Sites:
  - CB1: ['NULL', 'ALA A  14', 'ALA A  15', 'LYS A  16', 'GLU A  17', 'GLY A  18', 'ASP A  19', 'PRO A  20', 'ASN A  21', 'GLN A  22', 'LEU A  23', 'SER A  24', 'LYS A  25', 'GLU A  26', 'GLU A  27']
  - CB2: ['NULL', 'ASP A  54', 'LYS A  55', 'ASN A  56', 'GLY A  57', 'ASP A  58', 'GLY A  59', 'GLU A  60', 'VAL A  61', 'SER A  62', 'PHE A  63', 'GLU A  64', 'GLU A  65']
  - AC1: ['BINDING SITE FOR RESIDUE CA A 76', 'ALA A  14', 'GLU A  17', 'ASP A  19', 'GLN A  22', 'GLU A  27', 'HOH A  80']
  - AC2: ['BINDING SITE FOR RESIDUE CA A 77', 'ASP A  54', 'ASN A  56', 'ASP A  58', 'GLU A  60', 'GLU A  65', 'HOH A  79']
  - AC3: ['BINDING SITE FOR RESIDUE SO4 A 78', 'GLU A  35', 'PHE A  36', 'PRO A  37', 'SER A  38', 'HOH A  84']

## MCCE.Step1
### Renamed:
  - "CA    CA A  76" to "CA   _CA A  76"
  - "CA    CA A  77" to "CA   _CA A  77"

### Termini:
 - <strong>NTR</strong>: "LYS A   1"
 - <strong>CTR</strong>: "GLN A  75"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
SO4: https://pubchem.ncbi.nlm.nih.gov/#query=SO4&tab=substance; 

### Free Cofactors:
  - Removed all 36 HOH in A. Removed all 1 SO4 in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 37.
  - Species and properties with assigned default values in debug.log:

  - _CA+2: ['VDW_RAD', 'VDW_EPS']

  - SO4BK: ['VDW_RAD', 'VDW_EPS']


### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   1" to " CB  LYS A   1"

</details>

