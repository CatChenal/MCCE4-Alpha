---
# 1MBC :: X-Ray Structure And Refinement Of Carbon-Monoxy (Fe Ii)-Myoglobin At 1.5 Angstroms Resolution;
## PDB.Structure
### Function: OXYGEN STORAGE
### First Release: 15-SEP-88
### Method: X-RAY DIFFRACTION
### Resolution: 1.50 ANGSTROMS.
### Molecule: MYOGLOBIN
### Seqres Species: Residues: A:153; Total: 153
### Cofactors:
  - SO4:
 'SULFATE ION', 1
  - HEM:
 'HEME', 1
  - CMO:
 'CARBON MONOXIDE', 1

### Total cofactors: 3
### Total waters: 137
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 17, ARG: 4, ASN: 1, ASP: 7, GLN: 5, GLU: 14, GLY: 11, HIS: 12, ILE: 9, LEU: 18, LYS: 19, MET: 2, PHE: 6, PRO: 4, SER: 6, THR: 5, TRP: 2, TYR: 3, VAL: 8', 'Total: 153', 'Ionizable: 59',
              'Ratio: 38.6%')

### Model 1 Free Cofactors & Waters:
  - A:
 'CMO': 1, 'HEM': 1, 'HOH': 137, 'PAA': 1, 'PDD': 1, 'SO4': 1

### Links:
  - NE2 HIS A 93 -- 2.19 Å --> FE  HEM A 155
  - FE  HEM A 155 -- 1.92 Å --> C  CMO A 201
  - FE  HEM A 155 -- 2.93 Å --> O ACMO A 201
  - FE  HEM A 155 -- 2.73 Å --> O BCMO A 201

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE SO4 A 154', 'SER A  58', 'GLU A  59', 'ASP A  60', 'HOH A 311']
  - AC2: ['BINDING SITE FOR RESIDUE HEM A 155', 'THR A  39', 'LYS A  42', 'PHE A  43', 'ARG A  45', 'HIS A  64', 'THR A  67', 'VAL A  68', 'ALA A  71', 'SER A  92', 'HIS A  93', 'HIS A  97', 'ILE A  99', 'TYR A 103', 'ILE A 107', 'CMO A 201', 'HOH A 306', 'HOH A 314', 'HOH A 321']
  - AC3: ['BINDING SITE FOR RESIDUE CMO A 201', 'PHE A  43', 'HIS A  64', 'VAL A  68', 'HEM A 155']

## MCCE.Step1
### Renamed:
  - " C   CMO A 201" to " C_  CMO A 201"
  - " CAA HEM A 155" to " CAA PAA A 155"
  - " CAD HEM A 155" to " CAD PDD A 155"
  - " CBA HEM A 155" to " CBA PAA A 155"
  - " CBD HEM A 155" to " CBD PDD A 155"
  - " CGA HEM A 155" to " CGA PAA A 155"
  - " CGD HEM A 155" to " CGD PDD A 155"
  - " O   CMO A 201" to " O_  CMO A 201"
  - " O1A HEM A 155" to " O1A PAA A 155"
  - " O1D HEM A 155" to " O1D PDD A 155"
  - " O2A HEM A 155" to " O2A PAA A 155"
  - " O2D HEM A 155" to " O2D PDD A 155"

### Termini:
 - <strong>NTR</strong>: "VAL A   1"
 - <strong>CTR</strong>: "GLY A 153"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
SO4: https://pubchem.ncbi.nlm.nih.gov/#query=SO4&tab=substance; CMO: https://pubchem.ncbi.nlm.nih.gov/#query=CMO&tab=substance; 

### Free Cofactors:
  - Removed all 137 HOH in A. Removed all 1 SO4 in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 138.
  - Species and properties with assigned default values in debug.log:

  - SO4BK: ['VDW_RAD', 'VDW_EPS']

  - CMOBK: ['VDW_RAD', 'VDW_EPS']

  - HEM01: ['TORSION']

  - Empty connection slot(s):

  - NE2: ['HIL 97 to atom  ?']

  - FE: ['HEM 155 to atom  SD']


### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.55: " CA  NTR A   1" to " CB  VAL A   1"
- d= 1.51: " C2A HEM A 155" to " CAA PAA A 155"
- d= 1.54: " C3D HEM A 155" to " CAD PDD A 155"
- d= 1.92: "FE   HEM A 155" to " C_  CMO A 201"

</details>

