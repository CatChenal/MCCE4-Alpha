---
# 1A6K :: Aquomet-Myoglobin, Atomic Resolution
## PDB.Structure
### Function: HEME PROTEIN
### First Release: 26-FEB-98
### Method: X-RAY DIFFRACTION
### Resolution: 1.10 ANGSTROMS.
### Molecule: MYOGLOBIN
### Seqres Species: Residues: A:151; Total: 151
### Cofactors:
  - SO4:
 'SULFATE ION', 3
  - HEM:
 'HEME', 1

### Total cofactors: 4
### Total waters: 185
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 17, ARG: 4, ASN: 1, ASP: 7, GLN: 4, GLU: 14, GLY: 10, HIS: 12, ILE: 9, LEU: 18, LYS: 19, MET: 2, PHE: 6, PRO: 4, SER: 6, THR: 5, TRP: 2, TYR: 3, VAL: 8', 'Total: 151', 'Ionizable: 59',
              'Ratio: 39.1%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HEM': 1, 'HOH': 185, 'PAA': 1, 'PDD': 1, 'SO4': 3

### Links:
  - NE2 HIS A 93 -- 2.14 Å --> FE  HEM A 154
  - FE  HEM A 154 -- 2.13 Å --> O  HOH A1001

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE SO4 A 155', 'ALA A  57', 'SER A  58', 'GLU A  59', 'ASP A  60', 'HOH A1007', 'HOH A1026', 'HOH A1038', 'HOH A1086']
  - AC2: ['BINDING SITE FOR RESIDUE SO4 A 156', 'GLN A  26', 'LYS A  62', 'HOH A1043', 'HOH A1045', 'HOH A1063']
  - AC3: ['BINDING SITE FOR RESIDUE SO4 A 1188', 'ARG A  45', 'HIS A  64', 'THR A  67', 'HIS A 116', 'HOH A1088', 'HOH A1123']
  - AC4: ['BINDING SITE FOR RESIDUE HEM A 154', 'THR A  39', 'LYS A  42', 'PHE A  43', 'ARG A  45', 'HIS A  64', 'THR A  67', 'VAL A  68', 'LEU A  89', 'SER A  92', 'HIS A  93', 'HIS A  97', 'ILE A  99', 'TYR A 103', 'LEU A 104', 'PHE A 138', 'HOH A1001', 'HOH A1024', 'HOH A1033', 'HOH A1059', 'HOH A1088', 'HOH A1092', 'HOH A1118', 'HOH A1149']

## MCCE.Step1
### Renamed:
  - " CAA HEM A 154" to " CAA PAA A 154"
  - " CAD HEM A 154" to " CAD PDD A 154"
  - " CBA HEM A 154" to " CBA PAA A 154"
  - " CBD HEM A 154" to " CBD PDD A 154"
  - " CGA HEM A 154" to " CGA PAA A 154"
  - " CGD HEM A 154" to " CGD PDD A 154"
  - " O1A HEM A 154" to " O1A PAA A 154"
  - " O1D HEM A 154" to " O1D PDD A 154"
  - " O2A HEM A 154" to " O2A PAA A 154"
  - " O2D HEM A 154" to " O2D PDD A 154"

### Termini:
 - <strong>NTR</strong>: "VAL A   1"
 - <strong>CTR</strong>: "TYR A 151"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
SO4: https://pubchem.ncbi.nlm.nih.gov/#query=SO4&tab=substance; 

### Free Cofactors:
  - Removed all 185 HOH in A. Removed all 3 SO4 in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 187.
  - Species and properties with assigned default values in debug.log:

  - SO4BK: ['VDW_RAD', 'VDW_EPS']

  - HEM01: ['TORSION']

  - Empty connection slot(s):

  - FE: ['HEM 154 to atom  SD']


### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   1" to " CB  VAL A   1"
- d= 1.53: " C2A HEM A 154" to " CAA PAA A 154"
- d= 1.54: " C3D HEM A 154" to " CAD PDD A 154"

</details>

