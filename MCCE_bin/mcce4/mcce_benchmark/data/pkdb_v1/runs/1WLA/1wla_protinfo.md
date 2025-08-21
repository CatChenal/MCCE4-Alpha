---
# 1WLA :: Myoglobin (Horse Heart) Recombinant Wild-Type
## PDB.Structure
### Function: OXYGEN TRANSPORT
### First Release: 24-SEP-97
### Method: X-RAY DIFFRACTION
### Resolution: 1.70 ANGSTROMS.
### Molecule: MYOGLOBIN
### Seqres Species: Residues: A:153; Total: 153
### Cofactors:
  - SO4:
 'SULFATE ION', 1
  - HEM:
 'HEME', 1

### Total cofactors: 2
### Total waters: 74
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 15, ARG: 2, ASN: 2, ASP: 8, GLN: 6, GLU: 13, GLY: 15, HIS: 11, ILE: 9, LEU: 17, LYS: 19, MET: 2, PHE: 7, PRO: 4, SER: 5, THR: 7, TRP: 2, TYR: 2, VAL: 7', 'Total: 153', 'Ionizable: 55',
              'Ratio: 35.9%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HEM': 1, 'HOH': 74, 'PAA': 1, 'PDD': 1, 'SO4': 1

### Links:
  - NE2 HIS A 93 -- 2.12 Å --> FE  HEM A 154
  - FE  HEM A 154 -- 2.17 Å --> O  HOH A 156

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE SO4 A 155', 'ALA A  57', 'SER A  58', 'GLU A  59', 'ASP A  60', 'HOH A 206']
  - AC2: ['BINDING SITE FOR RESIDUE HEM A 154', 'THR A  39', 'LYS A  42', 'PHE A  43', 'LYS A  45', 'VAL A  68', 'SER A  92', 'HIS A  93', 'HIS A  97', 'ILE A  99', 'TYR A 103', 'LEU A 104', 'HIS A 113', 'HIS A 116', 'GLN A 128', 'HOH A 156', 'HOH A 174', 'HOH A 203', 'HOH A 219']

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
 - <strong>NTG</strong>: "GLY A   1"
 - <strong>CTR</strong>: "GLY A 153"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
SO4: https://pubchem.ncbi.nlm.nih.gov/#query=SO4&tab=substance; 

### Free Cofactors:
  - Removed all 74 HOH in A. Removed all 1 SO4 in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 75.
  - Species and properties with assigned default values in debug.log:

  - SO4BK: ['VDW_RAD', 'VDW_EPS']

  - NTG01: ['TORSION']

  - HEM01: ['TORSION']

  - Empty connection slot(s):

  - FE: ['HEM 154 to atom  SD']


### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.55: " C2A HEM A 154" to " CAA PAA A 154"
- d= 1.55: " C3D HEM A 154" to " CAD PDD A 154"

</details>

