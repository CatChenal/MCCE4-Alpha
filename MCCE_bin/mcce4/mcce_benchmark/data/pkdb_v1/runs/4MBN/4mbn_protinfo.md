---
# 4MBN :: Refinement Of Myoglobin And Cytochrome C
## PDB.Structure
### Function: OXYGEN STORAGE
### First Release: 14-JAN-88
### Method: X-RAY DIFFRACTION
### Resolution: 2.00 ANGSTROMS.
### Molecule: MYOGLOBIN
### Seqres Species: Residues: A:153; Total: 153
### Cofactors:
  - SO4:
 'SULFATE ION', 2
  - HEM:
 'HEME', 1

### Total cofactors: 3
### Total waters: 121
### Missing:
  - LigandAtoms: ['A: SO4 216, SO4 224']

### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 17, ARG: 4, ASN: 1, ASP: 7, GLN: 5, GLU: 14, GLY: 11, HIS: 12, ILE: 9, LEU: 18, LYS: 19, MET: 2, PHE: 6, PRO: 4, SER: 6, THR: 5, TRP: 2, TYR: 3, VAL: 8', 'Total: 153', 'Ionizable: 59',
              'Ratio: 38.6%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HEM': 1, 'HOH': 121, 'PAA': 1, 'PDD': 1, 'SO4': 2

### Links:
  - NE2 HIS A 93 -- 2.17 Å --> FE  HEM A 225
  - FE  HEM A 225 -- 2.10 Å --> O  HOH A 226

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE SO4 A 224', 'THR A  67']
  - AC2: ['BINDING SITE FOR RESIDUE HEM A 225', 'THR A  39', 'LYS A  42', 'PHE A  43', 'ARG A  45', 'HIS A  64', 'THR A  67', 'VAL A  68', 'SER A  92', 'HIS A  93', 'HIS A  97', 'ILE A  99', 'TYR A 103', 'GLN A 128', 'HOH A 226', 'HOH A 234', 'HOH A 300', 'HOH A 314', 'HOH A 334']

## MCCE.Step1
### Renamed:
  - " CAA HEM A 225" to " CAA PAA A 225"
  - " CAD HEM A 225" to " CAD PDD A 225"
  - " CBA HEM A 225" to " CBA PAA A 225"
  - " CBD HEM A 225" to " CBD PDD A 225"
  - " CGA HEM A 225" to " CGA PAA A 225"
  - " CGD HEM A 225" to " CGD PDD A 225"
  - " O1A HEM A 225" to " O1A PAA A 225"
  - " O1D HEM A 225" to " O1D PDD A 225"
  - " O2A HEM A 225" to " O2A PAA A 225"
  - " O2D HEM A 225" to " O2D PDD A 225"

### Termini:
 - <strong>NTR</strong>: "VAL A   1"
 - <strong>CTR</strong>: "GLY A 153"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
SO4: https://pubchem.ncbi.nlm.nih.gov/#query=SO4&tab=substance; 

### Free Cofactors:
  - Removed all 121 HOH in A. Removed all 2 SO4 in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 123.
  - Species and properties with assigned default values in debug.log:

  - SO4BK: ['VDW_RAD', 'VDW_EPS']

  - HEM01: ['TORSION']

  - Empty connection slot(s):

  - NE2: ['HIL 97 to atom  ?']

  - FE: ['HEM 225 to atom  SD']


### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   1" to " CB  VAL A   1"
- d= 1.51: " C2A HEM A 225" to " CAA PAA A 225"
- d= 1.52: " C3D HEM A 225" to " CAD PDD A 225"

</details>

