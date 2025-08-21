---
# 1DWR :: Myoglobin (Horse Heart) Wild-Type Complexed With Co
## PDB.Structure
### Function: OXYGEN TRANSPORT
### First Release: 11-DEC-99
### Method: X-RAY DIFFRACTION
### Resolution: 1.45 ANGSTROMS.
### Molecule: MYOGLOBIN
### Seqres Species: Residues: A:153; Total: 153
### Cofactors:
  - HEM:
 'HEME', 1
  - CMO:
 'CARBON MONOXIDE', 1
  - SO4:
 'SULFATE ION', 2

### Total cofactors: 4
### Total waters: 132
### Missing:
  - Residues:
 'A': [('GLY', 153), ('Count', 1)]

### Biounits: : A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 15, ARG: 2, ASN: 2, ASP: 8, GLN: 6, GLU: 13, GLY: 14, HIS: 11, ILE: 9, LEU: 17, LYS: 19, MET: 2, PHE: 7, PRO: 4, SER: 5, THR: 7, TRP: 2, TYR: 2, VAL: 7', 'Total: 152', 'Ionizable: 55',
              'Ratio: 36.2%')

### Model 1 Free Cofactors & Waters:
  - A:
 'CMO': 1, 'HEM': 1, 'HOH': 132, 'PAA': 1, 'PDD': 1, 'SO4': 2

### Links:
  - NE2 HIS A 93 -- 2.13 Å --> FE  HEM A 154

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE SO4 A 156', 'SER A  58', 'GLU A  59', 'ASP A  60', 'HOH A2130', 'HOH A2131']
  - AC2: ['BINDING SITE FOR RESIDUE SO4 A 900', 'GLY A   1', 'HOH A2132']
  - AC3: ['BINDING SITE FOR RESIDUE HEM A 154', 'THR A  39', 'LYS A  42', 'PHE A  43', 'LYS A  45', 'HIS A  64', 'VAL A  68', 'LEU A  89', 'SER A  92', 'HIS A  93', 'HIS A  97', 'ILE A  99', 'TYR A 103', 'HIS A 116', 'GLN A 128', 'CMO A 155', 'HOH A2080', 'HOH A2125', 'HOH A2126', 'HOH A2127', 'HOH A2128']
  - AC4: ['BINDING SITE FOR RESIDUE CMO A 155', 'PHE A  43', 'HIS A  64', 'VAL A  68', 'HEM A 154']

## MCCE.Step1
### Renamed:
  - " C   CMO A 155" to " C_  CMO A 155"
  - " CAA HEM A 154" to " CAA PAA A 154"
  - " CAD HEM A 154" to " CAD PDD A 154"
  - " CBA HEM A 154" to " CBA PAA A 154"
  - " CBD HEM A 154" to " CBD PDD A 154"
  - " CGA HEM A 154" to " CGA PAA A 154"
  - " CGD HEM A 154" to " CGD PDD A 154"
  - " O   CMO A 155" to " O_  CMO A 155"
  - " O1A HEM A 154" to " O1A PAA A 154"
  - " O1D HEM A 154" to " O1D PDD A 154"
  - " O2A HEM A 154" to " O2A PAA A 154"
  - " O2D HEM A 154" to " O2D PDD A 154"

### Termini:
 - <strong>NTG</strong>: "GLY A   1"
 - <strong>CTR</strong>: "GLN A 152"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
CMO: https://pubchem.ncbi.nlm.nih.gov/#query=CMO&tab=substance; SO4: https://pubchem.ncbi.nlm.nih.gov/#query=SO4&tab=substance; 

### Free Cofactors:
  - Removed all 132 HOH in A. Removed all 2 SO4 in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 134.
  - Species and properties with assigned default values in debug.log:

  - CMOBK: ['VDW_RAD', 'VDW_EPS']

  - SO4BK: ['VDW_RAD', 'VDW_EPS']

  - NTG01: ['TORSION']

  - HEM01: ['TORSION']

  - Empty connection slot(s):

  - FE: ['HEM 154 to atom  SD']


### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 152":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.51: " C2A HEM A 154" to " CAA PAA A 154"
- d= 1.49: " C3D HEM A 154" to " CAD PDD A 154"
- d= 1.90: "FE   HEM A 154" to " C_  CMO A 155"

</details>

