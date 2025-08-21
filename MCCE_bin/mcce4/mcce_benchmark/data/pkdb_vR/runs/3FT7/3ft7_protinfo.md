---
# 3FT7 :: Crystal Structure Of An Extremely Stable Dimeric Protein From Sulfolobus Islandicus;
## PDB.Structure
### Function: DNA BINDING PROTEIN
### First Release: 12-JAN-09
### Method: X-RAY DIFFRACTION
### Resolution: 2.00 ANGSTROMS.
### Molecule: UNCHARACTERIZED PROTEIN ORF56
### Seqres Species: Residues: A:55, B:55; Total: 110
### Cofactors:
  - GOL:
 'GLYCERIN; PROPANE-1,2,3-TRIOL', 1

### Total cofactors: 1
### Total waters: 41
### Missing:
  - Residues:
 'A': [('GLY', 2), ('ARG', 3), ('PRO', 4), ('TYR', 5), ('LYS', 6), ('LEU', 7), ('LYS', 53), ('GLN', 54), ('LYS', 55), ('LYS', 56), ('Count', 10)],
 'B': [('GLY', 2), ('ARG', 3), ('PRO', 4), ('TYR', 5), ('LYS', 6), ('LYS', 53), ('GLN', 54), ('LYS', 55), ('LYS', 56), ('Count', 9)]
  - ResAtoms:
 'A': ['ASP_52=(CG,OD1,OD2)'], 'B': ['ASP_52=(CG,OD1,OD2)']

### Biounits: DIMERIC: A, B;
### Models: 1
### Chains: A, B
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 2, ARG: 2, ASN: 3, ASP: 4, CYS: 1, GLN: 1, GLU: 4, GLY: 2, HIS: 2, ILE: 4, LEU: 7, LYS: 4, MET: 1, PRO: 1, SER: 1, THR: 1, TRP: 1, TYR: 2, VAL: 2', 'Total: 45', 'Ionizable: 19',
              'Ratio: 42.2%')
  - B:
 'RESIDUES': ('ALA: 2, ARG: 2, ASN: 3, ASP: 4, CYS: 1, GLN: 1, GLU: 4, GLY: 2, HIS: 2, ILE: 4, LEU: 8, LYS: 4, MET: 1, PRO: 1, SER: 1, THR: 1, TRP: 1, TYR: 2, VAL: 2', 'Total: 46', 'Ionizable: 19',
              'Ratio: 41.3%')

### Model 1 Free Cofactors & Waters:
  - A:
 'GOL': 1, 'HOH': 24
  - B:
 'HOH': 17

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE GOL A 57', 'ASP A  38', 'ARG A  41']

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "LEU A   8", "LEU B   7"
 - <strong>CTR</strong>: "ASP A  52", "ASP B  52"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
GOL: https://pubchem.ncbi.nlm.nih.gov/#query=GOL&tab=substance; 

### Free Cofactors:
  - Removed all 24 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 41.
  - Species and properties with assigned default values in debug.log:

  - GOLBK: ['VDW_RAD', 'VDW_EPS']


### Missing Heavy Atoms:
  -    Missing heavy atoms for ASP01 in "ASP A  52":   CG ,  OD1,  OD2
  -    Missing heavy atoms for CTR01 in "CTR A  52":   OXT
  -    Missing heavy atoms for ASP01 in "ASP B  52":   CG ,  OD1,  OD2
  -    Missing heavy atoms for CTR01 in "CTR B  52":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.54: " CA  NTR A   8" to " CB  LEU A   8"
- d= 1.52: " CA  NTR B   7" to " CB  LEU B   7"

</details>

