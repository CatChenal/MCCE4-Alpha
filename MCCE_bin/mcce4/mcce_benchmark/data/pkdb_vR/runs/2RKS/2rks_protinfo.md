---
# 2RKS :: Crystal Structure Of Staphylococcal Nuclease Variant Phs L38K At Cryogenic Temperature;
## PDB.Structure
### Function: HYDROLASE
### First Release: 17-OCT-07
### Method: X-RAY DIFFRACTION
### Resolution: 2.01 ANGSTROMS.
### Molecule: THERMONUCLEASE
### Seqres Species: Residues: A:149; Total: 149
### Cofactors:
  - PO4:
 'PHOSPHATE ION', 2

### Total cofactors: 2
### Total waters: 86
### Missing:
  - Residues:
 'A': [('ALA', 1), ('THR', 2), ('SER', 3), ('THR', 4), ('LYS', 5), ('LYS', 6), ('GLU', 142), ('ASP', 143), ('ASN', 144), ('ALA', 145), ('ASP', 146), ('SER', 147), ('GLY', 148), ('GLN', 149),
       ('Count', 14)]

### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 13, ARG: 5, ASN: 5, ASP: 6, GLN: 5, GLU: 11, GLY: 10, HIS: 3, ILE: 5, LEU: 11, LYS: 22, MET: 4, PHE: 3, PRO: 5, SER: 2, THR: 8, TRP: 1, TYR: 7, VAL: 9', 'Total: 135', 'Ionizable: 54',
              'Ratio: 40.0%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 86, 'PO4': 2

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE PO4 A 251', 'LYS A  71', 'LYS A  84', 'TYR A  85']
  - AC2: ['BINDING SITE FOR RESIDUE PO4 A 252', 'ARG A  35', 'LYS A  71', 'ARG A  87', 'HOH A 273', 'HOH A 277', 'HOH A 330', 'HOH A 334', 'HOH A 335']

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "LEU A   7"
 - <strong>CTR</strong>: "SER A 141"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
PO4: https://pubchem.ncbi.nlm.nih.gov/#query=PO4&tab=substance; 

### Free Cofactors:
  - Removed all 86 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 86.
  - Species and properties with assigned default values in debug.log:

  - PO4BK: ['VDW_RAD', 'VDW_EPS']


### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 141":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.54: " CA  NTR A   7" to " CB  LEU A   7"

</details>

