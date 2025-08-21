---
# 3DMU :: Crystal Structure Of Staphylococcal Nuclease Variant Phs T62K At Cryogenic Temperature;
## PDB.Structure
### Function: HYDROLASE
### First Release: 01-JUL-08
### Method: X-RAY DIFFRACTION
### Resolution: 1.80 ANGSTROMS.
### Molecule: THERMONUCLEASE
### Seqres Species: Residues: A:149; Total: 149
### Cofactors:
  - PO4:
 'PHOSPHATE ION', 1
  - MPD:
 '(4S)-2-METHYL-2,4-PENTANEDIOL', 1

### Total cofactors: 2
### Total waters: 73
### Missing:
  - Residues:
 'A': [('ALA', 1), ('THR', 2), ('SER', 3), ('THR', 4), ('LYS', 5), ('LYS', 6), ('THR', 44), ('LYS', 45), ('HIS', 46), ('PRO', 47), ('LYS', 48), ('LYS', 49), ('GLY', 50), ('VAL', 51), ('GLU', 52), ('GLU', 142), ('ASP', 143), ('ASN', 144), ('ALA', 145), ('ASP', 146), ('SER', 147), ('GLY', 148), ('GLN', 149),
       ('Count', 23)]

### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 13, ARG: 5, ASN: 5, ASP: 6, GLN: 5, GLU: 10, GLY: 9, HIS: 2, ILE: 5, LEU: 12, LYS: 19, MET: 4, PHE: 3, PRO: 4, SER: 2, THR: 6, TRP: 1, TYR: 7, VAL: 8', 'Total: 126', 'Ionizable: 49',
              'Ratio: 38.9%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 73, 'MPD': 1, 'PO4': 1

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE PO4 A 150', 'TYR A  27', 'LYS A  28', 'LYS A 134', 'LYS A 136']
  - AC2: ['BINDING SITE FOR RESIDUE MPD A 151', 'ASP A  21', 'ARG A  35', 'LYS A  70', 'LYS A  71', 'ARG A  87']

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "LEU A   7", "LYS A  53"
 - <strong>CTR</strong>: "GLU A  43", "SER A 141"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
PO4: https://pubchem.ncbi.nlm.nih.gov/#query=PO4&tab=substance; MPD: https://pubchem.ncbi.nlm.nih.gov/#query=MPD&tab=substance; 

### Free Cofactors:
  - Removed all 73 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 73.
  - Species and properties with assigned default values in debug.log:

  - PO4BK: ['VDW_RAD', 'VDW_EPS']

  - MPDBK: ['VDW_RAD', 'VDW_EPS']


### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A  43":   OXT
  -    Missing heavy atoms for CTR01 in "CTR A 141":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   7" to " CB  LEU A   7"
- d= 1.54: " CA  NTR A  53" to " CB  LYS A  53"

</details>

