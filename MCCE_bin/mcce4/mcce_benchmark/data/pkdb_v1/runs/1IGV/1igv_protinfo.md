---
# 1IGV :: Bovine Calbindin D9K Binding Mn2+
## PDB.Structure
### Function: METAL BINDING PROTEIN
### First Release: 18-APR-01
### Method: X-RAY DIFFRACTION
### Resolution: 1.85 ANGSTROMS.
### Molecule: VITAMIN D-DEPENDENT CALCIUM-BINDING PROTEIN, INTESTINAL
### Seqres Species: Residues: A:75; Total: 75
### Cofactors:
  -  MN:
 'MANGANESE (II) ION', 1

### Total cofactors: 1
### Total waters: 39
### Missing:
  - ResAtoms:
 'A': ['GLU_51=(CG,CD,OE1,OE2)']

### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 2, ASN: 2, ASP: 4, GLN: 4, GLU: 13, GLY: 5, ILE: 2, LEU: 12, LYS: 10, PHE: 5, PRO: 4, SER: 6, THR: 2, TYR: 1, VAL: 3', 'Total: 75', 'Ionizable: 28',
              'Ratio: 37.3%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 39, '_MN': 1

### Links:
  - OD1 ASP A 54 -- 2.02 Å --> MN  MN A 76
  - OD1 ASN A 56 -- 2.23 Å --> MN  MN A 76
  - OD1 ASP A 58 -- 2.15 Å --> MN  MN A 76
  - O  GLU A 60 -- 2.00 Å --> MN  MN A 76
  - MN  MN A 76 -- 2.12 Å --> O  HOH A 106
  - MN  MN A 76 -- 2.12 Å --> O  HOH A 108

### Sites:
  - MSE: ['MN BINDING EF-HAND LOOP.', 'ASP A  54', 'ASN A  56', 'ASP A  58', 'GLU A  60']
  - AC1: ['BINDING SITE FOR RESIDUE MN A 76', 'ASP A  54', 'ASN A  56', 'ASP A  58', 'GLU A  60', 'HOH A 106', 'HOH A 108']

## MCCE.Step1
### Renamed:
  - "MN    MN A  76" to "MN   _MN A  76"

### Termini:
 - <strong>NTR</strong>: "LYS A   1"
 - <strong>CTR</strong>: "GLN A  75"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
_MN: https://pubchem.ncbi.nlm.nih.gov/#query=MN&tab=substance; 

### Free Cofactors:
  - Removed all 39 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 39.
  - Species and properties with assigned default values in debug.log:

  - _MNBK: ['VDW_RAD', 'VDW_EPS']


### Missing Heavy Atoms:
  -    Missing heavy atoms for GLU01 in "GLU A  51":   CG ,  CD ,  OE1,  OE2

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.54: " CA  NTR A   1" to " CB  LYS A   1"
- d= 2.00: " O   GLU A  60" to "MN   _MN A  76"

</details>

