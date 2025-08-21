---
# 1U2U :: Nmr Solution Structure Of A Designed Heterodimeric Leucine Zipper
## PDB.Structure
### Function: TRANSCRIPTION
### First Release: 20-JUL-04
### Method: SOLUTION NMR
### Resolution: NOT APPLICABLE.
### Molecule: GENERAL CONTROL PROTEIN GCN4
### Seqres Species: Residues: A:32, B:32; Total: 64
### Cofactors:
  - ACE:
 'ACETYL GROUP', 2
  - NH2:
 'AMINO GROUP', 2

### Total cofactors: 4
### Biounits: DIMERIC: A, B;
### Models: 1
### Chains: A, B
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 4, ASN: 1, GLN: 5, GLU: 9, GLY: 1, HIS: 1, LEU: 4, LYS: 1, TYR: 1, VAL: 3', 'Total: 30', 'Ionizable: 12',
              'Ratio: 40.0%')
  - B:
 'RESIDUES': ('ALA: 5, ARG: 3, ASN: 1, GLN: 4, GLU: 1, GLY: 1, HIS: 1, LEU: 4, LYS: 6, TYR: 1, VAL: 3', 'Total: 30', 'Ionizable: 12',
              'Ratio: 40.0%')

### Model 1 Free Cofactors & Waters:
  - A:
 'ACE': 1, 'NH2': 1
  - B:
 'ACE': 1, 'NH2': 1

### Links:
  - C  ACE A  0 -- 1.33 Å --> N  GLU A  1
  - C  GLY A 30 -- 1.33 Å --> N  NH2 A 31
  - C  ACE B  0 -- 1.34 Å --> N  GLU B  1
  - C  GLY B 30 -- 1.33 Å --> N  NH2 B 31

### Sites:
  - AC3: ['BINDING SITE FOR RESIDUE NH2 A 31', 'GLU A  29', 'GLY A  30']
  - AC4: ['BINDING SITE FOR RESIDUE NH2 B 31', 'GLY B  30']

## MCCE.Step1
### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
ACE: https://pubchem.ncbi.nlm.nih.gov/#query=ACE&tab=substance; NH2: https://pubchem.ncbi.nlm.nih.gov/#query=NH2&tab=substance; 

### Load Structure:
  -    Warning: load_pdb(): "Duplicate backbone atom ignored, " CH3 ACE A   0".
  -    Warning: load_pdb(): "Duplicate backbone atom ignored, " C   ACE A   0".
  -    Warning: load_pdb(): "Duplicate backbone atom ignored, " O   ACE A   0".

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Species and properties with assigned default values in debug.log:

  - ACEBK: ['VDW_RAD', 'VDW_EPS']

  - NH2BK: ['VDW_RAD', 'VDW_EPS']


### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 2.00: "HG21 VAL A   9" to "HG21 VAL B   9"
- d= 1.97: "HD11 LEU A  12" to "HD13 LEU B  12"
- d= 1.80: "HD21 ASN A  16" to " O   LEU B  12"
- d= 1.71: " O   VAL B   2" to " N   ALA B   4"

</details>

