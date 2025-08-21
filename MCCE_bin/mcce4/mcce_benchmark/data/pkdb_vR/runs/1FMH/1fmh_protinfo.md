---
# 1FMH :: Nmr Solution Structure Of A Designed Heterodimeric Leucine Zipper
## PDB.Structure
### Function: TRANSCRIPTION
### First Release: 17-AUG-00
### Method: SOLUTION NMR
### Resolution: NOT APPLICABLE.
### Molecule: GENERAL CONTROL PROTEIN GCN4
### Seqres Species: Residues: A:33, B:33; Total: 66
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
 'RESIDUES': ('ALA: 5, ASN: 1, CYS: 1, GLN: 5, GLU: 9, GLY: 1, HIS: 1, LEU: 3, LYS: 1, TYR: 1, VAL: 3', 'Total: 31', 'Ionizable: 13',
              'Ratio: 41.9%')
  - B:
 'RESIDUES': ('ALA: 6, ARG: 3, ASN: 1, CYS: 1, GLN: 4, GLU: 1, GLY: 1, HIS: 1, LEU: 3, LYS: 6, TYR: 1, VAL: 3', 'Total: 31', 'Ionizable: 13',
              'Ratio: 41.9%')

### Model 1 Free Cofactors & Waters:
  - A:
 'ACE': 1, 'NH2': 1
  - B:
 'ACE': 1, 'NH2': 1

### Disulfides:
  - CYS A  30 -- 2.03 Å --> CYS B  30

### Links:
  - C  ACE A  0 -- 1.33 Å --> N  GLU A  1
  - C  GLY A 31 -- 1.33 Å --> N  NH2 A 32
  - C  ACE B  0 -- 1.32 Å --> N  GLU B  1
  - C  GLY B 31 -- 1.33 Å --> N  NH2 B 32

### Sites:
  - AC3: ['BINDING SITE FOR RESIDUE NH2 A 32', 'CYS A  30', 'GLY A  31']
  - AC4: ['BINDING SITE FOR RESIDUE NH2 B 32', 'CYS B  30', 'GLY B  31']

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

- d= 1.75: "HG11 VAL A   2" to "HD13 LEU B   5"
- d= 1.69: "HD13 LEU A   5" to "HG12 VAL B   2"
- d= 1.92: "HD21 ASN A  16" to " OD1 ASN B  16"

</details>

