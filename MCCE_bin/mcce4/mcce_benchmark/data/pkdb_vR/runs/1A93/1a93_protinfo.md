---
# 1A93 :: Nmr Solution Structure Of The C-Myc-Max Heterodimeric Leucine Zipper, Nmr, Minimized Average Structure;
## PDB.Structure
### Function: LEUCINE ZIPPER
### First Release: 15-APR-98
### Method: SOLUTION NMR
### Resolution: NOT APPLICABLE.
### Molecule: MYC PROTO-ONCOGENE PROTEIN, BOTH C-MYC (CHAIN A) AND MAX (CHAIN B) HAVE A NON-
### Seqres Species: Residues: A:34, B:34; Total: 68
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
 'RESIDUES': ('ALA: 1, ARG: 3, ASP: 1, CYS: 1, GLN: 4, GLU: 6, GLY: 2, HIS: 1, ILE: 1, LEU: 6, LYS: 4, SER: 1, VAL: 1', 'Total: 32', 'Ionizable: 16',
              'Ratio: 50.0%')
  - B:
 'RESIDUES': ('ALA: 2, ARG: 4, ASN: 2, ASP: 4, CYS: 1, GLN: 5, GLU: 1, GLY: 2, HIS: 1, ILE: 1, LEU: 4, LYS: 2, MET: 1, THR: 1, VAL: 1', 'Total: 32', 'Ionizable: 13',
              'Ratio: 40.6%')

### Model 1 Free Cofactors & Waters:
  - A:
 'ACE': 1, 'NH2': 1
  - B:
 'ACE': 1, 'NH2': 1

### Disulfides:
  - CYS A  3 -- 2.02 Å --> CYS B  3

### Links:
  - C  ACE A  2 -- 1.33 Å --> N  CYS A  3
  - C  LEU A 34 -- 1.31 Å --> N  NH2 A 35
  - C  ACE B  2 -- 1.33 Å --> N  CYS B  3
  - C  LEU B 34 -- 1.31 Å --> N  NH2 B 35

### Sites:
  - AC3: ['BINDING SITE FOR RESIDUE NH2 A 35', 'LEU A  34']
  - AC4: ['BINDING SITE FOR RESIDUE NH2 B 35', 'LEU B  34']

## MCCE.Step1
### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
ACE: https://pubchem.ncbi.nlm.nih.gov/#query=ACE&tab=substance; NH2: https://pubchem.ncbi.nlm.nih.gov/#query=NH2&tab=substance; 

### Load Structure:
  -    Warning: load_pdb(): "Duplicate backbone atom ignored, " CH3 ACE A   2".
  -    Warning: load_pdb(): "Duplicate backbone atom ignored, " C   ACE A   2".
  -    Warning: load_pdb(): "Duplicate backbone atom ignored, " O   ACE A   2".

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Species and properties with assigned default values in debug.log:

  - ACEBK: ['VDW_RAD', 'VDW_EPS']

  - NH2BK: ['VDW_RAD', 'VDW_EPS']


### Distance Clashes:
<details><summary>Clashes found</summary>

- No clash found.

</details>

