---
# 2ZTA :: X-Ray Structure Of The Gcn4 Leucine Zipper, A Two-Stranded, Parallel Coiled Coil;
## PDB.Structure
### Function: LEUCINE ZIPPER
### First Release: 05-JUL-91
### Method: X-RAY DIFFRACTION
### Resolution: 1.80 ANGSTROMS.
### Molecule: GCN4 LEUCINE ZIPPER
### Seqres Species: Residues: A:34, B:34; Total: 68
### Cofactors:
  - ACE:
 'ACETYL GROUP', 2

### Total cofactors: 2
### Total waters: 52
### Missing:
  - Residues:
 'A': [('GLU', 32), ('ARG', 33), ('Count', 2)], 'B': [('GLU', 32), ('ARG', 33), ('Count', 2)]

### Biounits: DIMERIC: A, B;
### Models: 1
### Chains: A, B
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 1, ARG: 2, ASN: 2, ASP: 1, GLN: 1, GLU: 5, GLY: 1, HIS: 1, LEU: 6, LYS: 5, MET: 1, SER: 1, TYR: 1, VAL: 3', 'Total: 31', 'Ionizable: 15',
              'Ratio: 48.4%')
  - B:
 'RESIDUES': ('ALA: 1, ARG: 2, ASN: 2, ASP: 1, GLN: 1, GLU: 5, GLY: 1, HIS: 1, LEU: 6, LYS: 5, MET: 1, SER: 1, TYR: 1, VAL: 3', 'Total: 31', 'Ionizable: 15',
              'Ratio: 48.4%')

### Model 1 Free Cofactors & Waters:
  - A:
 'ACE': 1, 'HOH': 23
  - B:
 'ACE': 1, 'HOH': 29

### Links:
  - C  ACE A  0 -- 1.35 Å --> N  ARG A  1
  - C  ACE B  0 -- 1.31 Å --> N  ARG B  1

## MCCE.Step1
### Termini:
 - <strong>CTR</strong>: "GLY A  31", "GLY B  31"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
ACE: https://pubchem.ncbi.nlm.nih.gov/#query=ACE&tab=substance; 

### Load Structure:
  -    Warning: load_pdb(): "Duplicate backbone atom ignored, " CH3 ACE A   0".
  -    Warning: load_pdb(): "Duplicate backbone atom ignored, " C   ACE A   0".
  -    Warning: load_pdb(): "Duplicate backbone atom ignored, " O   ACE A   0".

### Free Cofactors:
  - Removed all 23 HOH in A. Removed all 29 HOH in B.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 52.
  - Species and properties with assigned default values in debug.log:

  - ACEBK: ['VDW_RAD', 'VDW_EPS']


### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A  31":   OXT
  -    Missing heavy atoms for CTR01 in "CTR B  31":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- No clash found.

</details>

