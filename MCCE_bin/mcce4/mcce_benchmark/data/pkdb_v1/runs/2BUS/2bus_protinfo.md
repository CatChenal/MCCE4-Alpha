---
# 2BUS :: Solution Conformation Of Proteinase Inhibitor Iia From Bull Seminal Plasma By 1H Nuclear Magnetic Resonance And Distance Geometry;
## PDB.Structure
### Function: PROTEINASE INHIBITOR
### First Release: 14-MAY-90
### Method: SOLUTION NMR
### Resolution: NOT APPLICABLE.
### Molecule: PROTEINASE INHIBITOR IIA
### Seqres Species: Residues: A:57; Total: 57
### Cofactors:
  - PCA:
 'PYROGLUTAMIC ACID', 1

### Total cofactors: 1
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 4, ARG: 2, ASN: 4, ASP: 2, CYS: 6, GLN: 1, GLU: 3, GLY: 7, HIS: 2, ILE: 1, LEU: 1, LYS: 8, MET: 1, PHE: 2, PRO: 2, SER: 3, THR: 2, TYR: 2, VAL: 3', 'Total: 56', 'Ionizable: 25',
              'Ratio: 44.6%')

### Model 1 Free Cofactors & Waters:
  - A:
 'PCA': 1

### Disulfides:
  - CYS A  7 -- 2.04 Å --> CYS A  39
  - CYS A  17 -- 2.08 Å --> CYS A  36
  - CYS A  25 -- 1.92 Å --> CYS A  57

### Links:
  - C  PCA A  1 -- 1.32 Å --> N  GLY A  2

## MCCE.Step1
### Renamed:
  - " C   PCA A   1" to " C_  PCA A   1"
  - " C   PCA A   1" to " C_  PCA A   1"
  - " N   PCA A   1" to " N_  PCA A   1"
  - " N   PCA A   1" to " N_  PCA A   1"
  - " O   PCA A   1" to " O_  PCA A   1"
  - " O   PCA A   1" to " O_  PCA A   1"

### Termini:
 - <strong>NTG</strong>: "GLY A   2"
 - <strong>CTR</strong>: "CYS A  57"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
PCA: https://pubchem.ncbi.nlm.nih.gov/#query=PCA&tab=substance; 

### Load Structure:
  -    Warning: load_pdb(): "Duplicate backbone atom ignored, " N_  PCA A   1".
  -    Warning: load_pdb(): "Duplicate backbone atom ignored, " CA  PCA A   1".
  -    Warning: load_pdb(): "Duplicate backbone atom ignored, " CB  PCA A   1".
  -    Warning: load_pdb(): "Duplicate backbone atom ignored, " CG  PCA A   1".
  -    Warning: load_pdb(): "Duplicate backbone atom ignored, " CD  PCA A   1".
  -    Warning: load_pdb(): "Duplicate backbone atom ignored, " OE  PCA A   1".
  -    Warning: load_pdb(): "Duplicate backbone atom ignored, " C_  PCA A   1".
  -    Warning: load_pdb(): "Duplicate backbone atom ignored, " O_  PCA A   1".

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Species and properties with assigned default values in debug.log:

  - PCABK: ['VDW_RAD', 'VDW_EPS']

  - NTG01: ['TORSION']


### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A  57":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.32: " C_  PCA A   1" to " N   NTG A   2"
- d= 1.92: " SG  CYS A  25" to " SG  CYS A  57"

</details>

