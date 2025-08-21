---
# 1QH7 :: Catalysis And Specificity In Enzymatic Glycoside Hydrolases: A 2,5B Conformation For The Glycosyl-Enzyme Intermidiate Revealed By The; Structure Of The Bacillus Agaradhaerens Family 11 Xylanase;
## PDB.Structure
### Function: HYDROLASE
### First Release: 11-MAY-99
### Method: X-RAY DIFFRACTION
### Resolution: 1.78 ANGSTROMS.
### Molecule: XYLANASE, B-D-XYLANOPYRANOSIDE PRESENT IN THE ACTIVE SITE
### Seqres Species: Residues: A:207, B:207; Total: 414
### Cofactors:
  - PCA:
 'PYROGLUTAMIC ACID', 2
  - XYP:
 'BETA-D-XYLOSE; D-XYLOSE; XYLOSE', 2

### Total cofactors: 4
### Total waters: 611
### Biounits: MONOMERIC: A,B;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 8, ARG: 8, ASN: 22, ASP: 8, CYS: 1, GLN: 8, GLU: 7, GLY: 24, HIS: 4, ILE: 12, LEU: 9, LYS: 9, MET: 5, PHE: 7, PRO: 7, SER: 17, THR: 16, TRP: 7, TYR: 13, VAL: 14', 'Total: 206', 'Ionizable: 50',
              'Ratio: 24.3%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 318, 'PCA': 1, 'XYP': 1

### Links:
  - C  PCA A  1 -- 1.41 Å --> N  ILE A  2
  - C  PCA B  1 -- 1.39 Å --> N  ILE B  2

### Sites:
  - NUA: ['CATALYTIC NUCLEOPHILE', 'GLU A  94']
  - NUB: ['CATALYTIC NUCLEOPHILE', 'GLU B  94']
  - ACA: ['CATALYTIC ACID / BASE', 'GLU A 184']
  - ACB: ['CATALYTIC ACID / BASE', 'GLU B 184']

## MCCE.Step1
### Renamed:
  - " C   PCA A   1" to " C_  PCA A   1"
  - " C   PCA A   1" to " C_  PCA A   1"
  - " N   PCA A   1" to " N_  PCA A   1"
  - " N   PCA A   1" to " N_  PCA A   1"
  - " O   PCA A   1" to " O_  PCA A   1"
  - " O   PCA A   1" to " O_  PCA A   1"

### Termini:
 - <strong>NTR</strong>: "ILE A   2"
 - <strong>CTR</strong>: "SER A 207"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
PCA: https://pubchem.ncbi.nlm.nih.gov/#query=PCA&tab=substance; XYP: https://pubchem.ncbi.nlm.nih.gov/#query=XYP&tab=substance; 

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
  - Removed all 318 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 318.
  - Species and properties with assigned default values in debug.log:

  - PCABK: ['VDW_RAD', 'VDW_EPS']

  - XYPBK: ['VDW_RAD', 'VDW_EPS']


### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.41: " C_  PCA A   1" to " N   NTR A   2"
- d= 1.52: " CA  NTR A   2" to " CB  ILE A   2"

</details>

