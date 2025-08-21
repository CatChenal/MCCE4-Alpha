---
# 1IG5 :: Bovine Calbindin D9K Binding Mg2+
## PDB.Structure
### Function: METAL BINDING PROTEIN
### First Release: 17-APR-01
### Method: X-RAY DIFFRACTION
### Resolution: 1.50 ANGSTROMS.
### Molecule: VITAMIN D-DEPENDENT CALCIUM-BINDING PROTEIN, INTESTINAL
### Seqres Species: Residues: A:75; Total: 75
### Cofactors:
  -  MG:
 'MAGNESIUM ION', 1

### Total cofactors: 1
### Total waters: 50
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
 'HOH': 50, '_MG': 1

### Links:
  - OD1 ASP A 54 -- 2.06 Å --> MG  MG A 78
  - OD1 ASN A 56 -- 2.05 Å --> MG  MG A 78
  - OD1 ASP A 58 -- 2.11 Å --> MG  MG A 78
  - O  GLU A 60 -- 1.99 Å --> MG  MG A 78
  - MG  MG A 78 -- 1.96 Å --> O  HOH A 110
  - MG  MG A 78 -- 2.10 Å --> O  HOH A 112

### Sites:
  - MUM: ['MG BINDING EF-HAND LOOP.', 'ASP A  54', 'ASN A  56', 'ASP A  58', 'GLU A  60']
  - AC1: ['BINDING SITE FOR RESIDUE MG A 78', 'ASP A  54', 'ASN A  56', 'ASP A  58', 'GLU A  60', 'GLN A  75', 'HOH A 110', 'HOH A 112']

## MCCE.Step1
### Renamed:
  - "MG    MG A  78" to "MG   _MG A  78"

### Termini:
 - <strong>NTR</strong>: "LYS A   1"
 - <strong>CTR</strong>: "GLN A  75"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
_MG: https://pubchem.ncbi.nlm.nih.gov/#query=MG&tab=substance; 

### Free Cofactors:
  - Removed all 50 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 50.
  - Species and properties with assigned default values in debug.log:

  - _MGBK: ['VDW_RAD', 'VDW_EPS']


### Missing Heavy Atoms:
  -    Missing heavy atoms for GLU01 in "GLU A  51":   CG ,  CD ,  OE1,  OE2

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.55: " CA  NTR A   1" to " CB  LYS A   1"
- d= 1.99: " O   GLU A  60" to "MG   _MG A  78"

</details>

