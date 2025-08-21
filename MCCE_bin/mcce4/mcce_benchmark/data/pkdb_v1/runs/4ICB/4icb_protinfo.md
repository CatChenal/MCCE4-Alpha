---
# 4ICB :: Proline Cis-Trans Isomers In Calbindin D9K Observed By X-Ray Crystallography;
## PDB.Structure
### Function: CALCIUM-BINDING PROTEIN
### First Release: 27-AUG-91
### Method: X-RAY DIFFRACTION
### Resolution: 1.60 ANGSTROMS.
### Molecule: CALBINDIN D9K
### Seqres Species: Residues: A:76; Total: 76
### Cofactors:
  -  CA:
 'CALCIUM ION', 2

### Total cofactors: 2
### Total waters: 57
### Missing:
  - ResAtoms:
 'A': ['MET_0=(CG,SD,CE)']

### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 2, ASN: 2, ASP: 4, GLN: 4, GLU: 13, GLY: 5, ILE: 2, LEU: 12, LYS: 10, MET: 1, PHE: 5, PRO: 4, SER: 6, THR: 2, TYR: 1, VAL: 3', 'Total: 76', 'Ionizable: 28',
              'Ratio: 36.8%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 57, '_CA': 2

### Links:
  - O  ALA A 14 -- 2.33 Å --> CA  CA A 76
  - O  GLU A 17 -- 2.47 Å --> CA  CA A 76
  - O  ASP A 19 -- 2.25 Å --> CA  CA A 76
  - O  GLN A 22 -- 2.36 Å --> CA  CA A 76
  - OE1 GLU A 27 -- 2.42 Å --> CA  CA A 76
  - OE2 GLU A 27 -- 2.60 Å --> CA  CA A 76
  - OD1 ASP A 54 -- 2.41 Å --> CA  CA A 77
  - OD1 ASN A 56 -- 2.34 Å --> CA  CA A 77
  - OD1 ASP A 58 -- 2.39 Å --> CA  CA A 77
  - O  GLU A 60 -- 2.39 Å --> CA  CA A 77
  - OE1 GLU A 65 -- 2.54 Å --> CA  CA A 77
  - OE2 GLU A 65 -- 2.54 Å --> CA  CA A 77
  - CA  CA A 76 -- 2.43 Å --> O  HOH A 78
  - CA  CA A 77 -- 2.54 Å --> O  HOH A 79

### Sites:
  - CB1: ['CALCIUM BINDING SITE', 'ALA A  14', 'ALA A  15', 'LYS A  16', 'GLU A  17', 'GLY A  18', 'ASP A  19', 'PRO A  20', 'ASN A  21', 'GLN A  22', 'LEU A  23', 'SER A  24', 'LYS A  25', 'GLU A  26', 'GLU A  27']
  - CB2: ['CALCIUM BINDING SITE', 'ASP A  54', 'LYS A  55', 'ASN A  56', 'GLY A  57', 'ASP A  58', 'GLY A  59', 'GLU A  60', 'VAL A  61', 'SER A  62', 'PHE A  63', 'GLU A  64', 'GLU A  65']
  - AC1: ['BINDING SITE FOR RESIDUE CA A 76', 'ALA A  14', 'GLU A  17', 'ASP A  19', 'GLN A  22', 'GLU A  27', 'HOH A  78']
  - AC2: ['BINDING SITE FOR RESIDUE CA A 77', 'ASP A  54', 'ASN A  56', 'ASP A  58', 'GLU A  60', 'GLU A  65', 'HOH A  79']

## MCCE.Step1
### Renamed:
  - "CA    CA A  76" to "CA   _CA A  76"
  - "CA    CA A  77" to "CA   _CA A  77"

### Termini:
 - <strong>NTR</strong>: "MET A   0"
 - <strong>CTR</strong>: "GLN A  75"

### Free Cofactors:
  - Removed all 57 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 57.
  - Species and properties with assigned default values in debug.log:

  - _CA+2: ['VDW_RAD', 'VDW_EPS']


### Missing Heavy Atoms:
  -    Missing heavy atoms for MET01 in "MET A   0":   CG ,  SD ,  CE 

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.55: " CA  NTR A   0" to " CB  MET A   0"

</details>

