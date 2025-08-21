---
# 2LZT :: Refinement Of Triclinic Lysozyme. Ii. The Method Of Stereochemically Restrained Least-Squares;
## PDB.Structure
### Function: HYDROLASE(O-GLYCOSYL)
### First Release: 21-AUG-89
### Method: X-RAY DIFFRACTION
### Resolution: 1.97 ANGSTROMS.
### Molecule: HEN EGG WHITE LYSOZYME
### Seqres Species: Residues: A:129; Total: 129
### Cofactors:
  - NO3:
 'NITRATE ION', 5

### Total cofactors: 5
### Total waters: 249
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 12, ARG: 11, ASN: 14, ASP: 7, CYS: 8, GLN: 3, GLU: 2, GLY: 12, HIS: 1, ILE: 6, LEU: 8, LYS: 6, MET: 2, PHE: 3, PRO: 2, SER: 10, THR: 7, TRP: 6, TYR: 3, VAL: 6', 'Total: 129', 'Ionizable: 38',
              'Ratio: 29.5%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 249, 'NO3': 5

### Disulfides:
  - CYS A  6 -- 2.06 Å --> CYS A 127
  - CYS A  30 -- 2.04 Å --> CYS A 115
  - CYS A  64 -- 2.00 Å --> CYS A  80
  - CYS A  76 -- 1.94 Å --> CYS A  94

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE NO3 A 130', 'ASN A  65', 'ASN A  74', 'ASN A  77', 'ILE A  78', 'PRO A  79', 'ARG A 112', 'LYS A 116']
  - AC2: ['BINDING SITE FOR RESIDUE NO3 A 131', 'ARG A  21', 'ASP A  66', 'CYS A  80', 'SER A  81', 'HOH A 175', 'HOH A 187']
  - AC3: ['BINDING SITE FOR RESIDUE NO3 A 132', 'LYS A  33', 'PHE A  38', 'ARG A  73', 'HOH A 305', 'HOH A 319', 'HOH A 340']
  - AC4: ['BINDING SITE FOR RESIDUE NO3 A 133', 'ALA A  11', 'HIS A  15', 'ASP A  87', 'ILE A  88', 'HOH A 317', 'HOH A 349']
  - AC5: ['BINDING SITE FOR RESIDUE NO3 A 134', 'SER A  24', 'LEU A  25', 'GLY A  26', 'GLN A  41', 'GLN A 121', 'ILE A 124', 'HOH A 193']

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "LYS A   1"
 - <strong>CTR</strong>: "LEU A 129"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 254.
  - Species and properties with assigned default values in debug.log:

  - NO3-1: ['VDW_RAD', 'VDW_EPS']


### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.54: " CA  NTR A   1" to " CB  LYS A   1"
- d= 2.00: " SG  CYS A  64" to " SG  CYS A  80"
- d= 1.94: " SG  CYS A  76" to " SG  CYS A  94"

</details>

