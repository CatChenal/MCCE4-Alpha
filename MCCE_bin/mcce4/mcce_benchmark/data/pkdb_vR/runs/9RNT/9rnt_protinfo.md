---
# 9RNT :: Ribonuclease T1 With Free Recognition And Catalytic Site: Crystal Structure Analysis At 1.5 Angstroms Resolution;
## PDB.Structure
### Function: HYDROLASE(ENDORIBONUCLEASE)
### First Release: 25-SEP-91
### Method: X-RAY DIFFRACTION
### Resolution: 1.50 ANGSTROMS.
### Molecule: RIBONUCLEASE T1
### Seqres Species: Residues: A:104; Total: 104
### Cofactors:
  -  CA:
 'CALCIUM ION', 1

### Total cofactors: 1
### Total waters: 121
### Missing:
  - ResAtoms:
 'A': ['ASN_98=(OD1,ND2)']

### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 7, ARG: 1, ASN: 9, ASP: 6, CYS: 4, GLN: 2, GLU: 6, GLY: 12, HIS: 3, ILE: 2, LEU: 3, LYS: 2, PHE: 4, PRO: 4, SER: 15, THR: 6, TRP: 1, TYR: 9, VAL: 8', 'Total: 104', 'Ionizable: 31',
              'Ratio: 29.8%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 121, '_CA': 1

### Disulfides:
  - CYS A  2 -- 2.03 Å --> CYS A  10
  - CYS A  6 -- 1.99 Å --> CYS A 103

### Links:
  - OD2 ASP A 15 -- 2.54 Å --> CA  CA A 105
  - OD1 ASP A 15 -- 2.56 Å --> CA  CA A 105
  - CA  CA A 105 -- 2.40 Å --> O  HOH A 109
  - CA  CA A 105 -- 2.49 Å --> O  HOH A 114
  - CA  CA A 105 -- 2.36 Å --> O  HOH A 118
  - CA  CA A 105 -- 2.47 Å --> O  HOH A 122
  - CA  CA A 105 -- 2.42 Å --> O  HOH A 132
  - CA  CA A 105 -- 2.62 Å --> O  HOH A 154

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE CA A 105', 'ASP A  15', 'HOH A 109', 'HOH A 114', 'HOH A 118', 'HOH A 122', 'HOH A 132', 'HOH A 154']

## MCCE.Step1
### Renamed:
  - "CA    CA A 105" to "CA   _CA A 105"

### Termini:
 - <strong>NTR</strong>: "ALA A   1"
 - <strong>CTR</strong>: "THR A 104"

### Free Cofactors:
  - Removed all 121 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 121.
  - Species and properties with assigned default values in debug.log:

  - _CA+2: ['VDW_RAD', 'VDW_EPS']


### Missing Heavy Atoms:
  -    Missing heavy atoms for ASN01 in "ASN A  98":   OD1,  ND2

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   1" to " CB  ALA A   1"
- d= 1.99: " SG  CYS A   6" to " SG  CYS A 103"

</details>

