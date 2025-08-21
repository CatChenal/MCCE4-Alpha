---
# 2TGA :: On The Disordered Activation Domain In Trypsinogen. Chemical Labelling And Low-Temperature Crystallography;
## PDB.Structure
### Function: HYDROLASE ZYMOGEN (SERINE PROTEINASE)
### First Release: 26-OCT-81
### Method: X-RAY DIFFRACTION
### Resolution: 1.80 ANGSTROMS.
### Molecule: TRYPSINOGEN
### Seqres Species: Residues: A:229; Total: 229
### Cofactors:
  -  CA:
 'CALCIUM ION', 1

### Total cofactors: 1
### Total waters: 93
### Missing:
  - Residues:
 'A': [('VAL', 10), ('ASP', 11), ('ASP', 12), ('ASP', 13), ('ASP', 14), ('LYS', 15), ('Count', 6)]

### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 14, ARG: 2, ASN: 16, ASP: 6, CYS: 12, GLN: 10, GLU: 4, GLY: 25, HIS: 3, ILE: 15, LEU: 14, LYS: 14, MET: 2, PHE: 3, PRO: 8, SER: 34, THR: 10, TRP: 4, TYR: 10, VAL: 17', 'Total: 223', 'Ionizable: 51',
              'Ratio: 22.9%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 93, '_CA': 1

### Disulfides:
  - CYS A  22 -- 2.06 Å --> CYS A 157
  - CYS A  42 -- 2.05 Å --> CYS A  58
  - CYS A 128 -- 2.05 Å --> CYS A 232
  - CYS A 136 -- 2.05 Å --> CYS A 201
  - CYS A 168 -- 2.05 Å --> CYS A 182
  - CYS A 191 -- 2.06 Å --> CYS A 220

### Links:
  - OE1 GLU A 70 -- 2.23 Å --> CA  CA A 480
  - O  ASN A 72 -- 2.24 Å --> CA  CA A 480
  - O  VAL A 75 -- 2.25 Å --> CA  CA A 480
  - OE2 GLU A 80 -- 2.23 Å --> CA  CA A 480
  - CA  CA A 480 -- 2.33 Å --> O  HOH A 711
  - CA  CA A 480 -- 2.27 Å --> O  HOH A 714

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE CA A 480', 'GLU A  70', 'ASN A  72', 'VAL A  75', 'GLU A  80', 'HOH A 711', 'HOH A 714']

## MCCE.Step1
### Renamed:
  - "CA    CA A 480" to "CA   _CA A 480"

### Termini:
 - <strong>NTR</strong>: "ILE A  16"
 - <strong>CTR</strong>: "ASN A 245"

### Free Cofactors:
  - Removed all 93 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 93.
  - Species and properties with assigned default values in debug.log:

  - _CA+2: ['VDW_RAD', 'VDW_EPS']


### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.55: " CA  NTR A  16" to " CB  ILE A  16"

</details>

