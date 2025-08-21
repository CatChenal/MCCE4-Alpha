---
# 1PPO :: Determination Of The Structure Of Papaya Protease Omega
## PDB.Structure
### Function: HYDROLASE(THIOL PROTEASE)
### First Release: 12-JUL-91
### Method: X-RAY DIFFRACTION
### Resolution: 1.80 ANGSTROMS.
### Molecule: PROTEASE OMEGA
### Seqres Species: Residues: A:216; Total: 216
### Cofactors:
  -  HG:
 'MERCURY (II) ION', 1

### Total cofactors: 1
### Total waters: 133
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 14, ARG: 11, ASN: 10, ASP: 3, CYS: 7, GLN: 7, GLU: 11, GLY: 31, HIS: 4, ILE: 10, LEU: 11, LYS: 22, PHE: 3, PRO: 13, SER: 14, THR: 9, TRP: 4, TYR: 12, VAL: 20', 'Total: 216', 'Ionizable: 70',
              'Ratio: 32.4%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 133, '_HG': 1

### Disulfides:
  - CYS A  22 -- 1.99 Å --> CYS A  63
  - CYS A  56 -- 2.01 Å --> CYS A  95
  - CYS A 153 -- 2.05 Å --> CYS A 204

### Links:
  - SG CYS A 25 -- 2.42 Å --> HG  HG A 217
  - NZ LYS A 64 -- 2.89 Å --> HG  HG A 217
  - ND1 HIS A 159 -- 2.55 Å --> HG  HG A 217
  - HG  HG A 217 -- 2.75 Å --> O  HOH A 229

### Sites:
  - ACT: ['ACTIVE SITE', 'CYS A  25', 'HIS A 159', 'ASN A 179']
  - AC1: ['BINDING SITE FOR RESIDUE HG A 217', 'CYS A  25', 'LYS A  64', 'ASP A 158', 'HIS A 159', 'HOH A 229']

## MCCE.Step1
### Renamed:
  - "HG    HG A 217" to "HG   _HG A 217"

### Termini:
 - <strong>NTR</strong>: "LEU A   1"
 - <strong>CTR</strong>: "ASN A 216"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
_HG: https://pubchem.ncbi.nlm.nih.gov/#query=HG&tab=substance; 

### Free Cofactors:
  - Removed all 133 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 133.
  - Species and properties with assigned default values in debug.log:

  - _HGBK: ['VDW_RAD', 'VDW_EPS']


### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   1" to " CB  LEU A   1"
- d= 1.99: " SG  CYS A  22" to " SG  CYS A  63"
- d= 1.98: " C   GLY A 102" to " CD  PRO A 103"

</details>

