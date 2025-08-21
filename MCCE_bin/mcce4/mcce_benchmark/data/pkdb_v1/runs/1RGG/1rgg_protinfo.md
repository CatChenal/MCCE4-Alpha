---
# 1RGG :: Hydrolase, Guanyloribonuclease
## PDB.Structure
### Function: HYDROLASE (GUANYLORIBONUCLEASE)
### First Release: 05-JUN-95
### Method: X-RAY DIFFRACTION
### Resolution: 1.20 ANGSTROMS.
### Molecule: RIBONUCLEASE
### Seqres Species: Residues: A:96, B:96; Total: 192
### Cofactors:
  - SO4:
 'SULFATE ION', 2

### Total cofactors: 2
### Total waters: 338
### Biounits: DIMERIC: A, B;
### Models: 1
### Chains: A, B
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 6, ARG: 5, ASN: 2, ASP: 7, CYS: 2, GLN: 5, GLU: 5, GLY: 8, HIS: 2, ILE: 5, LEU: 6, PHE: 3, PRO: 6, SER: 7, THR: 13, TYR: 8, VAL: 6', 'Total: 96', 'Ionizable: 29',
              'Ratio: 30.2%')
  - B:
 'RESIDUES': ('ALA: 6, ARG: 5, ASN: 2, ASP: 7, CYS: 2, GLN: 5, GLU: 5, GLY: 8, HIS: 2, ILE: 5, LEU: 6, PHE: 3, PRO: 6, SER: 7, THR: 13, TYR: 8, VAL: 6', 'Total: 96', 'Ionizable: 29',
              'Ratio: 30.2%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 180, 'SO4': 2
  - B:
 'HOH': 158

### Disulfides:
  - CYS A  7 -- 2.03 Å --> CYS A  96
  - CYS B  7 -- 2.04 Å --> CYS B  96

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE SO4 A 97', 'ALA A  62', 'ARG A  63', 'THR A  64', 'HOH A 208', 'HOH A 211', 'PRO B  12', 'ARG B  68']
  - AC2: ['BINDING SITE FOR RESIDUE SO4 A 98', 'GLU A  54', 'ARG A  65', 'ARG A  69', 'HIS A  85', 'TYR A  86', 'HOH A 175', 'HOH A 222', 'HOH A 247']

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "ASP A   1", "ASP B   1"
 - <strong>CTR</strong>: "CYS A  96", "CYS B  96"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
SO4: https://pubchem.ncbi.nlm.nih.gov/#query=SO4&tab=substance; 

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 340.
  - Species and properties with assigned default values in debug.log:

  - SO4BK: ['VDW_RAD', 'VDW_EPS']


### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.57: " CA  NTR A   1" to " CB  ASP A   1"
- d= 1.54: " CA  NTR B   1" to " CB  ASP B   1"

</details>

