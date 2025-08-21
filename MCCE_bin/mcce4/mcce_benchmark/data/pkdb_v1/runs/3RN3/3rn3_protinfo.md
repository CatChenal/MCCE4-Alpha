---
# 3RN3 :: Segmented Anisotropic Refinement Of Bovine Ribonuclease A By The Application Of The Rigid-Body Tls Model;
## PDB.Structure
### Function: HYDROLASE (NUCLEIC ACID,RNA)
### First Release: 30-OCT-91
### Method: X-RAY DIFFRACTION
### Resolution: 1.45 ANGSTROMS.
### Molecule: RIBONUCLEASE A
### Seqres Species: Residues: A:124; Total: 124
### Cofactors:
  - SO4:
 'SULFATE ION', 1

### Total cofactors: 1
### Total waters: 107
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 12, ARG: 4, ASN: 10, ASP: 5, CYS: 8, GLN: 7, GLU: 5, GLY: 3, HIS: 4, ILE: 3, LEU: 2, LYS: 10, MET: 4, PHE: 3, PRO: 4, SER: 15, THR: 10, TYR: 6, VAL: 9', 'Total: 124', 'Ionizable: 42',
              'Ratio: 33.9%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 107, 'SO4': 1

### Disulfides:
  - CYS A  26 -- 2.06 Å --> CYS A  84
  - CYS A  40 -- 2.09 Å --> CYS A  95
  - CYS A  58 -- 2.02 Å --> CYS A 110
  - CYS A  65 -- 2.07 Å --> CYS A  72

### Sites:
  - CAT: ['NULL', 'HIS A  12', 'LYS A  41', 'THR A  45', 'HIS A 119', 'PHE A 120']
  - AC1: ['BINDING SITE FOR RESIDUE SO4 A 155', 'GLN A  11', 'HIS A  12', 'HIS A 119', 'PHE A 120', 'HOH A 202', 'HOH A 209', 'HOH A 282']

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "LYS A   1"
 - <strong>CTR</strong>: "VAL A 124"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
SO4: https://pubchem.ncbi.nlm.nih.gov/#query=SO4&tab=substance; 

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 108.
  - Species and properties with assigned default values in debug.log:

  - SO4BK: ['VDW_RAD', 'VDW_EPS']


### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.50: " CA  NTR A   1" to " CB  LYS A   1"

</details>

