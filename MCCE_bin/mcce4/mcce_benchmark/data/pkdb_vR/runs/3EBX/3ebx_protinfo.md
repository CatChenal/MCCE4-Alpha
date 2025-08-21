---
# 3EBX :: Refinement At 1.4 Angstroms Resolution Of A Model Of Erabutoxin B. Treatment Of Ordered Solvent And Discrete Disorder;
## PDB.Structure
### Function: TOXIN
### First Release: 15-JAN-88
### Method: X-RAY DIFFRACTION
### Resolution: 1.40 ANGSTROMS.
### Molecule: ERABUTOXIN B
### Seqres Species: Residues: A:62; Total: 62
### Cofactors:
  - SO4:
 'SULFATE ION', 1

### Total cofactors: 1
### Total waters: 111
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ARG: 3, ASN: 3, ASP: 1, CYS: 8, GLN: 4, GLU: 4, GLY: 5, HIS: 2, ILE: 4, LEU: 1, LYS: 4, PHE: 2, PRO: 4, SER: 8, THR: 5, TRP: 1, TYR: 1, VAL: 2', 'Total: 62', 'Ionizable: 23',
              'Ratio: 37.1%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 111, 'SO4': 1

### Disulfides:
  - CYS A  3 -- 2.04 Å --> CYS A  24
  - CYS A  17 -- 2.00 Å --> CYS A  41
  - CYS A  43 -- 2.06 Å --> CYS A  54
  - CYS A  55 -- 2.04 Å --> CYS A  60

### Sites:
  - RCT: ['SERIES INVARIANT RESIDUES OF THE REACTIVE site', 'TYR A  25', 'LYS A  27', 'TRP A  29', 'ASP A  31', 'PHE A  32', 'ARG A  33', 'GLY A  34', 'ILE A  36', 'GLU A  38', 'GLY A  40', 'CYS A  41', 'GLY A  42', 'CYS A  43', 'PRO A  44', 'VAL A  46', 'LYS A  47', 'GLY A  49', 'ILE A  50', 'LEU A  52', 'CYS A  54']
  - FNR: ['RESIDUES IN THE REACTIVE SITE SHOWN CHEMICALLY', 'LYS A  27', 'TRP A  29', 'ARG A  33', 'LYS A  47']
  - CMR: ['RESIDUES, INCLUDING FOUR CYSTINE LINKAGES, WHICH', 'CYS A   3', 'PHE A   4', 'CYS A  17', 'CYS A  24', 'TYR A  25', 'GLY A  40', 'CYS A  41', 'GLY A  42', 'CYS A  43', 'CYS A  54', 'CYS A  55', 'CYS A  60', 'ASN A  61']
  - AC1: ['BINDING SITE FOR RESIDUE SO4 A 63', 'ARG A   1', 'ASN A   5', 'LYS A  15', 'HOH A  74', 'HOH A 146', 'HOH A 150', 'HOH A 163']

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "ARG A   1"
 - <strong>CTR</strong>: "ASN A  62"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
SO4: https://pubchem.ncbi.nlm.nih.gov/#query=SO4&tab=substance; 

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 112.
  - Species and properties with assigned default values in debug.log:

  - SO4BK: ['VDW_RAD', 'VDW_EPS']


### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.52: " CA  NTR A   1" to " CB  ARG A   1"
- d= 2.00: " SG  CYD A  17" to " SG  CYS A  41"

</details>

