---
# 1PPN :: Structure Of Monoclinic Papain At 1.60 Angstroms Resolution
## PDB.Structure
### Function: HYDROLASE(SULFHYDRYL PROTEINASE)
### First Release: 25-OCT-91
### Method: X-RAY DIFFRACTION
### Resolution: 1.60 ANGSTROMS.
### Molecule: PAPAIN
### Seqres Species: Residues: A:212; Total: 212
### Cofactors:
  - UNL:
 'UNKNOWN LIGAND',
  - MOH:
 'METHANOL', 1

### Total cofactors: 1
### Total waters: 226
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 14, ARG: 12, ASN: 13, ASP: 6, CYS: 7, GLN: 10, GLU: 10, GLY: 28, HIS: 2, ILE: 12, LEU: 11, LYS: 10, PHE: 4, PRO: 10, SER: 13, THR: 8, TRP: 5, TYR: 19, VAL: 18', 'Total: 212', 'Ionizable: 66',
              'Ratio: 31.1%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 226, 'MOH': 1, 'UNL': 1

### Disulfides:
  - CYS A  22 -- 2.01 Å --> CYS A  63
  - CYS A  56 -- 2.05 Å --> CYS A  95
  - CYS A 153 -- 2.09 Å --> CYS A 200

### Links:
  - SG CYS A 25 -- 1.70 Å --> S  UNL A 301

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "ILE A   1"
 - <strong>CTR</strong>: "ASN A 212", "MOH A 302"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
UNL: https://pubchem.ncbi.nlm.nih.gov/#query=UNL&tab=substance; 

### Free Cofactors:
  - Removed all 226 HOH in A. Removed all 1 MOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 226.
  - Species and properties with assigned default values in debug.log:

  - UNLBK: ['VDW_RAD', 'VDW_EPS']

  - get_connect12():: ['Error!']

  - atom: ['get_connect12():']

  - CTR01: ['TORSION']


### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 302":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.49: " CA  NTR A   1" to " CB  ILE A   1"
- d= 1.70: " SG  CYS A  25" to " S   UNL A 301"

</details>

