---
# 3SRN :: Structural Changes That Accompany The Reduced Catalytic Efficiency Of Two Semisynthetic Ribonuclease Analogs;
## PDB.Structure
### Function: HYDROLASE(NUCLEIC ACID,RNA)
### First Release: 20-MAY-91
### Method: X-RAY DIFFRACTION
### Resolution: 2.00 ANGSTROMS.
### Molecule: RIBONUCLEASE A
### Seqres Species: Residues: A:113, B:11; Total: 124
### Cofactors:
  - SO4:
 'SULFATE ION', 1

### Total cofactors: 1
### Total waters: 133
### Biounits: DIMERIC: A, B;
### Models: 1
### Chains: A, B
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 11, ARG: 4, ASN: 10, ASP: 4, CYS: 8, GLN: 7, GLU: 5, GLY: 3, HIS: 3, ILE: 3, LEU: 2, LYS: 10, MET: 4, PHE: 2, PRO: 2, SER: 14, THR: 10, TYR: 5, VAL: 6', 'Total: 113', 'Ionizable: 39',
              'Ratio: 34.5%')
  - B:
 'RESIDUES': ('ALA: 1, ASN: 1, HIS: 1, PHE: 1, PRO: 2, SER: 1, TYR: 1, VAL: 3', 'Total: 11', 'Ionizable: 2',
              'Ratio: 18.2%')

### Model 1 Free Cofactors & Waters:
  - B:
 'HOH': 9, 'SO4': 1
  - A:
 'HOH': 124

### Disulfides:
  - CYS A  26 -- 2.00 Å --> CYS A  84
  - CYS A  40 -- 2.01 Å --> CYS A  95
  - CYS A  58 -- 2.08 Å --> CYS A 110
  - CYS A  65 -- 2.08 Å --> CYS A  72

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE SO4 B 125', 'GLN A  11', 'HIS A  12', 'HOH A 153', 'HIS B 119', 'PHE B 120', 'HOH B 223']

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "LYS A   1", "PRO B 114"
 - <strong>CTR</strong>: "ASN A 113", "VAL B 124"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
SO4: https://pubchem.ncbi.nlm.nih.gov/#query=SO4&tab=substance; 

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 134.
  - Species and properties with assigned default values in debug.log:

  - SO4BK: ['VDW_RAD', 'VDW_EPS']

  - get_connect12():: ['Error!']

  - atom: ['get_connect12():']

  - PRO01: ['TORSION']

  - PROBK: ['TORSION']


### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 113":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.55: " CA  NTR A   1" to " CB  LYS A   1"
- d= 1.29: " O   LYS A   1" to " N   GLU A   2"
- d= 2.00: " SG  CYS A  26" to " SG  CYS A  84"
- d= 1.53: " CA  NTR B 114" to " CB  PRO B 114"
- d= 1.53: " N   NTR B 114" to " CD  PRO B 114"

</details>

