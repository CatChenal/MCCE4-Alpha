---
# 7RSA :: Structure Of Phosphate-Free Ribonuclease A Refined At 1.26 Angstroms
## PDB.Structure
### Function: HYDROLASE (PHOSPHORIC DIESTER)
### First Release: 10-JUN-88
### Method: X-RAY DIFFRACTION
### Resolution: 1.26 ANGSTROMS.
### Molecule: RIBONUCLEASE A
### Seqres Species: Residues: A:124; Total: 124
### Cofactors:
  - TBU:
 '2-METHYL-2-PROPANOL', 1

### Total cofactors: 1
### Total waters: 188
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 12, ARG: 4, ASN: 10, ASP: 5, CYS: 8, GLN: 7, GLU: 5, GLY: 3, HIS: 4, ILE: 3, LEU: 2, LYS: 10, MET: 4, PHE: 3, PRO: 4, SER: 15, THR: 10, TYR: 6, VAL: 9', 'Total: 124', 'Ionizable: 42',
              'Ratio: 33.9%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 188, 'TBU': 1

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "LYS A   1"
 - <strong>CTR</strong>: "VAL A 124", "TBU A 125"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
TBU: https://pubchem.ncbi.nlm.nih.gov/#query=TBU&tab=substance; 

### Free Cofactors:
  - Removed all 188 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 188.
  - Species and properties with assigned default values in debug.log:

  - TBUBK: ['VDW_RAD', 'VDW_EPS']

  - get_connect12():: ['Error!']

  - atom: ['get_connect12():']

  - CTR01: ['TORSION']


### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 125":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.47: " CA  NTR A   1" to " CB  LYS A   1"
- d= 1.83: " O   GLN A  11" to "HD21 ASN A  44"
- d= 1.98: " SG  CYS A  58" to " SG  CYS A 110"
- d= 1.80: "HD22 ASN A  62" to " O   THR A  70"
- d= 1.97: " SG  CYS A  65" to " SG  CYS A  72"
- d= 1.93: "HD22 ASN A  67" to " OE1 GLN A  69"
- d= 2.00: "HE21 GLN A  69" to " OD1 ASN A  71"
- d= 1.96: "HD21 ASN A  71" to " O   CYD A 110"
- d= 1.52: " C1  TBU A 125" to " C   CTR A 125"
- d= 1.51: " C2  TBU A 125" to " C   CTR A 125"
- d= 1.53: " C3  TBU A 125" to " C   CTR A 125"

</details>

