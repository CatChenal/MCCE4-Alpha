---
# 1L98 :: Perturbation Of Trp 138 In T4 Lysozyme By Mutations At Gln 105 Used To Correlate Changes In Structure, Stability, Solvation, And; Spectroscopic Properties;
## PDB.Structure
### Function: HYDROLASE(O-GLYCOSYL)
### First Release: 13-JUL-92
### Method: X-RAY DIFFRACTION
### Resolution: 1.80 ANGSTROMS.
### Molecule: T4 LYSOZYME
### Seqres Species: Residues: A:164; Total: 164
### Cofactors:
  - BME:
 'BETA-MERCAPTOETHANOL', 1

### Total cofactors: 1
### Total waters: 165
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 15, ARG: 13, ASN: 12, ASP: 10, CYS: 2, GLN: 4, GLU: 9, GLY: 11, HIS: 1, ILE: 10, LEU: 16, LYS: 13, MET: 5, PHE: 5, PRO: 3, SER: 6, THR: 11, TRP: 3, TYR: 6, VAL: 9', 'Total: 164', 'Ionizable: 54',
              'Ratio: 32.9%')

### Model 1 Free Cofactors & Waters:
  - A:
 'BME': 1, 'HOH': 165

### Links:
  - SG CYS A 97 -- 2.06 Å --> S2 BME A 165

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE BME A 165', 'ILE A   3', 'ASN A  68', 'ALA A  93', 'ARG A  96', 'CYS A  97', 'ILE A 100', 'HOH A 317']

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "MET A   1"
 - <strong>CTR</strong>: "LEU A 164"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
BME: https://pubchem.ncbi.nlm.nih.gov/#query=BME&tab=substance; 

### Free Cofactors:
  - Removed all 165 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 165.
  - Species and properties with assigned default values in debug.log:

  - BMEBK: ['VDW_RAD', 'VDW_EPS']


### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.55: " CA  NTR A   1" to " CB  MET A   1"

</details>

