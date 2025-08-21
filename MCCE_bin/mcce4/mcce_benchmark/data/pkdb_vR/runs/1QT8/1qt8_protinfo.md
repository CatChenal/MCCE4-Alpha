---
# 1QT8 :: T26H Mutant Of T4 Lysozyme
## PDB.Structure
### Function: HYDROLASE
### First Release: 30-JUN-99
### Method: X-RAY DIFFRACTION
### Resolution: 1.90 ANGSTROMS.
### Molecule: PROTEIN (T4 LYSOZYME)
### Seqres Species: Residues: A:164; Total: 164
### Cofactors:
  - HED:
 '2-HYDROXYETHYL DISULFIDE', 1

### Total cofactors: 1
### Total waters: 140
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 16, ARG: 13, ASN: 12, ASP: 10, GLN: 5, GLU: 8, GLY: 11, HIS: 2, ILE: 10, LEU: 16, LYS: 13, MET: 5, PHE: 5, PRO: 3, SER: 6, THR: 11, TRP: 3, TYR: 6, VAL: 9', 'Total: 164', 'Ionizable: 52',
              'Ratio: 31.7%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HED': 1, 'HOH': 140

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE HED A 170', 'ASP A  72', 'HOH A 203']

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "MET A   1"
 - <strong>CTR</strong>: "LEU A 164"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
HED: https://pubchem.ncbi.nlm.nih.gov/#query=HED&tab=substance; 

### Free Cofactors:
  - Removed all 140 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 140.
  - Species and properties with assigned default values in debug.log:

  - HEDBK: ['VDW_RAD', 'VDW_EPS']


### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.55: " CA  NTR A   1" to " CB  MET A   1"

</details>

