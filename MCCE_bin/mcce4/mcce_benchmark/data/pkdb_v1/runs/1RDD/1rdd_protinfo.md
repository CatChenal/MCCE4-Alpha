---
# 1RDD :: Crystal Structure Of Escherichia Coli Rnase Hi In Complex With Mg2+ At 2.8 Angstroms Resolution: Proof For A Single Mg2+ Site;
## PDB.Structure
### Function: HYDROLASE(ENDORIBONUCLEASE)
### First Release: 23-JUN-93
### Method: X-RAY DIFFRACTION
### Resolution: 2.80 ANGSTROMS.
### Molecule: RIBONUCLEASE H
### Seqres Species: Residues: A:155; Total: 155
### Cofactors:
  -  MG:
 'MAGNESIUM ION', 1

### Total cofactors: 1
### Total waters: 35
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 14, ARG: 10, ASN: 7, ASP: 7, CYS: 3, GLN: 8, GLU: 12, GLY: 14, HIS: 5, ILE: 7, LEU: 12, LYS: 11, MET: 4, PHE: 2, PRO: 5, SER: 4, THR: 10, TRP: 6, TYR: 5, VAL: 9', 'Total: 155', 'Ionizable: 53',
              'Ratio: 34.2%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 35, '_MG': 1

### Links:
  - OD1 ASP A 10 -- 1.97 Å --> MG  MG A 300
  - O  GLY A 11 -- 2.44 Å --> MG  MG A 300
  - OE2 GLU A 48 -- 2.45 Å --> MG  MG A 300

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE MG A 300', 'ASP A  10', 'GLY A  11', 'ASN A  44', 'GLU A  48']

## MCCE.Step1
### Renamed:
  - "MG    MG A 300" to "MG   _MG A 300"

### Termini:
 - <strong>NTR</strong>: "MET A   1"
 - <strong>CTR</strong>: "VAL A 155"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
_MG: https://pubchem.ncbi.nlm.nih.gov/#query=MG&tab=substance; 

### Free Cofactors:
  - Removed all 35 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 35.
  - Species and properties with assigned default values in debug.log:

  - _MGBK: ['VDW_RAD', 'VDW_EPS']


### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.54: " CA  NTR A   1" to " CB  MET A   1"
- d= 1.97: " OD1 ASP A  10" to "MG   _MG A 300"

</details>

