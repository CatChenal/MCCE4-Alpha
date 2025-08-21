---
# 1DG9 :: Crystal Structure Of Bovine Low Molecular Weight Ptpase Complexed With Hepes;
## PDB.Structure
### Function: HYDROLASE
### First Release: 23-NOV-99
### Method: X-RAY DIFFRACTION
### Resolution: 1.90 ANGSTROMS.
### Molecule: TYROSINE PHOSPHATASE
### Seqres Species: Residues: A:157; Total: 157
### Cofactors:
  - EPE:
 'HEPES', 1

### Total cofactors: 1
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 10, ARG: 12, ASN: 11, ASP: 13, CYS: 8, GLN: 8, GLU: 8, GLY: 6, HIS: 2, ILE: 9, LEU: 11, LYS: 9, MET: 1, PHE: 6, PRO: 5, SER: 10, THR: 6, TRP: 2, TYR: 5, VAL: 15', 'Total: 157', 'Ionizable: 57',
              'Ratio: 36.3%')

### Model 1 Free Cofactors & Waters:
  - A:
 'EPE': 1

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE EPE A 201', 'CYS A  12', 'LEU A  13', 'GLY A  14', 'ASN A  15', 'ILE A  16', 'CYS A  17', 'ARG A  18', 'ASP A 129', 'TYR A 131', 'TYR A 132']

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "ALA A   1"
 - <strong>CTR</strong>: "ARG A 157"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
EPE: https://pubchem.ncbi.nlm.nih.gov/#query=EPE&tab=substance; 

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Species and properties with assigned default values in debug.log:

  - EPEBK: ['VDW_RAD', 'VDW_EPS']


### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   1" to " CB  ALA A   1"

</details>

