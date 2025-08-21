---
# 1Z12 :: Crystal Structure Of Bovine Low Molecular Weight Ptpase Complexed With Vanadate;
## PDB.Structure
### Function: HYDROLASE
### First Release: 03-MAR-05
### Method: X-RAY DIFFRACTION
### Resolution: 2.20 ANGSTROMS.
### Molecule: LOW MOLECULAR WEIGHT PHOSPHOTYROSINE PROTEIN PHOSPHATASE
### Seqres Species: Residues: A:157; Total: 157
### Cofactors:
  - VO4:
 'VANADATE ION', 1

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
 'VO4': 1

### Links:
  - SG CYS A 12 -- 2.16 Å --> V  VO4 A 158

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE VO4 A 158', 'CYS A  12', 'LEU A  13', 'GLY A  14', 'ASN A  15', 'ILE A  16', 'CYS A  17', 'ARG A  18']

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "ALA A   1"
 - <strong>CTR</strong>: "ARG A 157"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
VO4: https://pubchem.ncbi.nlm.nih.gov/#query=VO4&tab=substance; 

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 1.
  - Species and properties with assigned default values in debug.log:

  - VO4BK: ['VDW_RAD', 'VDW_EPS']


### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 157":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.52: " CA  NTR A   1" to " CB  ALA A   1"

</details>

