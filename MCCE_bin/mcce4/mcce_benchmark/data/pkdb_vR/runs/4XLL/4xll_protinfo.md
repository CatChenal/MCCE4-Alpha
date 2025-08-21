---
# 4XLL :: Toxoplasma Gondii Dj-1, Oxidized
## PDB.Structure
### Function: UNKNOWN FUNCTION
### First Release: 13-JAN-15
### Method: X-RAY DIFFRACTION
### Resolution: 2.08 ANGSTROMS.
### Molecule: DJ-1 FAMILY PROTEIN
### Seqres Species: Residues: A:188, B:188; Total: 376
### Cofactors:
  - CSD:
 'S-CYSTEINESULFINIC ACID; S-SULFINOCYSTEINE', 2

### Total cofactors: 2
### Total waters: 245
### Missing:
  - Residues:
 'A': [('GLY', -2), ('SER', -1), ('HIS', 0), ('MET', 1), ('Count', 4)],
 'B': [('GLY', -2), ('SER', -1), ('HIS', 0), ('MET', 1), ('Count', 4)]

### Biounits: DIMERIC: A, B;
### Models: 1
### Chains: A, B
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 27, ARG: 8, ASN: 3, ASP: 8, CYS: 5, GLN: 6, GLU: 14, GLY: 11, HIS: 3, ILE: 13, LEU: 14, LYS: 12, MET: 5, PHE: 3, PRO: 7, SER: 11, THR: 7, TYR: 5, VAL: 21', 'Total: 183', 'Ionizable: 55',
              'Ratio: 30.1%')
  - B:
 'RESIDUES': ('ALA: 27, ARG: 8, ASN: 3, ASP: 8, CYS: 5, GLN: 6, GLU: 14, GLY: 11, HIS: 3, ILE: 13, LEU: 14, LYS: 12, MET: 5, PHE: 3, PRO: 7, SER: 11, THR: 7, TYR: 5, VAL: 21', 'Total: 183', 'Ionizable: 55',
              'Ratio: 30.1%')

### Model 1 Free Cofactors & Waters:
  - A:
 'CSD': 1, 'HOH': 118
  - B:
 'CSD': 1, 'HOH': 127

### Links:
  - C  ILE A 103 -- 1.34 Å --> N  CSD A 104
  - C  CSD A 104 -- 1.33 Å --> N  ALA A 105
  - C  ILE B 103 -- 1.34 Å --> N  CSD B 104
  - C  CSD B 104 -- 1.34 Å --> N  ALA B 105

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "ALA A   2", "ALA B   2"
 - <strong>CTR</strong>: "TYR A 185", "TYR B 185"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
CSD: https://pubchem.ncbi.nlm.nih.gov/#query=CSD&tab=substance; 

### Free Cofactors:
  - Removed all 118 HOH in A. Removed all 127 HOH in B.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 245.
  - Species and properties with assigned default values in debug.log:

  - CSDBK: ['VDW_RAD', 'VDW_EPS']


### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 185":   OXT
  -    Missing heavy atoms for CTR01 in "CTR B 185":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   2" to " CB  ALA A   2"
- d= 1.53: " CA  NTR B   2" to " CB  ALA B   2"

</details>

