---
# 2L90 :: Solution Structure Of Murine Myristoylated Msra
## PDB.Structure
### Function: OXIDOREDUCTASE
### First Release: 27-JAN-11
### Method: SOLUTION NMR
### Resolution: NOT APPLICABLE.
### Molecule: PEPTIDE METHIONINE SULFOXIDE REDUCTASE
### Seqres Species: Residues: A:212; Total: 212
### Cofactors:
  - MYR:
 'MYRISTIC ACID', 1

### Total cofactors: 1
### Missing:
  - ResAtoms:
 'A': ['LYS_233=(CA,C,O,CB,CG,CD,CE)', 'LYS_233=(NZ)']

### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 14, ARG: 10, ASN: 6, ASP: 6, CYS: 4, GLN: 10, GLU: 18, GLY: 22, HIS: 8, ILE: 6, LEU: 8, LYS: 13, MET: 5, PHE: 10, PRO: 12, SER: 13, THR: 14, TRP: 3, TYR: 11, VAL: 19', 'Total: 212', 'Ionizable: 70',
              'Ratio: 33.0%')

### Model 1 Free Cofactors & Waters:
  - A:
 'MYR': 1

### Links:
  - C1 MYR A 21 -- 1.35 Å --> N  GLY A 22

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE MYR A 21', 'GLY A  22', 'SER A  24', 'LYS A  27', 'GLY A 112', 'ARG A 156']

## MCCE.Step1
### Termini:
 - <strong>NTG</strong>: "GLY A  22"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
MYR: https://pubchem.ncbi.nlm.nih.gov/#query=MYR&tab=substance; 

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Species and properties with assigned default values in debug.log:

  - MYRBK: ['VDW_RAD', 'VDW_EPS']

  - get_connect12():: ['Error!']

  - parameter: ['Check']


### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.35: " N   NTG A  22" to " C1  MYR A  21"
- d= 2.00: " O   ARG A  38" to "HD21 ASN A 100"
- d= 1.77: "HD22 ASN A  53" to " O   TYR A  87"
- d= 1.58: " OE2 GLU A  77" to "HE22 GLN A 208"
- d= 1.72: " O   VAL A  91" to "HE22 GLN A 207"
- d= 1.83: "HG13 VAL A 134" to "HD22 ASN A 138"
- d= 1.57: " O   ASP A 140" to "HE22 GLN A 147"
- d= 1.78: " O   GLN A 143" to "HE21 GLN A 147"

</details>

### Connectivity:
  -    Error! get_connect12(): Error in CONNECT parameter, refer to debug log file: debug.log for detail
  -    Errors were detected when making connectivity.

