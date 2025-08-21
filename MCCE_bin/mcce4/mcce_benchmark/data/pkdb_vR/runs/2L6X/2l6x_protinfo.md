---
# 2L6X :: Solution Nmr Structure Of Proteorhodopsin.
## PDB.Structure
### Function: PROTON TRANSPORT
### First Release: 29-NOV-10
### Method: SOLUTION NMR
### Resolution: NOT APPLICABLE.
### Molecule: GREEN-LIGHT ABSORBING PROTEORHODOPSIN
### Seqres Species: Residues: A:243; Total: 243
### Cofactors:
  - RET:
 'RETINAL', 1

### Total cofactors: 1
### Missing:
  - Residues:
 'A': [('GLY', 253), ('SER', 254), ('HIS', 255), ('HIS', 256), ('HIS', 257), ('HIS', 258), ('HIS', 259), ('HIS', 260),
       ('Count', 8)]

### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 25, ARG: 4, ASN: 9, ASP: 8, CYS: 3, GLN: 1, GLU: 7, GLY: 25, HIS: 1, ILE: 17, LEU: 26, LYS: 7, MET: 10, PHE: 14, PRO: 6, SER: 15, THR: 13, TRP: 10, TYR: 14, VAL: 20', 'Total: 235', 'Ionizable: 44',
              'Ratio: 18.7%')

### Model 1 Free Cofactors & Waters:
  - A:
 'RET': 1

### Links:
  - NZ LYS A 231 -- 1.33 Å --> C15 RET A 301

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE RET A 301', 'TRP A  98', 'VAL A 133', 'PHE A 137', 'ALA A 151', 'ILE A 154', 'TRP A 197', 'TYR A 200', 'ASP A 227', 'ASN A 230', 'LYS A 231']

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "MET A  18"
 - <strong>CTR</strong>: "GLY A 252"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
RET: https://pubchem.ncbi.nlm.nih.gov/#query=RET&tab=substance; 

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Species and properties with assigned default values in debug.log:

  - RETBK: ['VDW_RAD', 'VDW_EPS']

  - get_connect12():: ['Error!']

  - parameter: ['Check']


### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 252":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A  18" to " CB  MET A  18"
- d= 1.97: "HG21 VAL A  54" to "HG22 VAL A 243"
- d= 2.00: "HD13 LEU A 105" to "HG12 ILE A 106"
- d= 1.96: "HG11 VAL A 133" to "H161 RET A 301"
- d= 1.99: "HG12 ILE A 163" to "HD11 ILE A 194"
- d= 1.96: "HD11 ILE A 163" to "HG12 ILE A 194"
- d= 1.97: "HG13 ILE A 192" to "HG13 ILE A 237"
- d= 1.33: " NZ  LYS A 231" to " C15 RET A 301"

</details>

### Connectivity:
  -    Error! get_connect12(): Error in CONNECT parameter, refer to debug log file: debug.log for detail
  -    Errors were detected when making connectivity.

