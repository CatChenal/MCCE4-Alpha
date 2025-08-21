---
# 3FX5 :: Structure Of Hiv-1 Protease In Complex With Potent Inhibitor Kni-272 Determined By High Resolution X-Ray Crystallography;
## PDB.Structure
### Function: HYDROLASE/HYDROLASE INHIBITOR
### First Release: 20-JAN-09
### Method: X-RAY DIFFRACTION
### Resolution: 0.93 ANGSTROMS.
### Molecule: PROTEASE
### Seqres Species: Residues: A:99, B:99; Total: 198
### Cofactors:
  - GOL:
 'GLYCERIN; PROPANE-1,2,3-TRIOL', 3
  - KNI:
 'ALLOPHENYLNORSTATINE', 1

### Total cofactors: 4
### Total waters: 470
### Biounits: DIMERIC: A, B;
### Models: 1
### Chains: A, B
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 5, ARG: 4, ASN: 3, ASP: 4, GLN: 5, GLU: 4, GLY: 13, HIS: 1, ILE: 15, LEU: 10, LYS: 7, MET: 2, PHE: 2, PRO: 6, SER: 1, THR: 8, TRP: 2, TYR: 1, VAL: 6', 'Total: 99', 'Ionizable: 21',
              'Ratio: 21.2%')
  - B:
 'RESIDUES': ('ALA: 5, ARG: 4, ASN: 3, ASP: 4, GLN: 5, GLU: 4, GLY: 13, HIS: 1, ILE: 15, LEU: 10, LYS: 7, MET: 2, PHE: 2, PRO: 6, SER: 1, THR: 8, TRP: 2, TYR: 1, VAL: 6', 'Total: 99', 'Ionizable: 21',
              'Ratio: 21.2%')

### Model 1 Free Cofactors & Waters:
  - A:
 'GOL': 1, 'HOH': 240
  - B:
 'GOL': 2, 'HOH': 230, 'KNI': 1

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE GOL A 901', 'LYS A  45', 'MET A  46', 'PHE A  53', 'HOH A 324', 'GLN B 161', 'HOH B 306', 'HOH B 313', 'HOH B 509']
  - AC2: ['BINDING SITE FOR RESIDUE KNI B 900', 'ARG A   8', 'ASP A  25', 'GLY A  27', 'ALA A  28', 'GLY A  48', 'GLY A  49', 'ILE A  50', 'PRO A  81', 'ILE A  84', 'HOH A 301', 'HOH A 608', 'ASP B 125', 'GLY B 127', 'ALA B 128', 'ASP B 129', 'ASP B 130', 'VAL B 132', 'GLY B 148', 'GLY B 149', 'ILE B 150', 'PRO B 181', 'ILE B 184', 'HOH B 560', 'HOH B 566', 'HOH B 607']
  - AC3: ['BINDING SITE FOR RESIDUE GOL B 902', 'GLN A  18', 'MET A  36', 'SER A  37', 'THR B 112', 'GLU B 165', 'ILE B 166', 'ALA B 167', 'GLY B 168', 'HOH B 559']
  - AC4: ['BINDING SITE FOR RESIDUE GOL B 903', 'HOH A 771', 'GLU B 135', 'LYS B 155', 'ARG B 157', 'VAL B 177', 'GLY B 178', 'HOH B 528', 'HOH B 675']

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "PRO A   1", "PRO B 101"
 - <strong>CTR</strong>: "PHE A  99", "PHE B 199"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
GOL: https://pubchem.ncbi.nlm.nih.gov/#query=GOL&tab=substance; KNI: https://pubchem.ncbi.nlm.nih.gov/#query=KNI&tab=substance; 

### Free Cofactors:
  - Removed all 240 HOH in A. Removed all 230 HOH in B.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 457.
  - Species and properties with assigned default values in debug.log:

  - GOLBK: ['VDW_RAD', 'VDW_EPS']

  - KNIBK: ['VDW_RAD', 'VDW_EPS']

  - get_connect12():: ['Error!']

  - parameter: ['Check']

  - atom: ['get_connect12():']

  - PRO01: ['TORSION']


### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   1" to " CB  PRO A   1"
- d= 1.45: " N   NTR A   1" to " CD  PRO A   1"
- d= 1.53: " CA  NTR B 101" to " CB  PRO B 101"
- d= 1.49: " N   NTR B 101" to " CD  PRO B 101"

</details>

### Connectivity:
  -    Error! get_connect12(): Error in CONNECT parameter, refer to debug log file: debug.log for detail
  -    Errors were detected when making connectivity.

