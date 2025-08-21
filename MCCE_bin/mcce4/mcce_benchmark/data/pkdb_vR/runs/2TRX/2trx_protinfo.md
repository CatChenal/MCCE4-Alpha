---
# 2TRX :: Crystal Structure Of Thioredoxin From Escherichia Coli At 1.68 Angstroms Resolution;
## PDB.Structure
### Function: ELECTRON TRANSPORT
### First Release: 19-MAR-90
### Method: X-RAY DIFFRACTION
### Resolution: 1.68 ANGSTROMS.
### Molecule: THIOREDOXIN
### Seqres Species: Residues: A:108, B:108; Total: 216
### Cofactors:
  -  CU:
 'COPPER (II) ION', 2
  - MPD:
 '(4S)-2-METHYL-2,4-PENTANEDIOL', 7

### Total cofactors: 9
### Total waters: 140
### Biounits: MONOMERIC: A,B; : A,B;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 12, ARG: 1, ASN: 4, ASP: 11, CYS: 2, GLN: 3, GLU: 5, GLY: 9, HIS: 1, ILE: 9, LEU: 13, LYS: 10, MET: 1, PHE: 4, PRO: 5, SER: 3, THR: 6, TRP: 2, TYR: 2, VAL: 5', 'Total: 108', 'Ionizable: 32',
              'Ratio: 29.6%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 75, 'MPD': 4, '_CU': 1

### Disulfides:
  - CYS A  32 -- 2.09 Å --> CYS A  35
  - CYS B  32 -- 2.05 Å --> CYS B  35

### Links:
  - N  SER A  1 -- 2.05 Å --> CU  CU A 109
  - N  ASP A  2 -- 2.06 Å --> CU  CU A 109
  - OD1 ASP A  2 -- 2.00 Å --> CU  CU A 109
  - OD1 ASP A 10 -- 1.97 Å --> CU  CU A 109
  - OD2 ASP A 10 -- 2.62 Å --> CU  CU A 109
  - CU  CU A 109 -- 2.65 Å --> O  HOH A 405
  - N  SER B  1 -- 2.09 Å --> CU  CU B 109
  - N  ASP B  2 -- 2.05 Å --> CU  CU B 109
  - OD1 ASP B  2 -- 2.06 Å --> CU  CU B 109
  - OD1 ASP B 10 -- 2.08 Å --> CU  CU B 109
  - OD2 ASP B 10 -- 2.54 Å --> CU  CU B 109
  - CU  CU B 109 -- 2.63 Å --> O  HOH B 478

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE CU A 109', 'SER A   1', 'ASP A   2', 'LYS A   3', 'ASP A  10', 'HOH A 405']
  - AC2: ['BINDING SITE FOR RESIDUE CU B 109', 'SER B   1', 'ASP B   2', 'LYS B   3', 'ASP B  10', 'HOH B 478']
  - AC3: ['BINDING SITE FOR RESIDUE MPD A 601', 'ASP A  10', 'ASP A  43', 'GLU A  44', 'HOH A 442']
  - AC4: ['BINDING SITE FOR RESIDUE MPD B 602', 'GLU A  44', 'HOH A 524', 'GLU B  30', 'TRP B  31', 'GLY B  33', 'LYS B  36']
  - AC5: ['BINDING SITE FOR RESIDUE MPD B 603', 'TYR B  70', 'ILE B  72', 'THR B  77', 'THR B  89', 'VAL B  91']
  - AC6: ['BINDING SITE FOR RESIDUE MPD B 604', 'ILE B  60', 'ALA B  67', 'ILE B  72']
  - AC7: ['BINDING SITE FOR RESIDUE MPD A 605', 'MET A  37', 'ILE A  38', 'ALA A  93', 'LEU A  94']
  - AC8: ['BINDING SITE FOR RESIDUE MPD A 606', 'TYR A  70', 'GLY A  71', 'THR A  89', 'VAL A  91']
  - AC9: ['BINDING SITE FOR RESIDUE MPD A 607', 'ILE A  60', 'ALA A  67', 'ILE A  72', 'ARG A  73', 'GLY A  74', 'ILE A  75', 'HOH A 494', 'HOH B 528']

## MCCE.Step1
### Renamed:
  - "CU    CU A 109" to "CU   _CU A 109"

### Termini:
 - <strong>NTR</strong>: "SER A   1"
 - <strong>CTR</strong>: "ALA A 108"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
_CU: https://pubchem.ncbi.nlm.nih.gov/#query=CU&tab=substance; MPD: https://pubchem.ncbi.nlm.nih.gov/#query=MPD&tab=substance; 

### Free Cofactors:
  - Removed all 75 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 75.
  - Species and properties with assigned default values in debug.log:

  - _CUBK: ['VDW_RAD', 'VDW_EPS']

  - MPDBK: ['VDW_RAD', 'VDW_EPS']


### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.57: " CA  NTR A   1" to " CB  SER A   1"
- d= 1.99: " OD1 ASP A   2" to "CU   _CU A 109"

</details>

