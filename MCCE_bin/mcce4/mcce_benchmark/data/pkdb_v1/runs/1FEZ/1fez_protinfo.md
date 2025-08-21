---
# 1FEZ :: The Crystal Structure Of Bacillus Cereus Phosphonoacetaldehyde Hydrolase Complexed With Tungstate, A Product Analog;
## PDB.Structure
### Function: HYDROLASE
### First Release: 24-JUL-00
### Method: X-RAY DIFFRACTION
### Resolution: 3.00 ANGSTROMS.
### Molecule: PHOSPHONOACETALDEHYDE HYDROLASE
### Seqres Species: Residues: A:256, B:256, C:256, D:256; Total: 1024
### Cofactors:
  -  MG:
 'MAGNESIUM ION', 4
  - WO4:
 'TUNGSTATE(VI)ION', 2

### Total cofactors: 6
### Total waters: 16
### Biounits: DIMERIC: A, B,C, D;
### Models: 1
### Chains: A, B
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 20, ARG: 15, ASN: 8, ASP: 11, CYS: 2, GLN: 4, GLU: 33, GLY: 17, HIS: 5, ILE: 19, LEU: 16, LYS: 12, MET: 15, PHE: 10, PRO: 12, SER: 9, THR: 14, TRP: 4, TYR: 8, VAL: 22', 'Total: 256', 'Ionizable: 86',
              'Ratio: 33.6%')
  - B:
 'RESIDUES': ('ALA: 20, ARG: 15, ASN: 8, ASP: 11, CYS: 2, GLN: 4, GLU: 33, GLY: 17, HIS: 5, ILE: 19, LEU: 16, LYS: 12, MET: 15, PHE: 10, PRO: 12, SER: 9, THR: 14, TRP: 4, TYR: 8, VAL: 22', 'Total: 256', 'Ionizable: 86',
              'Ratio: 33.6%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 4, '_MG': 1
  - B:
 'HOH': 4, 'WO4': 1, '_MG': 1

### Links:
  - OD1 ASP A 12 -- 2.92 Å --> MG  MG A 801
  - O  ALA A 14 -- 2.53 Å --> MG  MG A 801
  - OD1 ASP A 186 -- 2.52 Å --> MG  MG A 801
  - OD2 ASP A 186 -- 2.70 Å --> MG  MG A 801
  - O  ALA B 14 -- 2.31 Å --> MG  MG B 802
  - OD1 ASP B 186 -- 2.69 Å --> MG  MG B 802
  - OD2 ASP B 186 -- 2.84 Å --> MG  MG B 802
  - OG SER B 209 -- 2.99 Å --> MG  MG B 802
  - O  HOH B 405 -- 2.85 Å --> MG  MG B 802
  - O3 WO4 B 800 -- 2.88 Å --> MG  MG B 802
  - OD1 ASP C 12 -- 2.87 Å --> MG  MG C 804
  - O  ALA C 14 -- 2.56 Å --> MG  MG C 804
  - OD1 ASP C 186 -- 2.44 Å --> MG  MG C 804
  - OD2 ASP C 186 -- 2.65 Å --> MG  MG C 804
  - O  ALA D 14 -- 2.47 Å --> MG  MG D 805
  - OD1 ASP D 186 -- 2.70 Å --> MG  MG D 805
  - OD2 ASP D 186 -- 2.81 Å --> MG  MG D 805
  - OG SER D 209 -- 2.88 Å --> MG  MG D 805
  - O  HOH D 413 -- 2.79 Å --> MG  MG D 805
  - O2 WO4 D 803 -- 2.90 Å --> MG  MG D 805

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE WO4 B 800', 'ASP B  12', 'TRP B  13', 'ALA B  14', 'THR B 126', 'GLY B 127', 'TYR B 128', 'ARG B 160', 'HOH B 405', 'HOH B 407', 'MG B 802']
  - AC2: ['BINDING SITE FOR RESIDUE MG A 801', 'ASP A  12', 'ALA A  14', 'ASP A 186']
  - AC3: ['BINDING SITE FOR RESIDUE MG B 802', 'ASP B  12', 'ALA B  14', 'GLY B  15', 'ASP B 186', 'THR B 187', 'SER B 209', 'HOH B 405', 'WO4 B 800']
  - AC4: ['BINDING SITE FOR RESIDUE WO4 D 803', 'ASP D  12', 'TRP D  13', 'ALA D  14', 'GLY D  50', 'THR D 126', 'GLY D 127', 'TYR D 128', 'ARG D 160', 'HOH D 413', 'HOH D 415', 'MG D 805']
  - AC5: ['BINDING SITE FOR RESIDUE MG C 804', 'ASP C  12', 'ALA C  14', 'ASP C 186']
  - AC6: ['BINDING SITE FOR RESIDUE MG D 805', 'ASP D  12', 'ALA D  14', 'GLY D  15', 'ASP D 186', 'THR D 187', 'SER D 209', 'HOH D 413', 'WO4 D 803']

## MCCE.Step1
### Renamed:
  - "MG    MG A 801" to "MG   _MG A 801"
  - "MG    MG B 802" to "MG   _MG B 802"

### Termini:
 - <strong>NTR</strong>: "LYS A   5", "LYS B   5"
 - <strong>CTR</strong>: "GLU A 260", "GLU B 260"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
_MG: https://pubchem.ncbi.nlm.nih.gov/#query=MG&tab=substance; WO4: https://pubchem.ncbi.nlm.nih.gov/#query=WO4&tab=substance; 

### Free Cofactors:
  - Removed all 4 HOH in A. Removed all 4 HOH in B.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 8.
  - Species and properties with assigned default values in debug.log:

  - _MGBK: ['VDW_RAD', 'VDW_EPS']

  - WO4BK: ['VDW_RAD', 'VDW_EPS']


### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.52: " CA  NTR A   5" to " CB  LYS A   5"
- d= 1.54: " CA  NTR B   5" to " CB  LYS B   5"

</details>

