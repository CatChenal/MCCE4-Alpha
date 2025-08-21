---
# 1D0D :: Crystal Structure Of Tick Anticoagulant Protein Complexed With Bovine Pancreatic Trypsin Inhibitor;
## PDB.Structure
### Function: BLOOD CLOTTING INHIBITOR
### First Release: 09-SEP-99
### Method: X-RAY DIFFRACTION
### Resolution: 1.62 ANGSTROMS.
### Molecule: ANTICOAGULANT PROTEIN, PANCREATIC TRYPSIN INHIBITOR
### Seqres Species: Residues: A:60, B:58; Total: 118
### Cofactors:
  - SO4:
 'SULFATE ION', 3

### Total cofactors: 3
### Total waters: 124
### Biounits: DIMERIC: A, B;
### Models: 1
### Chains: A, B
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 3, ARG: 5, ASN: 4, ASP: 7, CYS: 6, GLU: 4, GLY: 6, HIS: 1, ILE: 4, LEU: 1, LYS: 2, PHE: 3, PRO: 2, SER: 4, THR: 1, TRP: 2, TYR: 5', 'Total: 60', 'Ionizable: 30',
              'Ratio: 50.0%')
  - B:
 'RESIDUES': ('ALA: 6, ARG: 6, ASN: 3, ASP: 2, CYS: 6, GLN: 1, GLU: 2, GLY: 6, ILE: 2, LEU: 2, LYS: 4, MET: 1, PHE: 4, PRO: 4, SER: 1, THR: 3, TYR: 4, VAL: 1', 'Total: 58', 'Ionizable: 24',
              'Ratio: 41.4%')

### Model 1 Free Cofactors & Waters:
  - B:
 'HOH': 50, 'SO4': 3
  - A:
 'HOH': 74

### Disulfides:
  - CYS A  5 -- 2.03 Å --> CYS A  59
  - CYS A  15 -- 2.03 Å --> CYS A  39
  - CYS A  33 -- 2.03 Å --> CYS A  55
  - CYS B  5 -- 2.03 Å --> CYS B  55
  - CYS B  14 -- 2.04 Å --> CYS B  38
  - CYS B  30 -- 2.04 Å --> CYS B  51

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE SO4 B 59', 'ASN A   2', 'TYR A  49', 'ARG B  42', 'HOH B 121', 'HOH B 253', 'HOH B 279']
  - AC2: ['BINDING SITE FOR RESIDUE SO4 B 60', 'LYS B  15', 'ALA B  16', 'ARG B  17', 'HOH B 217', 'HOH B 254']
  - AC3: ['BINDING SITE FOR RESIDUE SO4 B 61', 'ARG B  20', 'TYR B  35', 'ALA B  40']

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "TYR A   1", "ARG B   1"
 - <strong>CTR</strong>: "ILE A  60", "ALA B  58"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
SO4: https://pubchem.ncbi.nlm.nih.gov/#query=SO4&tab=substance; 

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 127.
  - Species and properties with assigned default values in debug.log:

  - SO4BK: ['VDW_RAD', 'VDW_EPS']


### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.52: " CA  NTR A   1" to " CB  TYR A   1"
- d= 1.53: " CA  NTR B   1" to " CB  ARG B   1"

</details>

