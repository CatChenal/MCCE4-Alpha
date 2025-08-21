---
# 1POH :: The 2.0 Angstroms Resolution Structure Of Escherichia Coli Histidine- Containing Phosphocarrier Protein Hpr: A Redetermination;
## PDB.Structure
### Function: PHOSPHOTRANSFERASE
### First Release: 19-OCT-93
### Method: X-RAY DIFFRACTION
### Resolution: 2.00 ANGSTROMS.
### Molecule: HISTIDINE-CONTAINING PHOSPHOCARRIER PROTEIN HPR
### Seqres Species: Residues: A:85; Total: 85
### Cofactors:
  - SO4:
 'SULFATE ION', 2

### Total cofactors: 2
### Total waters: 89
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 9, ARG: 1, ASN: 2, ASP: 1, GLN: 6, GLU: 9, GLY: 6, HIS: 2, ILE: 3, LEU: 8, LYS: 7, MET: 2, PHE: 4, PRO: 2, SER: 6, THR: 10, VAL: 7', 'Total: 85', 'Ionizable: 20',
              'Ratio: 23.5%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 89, 'SO4': 2

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE SO4 A 90', 'HIS A  15', 'THR A  16', 'ARG A  17', 'LYS A  40', 'SER A  41', 'HOH A 113', 'HOH A 136']
  - AC2: ['BINDING SITE FOR RESIDUE SO4 A 91', 'LYS A  27', 'ASN A  38', 'LYS A  40', 'HOH A 130', 'HOH A 168']

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "MET A   1"
 - <strong>CTR</strong>: "GLU A  85"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
SO4: https://pubchem.ncbi.nlm.nih.gov/#query=SO4&tab=substance; 

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 91.
  - Species and properties with assigned default values in debug.log:

  - SO4BK: ['VDW_RAD', 'VDW_EPS']


### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   1" to " CB  MET A   1"

</details>

