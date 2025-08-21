---
# 1BPI :: The Structure Of Bovine Pancreatic Trypsin Inhibitor At 125K: Definition Of Carboxyl-Terminal Residues Glycine-57 And Alanine-58;
## PDB.Structure
### Function: PROTEINASE INHIBITOR (TRYPSIN)
### First Release: 18-FEB-95
### Method: X-RAY DIFFRACTION
### Resolution: 1.09 ANGSTROMS.
### Molecule: BOVINE PANCREATIC TRYPSIN INHIBITOR
### Seqres Species: Residues: A:58; Total: 58
### Cofactors:
  - PO4:
 'PHOSPHATE ION', 1

### Total cofactors: 1
### Total waters: 167
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 6, ARG: 6, ASN: 3, ASP: 2, CYS: 6, GLN: 1, GLU: 2, GLY: 6, ILE: 2, LEU: 2, LYS: 4, MET: 1, PHE: 4, PRO: 4, SER: 1, THR: 3, TYR: 4, VAL: 1', 'Total: 58', 'Ionizable: 24',
              'Ratio: 41.4%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 167, 'PO4': 1

### Disulfides:
  - CYS A  5 -- 2.01 Å --> CYS A  55
  - CYS A  14 -- 2.03 Å --> CYS A  38
  - CYS A  30 -- 2.02 Å --> CYS A  51

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE PO4 A 59', 'ARG A   1', 'LYS A  15', 'ARG A  20', 'TYR A  35', 'LYS A  41', 'LYS A  46', 'ALA A  58', 'HOH A 135', 'HOH A 142', 'HOH A 146', 'HOH A 151', 'HOH A 167', 'HOH A 192', 'HOH A 194', 'HOH A 195', 'HOH A 196', 'HOH A 197', 'HOH A 214', 'HOH A 215']

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "ARG A   1"
 - <strong>CTR</strong>: "ALA A  58"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
PO4: https://pubchem.ncbi.nlm.nih.gov/#query=PO4&tab=substance; 

### Free Cofactors:
  - Removed all 167 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 167.
  - Species and properties with assigned default values in debug.log:

  - PO4BK: ['VDW_RAD', 'VDW_EPS']


### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.52: " CA  NTR A   1" to " CB  ARG A   1"

</details>

