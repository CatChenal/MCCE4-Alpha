---
# 2QMT :: Crystal Polymorphism Of Protein Gb1 Examined By Solid-State Nmr And X- Ray Diffraction;
## PDB.Structure
### Function: IMMUNE SYSTEM
### First Release: 16-JUL-07
### Method: X-RAY DIFFRACTION
### Resolution: 1.05 ANGSTROMS.
### Molecule: IMMUNOGLOBULIN G-BINDING PROTEIN G
### Seqres Species: Residues: A:56; Total: 56
### Cofactors:
  - PO4:
 'PHOSPHATE ION', 1
  - MRD:
 '(4R)-2-METHYLPENTANE-2,4-DIOL', 1
  - IPA:
 '2-PROPANOL', 2

### Total cofactors: 4
### Total waters: 135
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 6, ASN: 3, ASP: 5, GLN: 2, GLU: 5, GLY: 4, ILE: 1, LEU: 3, LYS: 6, MET: 1, PHE: 2, THR: 10, TRP: 1, TYR: 3, VAL: 4', 'Total: 56', 'Ionizable: 19',
              'Ratio: 33.9%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 135, 'IPA': 2, 'MRD': 1, 'PO4': 1

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE PO4 A 101', 'TYR A   3', 'ASP A  22', 'ALA A  23', 'ASP A  47', 'LYS A  50', 'HOH A1058', 'HOH A1073', 'HOH A1077']
  - AC2: ['BINDING SITE FOR RESIDUE MRD A 102', 'ALA A  20', 'VAL A  21', 'GLU A  27', 'LYS A  31', 'TRP A  43', 'HOH A1021', 'HOH A1089']
  - AC3: ['BINDING SITE FOR RESIDUE IPA A 203', 'ALA A  24', 'GLU A  27', 'TYR A  45', 'PHE A  52', 'HOH A1030']
  - AC4: ['BINDING SITE FOR RESIDUE IPA A 204', 'ASN A   8', 'THR A  55', 'HOH A1094', 'HOH A1147']

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "MET A   1"
 - <strong>CTR</strong>: "GLU A  56"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
PO4: https://pubchem.ncbi.nlm.nih.gov/#query=PO4&tab=substance; MRD: https://pubchem.ncbi.nlm.nih.gov/#query=MRD&tab=substance; IPA: https://pubchem.ncbi.nlm.nih.gov/#query=IPA&tab=substance; 

### Free Cofactors:
  - Removed all 135 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 135.
  - Species and properties with assigned default values in debug.log:

  - PO4BK: ['VDW_RAD', 'VDW_EPS']

  - MRDBK: ['VDW_RAD', 'VDW_EPS']

  - IPABK: ['VDW_RAD', 'VDW_EPS']


### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.52: " CA  NTR A   1" to " CB  MET A   1"

</details>

