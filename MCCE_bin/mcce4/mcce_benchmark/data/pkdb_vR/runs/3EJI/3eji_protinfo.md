---
# 3EJI :: Crystal Structure Of Staphylococcal Nuclease Variant Delta+Phs L36K At Cryogenic Temperature;
## PDB.Structure
### Function: HYDROLASE
### First Release: 18-SEP-08
### Method: X-RAY DIFFRACTION
### Resolution: 1.90 ANGSTROMS.
### Molecule: THERMONUCLEASE
### Seqres Species: Residues: A:143; Total: 143
### Cofactors:
  - PO4:
 'PHOSPHATE ION', 1
  - MPD:
 '(4S)-2-METHYL-2,4-PENTANEDIOL', 1

### Total cofactors: 2
### Total waters: 71
### Missing:
  - Residues:
 'A': [('ALA', 1), ('THR', 2), ('SER', 3), ('THR', 4), ('LYS', 5), ('GLU', 142), ('ASP', 143), ('ASN', 144), ('ALA', 145), ('ASP', 146), ('SER', 147), ('GLY', 148), ('GLN', 149),
       ('Count', 13)]

### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 13, ARG: 5, ASN: 6, ASP: 6, GLN: 5, GLU: 11, GLY: 9, HIS: 2, ILE: 5, LEU: 11, LYS: 20, MET: 4, PHE: 4, PRO: 4, SER: 2, THR: 7, TRP: 1, TYR: 7, VAL: 8', 'Total: 130', 'Ionizable: 51',
              'Ratio: 39.2%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 71, 'MPD': 1, 'PO4': 1

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE PO4 A 201', 'ARG A  35', 'LYS A  71', 'ARG A  87', 'HOH A 356']
  - AC2: ['BINDING SITE FOR RESIDUE MPD A 301', 'GLU A  10', 'TYR A  27', 'LYS A  28', 'LYS A 134', 'LYS A 136']

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "LYS A   6"
 - <strong>CTR</strong>: "SER A 141"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
PO4: https://pubchem.ncbi.nlm.nih.gov/#query=PO4&tab=substance; MPD: https://pubchem.ncbi.nlm.nih.gov/#query=MPD&tab=substance; 

### Free Cofactors:
  - Removed all 71 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 71.
  - Species and properties with assigned default values in debug.log:

  - PO4BK: ['VDW_RAD', 'VDW_EPS']

  - MPDBK: ['VDW_RAD', 'VDW_EPS']


### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 141":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   6" to " CB  LYS A   6"

</details>

