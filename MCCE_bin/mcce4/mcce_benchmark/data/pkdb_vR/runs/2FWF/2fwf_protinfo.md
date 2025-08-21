---
# 2FWF :: High Resolution Crystal Structure Of The C-Terminal Domain Of The Electron Transfer Catalyst Dsbd (Reduced Form);
## PDB.Structure
### Function: OXIDOREDUCTASE
### First Release: 02-FEB-06
### Method: X-RAY DIFFRACTION
### Resolution: 1.30 ANGSTROMS.
### Molecule: DISULFIDE INTERCHANGE PROTEIN DSBD
### Seqres Species: Residues: A:134; Total: 134
### Cofactors:
  - IOD:
 'IODIDE ION', 7
  -  NA:
 'SODIUM ION', 2

### Total cofactors: 9
### Total waters: 164
### Missing:
  - Residues:
 'A': [('ALA', 419), ('THR', 420), ('HIS', 421), ('THR', 422), ('ALA', 423), ('GLN', 424), ('THR', 425), ('GLN', 426), ('THR', 427), ('HIS', 551), ('HIS', 552),
       ('Count', 11)]

### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 13, ARG: 3, ASN: 5, ASP: 10, CYS: 2, GLN: 10, GLU: 6, GLY: 5, HIS: 8, ILE: 2, LEU: 15, LYS: 7, MET: 2, PHE: 7, PRO: 5, SER: 2, THR: 8, TRP: 1, TYR: 2, VAL: 10', 'Total: 123', 'Ionizable: 38',
              'Ratio: 30.9%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 164, 'IOD': 7, '_NA': 2

### Links:
  - O  HOH A 37 -- 2.85 Å --> NA  NA A1102
  - O  HOH A 79 -- 2.74 Å --> NA A NA A1101
  - O  HOH A 87 -- 2.80 Å --> NA A NA A1101
  - O  HOH A 131 -- 1.72 Å --> NA  NA A1102
  - OE1 GLN A 442 -- 2.73 Å --> NA  NA A1102

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE IOD A 1001', 'HOH A  33', 'LEU A 429', 'GLU A 468', 'GLN A 488']
  - AC2: ['BINDING SITE FOR RESIDUE IOD A 1002', 'HOH A  68', 'HIS A 523', 'NA A1101']
  - AC3: ['BINDING SITE FOR RESIDUE IOD A 1003', 'HOH A 124', 'GLY A 518', 'GLN A 545', 'HIS A 547']
  - AC4: ['BINDING SITE FOR RESIDUE IOD A 1005', 'HOH A 103', 'ASP A 438']
  - AC5: ['BINDING SITE FOR RESIDUE IOD A 1006', 'HOH A 150', 'GLU A 466']
  - AC6: ['BINDING SITE FOR RESIDUE IOD A 1010', 'HOH A 136', 'HOH A 142', 'GLY A 509', 'LEU A 510']
  - AC7: ['BINDING SITE FOR RESIDUE IOD A 1013', 'HOH A   8', 'TRP A 460', 'GLU A 535']
  - AC8: ['BINDING SITE FOR RESIDUE NA A 1101', 'HOH A  79', 'HOH A  87', 'HIS A 523', 'ARG A 544', 'GLN A 545', 'IOD A1002']
  - AC9: ['BINDING SITE FOR RESIDUE NA A 1102', 'HOH A  37', 'HOH A 131', 'GLN A 442']

## MCCE.Step1
### Renamed:
  - "NA    NA A1101" to "NA   _NA A1101"
  - "NA    NA A1102" to "NA   _NA A1102"

### Termini:
 - <strong>NTR</strong>: "HIS A 428"
 - <strong>CTR</strong>: "HIS A 550"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
IOD: https://pubchem.ncbi.nlm.nih.gov/#query=IOD&tab=substance; 

### Free Cofactors:
  - Removed all 164 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 164.
  - Species and properties with assigned default values in debug.log:

  - IODBK: ['VDW_RAD', 'VDW_EPS']

  - _NA+1: ['VDW_RAD', 'VDW_EPS']


### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 550":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A 428" to " CB  HIS A 428"

</details>

