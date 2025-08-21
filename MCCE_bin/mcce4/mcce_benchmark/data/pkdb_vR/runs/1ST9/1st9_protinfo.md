---
# 1ST9 :: Crystal Structure Of A Soluble Domain Of Resa In The Oxidised Form
## PDB.Structure
### Function: OXIDOREDUCTASE
### First Release: 25-MAR-04
### Method: X-RAY DIFFRACTION
### Resolution: 1.50 ANGSTROMS.
### Molecule: THIOL-DISULFIDE OXIDOREDUCTASE RESA
### Seqres Species: Residues: A:143, B:143; Total: 286
### Cofactors:
  - EDO:
 'ETHYLENE GLYCOL', 4

### Total cofactors: 4
### Total waters: 434
### Missing:
  - Residues:
 'A': [('SER', 36), ('GLY', 174), ('GLU', 175), ('THR', 176), ('SER', 177), ('GLY', 178), ('Count', 6)],
 'B': [('SER', 36), ('GLU', 37), ('GLY', 38), ('GLY', 174), ('GLU', 175), ('THR', 176), ('SER', 177), ('GLY', 178), ('Count', 8)]
  - ResAtoms:
 'A': ['GLU_37=(CG,CD,OE1,OE2)']

### Biounits: MONOMERIC: A,B;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 5, ARG: 2, ASN: 9, ASP: 8, CYS: 2, GLN: 3, GLU: 9, GLY: 10, HIS: 3, ILE: 6, LEU: 9, LYS: 12, MET: 5, PHE: 8, PRO: 8, SER: 7, THR: 8, TRP: 2, TYR: 5, VAL: 16', 'Total: 137', 'Ionizable: 41',
              'Ratio: 29.9%')

### Model 1 Free Cofactors & Waters:
  - A:
 'EDO': 3, 'HOH': 218

### Disulfides:
  - CYS A  73 -- 2.17 Å --> CYS A  76
  - CYS B  73 -- 2.15 Å --> CYS B  76

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE EDO A 501', 'ALA A  41', 'SER A  57', 'LYS A  60', 'GLY A 149', 'HOH A 515', 'HOH A 582']
  - AC2: ['BINDING SITE FOR RESIDUE EDO A 502', 'LYS A 153', 'HOH A 529', 'HOH A 595', 'TYR B 134', 'ASP B 135', 'VAL B 154', 'HOH B 715']
  - AC3: ['BINDING SITE FOR RESIDUE EDO B 503', 'ALA B  41', 'ASN B  43', 'LEU B  56', 'SER B  57', 'LYS B  60', 'GLY B 149', 'HOH B 589', 'HOH B 661']
  - AC4: ['BINDING SITE FOR RESIDUE EDO A 504', 'ALA A 133', 'TYR A 134', 'VAL A 151', 'HOH A 513', 'HOH A 542', 'HOH A 696', 'HOH A 708']

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "GLU A  37"
 - <strong>CTR</strong>: "PRO A 173"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
EDO: https://pubchem.ncbi.nlm.nih.gov/#query=EDO&tab=substance; 

### Free Cofactors:
  - Removed all 218 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 218.
  - Species and properties with assigned default values in debug.log:

  - EDOBK: ['VDW_RAD', 'VDW_EPS']


### Missing Heavy Atoms:
  -    Missing heavy atoms for GLU01 in "GLU A  37":   CG ,  CD ,  OE1,  OE2
  -    Missing heavy atoms for CTR01 in "CTR A 173":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.50: " CA  NTR A  37" to " CB  GLU A  37"

</details>

