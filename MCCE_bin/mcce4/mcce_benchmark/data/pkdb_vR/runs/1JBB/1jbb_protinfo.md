---
# 1JBB :: Ubiquitin Conjugating Enzyme, Ubc13
## PDB.Structure
### Function: LIGASE
### First Release: 03-JUN-01
### Method: X-RAY DIFFRACTION
### Resolution: 2.00 ANGSTROMS.
### Molecule: UBIQUITIN CONJUGATING ENZYME E2-17.5 KDA
### Seqres Species: Residues: A:153, B:153; Total: 306
### Total waters: 131
### Missing:
  - Residues:
 'A': [('MET', 1), ('ALA', 2), ('LYS', 151), ('PRO', 152), ('GLU', 153), ('Count', 5)],
 'B': [('MET', 1), ('ALA', 2), ('LYS', 150), ('LYS', 151), ('PRO', 152), ('GLU', 153), ('Count', 6)]
  - ResAtoms:
 'A': ['LYS_6=(CG,CD,CE,NZ)', 'ASN_31=(CG,OD1,ND2)', 'GLN_36=(CG,CD,OE1,NE2)', 'ARG_82=(CG,CD,NE,CZ,NH1,NH2)', 'ASP_128=(CG,OD1,OD2)', 'LYS_150=(CG,CD,CE,NZ)'],
 'B': ['LYS_6=(CG,CD,CE,NZ)', 'LYS_10=(CG,CD,CE,NZ)', 'LEU_121=(CG,CD1,CD2)']

### Biounits: MONOMERIC: A,B;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 11, ARG: 7, ASN: 7, ASP: 11, CYS: 1, GLN: 5, GLU: 12, GLY: 5, HIS: 2, ILE: 11, LEU: 17, LYS: 12, MET: 1, PHE: 3, PRO: 14, SER: 6, THR: 7, TRP: 3, TYR: 6, VAL: 7', 'Total: 148', 'Ionizable: 51',
              'Ratio: 34.5%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 63

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "SER A   3"
 - <strong>CTR</strong>: "LYS A 150"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 63.

### Missing Heavy Atoms:
  -    Missing heavy atoms for LYS01 in "LYS A   6":   CG ,  CD ,  CE ,  NZ 
  -    Missing heavy atoms for ASN01 in "ASN A  31":   CG ,  OD1,  ND2
  -    Missing heavy atoms for GLN01 in "GLN A  36":   CG ,  CD ,  OE1,  NE2
  -    Missing heavy atoms for ARG01 in "ARG A  82":   CG ,  CD ,  NE ,  CZ ,  NH1,  NH2
  -    Missing heavy atoms for ASP01 in "ASP A 128":   CG ,  OD1,  OD2
  -    Missing heavy atoms for LYS01 in "LYS A 150":   CG ,  CD ,  CE ,  NZ 
  -    Missing heavy atoms for CTR01 in "CTR A 150":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   3" to " CB  SER A   3"

</details>

