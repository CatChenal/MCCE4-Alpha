---
# 1BNJ :: Barnase Wildtype Structure At Ph 9.0
## PDB.Structure
### Function: MICROBIAL RIBONUCLEASE
### First Release: 17-MAY-95
### Method: X-RAY DIFFRACTION
### Resolution: 2.10 ANGSTROMS.
### Molecule: BARNASE
### Seqres Species: Residues: A:110, B:110, C:110; Total: 330
### Total waters: 219
### Missing:
  - Residues:
 'A': [('ALA', 1), ('Count', 1)],
 'B': [('ALA', 1), ('GLN', 2), ('Count', 2)],
 'C': [('ALA', 1), ('GLN', 2), ('VAL', 3), ('Count', 3)]
  - ResAtoms:
 'A': ['GLN_2=(N,CA,O,CB,CG,CD,OE1)', 'GLN_2=(NE2)', 'GLU_60=(CG,CD,OE1,OE2)'],
 'B': ['VAL_3=(N,CB,CG1,CG2)', 'ARG_59=(CD,NE,CZ,NH1,NH2)']

### Biounits: MONOMERIC: A,B,C;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 7, ARG: 6, ASN: 6, ASP: 9, GLN: 4, GLU: 3, GLY: 10, HIS: 2, ILE: 8, LEU: 7, LYS: 8, PHE: 4, PRO: 3, SER: 9, THR: 9, TRP: 3, TYR: 7, VAL: 4', 'Total: 109', 'Ionizable: 35',
              'Ratio: 32.1%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 61

## MCCE.Step1
### Termini:
 - <strong>CTR</strong>: "ARG A 110"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 61.

### Missing Heavy Atoms:
  -    Missing heavy atoms for GLU01 in "GLU A  60":   CG ,  CD ,  OE1,  OE2
  -    Missing heavy atoms for CTR01 in "CTR A 110":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- No clash found.

</details>

