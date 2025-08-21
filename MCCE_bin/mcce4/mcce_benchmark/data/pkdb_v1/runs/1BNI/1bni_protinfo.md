---
# 1BNI :: Barnase Wildtype Structure At Ph 6.0
## PDB.Structure
### Function: MICROBIAL RIBONUCLEASE
### First Release: 17-MAY-95
### Method: X-RAY DIFFRACTION
### Resolution: 2.10 ANGSTROMS.
### Molecule: BARNASE
### Seqres Species: Residues: A:110, B:110, C:110; Total: 330
### Total waters: 216
### Missing:
  - Residues:
 'A': [('ALA', 1), ('GLN', 2), ('Count', 2)],
 'B': [('ALA', 1), ('GLN', 2), ('Count', 2)],
 'C': [('ALA', 1), ('GLN', 2), ('Count', 2)]
  - ResAtoms:
 'A': ['LYS_19=(CD,CE,NZ)', 'GLU_60=(CG,CD,OE1,OE2)', 'LYS_62=(CE,NZ)', 'GLN_104=(CD,OE1,NE2)'],
 'B': ['ARG_59=(CG,CD,NE,CZ,NH1,NH2)', 'GLU_60=(CG,CD,OE1,OE2)', 'LYS_62=(CG,CD,CE,NZ)', 'LYS_66=(NZ)', 'GLN_104=(CG,CD,OE1,NE2)'],
 'C': ['VAL_3=(N,CB,CG1,CG2)', 'LYS_39=(CE,NZ)', 'ARG_59=(CD,NE,CZ,NH1,NH2)']

### Biounits: MONOMERIC: A,B,C;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 7, ARG: 6, ASN: 6, ASP: 9, GLN: 3, GLU: 3, GLY: 10, HIS: 2, ILE: 8, LEU: 7, LYS: 8, PHE: 4, PRO: 3, SER: 9, THR: 9, TRP: 3, TYR: 7, VAL: 4', 'Total: 108', 'Ionizable: 35',
              'Ratio: 32.4%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 63

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "VAL A   3"
 - <strong>CTR</strong>: "ARG A 110"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 63.

### Missing Heavy Atoms:
  -    Missing heavy atoms for LYS01 in "LYS A  19":   CD ,  CE ,  NZ 
  -    Missing heavy atoms for GLU01 in "GLU A  60":   CG ,  CD ,  OE1,  OE2
  -    Missing heavy atoms for LYS01 in "LYS A  62":   CE ,  NZ 
  -    Missing heavy atoms for GLN01 in "GLN A 104":   CD ,  OE1,  NE2
  -    Missing heavy atoms for CTR01 in "CTR A 110":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.57: " CA  NTR A   3" to " CB  VAL A   3"

</details>

