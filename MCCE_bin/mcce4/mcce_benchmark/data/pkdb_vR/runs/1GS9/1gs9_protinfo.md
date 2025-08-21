---
# 1GS9 :: Apolipoprotein E4, 22K Domain
## PDB.Structure
### Function: LIPID BINDING PROTEIN
### First Release: 02-JAN-02
### Method: X-RAY DIFFRACTION
### Resolution: 1.70 ANGSTROMS.
### Molecule: APOLIPOPROTEIN E
### Seqres Species: Residues: A:165; Total: 165
### Total waters: 223
### Missing:
  - Residues:
 'A': [('LYS', 1), ('VAL', 2), ('GLU', 3), ('GLN', 4), ('ALA', 5), ('VAL', 6), ('GLU', 7), ('THR', 8), ('GLU', 9), ('PRO', 10), ('GLU', 11), ('PRO', 12), ('GLU', 13), ('LEU', 14), ('ARG', 15), ('GLN', 16), ('GLN', 17), ('THR', 18), ('GLU', 19), ('TRP', 20), ('GLN', 21),
       ('Count', 21)]

### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 14, ARG: 17, ASP: 7, GLN: 14, GLU: 17, GLY: 7, HIS: 1, LEU: 24, LYS: 7, MET: 4, PHE: 1, PRO: 1, SER: 8, THR: 6, TRP: 3, TYR: 4, VAL: 9', 'Total: 144', 'Ionizable: 53',
              'Ratio: 36.8%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 223

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "SER A  22"
 - <strong>CTR</strong>: "GLY A 165"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 223.

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.54: " CA  NTR A  22" to " CB  SER A  22"

</details>

