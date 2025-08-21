---
# 1SU9 :: Reduced Structure Of The Soluble Domain Of Resa
## PDB.Structure
### Function: OXIDOREDUCTASE
### First Release: 26-MAR-04
### Method: X-RAY DIFFRACTION
### Resolution: 1.95 ANGSTROMS.
### Molecule: THIOL-DISULFIDE OXIDOREDUCTASE RESA
### Seqres Species: Residues: A:143, B:143; Total: 286
### Total waters: 300
### Missing:
  - Residues:
 'A': [('SER', 36), ('GLU', 37), ('GLY', 38), ('THR', 176), ('SER', 177), ('GLY', 178), ('Count', 6)],
 'B': [('SER', 36), ('GLU', 37), ('GLY', 38), ('GLU', 175), ('THR', 176), ('SER', 177), ('GLY', 178), ('Count', 7)]

### Biounits: MONOMERIC: A,B;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 5, ARG: 2, ASN: 9, ASP: 8, CYS: 2, GLN: 3, GLU: 9, GLY: 10, HIS: 3, ILE: 6, LEU: 9, LYS: 12, MET: 5, PHE: 8, PRO: 8, SER: 7, THR: 8, TRP: 2, TYR: 5, VAL: 16', 'Total: 137', 'Ionizable: 41',
              'Ratio: 29.9%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 146

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "SER A  39"
 - <strong>CTR</strong>: "GLU A 175"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 146.

### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 175":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A  39" to " CB  SER A  39"

</details>

