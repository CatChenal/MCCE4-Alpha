---
# 1STN :: The Crystal Structure Of Staphylococcal Nuclease Refined At 1.7 Angstroms Resolution;
## PDB.Structure
### Function: HYDROLASE(PHOSPHORIC DIESTER)
### First Release: 17-FEB-93
### Method: X-RAY DIFFRACTION
### Resolution: 1.70 ANGSTROMS.
### Molecule: STAPHYLOCOCCAL NUCLEASE
### Seqres Species: Residues: A:149; Total: 149
### Total waters: 83
### Missing:
  - Residues:
 'A': [('ALA', 1), ('THR', 2), ('SER', 3), ('THR', 4), ('LYS', 5), ('GLU', 142), ('ASP', 143), ('ASN', 144), ('ALA', 145), ('ASP', 146), ('SER', 147), ('GLY', 148), ('GLN', 149),
       ('Count', 13)]

### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 12, ARG: 5, ASN: 5, ASP: 6, GLN: 5, GLU: 11, GLY: 9, HIS: 4, ILE: 5, LEU: 11, LYS: 22, MET: 4, PHE: 3, PRO: 6, SER: 3, THR: 8, TRP: 1, TYR: 7, VAL: 9', 'Total: 136', 'Ionizable: 55',
              'Ratio: 40.4%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 83

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "LYS A   6"
 - <strong>CTR</strong>: "SER A 141"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 83.

### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 141":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   6" to " CB  LYS A   6"

</details>

