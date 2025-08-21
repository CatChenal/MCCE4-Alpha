---
# 1EH6 :: Human O6-Alkylguanine-Dna Alkyltransferase
## PDB.Structure
### Function: TRANSFERASE
### First Release: 18-FEB-00
### Method: X-RAY DIFFRACTION
### Resolution: 2.00 ANGSTROMS.
### Molecule: O6-ALKYLGUANINE-DNA ALKYLTRANSFERASE
### Seqres Species: Residues: A:207; Total: 207
### Cofactors:
  -  ZN:
 'ZINC ION', 1

### Total cofactors: 1
### Total waters: 153
### Missing:
  - Residues:
 'A': [('MET', 1), ('ASP', 2), ('LYS', 3), ('ASP', 4), ('LYS', 36), ('GLY', 37), ('THR', 38), ('SER', 39), ('ALA', 40), ('ALA', 41), ('ASP', 42), ('ALA', 43), ('VAL', 44), ('GLY', 182), ('GLY', 183), ('SER', 184), ('SER', 185), ('GLY', 186), ('LEU', 187), ('ALA', 188), ('GLY', 189), ('ALA', 190), ('TRP', 191), ('LEU', 192), ('LYS', 193), ('GLY', 194), ('ALA', 195), ('GLY', 196), ('ALA', 197), ('THR', 198), ('SER', 199), ('GLY', 200), ('SER', 201), ('PRO', 202), ('PRO', 203), ('ALA', 204), ('GLY', 205), ('ARG', 206), ('ASN', 207),
       ('Count', 39)]

### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 17, ARG: 6, ASN: 4, ASP: 1, CYS: 5, GLN: 8, GLU: 13, GLY: 18, HIS: 7, ILE: 5, LEU: 22, LYS: 9, MET: 3, PHE: 5, PRO: 14, SER: 7, THR: 4, TRP: 3, TYR: 3, VAL: 14', 'Total: 168', 'Ionizable: 44',
              'Ratio: 26.2%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 153, '_ZN': 1

### Links:
  - SG CYS A  5 -- 2.47 Å --> ZN  ZN A 208
  - SG CYS A 24 -- 2.29 Å --> ZN  ZN A 208
  - NE2 HIS A 29 -- 2.05 Å --> ZN  ZN A 208
  - ND1 HIS A 85 -- 2.05 Å --> ZN  ZN A 208

### Sites:
  - ACT: ['ACTIVE SITE CYSTEINE', 'CYS A 145']
  - AC1: ['BINDING SITE FOR RESIDUE ZN A 208', 'CYS A   5', 'CYS A  24', 'HIS A  29', 'HIS A  85']

## MCCE.Step1
### Renamed:
  - "ZN    ZN A 208" to "ZN   _ZN A 208"

### Termini:
 - <strong>NTR</strong>: "CYS A   5", "GLU A  45"
 - <strong>CTR</strong>: "GLY A  35", "LEU A 181"

### Free Cofactors:
  - Removed all 153 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 153.

### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A  35":   OXT
  -    Missing heavy atoms for CTR01 in "CTR A 181":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.54: " CA  NTR A   5" to " CB  CYS A   5"
- d= 1.53: " CA  NTR A  45" to " CB  GLU A  45"

</details>

