---
# 1IGC :: Igg1 Fab Fragment (Mopc21) Complex With Domain Iii Of Protein G From Streptococcus;
## PDB.Structure
### Function: COMPLEX (ANTIBODY/BINDING PROTEIN)
### First Release: 05-AUG-94
### Method: X-RAY DIFFRACTION
### Resolution: 2.60 ANGSTROMS.
### Molecule: IGG1-KAPPA MOPC21 FAB (LIGHT CHAIN), IGG1-KAPPA MOPC21 FAB (HEAVY CHAIN)
### Seqres Species: Residues: L:213, H:222, A:61; Total: 496
### Total waters: 345
### Missing:
  - Residues:
 'A': [('MET', 1), ('THR', 2), ('PRO', 3), ('Count', 3)]

### Biounits: TRIMERIC: L, H, A;
### Models: 1
### Chains: A, H, L
### Model 1 Residues:
  - L:
 'RESIDUES': ('ALA: 10, ARG: 7, ASN: 12, ASP: 11, CYS: 5, GLN: 9, GLU: 11, GLY: 13, HIS: 3, ILE: 8, LEU: 11, LYS: 13, MET: 4, PHE: 7, PRO: 10, SER: 29, THR: 22, TRP: 3, TYR: 11, VAL: 14', 'Total: 213', 'Ionizable: 61',
              'Ratio: 28.6%')
  - H:
 'RESIDUES': ('ALA: 14, ARG: 8, ASN: 6, ASP: 8, CYS: 5, GLN: 7, GLU: 6, GLY: 19, HIS: 4, ILE: 3, LEU: 15, LYS: 9, MET: 5, PHE: 7, PRO: 16, SER: 32, THR: 21, TRP: 5, TYR: 11, VAL: 21', 'Total: 222', 'Ionizable: 51',
              'Ratio: 23.0%')
  - A:
 'RESIDUES': ('ALA: 7, ASN: 3, ASP: 5, GLN: 1, GLU: 4, GLY: 4, ILE: 1, LEU: 2, LYS: 7, PHE: 2, THR: 12, TRP: 1, TYR: 3, VAL: 6', 'Total: 58', 'Ionizable: 19',
              'Ratio: 32.8%')

### Model 1 Free Cofactors & Waters:
  - L:
 'HOH': 163
  - H:
 'HOH': 145
  - A:
 'HOH': 37

### Disulfides:
  - CYS L  23 -- 2.01 Å --> CYS L  88
  - CYS L 134 -- 2.08 Å --> CYS L 193
  - CYS L 213 -- 2.08 Å --> CYS H 222
  - CYS H  22 -- 2.04 Å --> CYS H  96
  - CYS H 147 -- 1.99 Å --> CYS H 202

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "ASN L   1", "ASP H   1", "ALA A   4"
 - <strong>CTR</strong>: "CYS L 213", "CYS H 222", "GLU A  61"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 345.

### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR L 213":   OXT
  -    Missing heavy atoms for CTR01 in "CTR H 222":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.50: " CA  NTR L   1" to " CB  ASN L   1"
- d= 1.89: " O   VAL L  30" to " O   SER L  67"
- d= 1.84: " O   ILE L 150" to " O   SER L 190"
- d= 1.80: " O   TYR L 191" to " O   SER L 207"
- d= 1.54: " CA  NTR H   1" to " CB  ASP H   1"
- d= 1.99: " SG  CYS H 147" to " SG  CYS H 202"
- d= 1.53: " CA  NTR A   4" to " CB  ALA A   4"

</details>

