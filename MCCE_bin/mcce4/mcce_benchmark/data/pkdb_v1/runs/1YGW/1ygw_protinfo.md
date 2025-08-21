---
# 1YGW :: Nmr Structure Of Ribonuclease T1, 34 Structures
## PDB.Structure
### Function: HYDROLASE
### First Release: 28-SEP-96
### Method: SOLUTION NMR
### Resolution: NOT APPLICABLE.
### Molecule: RIBONUCLEASE T1
### Seqres Species: Residues: A:104; Total: 104
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 7, ARG: 1, ASN: 9, ASP: 6, CYS: 4, GLN: 2, GLU: 6, GLY: 12, HIS: 3, ILE: 2, LEU: 3, LYS: 2, PHE: 4, PRO: 4, SER: 15, THR: 6, TRP: 1, TYR: 9, VAL: 8', 'Total: 104', 'Ionizable: 31',
              'Ratio: 29.8%')

### Disulfides:
  - CYS A  2 -- 2.10 Å --> CYS A  10
  - CYS A  6 -- 1.15 Å --> CYS A 103

### Sites:
  - REC: ['GUANINE RECOGNITION SITE.', 'TYR A  42', 'ASN A  43', 'ASN A  44', 'TYR A  45', 'GLU A  46', 'ASN A  98']
  - CAT: ['CATALYTIC SITE.', 'TYR A  38', 'HIS A  40', 'GLU A  58', 'ARG A  77', 'HIS A  92']

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "ALA A   1"
 - <strong>CTR</strong>: "THR A 104"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05

### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 104":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   1" to " CB  ALA A   1"
- d= 1.15: " SG  CYS A   6" to " SG  CYS A 103"
- d= 1.93: "HG11 VAL A  16" to "HG13 VAL A  78"
- d= 1.87: "HG12 VAL A  16" to "HG11 VAL A  78"
- d= 2.00: "HG12 VAL A  16" to "HD11 LEU A  86"
- d= 1.99: "HG13 VAL A  16" to "HD11 LEU A  86"
- d= 1.97: "HG23 VAL A  16" to "HG21 VAL A  89"
- d= 1.93: "HD21 LEU A  26" to "HG12 VAL A  33"
- d= 1.92: "HD23 LEU A  26" to "HG11 VAL A  33"
- d= 1.94: "HD23 LEU A  26" to "HG12 VAL A  33"
- d= 1.97: "HG22 VAL A  52" to " OD1 ASN A  81"
- d= 1.98: "HD11 ILE A  61" to "HG12 VAL A  78"
- d= 1.95: "HD13 ILE A  61" to "HG13 VAL A  78"
- d= 1.93: "HG21 VAL A  79" to "HD12 ILE A  90"
- d= 1.98: "HG22 VAL A  79" to "HD12 ILE A  90"
- d= 1.74: " ND2 ASN A  81" to "HD22 ASN A  83"
- d= 1.56: "HD21 ASN A  81" to "HD22 ASN A  83"
- d= 1.76: "HD21 ASN A  83" to "HE21 GLN A  85"

</details>

