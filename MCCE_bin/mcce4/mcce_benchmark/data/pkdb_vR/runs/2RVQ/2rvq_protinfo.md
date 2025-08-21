---
# 2RVQ :: Solution Structure Of The Isolated Histone H2A-H2B Heterodimer
## PDB.Structure
### Function: NUCLEAR PROTEIN/NUCLEAR PROTEIN
### First Release: 28-MAR-16
### Method: SOLUTION NMR
### Resolution: NOT APPLICABLE.
### Molecule: HISTONE H2A TYPE 1-B/E, HISTONE H2B TYPE 1-J
### Seqres Species: Residues: C:133, D:129; Total: 262
### Biounits: DIMERIC: C, D;
### Models: 1
### Chains: C, D
### Model 1 Residues:
  - C:
 'RESIDUES': ('ALA: 16, ARG: 13, ASN: 6, ASP: 2, GLN: 5, GLU: 7, GLY: 16, HIS: 4, ILE: 6, LEU: 16, LYS: 13, MET: 1, PHE: 1, PRO: 6, SER: 5, THR: 5, TYR: 3, VAL: 8', 'Total: 133', 'Ionizable: 42',
              'Ratio: 31.6%')
  - D:
 'RESIDUES': ('ALA: 14, ARG: 8, ASN: 3, ASP: 3, GLN: 3, GLU: 7, GLY: 9, HIS: 3, ILE: 7, LEU: 6, LYS: 20, MET: 3, PHE: 2, PRO: 7, SER: 13, THR: 8, TYR: 5, VAL: 8', 'Total: 129', 'Ionizable: 46',
              'Ratio: 35.7%')

## MCCE.Step1
### Termini:
 - <strong>NTG</strong>: "GLY C  -3", "GLY D  -3"
 - <strong>CTR</strong>: "LYS C 129", "LYS D 125"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.95: " O   GLY C   4" to "HE22 GLN C   6"
- d= 1.99: "HE22 GLN C  24" to " O   PHE C  25"
- d= 1.84: " OE1 GLU C  91" to "HD22 ASN C  94"
- d= 1.87: " O   ASN C 110" to "HE22 GLN C 112"
- d= 1.97: " O   LYS D  20" to "HE22 GLN D  22"
- d= 1.94: " OG  SER D  38" to "HD21 ASN D  63"
- d= 1.75: "HD22 ASN D  84" to " OE2 GLU D  93"
- d= 1.90: " O   SER D  91" to "HE22 GLN D  95"

</details>

