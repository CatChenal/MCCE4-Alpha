---
# 1BI6 :: Nmr Structure Of Bromelain Inhibitor Vi From Pineapple Stem
## PDB.Structure
### Function: CYSTEINE PROTEASE INHIBITOR
### First Release: 07-DEC-95
### Method: SOLUTION NMR
### Resolution: NOT APPLICABLE.
### Molecule: BROMELAIN INHIBITOR VI
### Seqres Species: Residues: L:11, H:41; Total: 52
### Biounits: DIMERIC: L, H;
### Models: 1
### Chains: H, L
### Model 1 Residues:
  - L:
 'RESIDUES': ('ALA: 1, ARG: 1, CYS: 3, GLU: 1, LEU: 1, PRO: 1, SER: 1, THR: 1, VAL: 1', 'Total: 11', 'Ionizable: 5',
              'Ratio: 45.5%')
  - H:
 'RESIDUES': ('ALA: 1, ASN: 1, ASP: 4, CYS: 7, GLU: 3, GLY: 2, ILE: 2, LEU: 2, LYS: 5, PHE: 2, PRO: 2, SER: 2, THR: 3, TYR: 4, VAL: 1', 'Total: 41', 'Ionizable: 23',
              'Ratio: 56.1%')

### Disulfides:
  - CYS L  3 -- 2.02 Å --> CYS H  7
  - CYS L  6 -- 2.02 Å --> CYS H  39
  - CYS L  8 -- 2.02 Å --> CYS H  5
  - CYS H  14 -- 2.02 Å --> CYS H  21
  - CYS H  18 -- 2.02 Å --> CYS H  30

### Sites:
  - B1: ['NULL']
  -  B1: ['LEU L  10', 'ARG L  11']

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "THR L   1", "GLU H   1"
 - <strong>CTR</strong>: "ARG L  11", "LYS H  41"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR L   1" to " CB  THR L   1"
- d= 1.98: "HG23 VAL L   7" to "HG22 THR H   8"
- d= 1.53: " CA  NTR H   1" to " CB  GLU H   1"

</details>

