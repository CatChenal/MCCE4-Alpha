---
# 1DUK :: Wild-Type Recombinant Sperm Whale Metaquomyoglobin
## PDB.Structure
### Function: OXYGEN STORAGE/TRANSPORT
### First Release: 17-JAN-00
### Method: X-RAY DIFFRACTION
### Resolution: 2.13 ANGSTROMS.
### Molecule: WILD-TYPE RECOMBINANT SPERM WHALE METAQUOMYOGLOBIN, FERRIC PROTOPORPHYRIN IX BOUND DURING BACTERIAL
### Seqres Species: Residues: A:153; Total: 153
### Cofactors:
  - HEM:
 'HEME', 1

### Total cofactors: 1
### Total waters: 12
### Missing:
  - Residues:
 'A': [('GLY', 153), ('Count', 1)]

### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 17, ARG: 4, ASN: 1, ASP: 7, GLN: 5, GLU: 14, GLY: 10, HIS: 12, ILE: 9, LEU: 18, LYS: 19, MET: 2, PHE: 6, PRO: 4, SER: 6, THR: 5, TRP: 2, TYR: 3, VAL: 8', 'Total: 152', 'Ionizable: 59',
              'Ratio: 38.8%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HEM': 1, 'HOH': 12, 'PAA': 1, 'PDD': 1

### Links:
  - NE2 HIS A 93 -- 2.29 Å --> FE  HEM A 155
  - FE  HEM A 155 -- 2.67 Å --> O  HOH A 156

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE HEM A 155', 'LYS A  42', 'PHE A  43', 'ARG A  45', 'VAL A  68', 'ALA A  71', 'SER A  92', 'HIS A  93', 'HIS A  97', 'ILE A  99', 'TYR A 103', 'PHE A 138', 'HOH A 156']

## MCCE.Step1
### Renamed:
  - " CAA HEM A 155" to " CAA PAA A 155"
  - " CAD HEM A 155" to " CAD PDD A 155"
  - " CBA HEM A 155" to " CBA PAA A 155"
  - " CBD HEM A 155" to " CBD PDD A 155"
  - " CGA HEM A 155" to " CGA PAA A 155"
  - " CGD HEM A 155" to " CGD PDD A 155"
  - " O1A HEM A 155" to " O1A PAA A 155"
  - " O1D HEM A 155" to " O1D PDD A 155"
  - " O2A HEM A 155" to " O2A PAA A 155"
  - " O2D HEM A 155" to " O2D PDD A 155"

### Termini:
 - <strong>NTR</strong>: "VAL A   1"
 - <strong>CTR</strong>: "GLN A 152"

### Free Cofactors:
  - Removed all 12 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 12.

### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 152":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.55: " CA  NTR A   1" to " CB  VAL A   1"
- d= 1.54: " C2A HEM A 155" to " CAA PAA A 155"
- d= 1.52: " C3D HEM A 155" to " CAD PDD A 155"

</details>

