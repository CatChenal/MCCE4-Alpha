---
# 1BEG :: Structure Of Fungal Elicitor, Nmr, 18 Structures
## PDB.Structure
### Function: SIGNAL
### First Release: 26-NOV-96
### Method: SOLUTION NMR
### Resolution: NOT APPLICABLE.
### Molecule: BETA-ELICITIN CRYPTOGEIN
### Seqres Species: Residues: A:98; Total: 98
### Missing:
  - ResAtoms:
 'A': ['LEU_98=(O)']

### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 11, ASN: 7, ASP: 3, CYS: 6, GLN: 4, GLY: 3, ILE: 3, LEU: 10, LYS: 6, MET: 3, PHE: 2, PRO: 4, SER: 12, THR: 14, TYR: 5, VAL: 5', 'Total: 98', 'Ionizable: 20',
              'Ratio: 20.4%')

### Disulfides:
  - CYS A  3 -- 2.02 Å --> CYS A  71
  - CYS A  27 -- 2.02 Å --> CYS A  56
  - CYS A  51 -- 2.02 Å --> CYS A  95

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "THR A   1"
 - <strong>CTR</strong>: "LEU A  98"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05

### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A  98":   O  ,  OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.54: " CA  NTR A   1" to " CB  THR A   1"
- d= 1.97: "HG23 THR A  14" to "HD11 ILE A  63"
- d= 1.95: "HG21 THR A  14" to "HD22 LEU A  66"
- d= 1.91: " O   CYS A  51" to "HD22 ASN A  57"
- d= 1.80: "HG23 THR A  74" to "HG23 VAL A  81"

</details>

