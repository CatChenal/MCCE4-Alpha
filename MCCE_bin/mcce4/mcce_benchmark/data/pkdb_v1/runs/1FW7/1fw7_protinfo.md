---
# 1FW7 :: Nmr Structure Of 15N-Labeled Barnase
## PDB.Structure
### Function: HYDROLASE
### First Release: 22-SEP-00
### Method: SOLUTION NMR
### Resolution: NOT APPLICABLE.
### Molecule: BARNASE
### Seqres Species: Residues: A:110; Total: 110
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 8, ARG: 6, ASN: 6, ASP: 9, GLN: 4, GLU: 3, GLY: 10, HIS: 2, ILE: 8, LEU: 7, LYS: 8, PHE: 4, PRO: 3, SER: 9, THR: 9, TRP: 3, TYR: 7, VAL: 4', 'Total: 110', 'Ionizable: 35',
              'Ratio: 31.8%')

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "ALA A   1"
 - <strong>CTR</strong>: "ARG A 110"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05

### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 110":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   1" to " CB  ALA A   1"
- d= 1.88: "HD22 ASN A  58" to " O   LEU A  63"
- d= 1.86: " O   ILE A  76" to "HD22 ASN A  84"

</details>

