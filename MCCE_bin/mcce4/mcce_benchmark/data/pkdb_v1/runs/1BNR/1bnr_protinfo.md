---
# 1BNR :: Barnase
## PDB.Structure
### Function: MICROBIAL RIBONUCLEASE
### First Release: 31-MAR-95
### Method: SOLUTION NMR
### Resolution: NOT APPLICABLE.
### Molecule: BARNASE (G SPECIFIC ENDONUCLEASE), NMR, 20 STRUCTURES
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

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   1" to " CB  ALA A   1"
- d= 1.97: "HG22 VAL A  10" to "HD11 ILE A  76"
- d= 1.99: "HD23 LEU A  42" to "HD13 ILE A  51"
- d= 1.88: "HD12 ILE A  76" to "HG23 ILE A  88"
- d= 1.94: "HD13 ILE A  88" to "HG23 ILE A  96"
- d= 1.86: "HD13 ILE A  88" to "HD13 ILE A 109"

</details>

