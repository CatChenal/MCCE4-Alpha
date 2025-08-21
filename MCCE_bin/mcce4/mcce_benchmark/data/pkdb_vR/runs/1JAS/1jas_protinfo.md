---
# 1JAS :: Hsubc2B
## PDB.Structure
### Function: LIGASE
### First Release: 31-MAY-01
### Method: SOLUTION NMR
### Resolution: NOT APPLICABLE.
### Molecule: UBIQUITIN-CONJUGATING ENZYME E2-17 KDA
### Seqres Species: Residues: A:152; Total: 152
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 8, ARG: 9, ASN: 11, ASP: 8, CYS: 1, GLN: 8, GLU: 12, GLY: 6, HIS: 1, ILE: 8, LEU: 10, LYS: 6, MET: 4, PHE: 7, PRO: 14, SER: 16, THR: 6, TRP: 3, TYR: 5, VAL: 9', 'Total: 152', 'Ionizable: 42',
              'Ratio: 27.6%')

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "MET A   1"
 - <strong>CTR</strong>: "SER A 152"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05

### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 152":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.54: " CA  NTR A   1" to " CB  MET A   1"
- d= 1.85: "HE21 GLN A  35" to " OE2 GLU A  58"
- d= 1.96: " O   TRP A  36" to "HD21 ASN A  37"
- d= 1.85: "HD21 ASN A  80" to " O   CYS A  88"
- d= 1.87: "HE21 GLN A 128" to " OE2 GLU A 132"
- d= 1.83: " O   ALA A 143" to "HE21 GLN A 147"

</details>

