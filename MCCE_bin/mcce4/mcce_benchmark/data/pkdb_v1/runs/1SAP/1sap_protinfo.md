---
# 1SAP :: Hyperthermophile Protein, Relaxation Matrix Refinement Structure
## PDB.Structure
### Function: DNA BINDING PROTEIN
### First Release: 25-APR-95
### Method: SOLUTION NMR
### Resolution: NOT APPLICABLE.
### Molecule: SAC7D
### Seqres Species: Residues: A:66; Total: 66
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 4, ARG: 4, ASN: 1, ASP: 5, GLU: 7, GLY: 5, ILE: 1, LEU: 3, LYS: 14, MET: 3, PHE: 2, PRO: 1, SER: 3, THR: 3, TRP: 1, TYR: 2, VAL: 7', 'Total: 66', 'Ionizable: 32',
              'Ratio: 48.5%')

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "MET A   1"
 - <strong>CTR</strong>: "LYS A  66"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05

### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A  66":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.54: " CA  NTR A   1" to " CB  MET A   1"

</details>

