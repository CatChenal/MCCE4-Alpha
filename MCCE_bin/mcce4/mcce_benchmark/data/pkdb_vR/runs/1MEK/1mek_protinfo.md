---
# 1MEK :: Human Protein Disulfide Isomerase, Nmr, 40 Structures
## PDB.Structure
### Function: ELECTRON TRANSPORT
### First Release: 16-APR-96
### Method: SOLUTION NMR
### Resolution: NOT APPLICABLE.
### Molecule: PROTEIN DISULFIDE ISOMERASE
### Seqres Species: Residues: A:120; Total: 120
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 20, ARG: 6, ASN: 3, ASP: 7, CYS: 2, GLN: 2, GLU: 12, GLY: 8, HIS: 3, ILE: 3, LEU: 10, LYS: 11, PHE: 4, PRO: 6, SER: 4, THR: 5, TRP: 2, TYR: 6, VAL: 6', 'Total: 120', 'Ionizable: 47',
              'Ratio: 39.2%')

### Disulfides:
  - CYS A  36 -- 2.02 Å --> CYS A  39

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "ASP A   1"
 - <strong>CTR</strong>: "ALA A 120"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   1" to " CB  ASP A   1"

</details>

