---
# 2FD3 :: Crystal Structure Of Thioredoxin Mutant P34H
## PDB.Structure
### Function: ELECTRON TRANSPORT
### First Release: 13-DEC-05
### Method: X-RAY DIFFRACTION
### Resolution: 2.45 ANGSTROMS.
### Molecule: THIOREDOXIN 1
### Seqres Species: Residues: A:108, B:108; Total: 216
### Total waters: 91
### Biounits: MONOMERIC: A,B;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 12, ARG: 1, ASN: 4, ASP: 11, CYS: 2, GLN: 3, GLU: 5, GLY: 9, HIS: 2, ILE: 9, LEU: 13, LYS: 10, MET: 1, PHE: 4, PRO: 4, SER: 3, THR: 6, TRP: 2, TYR: 2, VAL: 5', 'Total: 108', 'Ionizable: 33',
              'Ratio: 30.6%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 37

### Disulfides:
  - CYS A  32 -- 2.04 Å --> CYS A  35
  - CYS B  32 -- 2.04 Å --> CYS B  35

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "SER A   1"
 - <strong>CTR</strong>: "ALA A 108"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 37.

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   1" to " CB  SER A   1"

</details>

