---
# 1O73 :: Tryparedoxin From Trypanosoma Brucei
## PDB.Structure
### Function: ELECTRON TRANSPORT
### First Release: 23-OCT-02
### Method: X-RAY DIFFRACTION
### Resolution: 2.28 ANGSTROMS.
### Molecule: TRYPAREDOXIN
### Seqres Species: Residues: A:144; Total: 144
### Total waters: 63
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 10, ARG: 4, ASN: 6, ASP: 7, CYS: 2, GLN: 2, GLU: 9, GLY: 12, HIS: 3, ILE: 7, LEU: 13, LYS: 7, MET: 2, PHE: 9, PRO: 10, SER: 12, THR: 10, TRP: 4, TYR: 5, VAL: 10', 'Total: 144', 'Ionizable: 37',
              'Ratio: 25.7%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 63

### Disulfides:
  - CYS A  40 -- 2.03 Å --> CYS A  43

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "MET A   1"
 - <strong>CTR</strong>: "ASN A 144"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 63.

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.54: " CA  NTR A   1" to " CB  MET A   1"

</details>

