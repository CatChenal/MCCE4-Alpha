---
# 1YPI :: Structure Of Yeast Triosephosphate Isomerase At 1.9-Angstroms Resolution;
## PDB.Structure
### Function: ISOMERASE(INTRAMOLECULAR
### First Release: OXIDOREDUCTASE)12-JAN-90
### Method: X-RAY DIFFRACTION
### Resolution: 1.90 ANGSTROMS.
### Molecule: TRIOSEPHOSPHATE ISOMERASE
### Seqres Species: Residues: A:247, B:247; Total: 494
### Total waters: 119
### Biounits: DIMERIC: A, B;
### Models: 1
### Chains: A, B
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 25, ARG: 8, ASN: 12, ASP: 15, CYS: 2, GLN: 7, GLU: 17, GLY: 22, HIS: 3, ILE: 15, LEU: 19, LYS: 21, PHE: 11, PRO: 7, SER: 16, THR: 12, TRP: 3, TYR: 6, VAL: 26', 'Total: 247', 'Ionizable: 72',
              'Ratio: 29.1%')
  - B:
 'RESIDUES': ('ALA: 25, ARG: 8, ASN: 12, ASP: 15, CYS: 2, GLN: 7, GLU: 17, GLY: 22, HIS: 3, ILE: 15, LEU: 19, LYS: 21, PHE: 11, PRO: 7, SER: 16, THR: 12, TRP: 3, TYR: 6, VAL: 26', 'Total: 247', 'Ionizable: 72',
              'Ratio: 29.1%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 66
  - B:
 'HOH': 53

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "ALA A   2", "ALA B   2"
 - <strong>CTR</strong>: "ASN A 248", "ASN B 248"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 119.

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.48: " CA  NTR A   2" to " CB  ALA A   2"
- d= 1.53: " CA  NTR B   2" to " CB  ALA B   2"

</details>

