---
# 2RN2 :: Structural Details Of Ribonuclease H From Escherichia Coli As Refined To An Atomic Resolution;
## PDB.Structure
### Function: HYDROLASE(ENDORIBONUCLEASE)
### First Release: 15-APR-92
### Method: X-RAY DIFFRACTION
### Resolution: 1.48 ANGSTROMS.
### Molecule: RIBONUCLEASE H
### Seqres Species: Residues: A:155; Total: 155
### Total waters: 225
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 14, ARG: 10, ASN: 7, ASP: 7, CYS: 3, GLN: 8, GLU: 12, GLY: 14, HIS: 5, ILE: 7, LEU: 12, LYS: 11, MET: 4, PHE: 2, PRO: 5, SER: 4, THR: 10, TRP: 6, TYR: 5, VAL: 9', 'Total: 155', 'Ionizable: 53',
              'Ratio: 34.2%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 225

### Sites:
  - MG: ['NULL']
  -  MG: ['ASP A  10', 'GLU A  48', 'ASP A  70']

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "MET A   1"
 - <strong>CTR</strong>: "VAL A 155"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 225.

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.51: " CA  NTR A   1" to " CB  MET A   1"

</details>

