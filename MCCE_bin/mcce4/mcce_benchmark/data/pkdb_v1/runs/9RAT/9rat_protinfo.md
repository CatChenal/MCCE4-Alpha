---
# 9RAT :: Effects Of Temperature On Protein Structure And Dynamics: X-Ray Crystallographic Studies Of The Protein Ribonuclease-A At Nine; Different Temperatures From 98 To 320 K;
## PDB.Structure
### Function: HYDROLASE (NUCLEIC ACID,RNA)
### First Release: 13-AUG-91
### Method: X-RAY DIFFRACTION
### Resolution: 1.50 ANGSTROMS.
### Molecule: RIBONUCLEASE A
### Seqres Species: Residues: A:124; Total: 124
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 12, ARG: 4, ASN: 10, ASP: 5, CYS: 8, GLN: 7, GLU: 5, GLY: 3, HIS: 4, ILE: 3, LEU: 2, LYS: 10, MET: 4, PHE: 3, PRO: 4, SER: 15, THR: 10, TYR: 6, VAL: 9', 'Total: 124', 'Ionizable: 42',
              'Ratio: 33.9%')

### Disulfides:
  - CYS A  26 -- 2.13 Å --> CYS A  84
  - CYS A  40 -- 2.06 Å --> CYS A  95
  - CYS A  58 -- 2.14 Å --> CYS A 110
  - CYS A  65 -- 2.16 Å --> CYS A  72

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "LYS A   1"
 - <strong>CTR</strong>: "VAL A 124"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.52: " CA  NTR A   1" to " CB  LYS A   1"

</details>

