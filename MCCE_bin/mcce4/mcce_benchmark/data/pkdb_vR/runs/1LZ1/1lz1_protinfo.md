---
# 1LZ1 :: Refinement Of Human Lysozyme At 1.5 Angstroms Resolution. Analysis Of Non-Bonded And Hydrogen-Bond Interactions;
## PDB.Structure
### Function: HYDROLASE (O-GLYCOSYL)
### First Release: 12-OCT-84
### Method: X-RAY DIFFRACTION
### Resolution: 1.50 ANGSTROMS.
### Molecule: HUMAN LYSOZYME
### Seqres Species: Residues: A:130; Total: 130
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 14, ARG: 14, ASN: 10, ASP: 8, CYS: 8, GLN: 6, GLU: 3, GLY: 11, HIS: 1, ILE: 5, LEU: 8, LYS: 5, MET: 2, PHE: 2, PRO: 2, SER: 6, THR: 5, TRP: 5, TYR: 6, VAL: 9', 'Total: 130', 'Ionizable: 45',
              'Ratio: 34.6%')

### Disulfides:
  - CYS A  6 -- 2.07 Å --> CYS A 128
  - CYS A  30 -- 2.06 Å --> CYS A 116
  - CYS A  65 -- 2.08 Å --> CYS A  81
  - CYS A  77 -- 2.04 Å --> CYS A  95

### Sites:
  - 51: ['NULL']
  - 52: ['NULL']
  -  51: ['GLU A  35']
  -  52: ['ASP A  53']

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "LYS A   1"
 - <strong>CTR</strong>: "VAL A 130"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.51: " CA  NTR A   1" to " CB  LYS A   1"

</details>

