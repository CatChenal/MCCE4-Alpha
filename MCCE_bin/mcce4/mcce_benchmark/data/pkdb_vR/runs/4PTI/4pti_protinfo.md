---
# 4PTI :: The Geometry Of The Reactive Site And Of The Peptide Groups In Trypsin, Trypsinogen And Its Complexes With Inhibitors;
## PDB.Structure
### Function: PROTEINASE INHIBITOR (TRYPSIN)
### First Release: 27-SEP-82
### Method: X-RAY DIFFRACTION
### Resolution: 1.50 ANGSTROMS.
### Molecule: TRYPSIN INHIBITOR
### Seqres Species: Residues: A:58; Total: 58
### Total waters: 60
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 6, ARG: 6, ASN: 3, ASP: 2, CYS: 6, GLN: 1, GLU: 2, GLY: 6, ILE: 2, LEU: 2, LYS: 4, MET: 1, PHE: 4, PRO: 4, SER: 1, THR: 3, TYR: 4, VAL: 1', 'Total: 58', 'Ionizable: 24',
              'Ratio: 41.4%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 60

### Disulfides:
  - CYS A  5 -- 2.05 Å --> CYS A  55
  - CYS A  14 -- 2.09 Å --> CYS A  38
  - CYS A  30 -- 2.02 Å --> CYS A  51

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "ARG A   1"
 - <strong>CTR</strong>: "ALA A  58"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 60.

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.55: " CA  NTR A   1" to " CB  ARG A   1"

</details>

