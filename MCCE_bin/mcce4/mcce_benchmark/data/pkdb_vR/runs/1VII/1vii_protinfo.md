---
# 1VII :: Thermostable Subdomain From Chicken Villin Headpiece, Nmr, Minimized Average Structure;
## PDB.Structure
### Function: ACTIN BINDING
### First Release: 15-JAN-97
### Method: SOLUTION NMR
### Resolution: NOT APPLICABLE.
### Molecule: VILLIN
### Seqres Species: Residues: A:36; Total: 36
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 3, ARG: 1, ASN: 2, ASP: 2, GLN: 2, GLU: 2, GLY: 2, LEU: 5, LYS: 5, MET: 2, PHE: 4, PRO: 1, SER: 2, THR: 1, TRP: 1, VAL: 1', 'Total: 36', 'Ionizable: 10',
              'Ratio: 27.8%')

### Sites:
  - ABR: ['RESIDUES IMPLICATED IN ACTIN BINDING.', 'LYS A  65', 'LYS A  70', 'LYS A  71', 'GLU A  72', 'LYS A  73', 'LEU A  75', 'PHE A  76']

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "MET A  41"
 - <strong>CTR</strong>: "PHE A  76"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A  41" to " CB  MET A  41"
- d= 1.75: "HD22 LEU A  42" to "HG22 VAL A  50"
- d= 2.00: " O   PHE A  58" to "HD13 LEU A  61"

</details>

