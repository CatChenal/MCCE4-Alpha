---
# 1L54 :: The Structural And Thermodynamic Consequences Of Burying A Charged Residue Within The Hydrophobic Core Of T4 Lysozyme;
## PDB.Structure
### Function: HYDROLASE (O-GLYCOSYL)
### First Release: 28-JAN-91
### Method: X-RAY DIFFRACTION
### Resolution: 1.90 ANGSTROMS.
### Molecule: T4 LYSOZYME
### Seqres Species: Residues: A:164; Total: 164
### Total waters: 152
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 16, ARG: 13, ASN: 12, ASP: 10, GLN: 5, GLU: 8, GLY: 11, HIS: 1, ILE: 10, LEU: 16, LYS: 14, MET: 4, PHE: 5, PRO: 3, SER: 6, THR: 12, TRP: 3, TYR: 6, VAL: 9', 'Total: 164', 'Ionizable: 52',
              'Ratio: 31.7%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 152

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "MET A   1"
 - <strong>CTR</strong>: "LEU A 164"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 152.

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   1" to " CB  MET A   1"

</details>

