---
# 1CDC :: Cd2, N-Terminal Domain (1-99), Truncated Form
## PDB.Structure
### Function: IMMUNE SYSTEM PROTEIN, RECEPTOR
### First Release: 23-MAY-95
### Method: X-RAY DIFFRACTION
### Resolution: 2.00 ANGSTROMS.
### Molecule: CD2
### Seqres Species: Residues: B:99, A:99; Total: 198
### Total waters: 110
### Missing:
  - Residues:
 'A': [('ARG', 1), ('ASP', 2), ('SER', 3), ('Count', 3)],
 'B': [('ARG', 1), ('ASP', 2), ('SER', 3), ('Count', 3)]

### Biounits: DIMERIC: B, A;
### Models: 1
### Chains: A, B
### Model 1 Residues:
  - B:
 'RESIDUES': ('ALA: 5, ARG: 6, ASN: 7, ASP: 8, GLN: 1, GLU: 5, GLY: 9, HIS: 1, ILE: 7, LEU: 11, LYS: 7, MET: 2, PHE: 4, PRO: 2, SER: 4, THR: 8, TRP: 2, TYR: 2, VAL: 5', 'Total: 96', 'Ionizable: 29',
              'Ratio: 30.2%')
  - A:
 'RESIDUES': ('ALA: 5, ARG: 6, ASN: 7, ASP: 8, GLN: 1, GLU: 5, GLY: 9, HIS: 1, ILE: 7, LEU: 11, LYS: 7, MET: 2, PHE: 4, PRO: 2, SER: 4, THR: 8, TRP: 2, TYR: 2, VAL: 5', 'Total: 96', 'Ionizable: 29',
              'Ratio: 30.2%')

### Model 1 Free Cofactors & Waters:
  - B:
 'HOH': 54
  - A:
 'HOH': 56

## MCCE.Step1
### Termini:
 - <strong>NTG</strong>: "GLY B   4", "GLY A   4"
 - <strong>CTR</strong>: "GLU B  99", "GLU A  99"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 110.

### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR B  99":   OXT
  -    Missing heavy atoms for CTR01 in "CTR A  99":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.91: " CB  GLU A  99" to " C   CTR A  99"

</details>

