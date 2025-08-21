---
# 1DSB :: Crystal Structure Of The Dsba Protein Required For Disulphide Bond Formation In Vivo;
## PDB.Structure
### Function: DISULFIDE OXIDOREDUCTASE
### First Release: 24-MAY-93
### Method: X-RAY DIFFRACTION
### Resolution: 2.00 ANGSTROMS.
### Molecule: DSBA
### Seqres Species: Residues: A:189, B:189; Total: 378
### Total waters: 195
### Missing:
  - Residues:
 'A': [('LYS', 189), ('Count', 1)], 'B': [('LYS', 189), ('Count', 1)]

### Biounits: DIMERIC: A, B;
### Models: 1
### Chains: A, B
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 17, ARG: 3, ASN: 7, ASP: 12, CYS: 2, GLN: 14, GLU: 12, GLY: 13, HIS: 3, ILE: 5, LEU: 12, LYS: 16, MET: 6, PHE: 11, PRO: 7, SER: 8, THR: 9, TRP: 2, TYR: 8, VAL: 21', 'Total: 188', 'Ionizable: 56',
              'Ratio: 29.8%')
  - B:
 'RESIDUES': ('ALA: 17, ARG: 3, ASN: 7, ASP: 12, CYS: 2, GLN: 14, GLU: 12, GLY: 13, HIS: 3, ILE: 5, LEU: 12, LYS: 16, MET: 6, PHE: 11, PRO: 7, SER: 8, THR: 9, TRP: 2, TYR: 8, VAL: 21', 'Total: 188', 'Ionizable: 56',
              'Ratio: 29.8%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 107
  - B:
 'HOH': 88

### Disulfides:
  - CYS A  30 -- 2.00 Å --> CYS A  33
  - CYS B  30 -- 2.02 Å --> CYS B  33

### Sites:
  - CAA: ['NULL', 'CYS A  30', 'PRO A  31', 'HIS A  32', 'CYS A  33', 'PRO A 151']
  - CAB: ['NULL', 'CYS B  30', 'PRO B  31', 'HIS B  32', 'CYS B  33', 'PRO B 151']

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "ALA A   1", "ALA B   1"
 - <strong>CTR</strong>: "LYS A 188", "LYS B 188"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 195.

### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 188":   OXT
  -    Missing heavy atoms for CTR01 in "CTR B 188":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.52: " CA  NTR A   1" to " CB  ALA A   1"
- d= 1.47: " CA  NTR B   1" to " CB  ALA B   1"

</details>

