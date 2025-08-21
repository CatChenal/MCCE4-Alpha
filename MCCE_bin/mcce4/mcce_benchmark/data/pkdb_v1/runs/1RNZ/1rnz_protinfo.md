---
# 1RNZ :: Ribonuclease A Crystallized From 2.5M Sodium Chloride, 3.3M Sodium Formate;
## PDB.Structure
### Function: ENDONUCLEASE
### First Release: 08-NOV-96
### Method: X-RAY DIFFRACTION
### Resolution: 1.90 ANGSTROMS.
### Molecule: RIBONUCLEASE A
### Seqres Species: Residues: A:124; Total: 124
### Cofactors:
  -  CL:
 'CHLORIDE ION', 3

### Total cofactors: 3
### Total waters: 95
### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 12, ARG: 4, ASN: 10, ASP: 5, CYS: 8, GLN: 7, GLU: 5, GLY: 3, HIS: 4, ILE: 3, LEU: 2, LYS: 10, MET: 4, PHE: 3, PRO: 4, SER: 15, THR: 10, TYR: 6, VAL: 9', 'Total: 124', 'Ionizable: 42',
              'Ratio: 33.9%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 95, '_CL': 3

### Disulfides:
  - CYS A  26 -- 2.00 Å --> CYS A  84
  - CYS A  40 -- 1.97 Å --> CYS A  95
  - CYS A  58 -- 2.01 Å --> CYS A 110
  - CYS A  65 -- 2.00 Å --> CYS A  72

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE CL A 126', 'ARG A  10', 'ASN A  34', 'THR A  78']
  - AC2: ['BINDING SITE FOR RESIDUE CL A 127', 'THR A  45']
  - AC3: ['BINDING SITE FOR RESIDUE CL A 128', 'THR A   3', 'HOH A 219']

## MCCE.Step1
### Renamed:
  - "CL    CL A 126" to "CL   _CL A 126"
  - "CL    CL A 127" to "CL   _CL A 127"
  - "CL    CL A 128" to "CL   _CL A 128"

### Termini:
 - <strong>NTR</strong>: "LYS A   1"
 - <strong>CTR</strong>: "VAL A 124"

### Free Cofactors:
  - Removed all 95 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 95.

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.51: " CA  NTR A   1" to " CB  LYS A   1"
- d= 2.00: " SG  CYS A  26" to " SG  CYS A  84"
- d= 1.97: " SG  CYS A  40" to " SG  CYS A  95"
- d= 2.00: " SG  CYS A  65" to " SG  CYS A  72"

</details>

