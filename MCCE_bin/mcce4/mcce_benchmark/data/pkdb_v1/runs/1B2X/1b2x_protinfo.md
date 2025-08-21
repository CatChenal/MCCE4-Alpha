---
# 1B2X :: Barnase Wildtype Structure At Ph 7.5 From A Cryo_Cooled Crystal At 100K;
## PDB.Structure
### Function: HYDROLASE
### First Release: 03-DEC-98
### Method: X-RAY DIFFRACTION
### Resolution: 1.80 ANGSTROMS.
### Molecule: PROTEIN (BARNASE)
### Seqres Species: Residues: A:110, B:110, C:110; Total: 330
### Cofactors:
  -  ZN:
 'ZINC ION', 1

### Total cofactors: 1
### Total waters: 380
### Missing:
  - Residues:
 'A': [('ALA', 1), ('Count', 1)],
 'B': [('ALA', 1), ('GLN', 2), ('Count', 2)],
 'C': [('ALA', 1), ('GLN', 2), ('Count', 2)]
  - ResAtoms:
 'A': ['GLN_2=(N,CA,O,CB,CG,CD,OE1)', 'GLN_2=(NE2)']

### Biounits: MONOMERIC: A,B,C;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 7, ARG: 6, ASN: 6, ASP: 9, GLN: 4, GLU: 3, GLY: 10, HIS: 2, ILE: 8, LEU: 7, LYS: 8, PHE: 4, PRO: 3, SER: 9, THR: 9, TRP: 3, TYR: 7, VAL: 4', 'Total: 109', 'Ionizable: 35',
              'Ratio: 32.1%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 126

### Links:
  - O  HOH A 133 -- 2.13 Å --> ZN  ZN C 384
  - ND1 HIS C 18 -- 2.13 Å --> ZN  ZN C 384
  - OE2 GLU C 60 -- 2.19 Å --> ZN  ZN C 384
  - OE1 GLU C 60 -- 2.40 Å --> ZN  ZN C 384
  - NZ LYS C 62 -- 2.21 Å --> ZN  ZN C 384

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE ZN C 384', 'HOH A 133', 'HIS C  18', 'GLU C  60', 'LYS C  62']

## MCCE.Step1
### Termini:
 - <strong>CTR</strong>: "ARG A 110"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 126.

### Distance Clashes:
<details><summary>Clashes found</summary>

- No clash found.

</details>

