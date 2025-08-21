---
# 1PTD :: Phosphatidylinositol-Specific Phospholipase C
## PDB.Structure
### Function: HYDROLASE (PHOSPHORIC DIESTER)
### First Release: 24-MAY-95
### Method: X-RAY DIFFRACTION
### Resolution: 2.60 ANGSTROMS.
### Molecule: PHOSPHATIDYLINOSITOL-SPECIFIC PHOSPHOLIPASE C
### Seqres Species: Residues: A:298; Total: 298
### Total waters: 52
### Missing:
  - Residues:
 'A': [('LYS', 297), ('GLU', 298), ('Count', 2)]

### Biounits: MONOMERIC: A;
### Models: 1
### Chains: A
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 12, ARG: 9, ASN: 24, ASP: 19, GLN: 11, GLU: 19, GLY: 16, HIS: 6, ILE: 23, LEU: 20, LYS: 23, MET: 7, PHE: 12, PRO: 14, SER: 23, THR: 18, TRP: 7, TYR: 19, VAL: 14', 'Total: 296', 'Ionizable: 95',
              'Ratio: 32.1%')

### Model 1 Free Cofactors & Waters:
  - A:
 'HOH': 52

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "ALA A   1"
 - <strong>CTR</strong>: "ILE A 296"

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 52.

### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 296":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   1" to " CB  ALA A   1"

</details>

