---
# 1GYM :: Phosphatidylinositol-Specific Phospholipase C In Complex With Glucosamine-(Alpha-1-6)-Myo-Inositol;
## PDB.Structure
### Function: HYDROLASE (PHOSPHORIC DIESTER)
### First Release: 02-MAY-96
### Method: X-RAY DIFFRACTION
### Resolution: 2.20 ANGSTROMS.
### Molecule: PHOSPHATIDYLINOSITOL-SPECIFIC PHOSPHOLIPASE C
### Seqres Species: Residues: A:298; Total: 298
### Cofactors:
  - MYG:
 '2-AMINO-2-DEOXY-GLUCOSIDE', 1

### Total cofactors: 1
### Total waters: 127
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
 'HOH': 127, 'MYG': 1

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "ALA A   1"
 - <strong>CTR</strong>: "ILE A 296"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
MYG: https://pubchem.ncbi.nlm.nih.gov/#query=MYG&tab=substance; 

### Free Cofactors:
  - Removed all 127 HOH in A.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 127.
  - Species and properties with assigned default values in debug.log:

  - MYGBK: ['VDW_RAD', 'VDW_EPS']


### Missing Heavy Atoms:
  -    Missing heavy atoms for CTR01 in "CTR A 296":   OXT

### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   1" to " CB  ALA A   1"

</details>

