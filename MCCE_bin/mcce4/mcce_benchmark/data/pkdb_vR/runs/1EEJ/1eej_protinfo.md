---
# 1EEJ :: Crystal Structure Of The Protein Disulfide Bond Isomerase, Dsbc, From Escherichia Coli;
## PDB.Structure
### Function: ISOMERASE
### First Release: 31-JAN-00
### Method: X-RAY DIFFRACTION
### Resolution: 1.90 ANGSTROMS.
### Molecule: DISULFIDE INTERCHANGE PROTEIN
### Seqres Species: Residues: A:216, B:216; Total: 432
### Cofactors:
  - MES:
 '2-(N-MORPHOLINO)-ETHANESULFONIC ACID', 1

### Total cofactors: 1
### Total waters: 208
### Biounits: DIMERIC: A, B;
### Models: 1
### Chains: A, B
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 21, ARG: 2, ASN: 7, ASP: 17, CYS: 4, GLN: 11, GLU: 9, GLY: 16, HIS: 6, ILE: 12, LEU: 16, LYS: 20, MET: 10, PHE: 4, PRO: 11, SER: 10, THR: 13, TRP: 1, TYR: 8, VAL: 18', 'Total: 216', 'Ionizable: 66',
              'Ratio: 30.6%')
  - B:
 'RESIDUES': ('ALA: 21, ARG: 2, ASN: 7, ASP: 17, CYS: 4, GLN: 11, GLU: 9, GLY: 16, HIS: 6, ILE: 12, LEU: 16, LYS: 20, MET: 10, PHE: 4, PRO: 11, SER: 10, THR: 13, TRP: 1, TYR: 8, VAL: 18', 'Total: 216', 'Ionizable: 66',
              'Ratio: 30.6%')

### Model 1 Free Cofactors & Waters:
  - B:
 'HOH': 125, 'MES': 1
  - A:
 'HOH': 83

### Disulfides:
  - CYS A  98 -- 2.06 Å --> CYS A 101
  - CYS A 141 -- 2.03 Å --> CYS A 163
  - CYS B  98 -- 2.05 Å --> CYS B 101
  - CYS B 141 -- 2.03 Å --> CYS B 163

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE MES B 1001', 'SER B  17', 'LYS B  82', 'PRO B  84', 'GLN B  85', 'GLU B  86', 'ASP B 166', 'ASP B 169', 'HOH B1111']

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "ASP A   1", "ASP B   1"
 - <strong>CTR</strong>: "LYS A 216", "LYS B 216"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
MES: https://pubchem.ncbi.nlm.nih.gov/#query=MES&tab=substance; 

### Free Cofactors:
  - Removed all 125 HOH in B.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 208.
  - Species and properties with assigned default values in debug.log:

  - MESBK: ['VDW_RAD', 'VDW_EPS']


### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   1" to " CB  ASP A   1"
- d= 1.53: " CA  NTR B   1" to " CB  ASP B   1"

</details>

