---
# 8UOZ :: Emre Structure In The Tpp-Bound State (Wt/E14Q Heterodimer)
## PDB.Structure
### Function: MEMBRANE PROTEIN
### First Release: 20-OCT-23
### Method: SOLUTION NMR; SOLID-STATE NMR
### Resolution: NOT APPLICABLE.
### Molecule: SMR FAMILY MULTIDRUG EFFLUX PROTEIN EMRE
### Seqres Species: Residues: A:110, B:110; Total: 220
### Cofactors:
  - P4P:
 'TETRAPHENYLPHOSPHONIUM', 1

### Total cofactors: 1
### Biounits: DIMERIC: A, B;
### Models: 1
### Chains: A, B
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 9, ARG: 3, ASN: 2, ASP: 1, CYS: 3, GLN: 2, GLU: 2, GLY: 12, HIS: 1, ILE: 15, LEU: 16, LYS: 1, MET: 4, PHE: 5, PRO: 5, SER: 8, THR: 7, TRP: 4, TYR: 5, VAL: 5', 'Total: 110', 'Ionizable: 16',
              'Ratio: 14.5%')
  - B:
 'RESIDUES': ('ALA: 9, ARG: 3, ASN: 2, ASP: 1, CYS: 3, GLN: 3, GLU: 1, GLY: 12, HIS: 1, ILE: 15, LEU: 16, LYS: 1, MET: 4, PHE: 5, PRO: 5, SER: 8, THR: 7, TRP: 4, TYR: 5, VAL: 5', 'Total: 110', 'Ionizable: 15',
              'Ratio: 13.6%')

### Model 1 Free Cofactors & Waters:
  - A:
 'P4P': 1

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "MET A   1", "MET B   1"
 - <strong>CTR</strong>: "HIS A 110", "HIS B 110"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
P4P: https://pubchem.ncbi.nlm.nih.gov/#query=P4P&tab=substance; 

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Species and properties with assigned default values in debug.log:

  - P4PBK: ['VDW_RAD', 'VDW_EPS']


### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.56: " CA  NTR A   1" to " CB  MET A   1"
- d= 1.55: " CA  NTR B   1" to " CB  MET B   1"

</details>

