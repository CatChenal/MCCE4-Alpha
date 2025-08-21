---
# 6GST :: First-Sphere And Second-Sphere Electrostatic Effects In The Active Site Of A Class Mu Glutathione Transferase;
## PDB.Structure
### Function: TRANSFERASE
### First Release: 26-JAN-96
### Method: X-RAY DIFFRACTION
### Resolution: 2.20 ANGSTROMS.
### Molecule: MU CLASS GLUTATHIONE S-TRANSFERASE OF ISOENZYME 3-3
### Seqres Species: Residues: A:217, B:217; Total: 434
### Cofactors:
  - GSH:
 'GLUTATHIONE', 2

### Total cofactors: 2
### Total waters: 420
### Biounits: DIMERIC: A, B;
### Models: 1
### Chains: A, B
### Model 1 Residues:
  - A:
 'RESIDUES': ('ALA: 11, ARG: 13, ASN: 9, ASP: 15, CYS: 3, GLN: 7, GLU: 15, GLY: 9, HIS: 4, ILE: 13, LEU: 25, LYS: 19, MET: 8, PHE: 12, PRO: 12, SER: 12, THR: 7, TRP: 4, TYR: 14, VAL: 5', 'Total: 217', 'Ionizable: 83',
              'Ratio: 38.2%')
  - B:
 'RESIDUES': ('ALA: 11, ARG: 13, ASN: 9, ASP: 15, CYS: 3, GLN: 7, GLU: 15, GLY: 9, HIS: 4, ILE: 13, LEU: 25, LYS: 19, MET: 8, PHE: 12, PRO: 12, SER: 12, THR: 7, TRP: 4, TYR: 14, VAL: 5', 'Total: 217', 'Ionizable: 83',
              'Ratio: 38.2%')

### Model 1 Free Cofactors & Waters:
  - A:
 'GSH': 1, 'HOH': 229
  - B:
 'GSH': 1, 'HOH': 191

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE GSH A 218', 'TYR A   6', 'TRP A   7', 'TRP A  45', 'LYS A  49', 'ASN A  58', 'LEU A  59', 'GLN A  71', 'SER A  72', 'ASP B 105', 'HOH A 221', 'HOH A 223', 'HOH A 239', 'HOH A 277', 'HOH A 351', 'HOH A 413', 'HOH A 414']
  - AC2: ['BINDING SITE FOR RESIDUE GSH B 218', 'ASP A 105', 'TYR B   6', 'TRP B   7', 'TRP B  45', 'LYS B  49', 'ASN B  58', 'LEU B  59', 'GLN B  71', 'SER B  72', 'HOH B 221', 'HOH B 222', 'HOH B 261', 'HOH B 266', 'HOH B 278', 'HOH B 307', 'HOH B 321', 'HOH B 353', 'HOH B 386']

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "PRO A   1", "PRO B   1"
 - <strong>CTR</strong>: "LYS A 217", "LYS B 217"

### Labeling:
<strong><font color='red'>Generic topology file created for:</font></strong>  
GSH: https://pubchem.ncbi.nlm.nih.gov/#query=GSH&tab=substance; 

### Free Cofactors:
  - Removed all 229 HOH in A. Removed all 191 HOH in B.
  - NOTE: Include the '--wet' option at the command line to keep buried waters and cofactors. Alternatively, change the water SAS cutoff to a non-zero, positive number using the command line 'u' option:
  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05
  - Total deleted cofactors = 420.
  - Species and properties with assigned default values in debug.log:

  - GSHBK: ['VDW_RAD', 'VDW_EPS']

  - get_connect12():: ['Error!']

  - atom: ['get_connect12():']

  - PRO01: ['TORSION']


### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.54: " CA  NTR A   1" to " CB  PRO A   1"
- d= 1.49: " N   NTR A   1" to " CD  PRO A   1"
- d= 1.57: " CA  NTR B   1" to " CB  PRO B   1"
- d= 1.50: " N   NTR B   1" to " CD  PRO B   1"

</details>

