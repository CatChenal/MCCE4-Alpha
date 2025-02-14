# ProtInfo
This tool provides information about a protein pdb file on its structure, along with applied transformations and issues flagged in the log produced by MCCE's step1.py.  
This information is saved into a file named `ProtInfo.md`.

  * See sample markdown [reports](#Samples)

## USAGE:
```
Command line interface for the MCCE ProtInfo tool, which gathers:
 * Info about the input protein from mcce4.pdbio.Structure.
 * Info from MCCE step1 run.log & debug.log when step1 can be run.

Options:
 1. pdb (required): a pdb file name or pdbid.
 2. --fetch (False if not used): If 'pdb' is a pdbid and flag is used,
    the biological assembly is downloaded from rcsb.org.

 Step1 options (default value):
  --wet (False): Keep water molecules.
  --noter (False): Do not label terminal residues (for making ftpl).
  -d (4): Protein dielectric constant for delphi.
  -u (''): User selected, comma-separated KEY=var pairs from run.prm; e.g.:
           -u HOME_MCCE=/path/to/mcce_home,EXTRA=./extra.tpl.
  -e (mcce): mcce executable location.
  --fetch (False): Download the biological assembly of given pdb (if not a file).

Usage:
 > ProtInfo 1fat --fetch
 > ProtInfo 1fat.pdb --wet
 > ProtInfo 1fat.pdb --noter
 > ProtInfo  --help  # Show the tool help
```

## Basic Info:
Info recorded in the final report:
  * Number of models
  * Number of chains
  * Number of residues
  * Unknown residues
  * Residues with multiple alternate locations:
    - MCCE can only handle a single location and only considers the 'A' location (which may cause problems)
  * Number of water molecules
  * Number and identity of cofactors

## MCCE Info:
After running MCCE's step1.py, the run.log provides information about errors and repairs applied.  

Info recorded in the final report:
  * Terminal residues
  * Waters: count and buried list
  * Buried Cofactors:
    - Water: buried list & count
    - Other: buried list & count
    - Missing topology
  * __TODO__: Which repairs were applied  ::  Need more info.
  * Fatal errors

---
---

# Sample reports for 1A6K and 1AIG:

---
---
__ProtInfo__ - <u>User options:</u> __pdb__: '1a6k.pdb'; __fetch__: False; __d__: 4; __e__: 'mcce'; __u__: ''; __wet__: False; __noter__: False;

---
# 1A6K :: Aquomet-Myoglobin, Atomic Resolutio
## ParsedStructure
### Function: HEME PROTEIN
### First Release: 26-FEB-98
### Method: X-RAY DIFFRACTION
### Resolution: 1.10 ANGSTROMS.
### Molecule: MYOGLOBIN
### Models: 1
### Chains: A
### Seqres Species: Residues: A:151
### Hetero Species:
  - SO4 SULFATE ION
  - HEM PROTOPORPHYRIN IX CONTAINING FE

### Hetero_Aliases:
  - HEM HEME

### Links:
  - FE  HEM A 154 -- 2.14 Å --> NE2 HIS A 93
  - FE  HEM A 154 -- 2.13 Å --> O  HOH A1001

### Sites:
  - AC1: ['BINDING SITE FOR RESIDUE SO4 A 155', 'ALA A  57', 'SER A  58', 'GLU A  59', 'ASP A  60', 'HOH A1007', 'HOH A1026', 'HOH A1038', 'HOH A1086']
  - AC2: ['BINDING SITE FOR RESIDUE SO4 A 156', 'GLN A  26', 'LYS A  62', 'HOH A1043', 'HOH A1045', 'HOH A1063']
  - AC3: ['BINDING SITE FOR RESIDUE SO4 A 1188', 'ARG A  45', 'HIS A  64', 'THR A  67', 'HIS A 116', 'HOH A1088', 'HOH A1123']
  - AC4: ['BINDING SITE FOR RESIDUE HEM A 154', 'THR A  39', 'LYS A  42', 'PHE A  43', 'ARG A  45', 'HIS A  64', 'THR A  67', 'VAL A  68', 'LEU A  89', 'SER A  92', 'HIS A  93', 'HIS A  97', 'ILE A  99', 'TYR A 103', 'LEU A 104', 'PHE A 138', 'HOH A1001', 'HOH A1024', 'HOH A1033', 'HOH A1059', 'HOH A1088', 'HOH A1092', 'HOH A1118', 'HOH A1149']

## MCCE.Step1
### Renamed:
  - " CAA HEM A 154" to " CAA PAA A 154"
  - " CAD HEM A 154" to " CAD PDD A 154"
  - " CBA HEM A 154" to " CBA PAA A 154"
  - " CBD HEM A 154" to " CBD PDD A 154"
  - " CGA HEM A 154" to " CGA PAA A 154"
  - " CGD HEM A 154" to " CGD PDD A 154"
  - " O1A HEM A 154" to " O1A PAA A 154"
  - " O1D HEM A 154" to " O1D PDD A 154"
  - " O2A HEM A 154" to " O2A PAA A 154"
  - " O2D HEM A 154" to " O2D PDD A 154"

### Termini:
 - <strong>NTR</strong>: "VAL A   1"
 - <strong>CTR</strong>: "TYR A 151"

### Labeling:

<strong><font color='red'>Generic topology file created for</font></strong>:
SO4:  https://pubchem.ncbi.nlm.nih.gov/#query=SO4&tab=substance;

### Free Cofactors:
  - NOTE: Include the '--wet' option at the command line to keep waters and cofactors.
  - Total deleted cofactors = 187.
  - Species and properties with assigned default values in debug.log:

  - SO4BK: ['VDW_RAD', 'VDW_EPS']

  - HEM01: ['TORSION']

  - Empty connection slot(s):

  - FE: ['HEM 154 to atom  SD']


### Distance Clashes:
<details><summary>Clashes found</summary>

- d= 1.53: " CA  NTR A   1" to " CB  VAL A   1"
- d= 1.53: " C2A HEM A 154" to " CAA PAA A 154"
- d= 1.54: " C3D HEM A 154" to " CAD PDD A 154"

</details>
---

---
__ProtInfo__ - <u>User options:</u> __pdb__: '1aig.pdb'; __fetch__: False; __d__: 4; __e__: 'mcce'; __u__: ''; __wet__: False; __noter__: False;

---
# 1AIG :: Photosynthetic Reaction Center From Rhodobacter Sphaeroides In The D+Qb-Charge Separated State
## ParsedStructure
### Function: PHOTOSYNTHETIC REACTION CENTER
### First Release: 17-APR-97
### Method: X-RAY DIFFRACTION
### Models: 1
### Chains: H, L, M
### Seqres Species: Residues: L:281, M:307, H:260, N:281, O:307, P:260
### Hetero Species:
  - BCL BACTERIOCHLOROPHYLL A
  - BPH BACTERIOPHEOPHYTIN A
  - U10 UBIQUINONE-10
  - FE2 FE (II) ION

### Hetero_Aliases:
  - U10 COENZYME Q10

### Links:
  - NE2 HIS L 153 -- 2.57 Å --> MG  BCL L 283
  - NE2 HIS L 173 -- 2.37 Å --> MG  BCL L 282
  - NE2 HIS L 190 -- 2.13 Å --> FE  FE2 M 308
  - NE2 HIS L 230 -- 2.36 Å --> FE  FE2 M 308
  - NE2 HIS M 182 -- 2.42 Å --> MG  BCL M 309
  - NE2 HIS M 202 -- 2.34 Å --> MG  BCL M 310
  - NE2 HIS M 219 -- 2.14 Å --> FE  FE2 M 308
  - OE1 GLU M 234 -- 2.12 Å --> FE  FE2 M 308
  - OE2 GLU M 234 -- 2.15 Å --> FE  FE2 M 308
  - NE2 HIS M 266 -- 2.28 Å --> FE  FE2 M 308
  - NE2 HIS N 153 -- 2.68 Å --> MG  BCL N 284
  - NE2 HIS N 173 -- 2.22 Å --> MG  BCL N 283
  - NE2 HIS N 190 -- 2.16 Å --> FE  FE2 O 308
  - NE2 HIS N 230 -- 2.20 Å --> FE  FE2 O 308
  - MG  BCL N 282 -- 2.38 Å --> NE2 HIS O 182
  - MG  BCL N 283 -- 3.10 Å --> OBB BCL O 309
  - NE2 HIS O 202 -- 2.22 Å --> MG  BCL O 309
  - NE2 HIS O 219 -- 2.17 Å --> FE  FE2 O 308
  - OE1 GLU O 234 -- 2.18 Å --> FE  FE2 O 308
  - OE2 GLU O 234 -- 2.28 Å --> FE  FE2 O 308
  - NE2 HIS O 266 -- 2.27 Å --> FE  FE2 O 308

## MCCE.Step1
### Termini:
 - <strong>NTR</strong>: "ALA L   1", "TYR M   3", "ASP H  11"
 - <strong>CTR</strong>: "GLY L 281", "HIS M 301", "GLU H 258"

### Labeling:

<strong><font color='red'>Generic topology file created for</font></strong>:
BPH:  https://pubchem.ncbi.nlm.nih.gov/#query=BPH&tab=substance; U10:  https://pubchem.ncbi.nlm.nih.gov/#query=U10&tab=substance;

<strong><font color='red'>Unloadable topology</font></strong>:
Atoms of residue BCL (L 282, L 283, M 309, M 310), do not match the topology conformer BCL-1.

<strong><font color='red'>Likely cause</font></strong>: the renaming file (path: /home/cat/projects/MCCE4/name.txt) is missing entries for these species, resulting in unloadable topology files (path: /home/cat/projects/MCCE4/param).
---
