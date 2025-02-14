
# tcms: Tautomeric, TopN Charge Micro States

## USAGE:
```
This tool outputs:
  - A file listing the topN tautomeric charge microstates, along with
    their related properties: energy (E), net charge (sum_crg), count,
    and occupancy (occ).
  - The <topN> pdb files of each charge state.

Usage:
If called inside a mcce output folder (not a requirement) at pH7 (default) & n_top=5:
  tcms input.pdb

Otherwise:
  tcms path/to/mcce_dir/input.pdb -ph 4 n_top 10

positional arguments:
  inputpdb_filepath  The pdb filepath in a mcce output dir.

options:
  -h, --help         show this help message and exit
  -ph PH             pH point at which the charge microstates are retrieved; Default: 7.0.
  -n_top N_TOP       Number of most favorable charge microstates to return; Default: 5.
  --overwrite        Overwrite existing output files; Default: False.

Report issues here:
        https://github.com/GunnerLab/MCCE4/issues
```

### Note on the only required option:
The only required input to `tcms` is the filepath of a pdb that is in a mcce run dir.
It is not required to run the tool inside a mcce run_dir. However, the log file, `tcms.log`, __is created where the call is made, so it may be preferable to cd to the run directory before running `tcms`__.

### Note on the 'overwrite' flag:
When using the `tcms` tool in MCCE4/MCCE_bin, the `--overwrite` cli option is set interactively with a prompt to the user.
If the flag is used (programmatically), or set to True via the cli, the existing output folder is deleted and recreated.

### Output folder name:
To facilitate the retention of `tcms` calculations with different parameters, the default output folder
name now includes the ph and number of cms to return, i.e. the default folder name is: 'tcms_ph7.00_top5'.  

### Output files in pdb format have a `REMARK 250` section:
The charge state pdb files created contain a `REMARK 250` section, which retains the pdb name and the charge vector totals:
```
$ grep "REMARK" tcrgms_1.pdb
REMARK 250
REMARK 250 EXPERIMENTAL DETAILS
REMARK 250   EXPERIMENT TYPE               : MCCE simulation; MCCE4 tcms tool.
REMARK 250   DATE OF DATA COLLECTION       : 07-Jul-24
REMARK 250   REMARK: DATE OF DATA COLLECTION is the creation date of this pdb.
REMARK 250 EXPERIMENTAL CONDITIONS
REMARK 250   SIMULATION INPUT PDB          : 4lzt.pdb
REMARK 250   TEMPERATURE                   : 298.15 (K)
REMARK 250   PH                            : 7.00
REMARK 250 CHARGE MICROSTATE INFORMATION
REMARK 250   ENERGY                        : -181 (kcal/mol)
REMARK 250   NET CHARGE                    : 7
REMARK 250   COUNT                         : 1,669,043
REMARK 250   OCCUPANCY                     : 70.00%
REMARK 250 REMARK:
REMARK 250  This pdb was created from a tautomeric charge microstate vector
REMARK 250  extracted by tcms, a MCCE4 tool.
REMARK 250
```

## Using non-default inputs:
The pdb filepath is all you need if your pH is 7 and the 'top N' charge states to return is 5. 
Call example for non-default inputs:
```
 tcms <run_dir>/4lzt.pdb -ph 7.5 -n_top 3
```

## Files format:
 1. 'User file', topN_tcms.tsv:

 __Column names:__  
   * 'residues': residue ids or identifier for the `info` = `totals`rows ('E', 'sum_crg', 'count', 'occ');
   * columns with integers: 'top n' columns (1 - N), column '1' is the most favorable state;
   * 'info': for residues rows, possible values ('free' or 'fixed') denote their simulation states; for
   other rows, the value is 'totals';
   * 'chain': residues chain, empty for non-residue rows.

 __Column count:__
   * The number of columns in the user file is 3 + n_top if all requested cms pass the occ threshold (1%).

 2. 'Complete file', topN_master_tcms.tsv:
 Same as above, except that each 'top n' column is preceded by its 'conf_n' column holding the conformer index.
 __Column count:__ 3 + n_top x 2, with n_top being the effective number of returned cms as above.

## Final screen output:
When the processing is over, the ionizable residues in the user file are displayed. Note that the names of the 
'topN files' retains the requested number of cms to return, yet the files will have columns for only those states 
with occupancy beyond the 1% threshold as shown here:
```
Ionizable residues in '/path/to/4lzt4/tcms_ph7.00_top5/top5_tcms.tsv':

     residues        1        2      3      4    info chain
0    NTRA0001        0        0      1      0    free     A
1    LYSA0001        1        1      1      1    free     A
4    ARGA0005        1        1      1      1   fixed     A
6    GLUA0007       -1       -1     -1     -1    free     A
12   LYSA0013        1        1      1      1   fixed     A
13   ARGA0014        1        1      1      1   fixed     A
14   HISA0015      NE2        1    ND1    ND1    free     A
16   ASPA0018       -1       -1     -1     -1   fixed     A
18   TYRA0020        0        0      0      0    free     A
19   ARGA0021        1        1      1      1   fixed     A
20   TYRA0023        0        0      0      0    free     A
29   LYSA0033        1        1      1      1   fixed     A
31   GLUA0035       -1       -1     -1    OE1    free     A
41   ARGA0045        1        1      1      1   fixed     A
44   ASPA0048       -1       -1     -1     -1    free     A
47   ASPA0052       -1       -1     -1     -1    free     A
48   TYRA0053        0        0      0      0   fixed     A
55   ARGA0061        1        1      1      1   fixed     A
60   ASPA0066       -1       -1     -1     -1    free     A
61   ARGA0068        1        1      1      1    free     A
65   ARGA0073        1        1      1      1   fixed     A
79   ASPA0087       -1       -1     -1     -1   fixed     A
88   LYSA0096        1        1      1      1   fixed     A
89   LYSA0097        1        1      1      1   fixed     A
93   ASPA0101       -1       -1     -1     -1   fixed     A
102  ARGA0112        1        1      1      1    free     A
104  ARGA0114        1        1      1      1   fixed     A
106  LYSA0116        1        1      1      1    free     A
108  ASPA0119       -1       -1     -1     -1   fixed     A
114  ARGA0125        1        1      1      1   fixed     A
116  ARGA0128        1        1      1      1   fixed     A
118  CTRA0129       -1       -1     -1     -1   fixed     A
128         E     -174     -174   -172   -172  totals
129   sum_crg        7        8      8      8  totals
130     count  317,222  119,907  5,781  5,328  totals
131       occ   69.51%  26.27%  1.27%   1.17%  totals
```    

## All outputs:
The subfolder `tcms_ph7.00_top5` contains all output files, here with n_top=5 (and pdbs in s2 & standard format):
```
s2_tcms1_4lzt.pdb
s2_tcms2_4lzt.pdb 
s2_tcms3_4lzt.pdb
s2_tcms4_4lzt.pdb
tcms1_4lzt.pdb
tcms2_4lzt.pdb
tcms3_4lzt.pdb
tcms4_4lzt.pdb
top5_master_tcms.tsv
top5_tcms.tsv
```
The file `top5_master_tcms.tsv` retains the vectors of conformer indices that were used to create the pdbs.
