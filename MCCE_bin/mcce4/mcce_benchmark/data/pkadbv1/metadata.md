
## Experimental pKas data source
The original data comes from [Dr. Emil Axelov's pKa Database (1)](http://compbio.clemson.edu/lab/software/5/). The list for wild types and mutants alike was >
mainly remove membrane proteins and those containing nucleotides, and to select the biological unit(s). The resulting files: `WT_pkas.csv` and `MT_pkas.csv`, co>
However, only the WT data is currently used as most of pKas of the mutant proteins were measured by NMR, hence
yield pKa _ranges_ rather than values.

#### pKa file header:
**'PDB ID', 'Res Name', 'Chain', 'Res ID', 'Expt. pKa', 'Expt. Uncertainty', '%SASA', 'Expt. method', 'Expt. salt conc','Expt. pH', 'Expt. temp', 'Reference'**

### File 'proteins.txt':
Comments out the excluded pdbs and gives the reason.
### File 'metadata.md':
Experimental data source details (this file); to be kept in data folder.
### Folder clean_pdbs:
Holds the prepared pdb files which reside inside a folder with the same id in upper case.
This folder is re-used for every MCCE benchmarking job.

The "queue book" ('./clean_pdbs/book.txt' file), is used to monitor the submitted batch
of MCCE simulations, and the submit script 'run.sh' is used to run MCCE in each subfolder.
