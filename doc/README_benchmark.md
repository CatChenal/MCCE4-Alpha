<!-- Editing Guide: The pipe (|) position in this comment is 120:                                                       | -->

# Introduction to the MCCE4 benchmarking app
* Path to codebase: `MCCE4/MCCE_bin/mcce4/mcce_benchmark`
* Note: Command lines start with "> " to denote a generic prompt.
---
---

## Command line interface:
The functionality of the app is available via command line "entry points" which consist of:
 * A specific command and its options, e.g.: 
 ```
 > bench_analyze -bench_dir .
 ```
 * A specific command with a subcommand and its options, e.g.: 
 ```
 > bench_setup pkdb_pdbs -bench_dir .
 ```

## Main goal:
The benchmarking app sets up proteins subfolders so that a run script including the first 4 steps of  
an MCCE simulation can either be submitted in batches of 10, or be scheduled for submission every minute 
via the cron scheduler.

## Secondary goals:
 1. Analyze the results in a benchmark folder.
 2. Compare the results between two benchmark folders.

### Five main command choices in MCCE4 benchmarking app:

## 1. Output the list of pdb ids that belong to the pKaDBv1 database :: `bench_setup pdbids`
 * Purpose1: Provide a list of the available pdbs in file named "PDBIDS_WT" (with header: PDBID  res# method)
 * Purpose2: Amend  the file and passed to the '-pdbids_file' option (only the first column is considered).

## 2. Create a new benchmarking set (with a "launch right away" option) :: `bench_setup [pkdb_pdbs | user_pdbs]`
The main command is `bench_setup`, and the subcommand indicates the source of the pdbs to install.

## 3. Launch the mcce simulations over the set :: `bench_setup launch`

Finally, these are the two non-setup commands:  
## 4. Analyze one set of runs :: `bench_analyze`

## 5. Compare two sets of runs (after each set is analyzed) :: `bench_compare`

---
---

# Getting help
Each command & command subcommand combination can display their options information and usage via the `--help` flag, e.g.:

```
> bench_setup --help
> bench_setup user_pdbs --help
```

---
---

# Usage Examples
---

## 1. List the pdb ids
 * Command:
 ```
 > bench_setup pdbids
 ```
 * Outcome:
  The list is displayed and saved as 'PDBIDS_WT' at the call location. 

## 2. Create a new benchmarking set
Note: It is __strongly recommended__ to create the folder where you want the runs to be setup, then cd to it. This way,  
the benchmark.log will reside inside it and only pertain to this specific setup.


## 2.A Using the pdbs in the pKaDBv1 database


### Setup using all defaults
```
> mkdir allpdbs
> cd allpdbs
> bench_setup pkdb_pdbs -bench_dir .
```

#### Outcome
These files and the __runs__ folder are created:
```
> ls
benchmark.info
benchmark.log
bench_setup_options.txt  # also output when --help is run
cli_args.pickle
cli_args.txt
launch.sh
analyze.sh
prerun_report.md
runs/
```
__Executable files: `launch.sh`,  '`analyze.sh':__  
Provide preset commands to:
  * Launch the script over the pdb set via the scheduler __according to the initial setup__
  * Run the intra set analysis __according to the initial setup__

All the 118 pdbs are setup in their respective runs/PDBID folder.
The runs/ folder contains the default bash script ('default_run.sh') the will be used by the batcher or scheduler, 
along with a 'bookkeeping' file ('book.txt'), that will be updated with these status flags on each pdb folder line: 
'i': initial setup, 'r': running, 'c': completed, 'e': ended with error.

### Setup using all defaults, but subsetting the pKDB with a defined number of pdbs
Note: Subsetting the pkDB with `-n_pdbs N` will always output the same N first and smallest pdbs. 

```
> mkdir pdbs6
> cd pdbs6
> bench_setup pkdb_pdbs -bench_dir . -n_pdbs 6
```

#### Outcome
Same outcome except the runs folder will contain only 6 subfolders.

### Setup using all defaults, but subsetting the pKDB with specific pdbids:
 * Prior steps:
   * Run `bench_setup pdbids`
   * Delete or comment the lines of the ids you do not want -> slected_PDBIDS

Note: When this file is read, only the first column is used to obtain the pdbids.

```
> mkdir db_subset
> cd db_subset

> bench_setup pkdb_pdbs -bench_dir . -pdbids_file slected_PDBIDS
```

#### Outcome
Only your selected proteins are setup.


### Setup using non-defaults options
 * The non-default options pertain to each of the 4 mcce 'step commands', e.g. `step1.py`.
 * Changing ANY of the step parameter requires passing a value to the `-job_name` option.

### Setup with a non-default PBE solver

```
> mkdir pdbs_ngpb
> cd pdbs_ngpb
> bench_setup pkdb_pdbs -bench_dir . -s ngpb -job_name run_ngpb
```

#### Outcome

The benchmark.log will have recorded these new options:
```
  job_name: run_ngpb  # instead of default_run
  s: ngpb             # instead of delphi
```

The custom script ('run_ngpb.sh'), also listed in the log, wil be:
```
#!/bin/bash

step1.py prot.pdb
step2.py
step3.py -s ngpb
step4.py --xts

sleep 10
```

### Setup with a different topology folder

```
> mkdir pdbs_amber
> cd pdbs_amber

> bench_setup pkdb_pdbs -bench_dir . -ftpl amber_dir -job_name run_amber
```

#### Outcome

The benchmark.log will have recorded these new options:
```
  job_name: run_amber  # instead of default_run
  ftpl: amber_dir      # instead of "", which means mcce_exec/param/
```

The custom script ('run_amber.sh'), also listed in the log, will be: 
```
#!/bin/bash

step1.py prot.pdb
step2.py -ftpl amber_dir
step3.py -ftpl amber_dir
step4.py --xts

sleep 10
```

## 2.B Using the pdbs provided by the user
Typically, these pdbs do not have experimental pKas or they are not yet part of the pkDB.

### `user_pdbs` vs. `pkdb_pdbs` subcommands:
They have the same 'stepN.py' options, but differ on these three setup options:

<div style="display: block; float: left; padding-left:50px">
    
  |  pkdb_pdbs       | user_pdbs      | function |
  | :---             | :----          | :-----   |
  | _-n_pdbs_        | &#10006;       | Subset pKDB by count (N smallest) |
  | _-pdbids_file_   | &#10006;       | Subset pKDB by pdbids in file |
  | &#10006;         | _-pdbs_list_   | Provide source of pdbs as a pdb files folder or a file with pdbs fullpaths <br>|
</div>

<p><br><br><br><br><br><br></p>

__Example__:  
<div style="display: block; float: left;"></div>
<pre>
> mkdir pdbs_amber
> cd pdbs_amber <br>
> bench_setup user_pdbs -bench_dir . -pdbs_list PDBS_LIST
</pre>


## 3. Launch the mcce simulations over the set
The app has a 'batcher' and a 'scheduler', which uses the batcher on a cron schedule.  
The batcher is typically used for very small sets, or as a fallback command if there is a problem with  
the scheduler).
* Note: the sentinel_file option tells the program to mark a run has completed once this file is created.

### 3.1 Launching jobs with the 'batcher'
* Note: This batching option, is still 'manual': until the screen/log does not show the final count,
 the command has to be repeated.
```
> cd bench_dir

# with all default options:
> bench_batch -bench_dir .

# with non-default options:
> bench_batch -bench_dir . -job_name special -n_batch 12 -sentinel_file head3.lst
```
### Outcome
If you had a good estimate of the time needed to complete the runs, then re-running the command once more  
will likely produce the final counts, e.g.:
```
- Completed: 10
-    Failed: 2
-     Total: 12
```
If not, re-run the command until you do see it.

### 3.2 Launching jobs with the 'scheduler'
These commands will activate the cron scheduler, wich will launch the runs in batches of 10 (default):

### 3.2A Launch an all-default set of runs
```
> cd bench_dir
> bench_setup launch -bench_dir .
```

### 3.2B Launch an set of runs with customized options
```
> cd bench_dir
> bench_setup launch -bench_dir . -job_name special -n_batch 12
```

### 3.2C Monitoring of scheduled jobs
Because the app allows multiple benchmark sets to be scheduled, the crontab file will accumulate their scheduling definitions.  
The tail of the cron.log will display the final counts if the runs have completed, and will keep on outputing them until the  
corresponding crontab definition is removed.  
If the runs have completed in bench_dir X, then the cron job definition related to that folder in the crontab file can be deleted:
```
# cmd to edit contab:
contab -e
```

---
---
# Options passing in `step[1234].py` and `bench_setup`

## REMINDER: Expected usage for `stepx.py` command line options passing:
 * If using the `-load_runprm` option to provide a custom run.prm file:
   -  that runprm file should only contain the value/keys entries that need changing;
   -  if you only have a couple of entries in that file, consider using the `-u` option to pass them;
 * Several options are "shortcuts" that should only be set via their respective command line option identifer or flag, e.g. `--dry`, `level`;
 * Tip: try the glossary tool, to learn about each step keys, e.g. `glossary Step1`.

### Customizable input files: run.prm, name.txt, extra.tpl.
As with the 'legacy' `mcce` command, the `stepx.py` commands expect to find any of the customized input files in the run folder where 'prot.pdb' is.  
Additionally, the __names__ (not paths) given to any customized file in that folder must appear in the (customized) run.prm, i.e.: If using 'test_extra.tpl', then that name is the value of the `EXTRA` key in the local run.prm file.  

### Using customized input files in `bench_setup`
When using `bench_setup`, the program will need slightly different pieces of information as it will soft-link any customized input file into each of the proteins subfolder in the bench_dir/runs folder. 

#### Example:
Given this starting bench_dir folder structure, `bench_setup` will soft-link each into all protein run folders with these names:
```
.                                     PROT/
├── jts_extra.tpl                  -> extra.tpl
├── my_name.txt                    -> name.txt
└── long_wierd_name.prm            -> run.prm.user  # so it is not overwritten
```

#### Requirements:
 1. __Same as with `stepx.py`__: The local prm file MUST reference the EXTRA/RENAME files as they will be search for in the protein run folder, i.e. by their names only, e.g.: "extra.tpl", (__no path__).
 2. __New__: If needed, customized name.txt and extra.tpl files MUST be given via the `-u` option, i.e.:
    ```
    bench_setup pkdb_pdbs -bench_dir . -n_pdbs 2 -u EXTRA=jts_extra.tpl,RENAME_RULE=my_name.txt
    ```
    With that information, `bench_setup` will know what to link. it will then alter the -u option with the files' standard names before creating the run script.

---
---

## TODO: Complete Analyze

## TODO: Complete Compare
