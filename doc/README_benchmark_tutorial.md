<!---..................................................................................................................|:120 --->
# TUTORIAL:
  * Notation remarks:
    1. ">" denotes the command line prompt (whatever it may be in your profile)
    2. ">>" denotes the output of a command.
    3. On a command line, "\<x\>" is an example of an argument value.

  ### PRE-REQUISITES:
  1. MCCE_Benchmarking pip-installed into a dedicated environment, e.g. `mce`;
  ```
  (mce) >pip install git+https://github.com/GunnerLab/MCCE_Benchmarking.git#egg=mcce_benchmark

  ```
  2. The environment is activated: `>conda activate mce`

---
---

# USE CASE 1:
---
# Run steps 1 and 2 on a user-defined set of "md frames proteins" to add a membrane (Raihan's case)

## 1. `cd` to the folder of interest
```
(mce) >cd /home/muddin/raihan_work/complexI/complexI_E-channel/md_frames/mq_all_frames

# List its contents (truncated):
(mce) >ls -l
drwxrwxr-x  3 muddin muddin   31 Jul 30  2023 fr0/
drwxrwxr-x  7 muddin muddin 4.0K Aug  2  2023 frame0/
drwxrwxr-x  7 muddin muddin 4.0K Jul 30  2023 frame10/
drwxrwxr-x  6 muddin muddin 4.0K Jul 30  2023 frame20/
drwxrwxr-x  6 muddin muddin 4.0K Jul 30  2023 frame30/
drwxrwxr-x  6 muddin muddin 4.0K Jul 30  2023 frame40/
drwxrwxr-x  6 muddin muddin 4.0K Jul 30  2023 frame50/
drwxrwxr-x  6 muddin muddin 4.0K Mar 19 01:44 frame60/
#[...]
```

## 2. Create a folder that will contain the selected pdbs:
```
(mce) >mkdir ./pdblist
```

## 3. Copy the selected pdbs into ./pdblist with a unique name:
```
# examples:

(mce) >cp frame0/pqrs.pdb ./pdblist/frame0.pdb
(mce) >cp frame10/pqrs.pdb ./pdblist/frame10.pdb

```

## 4. Once your selection is complete, list the ./pdblist/ folder for inspection:
```
(mce) >ls -l ./pdblist/
rwxrwxrwx2 muddin muddin     25 Mar 20 16:53 frame0.pdb
rwxrwxrwx2 muddin muddin     25 Mar 20 16:53 frame10.pdb
#[...]
```

## 5. Now you are ready to use the benchmarking command: `bench_setup user_pdbs`
The 'setup' mode prepares the pdbs in their own folders where mcce will be run, and creates a custom run script if requested.  
The default run script is the following:  
```
#!/bin/bash

step1.py --dry prot.pdb
step2.py -d 4
step3.py -d 4
step4.py --xts

sleep 10
```

### Main inputs (and the only ones if you are using the default script):
  1. A directory where the runs will be setup (does not have to exists);
  2. The name of the folder (or file) listing the pdbs to use, i.e. 'pdblist' in our example;

### Parameters for each of the 4 steps if not default

**Note**: Changing even just one option, say "-d 4" to "-d 8" makes the script a custom one, in which case a "job name" must be passed to the `-job_name" option.

In our current use case, Raihan needs to run steps 1 and 2 to add a membrane (to be placed in another protein), so we know that the final script will look like:
```
#!/bin/bash

step1.py --dry prot.pdb [..]
step2.py -d 4
step3.py --norun
step4.py --norun

sleep 10
```

### Additional, required command line options
  1. `-job_name`
  2. `-sentinel_file`: it is "pK.out" if you are running step4, here it must be "step2_out.pdb"
  3.  The "off switch" for steps 3 and 4: `--s3_norun` and `--s4_norun`
  4. `-u`: Accepts all "not often used" parameters, e.g. to add a membrane

The fully specified command line for this task is the following (press Enter to run it):
```
(mce) >bench_setup user_pdbs -bench_dir ./S12M -pdbs_list ./pdblist -job_name s12m -sentinel_file step2_out.pdb --s3_norun --s4_norun -u IPECE_ADD_ME=t,IPECE_MEM_THICKNESS=28

```

### Inspect the setup folder `S12M`:
```

(mce) >ls -l S12M
drwxr-xr-x 5 muddin muddin 4.0K Mar 17 13:42 runs/


(mce) >ls -l S12M/runs
drwxr-xr-x 5 muddin muddin 4.0K Mar 17 13:42 runs/
-rwxr-xr-x 5 muddin muddin   90 Mar  3 10:34 default_run.sh*
-rwxr-xr-x 5 muddin muddin   87 Mar 17 13:42 s12m.sh*
drwxr-xr-x 5 muddin muddin 4.0K Mar 17 13:51 FRAME0/
drwxr-xr-x 5 muddin muddin 4.0K Mar 17 13:51 FRAME10/
drwxr-xr-x 5 muddin muddin 4.0K Mar 17 13:55 FRAME20/
# [...]
-rw-r--r-- 5 muddin muddin   27 Mar 17 13:56 book.txt

```

### Inspect the run script `s12m.sh`:
(You could also modify it.)  
```

(mce) >cat S12M/runs/s12m.sh
#!/bin/bash

step1.py prot.pdb --dry -u IPECE_ADD_ME=t,IPECE_MEM_THICKNESS=28
step2.py -u IPECE_ADD_ME=t,IPECE_MEM_THICKNESS=28
step3.py --norun -u IPECE_ADD_ME=t,IPECE_MEM_THICKNESS=28
step4.py --xts --norun -u IPECE_ADD_ME=t,IPECE_MEM_THICKNESS=28

sleep 10

```
<br>

**Dev TODO**: For `ibench`: Make `-u` step-specific => mapping of keys per step needed.


### Launch the mcce runs in batches using `bench_batch` (default batch size is 10):

```
# first batch:

(mce) >bench_batch -bench_dir ./S12M -job_name s12m -sentinel_file step2_out.pdb

# monitor the bookkeeping file; a state of 'r' (running), means mcce has been launched in that folder:
(mce) >cat S12M/runs/book.txt
FRAME0    r
FRAME10   r
FRAME20   r
#[...]

```
Wait a few minutes...  
```
# second batch:
(mce) >bench_batch -bench_dir ./S12M -job_name s12m -sentinel_file step2_out.pdb

# book file, some runs show a completed state ('c'):
(mce) >cat small/runs/book.txt
FRAME0    c
FRAME10   r
FRAME20   c
#[...]

```
Wait a few minutes...  
```
# third batch:
(mce) >bench_batch -bench_dir ./S12M -job_name s12m -sentinel_file step2_out.pdb

# book file:
(mce) >cat small/runs/book.txt
FRAME0    c
FRAME10   c
FRAME20   c
#[...]

```

If all the folders listed in the book file have a state of 'c' (complete) or 'e' (error): they all have been processed.  
You can list one of them to verify that no mcce output file for step3 or step4 was created.


# This concludes Use Case 1!
---
---

# Command line options
Here is a short survey of the available options.

## 1. To create a set using the pdbs from pKaDBv1, use the setup command & sub-command `bench_setup pkdb_pdbs`.
Note: `bench_setup` creates the run script with the command line options for the mcce steps; if none are given, the default run script is used.
  * Minimal input `-bench_dir`:
    ```sh

    (mce) >bench_setup pkdb_pdbs -bench_dir <some/dir>  # -bench_dir is a required arg

    ```

    You can limit the number of curated proteins (folder) to setup using `-n_pdbs`.
    Give a low number for testing (max is 120, default):
    ```sh

    (mce) >bench_setup pkdb_pdbs -bench_dir <some/dir> -n_pdbs <2>

    ```

  * `-job_name` option
   The default job name is 'default_run'; you can change it using `-job_name`.
    - If you use `-job_name` without amending any of the steps args, the default_run.sh file will be soft-linked as <job name>;
    - If you pass one or more non-default mcce steps parameters, then you MUST have a job name (& not 'default_run'):
    ```sh

    (mce) >bench_setup pkdb_pdbs -bench_dir <some/dir> -n_pdbs <2> -d 8 -job_name <foo_e8>

    ```

  * `-sentinel_file` option:
  The default sentinel_file is **pK.out** (means the run completed Step 4). It is part of script setup
  to ensure it is deleted prior to launching. You need to include it if you are not running all 4 steps. Example running only steps 1 & 2:
  ```sh

  (mce) >bench_setup pkdb_pdbs -bench_dir <some/dir> -s3_norun True -s4_norun True -job_name <up_to_s2> -sentinel_file step2_out.pdb

  # OR
  (mce) >bench_setup user_pdbs -bench_dir <some/dir> -pdbs_list ./pdbs -s3_norun True -s4_norun True -job_name <up_to_s2> -sentinel_file step2_out.pdb

  ```

### DO NOT USE `--launch` or `bench_setup launch`: the scheduling cannot run due to issues in `Stable-MCCE` packaging.
Instead, launch a batch of 10 runs every time this command is used:

```sh

(mce) >bench_batch -bench_dir <some/dir> -job_name <foo_e8> [-n_batch if not 10]   # 10 is the default batch size

```

  * `--launch` option (flag) means "create a crontab entry to scheduling the processing of all pdbs in the set".
  Launching the scheduled processing of the set with `bench_setup` is possible with the `--launch` flag:
  ```sh

  (mce) >bench_setup pkdb_pdbs -bench_dir <some/dir> -d 8 -job_name <foo_e8> --launch

  # OR
  (mce) >bench_setup user_pdbs -bench_dir <some/dir> -pdbs_list ./pdbs -d 8 -job_name <foo_e8> --launch

  ```

  In this case, you CANNOT review or amend the run script. You would not launch the job this way if, for instance,  
  you want to add Step 6.
  **TODO**: Step 6 will be part of the script once a hydrogen-bond network comparison procedure is finalized.

## 2. To create a set using a user-provided list of pdbs, use the setup command & sub-command `bench_setup user_pdbs`.

  * Minimal input `-bench_dir`, `-pdbs_list`:
  ```sh

  (mce) >bench_setup user_pdbs -bench_dir <some/dir>  -pdbs_list ./path/to/dir or ./path/to/file/of/pdbs filepaths

  ```

  * `-job_name` option: same usage as with `pkdb_pdbs` sub-command

  * `-sentinel_file` option: same usage as with `pkdb_pdbs` sub-command

  * `--launch` option (flag): same usage as with `pkdb_pdbs` sub-command

## 3. To launch batches of jobs, use `bench_setup launch`
   * Use the automated scheduling to process all the pdbs in batches of size n_batch (default 10):
   ```
   (mce) >bench_setup launch -bench_dir <some/dir> -job_name <foo> -n_batch <4> [-sentinel_file x]  # needed if used in setup

   ```

## 4. Analyze (one set):

The command line interface for analyzing a job setup via `bench_analyze [pkdb_pdbs or user_pdbs]`; requires 1 argument: `-bench_dir`.
The first step of the analysis checks that all the runs are completed. Incidentally, it's an easy way to check how far along the processing is.

  * Output files & sub-commands:
    - user_pdbs: Collated pK.out and sum_crg.out, res and confs count files & related figures
    - pkdb_pdbs: Analysis additionally includes res and confs stats viz experimental pKa values and related figures.


## 5. Compare two sets of runs, i.e. A/B testing (convention: B is reference):

To compare two (completed) sets of runs, use `bench_compare`:
  * Options: -di1, -dir2 -o (output folder), and two flags:
    * `--user_db_pdbs`      # absence means 'pkdb_pdbs'
    * `--dir2_is_refset`    # if used, --user_pdbs must NOT be present, e.g.:
    ```
    (mce) >bench_compare -dir1 <d1> dir2 parse.e4 --dir2_is_refset -o <d1/comp>

    ```


## 6. Give feedback in repo:
Please, open an issue with any problem installing or using the app, along with enhancement requests on the [MCCE_Benchmarking issues](https://github.com/GunnerLab/MCCE_Benchmarking/issues) page. Thanks.

<!---..................................................................................................................|:120 --->
