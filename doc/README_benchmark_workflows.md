# MCCE_Benchmarking app

## Four main choices in MCCE_Benchmarking:
### 1. Create a new benchmarking set (with a "launch right away" option)
### 2. Launch the automated runs for the set
### 3. Analyze one set of runs
### 4. Compare two sets of runs (after each set is analyzed)

__All these correspond to command line "entry points", which the flowcharts below overview.__  

## Flowchart for the current entry points:

### `bench_setup + [user_pdbs or pkdb_pdbs sub-command]` handles the data & script setup:
```mermaid
%%{init: {"flowchart": {"htmlLabels": false}} }%%
flowchart LR
    one(("`__1__`")) ~~~ create["`__Create a new benchmarking set__`"] --> pdbs{Use pdbs in pKaDB?}
    pdbs -->|N| l1{launch right away?}
    l1 -->|Y| mcbY("`__:: bench_setup user_pdbs__ <args> --launch`")
    l1 -->|N| mcbN("`__:: bench_setup user_pdbs__ <args>`")
    pdbs -->|Y| phDB{"`__pH?__`"}
    phDB -->|Y| l2{launch right away?}
    l2 -->|Y| explY("`__:: bench_setup pkdb_pdbs__ <args> --launch`")
    l2 -->|N| explN("`__:: bench_setup pkdb_pdbs__ <args>`")
    phDB -->|N| nogo["`#10006; __Invalid choice:__
    pKaDB pdbs have no cofactors`"]
```

### `bench_setup launch` is used if `--launch` was not used during data & script setup:
```mermaid
%%{init: {"flowchart": {"htmlLabels": false}} }%%
flowchart LR
    two(("`__2__`")) ~~~ launch["`__Launch automated processing of the set__`"] --> do("`:: __bench_setup launch__ <args>`")
```

### `bench_analyze` does an "intra set" analysis; output folder: `-bench_dir` path/analysis:

```mermaid
%%{init: {"flowchart": {"htmlLabels": false}} }%%

flowchart LR
    thre(("`__3__`")) ~~~ anlyze["`__Analyze 1 set__`"] --> orig{Set created
     by app?}
    orig -->|Y| pdbs{Set pdbs
    from pKaDB?}
    pdbs -->|N| mcb("`__:: bench_analyze user_pdbs__`")
    pdbs -->|Y| phDB{"`__pH?__`"}
    phDB -->|Y| expl("`__:: bench_analyze pkdb_pdbs__`")
    phDB -->|N| nogo1["`#10006; __Invalid choice__:
    pKaDB data is pKas only`"]
    orig -->|N| nogo2["`#10006; __Not possible__
    Expected structure:
    .'bench_dir'
    |--runs/
      |--PDB1/
         |--run.prm.record
      [...]
      --book.txt`"]
```

### `bench_compare` does an "inter sets" analysis; output folder: `-o` value:

```mermaid
%%{init: {"flowchart": {"htmlLabels": false}} }%%

flowchart LR
    four(("`__4__`")) ~~~ comp["`__Compare 2 sets__`"] --> has_analysis{Each set
    was analyzed?}
    has_analysis -->|N| stop("`__STOP__: Run _bench_analyze_
    on each set first.`")
    has_analysis -->|Y| orig{Sets created
     by app?}
    orig -->|Y| pdbs{Sets pdbs
    from pKaDB?}
    pdbs -->|N| user("`__:: bench_compare -dir1 A -dir2 B --user_pdbs__`")
    pdbs -->|Y| phDB{"`__pH?__`"}
    phDB -->|Y| ref{Is dir2 a
    reference set?
    i.e. parse.e4}
    ref -->|N| pkdbN("`__:: bench_compare -dir1 A -dir2 B__`")
    ref -->|Y| pkdbY("`__:: bench_compare -dir1 A -dir2 B --dir2_is_refset__`")
    phDB -->|N| nogo1["`#10006; __Invalid choice__:
    pKaDB data is pKas only`"]
    orig -->|N| nogo2["`#10006; __Not possible__
    Expected structure:
    .'bench_dir'
    |--runs/
      |--PDB1/
         |--run.prm.record
      [...]
      --book.txt`"]
    four ~~~ note["`__Important__
    In __A/B__ testing, B is the reference set, i.e:
    __dir1::A; dir2::B__
    A vs. B
    A - B
    New vs. existing
    Y vs. X`"]
```
  * Note: The flag `--dir2_is_refset` indicates that 'dir2' is the _name_ of a packaged reference dataset, currently __parse.e4__ (pH titrations using the pKaDBv1 pdbs); without it, dir2 is a path.


### Notes on processing:

In the flowcharts above, __`bench_setup launch`__ means starting the batch-processing of the entire set (via creation of a crontab entry).  
This can be done two ways:
  1) During `bench_setup [user_pdbs or pkdb_pdbs]` with the --launch flag;
  2) Using `bench_setup launch` + args;

In case there is a problem with the automated scheduling, the processing can still be done 'batch by batch' at the command line using `bench_bath`, the entry point for batching (which is what the crontab uses).  
After activating the conda env where MCCE_Benchmarking is installed, provide -job_name if your script setup was not default:

```
#  The batch size can be changed every time the cli is called:
> cd A
> bench_batch -bench_dir . -n_batch 15

```

#### Monitor the state of processing by listing the sentinel file:
```
# pattern:
# > ls -l A/runs/*/<sentinel_file> | wc -l  # how many runs have completed the step that creates <sentinel_file>?
  > ls -l A/runs/*/pK.out | wc -l           # how many runs have completed step4?
```
