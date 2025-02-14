# Welcome to MCCE

Multi-Conformational Continuum Electrostatics, or MCCE, is a computational biophysics tool that accepts Protein Data Bank (PDB) files and outputs useful statistics, including pKa values. In this manual, we will cover a few basic uses of these tools designed to help first users start using MCCE in earnest. For the purposes of this document, we assume that MCCE is installed. 

First, let's look at tools at the individual protein level. Copy and past the following code into a terminal:

```
mkdir 4PTI 
cd 4PTI
getpdb 4PTI
ls
```

"mkdir" creates the new folder/directory "4PTI" to work in, and "cd" changes directory to 4PTI. Then we use an executable from MCCE called "getpdb" to download a file representing a protein from RCSB.org, and "ls" lists the new files in our directory. You should now have two files, "4pti.pdb", and "protinfo.log" in 4PTI. Now run the command

```
ProtInfo 4pti.pdb
```

This will output information of interest before a protein run. MCCE runs may take significant computation time, so for future runs pay attention to the number of Residues in the terminal output. 4pti is a small protein with 58 residues, which we have chosen for ease of demonstration. Now, let's start a run with MCCE:

```
step1.py 4pti.pdb --dry
```

"step1.py 4pti.pdb" begins by reformatting the PDB file. The "--dry" flag removes waters before processing, which we have done here to speed calculation slightly. Next, enter

```
step2.py
```

"step2.py" automatically looks for a file named "step1_out.pdb" in the local directory, which is why each protein should have its own directory. Step 2 of MCCE make side chain conformers from Step 1's output.

```
step3.py
```

"step3.py" computes the energy lookup table for the protein, and is the most computationally expensive part of MCCE. 

```
step4.py --xts
ls
```

"step4.py" is the last step we will look at. It uses Monte Carlo sampling to compute pKa values based on the output of Step 3, head3.lst. You will now have about 24 files in your directory, but we're most interested in "pK.out".

```
awk '{ printf "%-12s %-12s %-12s %-12s\n", $1, $2, $3, $4 }' pK.out
```

This command prints the first four columns of pK.out, specifically, the residues, calculated pH values, slope, and chi-squared values for our original "4pti.pdb" file. 



Now that we understand how this works for one protein, we might ask "is there a way to do this for many proteins at once?", or "can I check if MCCE's pKa values are correct?" Yes we can, with MCCE's benchmarking tools. let's exit our 4PIT directory:

```
cd ../
```

The command "bench_setup" allows us to quickly set up a number of PDB directories:

```
bench_setup pkdb_pdbs -bench_dir basic_bench -n_pdbs 3
```

This command takes the three smallest PDB files from pkdb_pdbs (which contains useful information like experimental pKas), and creates a nested directory. Let's cd into it and check it out:

```
cd basic_bench
ls
```

The most important folder here is "runs". After we analyze the benchmarks, the analysis folder will be in this directory. But for now let's cd futher in runs:

```
cd runs
ls
```

The runs folder contains folders for each of the proteins we downloaded, as well a "book.txt" and "default_run.sh". "book.txt" keeps track of the status of runs, if they are running, completed, initialized, or an error has occured.

```
cat default_run.sh
```

This is the .sh, or shell script, that will be run on each protein file in the PDB directories. This can be edited directly, or a non-default .sh file can be chosen during setup with the -job_name flag. However, we will leave it as is for now. Let's try running "bench_batch" on our proteins. To avoid cluttering the terminal, we have used the command "nohup" so the output of bench_batch goes there, and the "&" at the end of the command so the command runs in the background. The "." refers to the local directory. It's possible to run bench_batch from outside the directory to run, but it's recommended to do it this way so the "nohup.out" file is always within its associated directory.

```
cd ../ # go back to the bench_batch directory
ls
nohup bench_batch -bench_dir . -n_batch 3 &
```

After two or minutes, the runs should be completed. Run the following to update "book.txt". 

```
bench_batch -bench_dir .
```

For our analysis to run, we need all entries changed to 'c'. You should get a message like the following:

```
WARNING! Some functionality may not work with python3 version minor other than 10.
INFO - batch_run: Read book entries
INFO - batch_run: Running jobs: 0
INFO - batch_run: Changed 1ANS: 'r' -> 'c'
INFO - batch_run: Changed 1BI6: 'r' -> 'c'
INFO - batch_run: Changed 1HIC: 'r' -> 'c'
INFO - launch_cli: All jobs processed. Final counts:

- Completed: 3
-    Failed: 0
-     Total: 3
```

If the runs have completed with no errors, we can run bench_analyze:

```
bench_analyze pkdb_pdbs -bench_dir .
```

This creates a new folder, analysis. Let's go inside:

```
cd analysis
ls
```

There's a variety of files here, but let's inspect pkas_fit.png. This file creates a helpful graph comparing experimentally verified pKa values with those computed by MCCE:

```
mimeopen pkas_fit.png # or the image opener of your choice, 'feh' can work well 
```

The same information in text form is available in "matched_pkas.txt", useful if you want to store data in a spreadsheet or other text form. Again, this information is available to MCCE because we are using proteins in the pkdb_pdbs category.

```
cat matched_pkas.txt
```


First, we need a text file containing the proteins we want to use. Use Vim to create an empty text file:

```
cd ../../ # exit analysis and basic_bench folder
vi base_pdbs.txt
```

Copy (Ctrl-C) the following text and paste it (Ctrl-V) into the empty text file created in the terminal. We'll be using the same proteins as our other folder so we can demonstrate bench_compare later.

```
PDBID
1ANS
1BI6
1HIC
```

Now, press the ESCAPE key, and type ":wq" to save and return to the previous directory. Then we'll set up a new directory with our custom set of proteins. This time, we'll use the "--launch" flag to immediately begin the runs, instead of having to run bench_batch seperately.

```
bench_setup pkdb_pdbs -bench_dir custom_bench -pdbids_file base_pdbs.txt --launch
```

It may take three or four minutes for the runs to process. Again, we use the "bench_batch" command to check:

```
bench_batch -bench_dir .
```

After they've all completed, let's analyze the results:

```
bench_analyze pkdb_pdbs -bench_dir .
```

Now, we'll cd out of the directory and try a bench_compare:

```
bench_compare -dir1 basic_bench -dir2 custom_bench -o .
```

Because the directories were made with the same proteins and instructions, variation between the two should be minimal, only whatever differences arising from random chance during MCCE's steps.

Another thing you might want to do is use your own custom PDB files. To do this, we'll use bench_setup with the "user_pdbs" flag. But first, we need to create a text file to direct the system to our custom pdbs. cd to a clean directory, and we'll use vim to create a new text file.

```
vi custom_pdbs.txt
```

We're giving MCCE PDB files it may not have in its system, so the text file will be formatted to give MCCE the paths to the files we want:

```
# Example of PDB file paths
/home/user/pdb_storage/custom_1.pdb
/home/user/pdb_storage/custom_2.pdb
```

Suppose the above file was named "custom_pdbs.txt". Then I could run bench_setup with those custom pdbs, using the following command:

```
bench_setup user_pdbs -bench_dir user_bench -pdbs_list custom_pdbs.txt
```

This establishes a directory "user_bench" with a runs folder containing our custom PDB files, which you can then run bench_batch as normal. Note that bench_analyze needs the "user_pdbs" flag:

```
bench_analyze user_pdbs -bench_dir .
```
