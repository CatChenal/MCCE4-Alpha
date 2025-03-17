# Welcome to MCCE tutorial

MCCE starts with atomic coordinates in PDB format.  The output includes:
ALWAYS
  The protonation states of all Asp, Glu, Arg, Lys, His, Tyr, Cys (and His tautomer) at a pH of interest. 
  The pKa of these residues
  Information about why pKas are shifted within the protein from the values of the isolated resides in solution. 

ADDITIONAL ANALYSIS 
  Analaysis of microstates generated in Monte Carlo sampling. There are tools to handel the microstate output.
  
WITH APPROPRIATE TOPOLOGY FILE (.tpl)
  Redox potententials at equilbriu with residue and cofactor protonation states.
  Ion binding.
------------------------------------------------------------------------
This tutorial will show you how to set up, run and interprete the results for the PDB file 4lzt (Hen Egg White Lysozyme).  MCCE must be installed. See ____ for instilation instructions). 

EACH MCCE RUN IS DESIGNED TO RUN IN ITS OWN DIRECTORY.  If you run a new protein or change conditions make a new directory.

(If you're new to LINUX - Copy and past the following code into a terminal:
 > mkdir 4PTI 
 > cd 4PTI
 > getpdb 4PTI
 > ls
"mkdir" creates the new folder/directory "4PTI" to work in, and "cd" changes directory to 4PTI. "getpdb" is a utility distributed with MCCE that downlads the PDB file from RCSB.org, and "ls" lists the new files in our directory. 
)
--------------------------------------------------------------------

FIND OUT HOW MCCE VIEWS YOUR PROTEIN BEFORE YOU RUN IT

>p_info.py 4lzt.pdb

p_info.py will tell you
(1) how big your protein is.  Then look at  timing data _ to see if it's minutes or hours;  
(2) if there are ligands that the program doesn't recognize. (What to do with unkown ligands ___).
(3) If MCCE is modifying your protein.  For example, it will make an the N- and C termini ionizable.  

RUN MCCE

>run_mcce 4lzt.pdb > run.log &

_________________________________________________________________________________________
KEY OUTPUT FILES

fort.38 > sum_crg.out > pK.out

pK.out has information about how the protein shifts the pK or each residue.  To start you can edit the file with 

To get the pK simply 
>cut -c1-48 pK.out

or 
>cut -c1-48 pK.out | sort
The sort will grop each residue type.  In  s hoof tun most residues should have pKs near their values in solution.  All Arg should have high pKs.  

If we look at pK.out Lys 1. This titration is close to ideal. 
    Output : [LYS+A0001_        9.545     0.988     0.002]
The pK is 9.545, 
The n value is 0.935 (1 is ideal). (the Henderson Henderson–Hasselbalch equation (10^(n(pK-pH))/[1+10^(n(pK-pH))]
The sum of the squared differences of the calculated points from an ideal titration curve is 0.002 (1000*chi2).  

MCCE calculated the pK using the data in sum_crg.out (I've deleted the NTR line)
MCCE equilbrated the protein from pH 0 to 14. This gives the carge as a function of the pH.  The probability of being ionized as a fundtion of pH --> pK
 pH           0     1     2     3     4     5     6     7     8     9    10    11    12    13    14
LYS+A0001_  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  0.97  0.77  0.27  0.04  0.00  0.00  0.00

If chi^2 is greater than 5 you should graph sum_crg.out vs pH for that residue.  The pK fit to a single titration may not be appropriate for this residue. So the pK in pK.out is no long a good measure of the MCCE titration.  

Sum_crg.out give the charge at each pH.  It is derived from fort.38 which gives the % of each conformer at each pH.  Here there is only 1 ionized Lys conformand 1 neutral Lys conformer.  the sum of all conformers of a residue =1.00.  

ph              0.0   1.0   2.0   3.0   4.0   5.0   6.0   7.0   8.0   9.0  10.0  11.0  12.0  13.0  14.0
LYS01A0001_001 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.004 0.031 0.228 0.734 0.961 0.995 1.000 1.000
LYS+1A0001_002 1.000 1.000 1.000 1.000 1.000 1.000 1.000 0.996 0.969 0.772 0.266 0.039 0.005 0.000 0.000

The pK.out file for 4lzt should be close to this.  The values will not be exactly the same.  
pKas >14.0 (or <0.0) are out of the range of the titration. The titration range can be modified.  
> cut -1-48 pK.out
  pH             pKa/Em  n(slope) 1000*chi2     
NTR+A0001_        6.415     1.040     0.007     
LYS+A0001_        9.545     0.988     0.002     
ARG+A0005_       12.677     0.855     0.012     
GLU-A0007_        3.200     0.965     0.004     
LYS+A0013_       11.328     0.960     0.027     
ARG+A0014_       13.325     0.934     0.003     
HIS+A0015_        6.869     0.967     0.022     
ASP-A0018_        2.040     0.922     0.050     
TYR-A0020_       13.293     0.672     0.498     
ARG+A0021_       13.314     0.781     0.058     
TYR-A0023_       10.347     0.856     0.022     
LYS+A0033_       10.494     0.965     0.031     
GLU-A0035_        5.080     0.955     0.457     
ARG+A0045_       13.026     0.836     1.177     
ASP-A0048_        1.389     0.939     0.032     
ASP-A0052_        2.854     0.965     0.672     
TYR-A0053_        >14.0                         
ARG+A0061_       13.733     0.829     0.003     
ASP-A0066_        1.271     0.965     0.055     
ARG+A0068_        >14.0                         
ARG+A0073_       12.825     0.897     0.237     
ASP-A0087_        1.638     1.022     0.027     
LYS+A0096_       10.607     0.812     0.480     
LYS+A0097_       10.818     0.913     0.017     
ASP-A0101_        3.951     1.019     0.134     
ARG+A0112_       12.899     0.894     0.006     
ARG+A0114_       13.233     0.919     0.004     
LYS+A0116_        9.541     0.933     0.092     
ASP-A0119_        3.583     1.015     0.100     
ARG+A0125_       13.096     0.902     0.103     
ARG+A0128_       13.248     0.925     0.010     
CTR-A0129_        2.450     0.906     0.032     



____________________________________________________________________________________

MCCE IS MODULAR AND CAN BE STOPPED AND MODIFIED AFTER EACH OF THE 4 STEPS.
All steps should be run in 1 directory.  Step(n) needs output of step(n-1).
Step1: read in and fix the file
Step2: make conformers using the sinstructions in head1.lst
Step3: Calculat self and pairwise interactions for all conformes in step2_out.pdb
Step4: Run MC (can control pH or eH limits of titration; can fix individual conformers on or off)

  BREAKDOWN OF 4 MCCE STEPS

```
>step1.py 4pti.pdb --dry
```

"step1.py 4lzt.pdb" begins by reformatting the PDB file. The "--dry" flag removes waters, which we suggest you do if you are calculating pKs, Kds or Ems.  If you are looking for hydrogen bond networks you need the waters. 

```
>step2.py
```

"step2.py" finds "step1_out.pdb" 
Step 2 of MCCE makes conformers for all residues and ligands.  Topology files indicate the making of different protonation and redox states as well as whether ligands can leave the protein.  The default sidechain rotomer making is 'isosteric' which does not rotate C-C boinds in side chains.  This changes hydroxyl positions, flips ASN and GLN amide dide chain oritentation and makes both neutral His tautomers.    .

```
>step3.py
```

"step3.py" computes the energy lookup table for the protein, and is the most computationally expensive part of MCCE. 

```
step4.py --xts
ls
```

"step4.py" is the last step we will look at. It uses Monte Carlo sampling to compute pKa values based on the output of Step 3, head3.lst. You will now have about 24 files in your directory, but we're most interested in "pK.out".



))))))))))))))))))))))))))))))))))))
This should move into a new file

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
