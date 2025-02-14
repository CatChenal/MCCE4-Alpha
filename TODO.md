# What needs to be done to have a working version (without step6)

## General
 1. [x] I’d like to have fort.38 in step 4 and sum_crg.out; pK.out etc in step 5 **already have it in step 4, script to do it again**
 2. [x] The corrected sum_crg (using head3.lst for charge)
 3. [x] the mfe analysis in pK.out (and an ability to change the mfe point from pK to a pH by retuning step 5 **mfe.py**
 4. [x] Bug, step1 can not handle NTR and CTR correctly when the NTR and CTR atoms are out of order.
 5. [ ] Change in chemical potential of a ligand chosen in run.prm, **temp fix script, automated**

## Error checking
 1. [x] `postrun` tool: write out for Asp, Glu, Arg, Lys: which ones are <90% ionized and the total number (if you have a lot of neutral Arg the run is likely to be bad) **A tool or step 4, print out abnoramlies, chi2, n etc**
 2. [ ] More than the maximum number of conformers **Why the error happens and not printed out clearly**
 3. [ ] The number of conformers in all input files (head3.lst; opp etc) don’t match **Check the consistency**
 4. [ ] The program stops and returns an error message **From users' feedback**

## Improvement
 1. [x] Do not pass 0 radius H to delphi to see if this fixes delphi surface error. **Implement in the next version**
 2. [x] Fix delphi run time error.
 3. [x] Strip down and reorganize run.prm.
 4. [x] Have full.prm (with all possibilities) on Wiki. **Good idea, group options that default to a pre-defined choices**
 5. [x] One run.prm with a toggle for run quick or default (rather than 2 different basic input files) **Questionaire for determining the run.prm**

## Open questions
 1. [x] Proton naming.
 2. [x] New tpl files, **take both tpl, ftpl**
 3. [ ] Python script to show titration.
 4. [ ] Phi Map instructions **Ask Dyvia? On a single delphi input file from microstate, or most occupied confomer structure**
