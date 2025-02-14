#! /bin/bash -x

## stop on errors
set -e


in_py="step1_mcce_to_namd_Chun.py"
in_pdb_path='../tests/test_MCCE_to_NAMD_Chun_ver/4lzt/4lzt.pdb'
to_mcce_path='./mcce'
out_dir='../tests/test_MCCE_to_NAMD_Chun_ver/4lzt'
number_of_mcce_runs="2"
pH_value="7.0"
pbs_solver_choice="delphi_pgf90"


#chmod +x $in_py


#-Junjun only need the below command
python $in_py $in_pdb_path $to_mcce_path $out_dir --run_mcce -n $number_of_mcce_runs -p $pH_value -s $pbs_solver_choice
