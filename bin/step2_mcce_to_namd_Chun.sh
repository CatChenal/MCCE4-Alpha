#! /bin/bash -x

## stop on errors
set -e


in_py="step2_mcce_to_namd_Chun.py"
in_step1_log="mcce_to_namd_step_1.log"
out_dir='../tests/test_MCCE_to_NAMD_Chun_ver/4lzt'
topology_files="../tests/test_MCCE_to_NAMD_Chun_ver/Topology_files/top_all36_prot.rtf"
num_of_ms_picked="5"


#chmod +x $in_py

#-Junjun only need the below command
python $in_py $in_step1_log $out_dir $topology_files -m $num_of_ms_picked

