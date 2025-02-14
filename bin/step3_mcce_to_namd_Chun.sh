#! /bin/bash -x

## stop on errors
set -e

in_py="step3_mcce_to_namd_Chun.py"
in_step1_log="mcce_to_namd_step_1.log"


#chmod +x $in_py

#-Junjun only need the below command
python $in_py $in_step1_log

