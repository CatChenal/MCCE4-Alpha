#! /bin/bash -x

## stop on errors
set -e


#--Paths designed to assume users will run the codes at bin
in_py="step4_mcce_to_namd_Chun.py"
out_dir_prefix="../tests/test_MCCE_to_NAMD_Chun_ver/4lzt"
in_psfgen_script_paths="${out_dir_prefix}_0/MCCE_case_0_mcce_to_namd.tcl,${out_dir_prefix}_1/MCCE_case_0_mcce_to_namd.tcl"
salt_level=0.15
map_dimensions="300,300,300"
vmd_path="/Common/linux/bin/vmd"
apbs_path="/Common/linux/bin/apbs"
run_psfgen="False"
run_apbs="False"

#chmod +x $in_py

#-Junjun only need the below command
python $in_py $in_psfgen_script_paths -s $salt_level -m $map_dimensions -v $vmd_path -a $apbs_path -r $run_psfgen -p $run_apbs

