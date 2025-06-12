#---This step focuses on preparing visualization scritps from previous MCCE runs.

import argparse
from MCCE_to_NAMD_Chun_ver.Vis_MCCE_on_VMD import Vis_MCCE_on_VMD as vmcce
from pathlib import Path
import os


def arg_parsing_for_step_3():
    #--Argument setup
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "in_step1_log", type=Path,
        help="the path to an the ouput log file of step1_mcce_to_namd_Chun, namely: 'mcce_to_namd_step_1.log'"
    )
    parser.add_argument(
        "-y", "--python_path", default='python',
        help="the alias command or the path to the python to be called (Default: python)"
    )
    args = parser.parse_args()
    return args


def parse_log_file(log_file_path):
    with open(log_file_path, 'r') as fin:
        lines = fin.readlines()
    csv_paths = [elm.strip() for elm in lines[:-1]]
    pdb_path = lines[-1].strip()
    return csv_paths, pdb_path


def main():

    args = arg_parsing_for_step_3()

    csv_paths, in_pdb_path = parse_log_file(args.in_step1_log)

    step_text = "step1_mcce_to_namd_Chun.py"
    print(f"Results from each MCCE run in {step_text} will be used to prepare VMD visualizations.")

    for ind, in_ms_csv in enumerate(csv_paths):

        out_dir = os.path.dirname(in_ms_csv)
        in_pK_out_path = os.path.join(out_dir, 'ms_out', "pK.out")
        in_res_pair_lst_path = os.path.join(out_dir, 'ms_out', "respair.lst")
        #print(in_pK_out_path)
        #print(in_res_pair_lst_path)
        vis_obj = vmcce(
            in_pdb_path=in_pdb_path,
            vis_option='all',
            in_ms_ana_csv_path=in_ms_csv,
            in_pK_out_path=in_pK_out_path,
            in_res_pair_lst_path=in_res_pair_lst_path,
            out_dir=out_dir,
        )
        vis_obj.gen_VMD_vis_script()



if __name__ == "__main__":
    main()

