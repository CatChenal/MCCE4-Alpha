#---This step focuses on running necessary MCCE steps for NAMD preparation.
#--1) This includes (i) running MCCE;
#--2) (ii) running subsequent analyses;

import argparse
from MCCE_to_NAMD_Chun_ver.MCCE_to_NAMD import MSAnaToNAMD
from MCCE_to_NAMD_Chun_ver.MCCE_to_NAMD import AggMCCEResults as amr
from pathlib import Path
import os


def arg_parsing_for_step_2():
    #--Argument setup
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "in_step1_log", type=Path,
        help="the path to an the ouput log file of step1_mcce_to_namd_Chun, namely: 'mcce_to_namd_step_1.log'"
    )
    parser.add_argument(
        "out_dir", type=Path,
        help="the directory to hold aggregated outputs for MCCE runs."
    )
    parser.add_argument(
        "topology_files", type=str,
        help="a list of paths, each to a sourced Charmm topology file"
    )
    parser.add_argument(
        "-y", "--python_path", default='python',
        help="the alias command or the path to the python to be called (Default: python)"
    )
    parser.add_argument(
        "-e", "--residue_of_interest", default='default',
        help="a list of residues selected for MCCE runs (default: all documented residues)"
    )
    parser.add_argument(
        "-m", "--num_of_ms_picked", type=int, default=int(5),
        help="the number of nth most popular microstates computed from MCCE runs (default: 5)"
    )
    args = parser.parse_args()
    args.topology_files = [elm.strip() for elm in args.topology_files.split(',')]
    return args


def parse_log_file(log_file_path):
    with open(log_file_path, 'r') as fin:
        lines = fin.readlines()
    csv_paths = [elm.strip() for elm in lines[:-1]]
    pdb_path = lines[-1].strip()
    return csv_paths, pdb_path


def main():

    args = arg_parsing_for_step_2()

    csv_paths, in_pdb_path = parse_log_file(args.in_step1_log)

    step_text = "step1_mcce_to_namd_Chun.py"
    print(f"Results from each MCCE run in {step_text} will be used to prepare NAMD.")

    for ind, in_ms_csv in enumerate(csv_paths):
        out_dir = os.path.dirname(in_ms_csv)
        psf_obj = MSAnaToNAMD(
            in_ms_csv=in_ms_csv,
            in_pdb=in_pdb_path,
            in_topologies=args.topology_files,
            protein_only=True
        )
        psf_obj.parse_info_from_ms_ana(
            number_of_chosen_cases=args.num_of_ms_picked
        )
        #psf_obj.print_action_messages()
        psf_obj.write_psfgen_file(out_dir=out_dir)
        print(f"NAMD preparation scripts for run {ind} are deposited at {out_dir}.")

    print(f"An aggregrated version for results from {step_text} will now be used to prepare NAMD.")
    amr_obj = amr(
        in_ms_csvs=csv_paths,
        out_dir=args.out_dir,
        number_of_chosen_cases=args.num_of_ms_picked,
    )
    _, out_paths = amr_obj.run()
    for ind, in_ms_csv in enumerate(out_paths):
        agg_out_dir = str(args.out_dir) + f"/Agg_result_{ind}"
        os.makedirs(agg_out_dir)
        psf_obj = MSAnaToNAMD(
            in_ms_csv=in_ms_csv,
            in_pdb=in_pdb_path,
            in_topologies=args.topology_files,
            protein_only=True
        )
        psf_obj.parse_info_from_ms_ana(
            number_of_chosen_cases=args.num_of_ms_picked
        )
        #psf_obj.print_action_messages()
        psf_obj.write_psfgen_file(out_dir=agg_out_dir)



if __name__ == "__main__":
    main()

