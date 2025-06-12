#---This step focuses on running necessary MCCE steps for NAMD preparation.
#--1) This includes (i) running MCCE;
#--2) (ii) running subsequent analyses;

import argparse
from MCCE_to_NAMD_Chun_ver.MCCE_to_NAMD import PDBtoMSAna as PMSana
from pathlib import Path
import os


def arg_parsing_for_step_1():
    #--Argument setup
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "in_pdb_path", type=Path,
        help="the path to the input pdb for MCCE runs"
    )
    parser.add_argument(
        "to_mcce_path", type=Path,
        help="the path to the mcce executible"
    )
    parser.add_argument(
        "out_dir", type=Path,
        help="the directory to hold relevnt outputs for MCCE runs (An index will be attached to the value for this argument to reflect the number of MCCE runs initiated.)"
    )
    parser.add_argument(
        "-p", "--pH_value", type=float, default=7.0,
        help="the pH-value to be accessed in MCCE calculations."
    )
    parser.add_argument(
        "-r", "--run_mcce",
        action="store_true", default=False,
        help="to run MCCE calculations or just to generate related scripts."
    )
    parser.add_argument(
        "-y", "--python_path", default='python',
        help="the alias command or the path to the python to be called."
    )
    parser.add_argument(
        "-n", "--number_of_mcce_runs", type=int, default=int(1),
        help="The number of MCCE runs to be performed. (Default: 1)"
    )
    parser.add_argument(
        "-e", "--residue_of_interest", default='default',
        help="a list of residues selected for MCCE runs (default: all documented residues)."
    )
    parser.add_argument(
        "-m", "--num_of_ms_picked", default=int(5),
        help="the number of nth most popular microstates computed from MCCE runs."
    )
    parser.add_argument(
        "-s", "--pbs_solver_choice", default="delphi",
        help="the possion-boltzman solver to be used."
    )

    args = parser.parse_args()
    return args


def main():

    args = arg_parsing_for_step_1()
    #create_default_values(args)

    #--Part 1: Running MCCEs and creating NAMD-ready strutcures
    print(f"Files and/or runs for {args.number_of_mcce_runs} MCCE calcualtions will be initiated.")
    out_dirs = []

    for run_id in range(args.number_of_mcce_runs):

        out_dir = str(args.out_dir) + f"_{run_id}"
        print(f"Run {run_id}. Output will be deposited at {out_dir}.")

        pms_obj = PMSana(
            in_pdb_path=args.in_pdb_path,
            to_mcce_path=args.to_mcce_path,
            out_dir=out_dir,
            python_path=args.python_path,
            pbs_solver_choice=args.pbs_solver_choice,
        )
        pms_obj.run_mcce(
            pH_value=args.pH_value,
            run_scripts=args.run_mcce,
        )
        pms_obj.run_MSana(
            residue_of_interest=args.residue_of_interest,
            num_of_ms_picked=args.num_of_ms_picked,
        )

        #out_dirs.append(out_dir)
        out_dirs.append(pms_obj.out_csv_for_namd)

    out_dirs.append(str(args.in_pdb_path))
    output_paths_record = os.path.join(os.path.dirname(args.to_mcce_path), "mcce_to_namd_step_1.log")
    with open(output_paths_record, 'w') as fout:
        fout.write('\n'.join(out_dirs))



if __name__ == "__main__":
    main()

