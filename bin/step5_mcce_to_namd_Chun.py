import argparse
from MCCE_to_NAMD_Chun_ver.MCCE_to_NAMD import PDBtoMSAna as PMSana
from MCCE_to_NAMD_Chun_ver.MCCE_to_NAMD import MSAnaToNAMD
from MCCE_to_NAMD_Chun_ver.MCCE_to_NAMD import AggMCCEResults as amr
from MCCE_to_NAMD_Chun_ver.Vis_MCCE_on_VMD import Vis_MCCE_on_VMD as vmcce
from APBS_for_MCCE_Chun_ver.APBS_for_MCCE import VisEMfromMCCEtoNAMD as vapbs
from pathlib import Path
import os


#--Things will be run locally.

def arg_parsing_for_step_5():
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
        "ms_base_dir", type=Path,
        help='''the directory to hold MS analyses outputs that follow each MCCE run.
Outputs for the i th MCCE run will be deposited in the subdirectory "i" in the provided directory here.'''
    )
    #--Optional arguments
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
        "-n", "--number_of_mcce_runs", type=int, default=int(1),
        help="The number of MCCE runs to be performed. (Default: 1)"
    )
    parser.add_argument(
        "-e", "--residue_of_interest", default='all',
        help="a list of residues selected for MCCE runs (default: all documented residues)."
    )
    parser.add_argument(
        "-m", "--num_of_ms_picked", default=int(5),
        help="the number of nth most popular microstates selected from each MCCE run."
    )
    parser.add_argument(
        "-s", "--pbs_solver_choice", default="delphi",
        help="the possion-boltzman solver to be used."
    )
    parser.add_argument(
        "-l", "--salt_level", type=float, default=0.15,
        help="the salinity level for the system (unit: M)."
    )
    parser.add_argument(
        "-d", "--map_dimensions", type=str, default='100,100,100',
        help="the dimension of the EM map to be computed (a list of float numbers in Å for x, y, z dimension.)."
    )
    parser.add_argument(
        "-t", "--topology_files", type=str, default="",
        help="a list of paths, each to a sourced Charmm topology file apart from the default Charmm36FF."
    )
    parser.add_argument(
        "-v", "--vmd_path", type=Path, default='/Common/linux/bin/vmd',
        help="the path to soruce VMD."
    )
    parser.add_argument(
        "-a", "--apbs_path", type=Path, default='/Common/linux/bin/apbs',
        help="the path to soruce APBS."
    )
    parser.add_argument(
        "-u", "--run_psfgen", default=False,
        help="to run psfgen in VMD to generate MD-ready structure."
    )
    parser.add_argument(
        "-b", "--run_apbs", default=False,
        help="to run APBS to get EM maps."
    )

    args = parser.parse_args()
    if not args.topology_files:
        args.topology_files = None
    else:
        args.topology_files = [elm.strip() for elm in args.topology_files.split(',')]
    args.map_dimensions = [float(elm.strip()) for elm in args.map_dimensions.split(',')]
    args.run_psfgen = bool(args.run_psfgen)
    args.run_apbs = bool(args.run_apbs)
    return args


def main():

    args = arg_parsing_for_step_5()

    #--Part 1: Running MCCEs and creating NAMD-ready strutcures.
    print(f"Files and/or runs for {args.number_of_mcce_runs} MCCE calcualtions will be initiated.")

    info_summary = []
    for run_id in range(args.number_of_mcce_runs):

        ms_base_dir = os.path.join(args.ms_base_dir, f"{run_id + 1}")
        print(f"Outputs for round {run_id + 1} of MCCE run will be deposited at {ms_base_dir}.")

        pms_obj = PMSana(
            in_pdb_path=args.in_pdb_path,
            to_mcce_path=args.to_mcce_path,
            ms_base_dir=ms_base_dir,
            pbs_solver_choice=args.pbs_solver_choice,
        )
        #--Part 2: Running MS analyses.
        if pms_obj:
            if args.run_mcce:
                print("Executing MCCE.")
            pms_obj.run_mcce(
                pH_value=args.pH_value,
                run_mcce=args.run_mcce,
            )
            signal = pms_obj.run_MSana(
                residue_of_interest=args.residue_of_interest,
                num_of_ms_picked=args.num_of_ms_picked,
            )
            #--Part 3: Gathering information as a log file
            if signal:
                info_summary.append(pms_obj.out_csv_for_namd)

    if info_summary:
        info_summary.append(str(args.in_pdb_path))
    info_summary_path = os.path.join(args.ms_base_dir, "mcce_to_namd_Chun.log")
    print(f'''
Log file at {info_summary_path}.
''')
    with open(info_summary_path, 'w') as fout:
        fout.write('\n'.join(info_summary))

    if info_summary:
        #--Part 4: Write psfgen file & vis file
        csv_paths = info_summary[:-1]
        in_pdb_path = info_summary[-1]
        default_topo_dir = os.path.join(os.path.dirname(args.to_mcce_path), 'param_Chun')
        for ind, in_ms_csv in enumerate(csv_paths):
            out_dir = os.path.dirname(in_ms_csv)
            psf_obj = MSAnaToNAMD(
                in_ms_csv=in_ms_csv,
                in_pdb=in_pdb_path,
                in_topologies=args.topology_files,
                protein_only=True,
                default_topo_dir=default_topo_dir
            )
            psf_obj.parse_info_from_ms_ana(
                number_of_chosen_cases=args.num_of_ms_picked
            )
            #psf_obj.print_action_messages()
            psf_obj.write_psfgen_file(out_dir=out_dir)
            print(f"NAMD preparation scripts for run {ind} are deposited at {out_dir}.")

            in_pK_out_path = os.path.join(out_dir, 'ms_out', "pK.out")
            in_res_pair_lst_path = os.path.join(out_dir, 'ms_out', "respair.lst")
            vis_obj = vmcce(
                in_pdb_path=in_pdb_path,
                vis_option='all',
                in_ms_ana_csv_path=in_ms_csv,
                in_pK_out_path=in_pK_out_path,
                in_res_pair_lst_path=in_res_pair_lst_path,
                out_dir=out_dir,
            )
            vis_obj.gen_VMD_vis_script()

        #--Part 5: Generate APBS or EM Map
        in_psfgen_script_paths = []
        for in_csv_path in info_summary[:-1]:
            in_dir = os.path.dirname(in_csv_path)
            files = os.listdir(in_dir)
            in_psfgen_scripts = [os.path.join(in_dir, elm) for elm in files if 'mcce_to_namd' in elm]
            in_psfgen_script_paths.extend(in_psfgen_scripts)

        for ind, in_psfgen in enumerate(in_psfgen_script_paths):
            print(in_psfgen)
            vEM_obj = vapbs(
                in_psfgen_path=in_psfgen,
                vmd_path=args.vmd_path,
                apbs_path=args.apbs_path,
                map_dimensions=args.map_dimensions,
                salinity=args.salt_level,
            )
            vEM_obj.get_APBS_EM_map_and_Vis()

        #--Part 6: Aggregated results (psfgen only)
        print(f"An aggregrated version for results from all MCCE runs will now be used to prepare NAMD.")
        amr_obj = amr(
            in_ms_csvs=csv_paths,
            out_dir=args.ms_base_dir,
            number_of_chosen_cases=args.num_of_ms_picked,
        )
        _, out_paths = amr_obj.run()
        for ind, in_ms_csv in enumerate(out_paths):
            agg_out_dir = os.path.join(args.ms_base_dir, f"Agg_result_{ind}")
            print(f"{ind + 1} aggregrated results will be deposted at {agg_out_dir}.")
            if not os.path.exists(agg_out_dir):
                os.makedirs(agg_out_dir)
            psf_obj = MSAnaToNAMD(
                in_ms_csv=in_ms_csv,
                in_pdb=in_pdb_path,
                in_topologies=args.topology_files,
                protein_only=True,
                default_topo_dir=default_topo_dir
            )
            psf_obj.parse_info_from_ms_ana(
                number_of_chosen_cases=args.num_of_ms_picked
            )
            #psf_obj.print_action_messages()
            psf_obj.write_psfgen_file(out_dir=agg_out_dir)



if __name__ == "__main__":
    main()

