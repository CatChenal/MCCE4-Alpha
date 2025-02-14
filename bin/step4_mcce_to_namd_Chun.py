#---This step focuses on preparing visualization scritps from previous MCCE runs.

import argparse
#from MCCE_to_NAMD_Chun_ver.Vis_MCCE_on_VMD import Vis_MCCE_on_VMD as vmcce
from APBS_for_MCCE_Chun_ver.APBS_for_MCCE import VisEMfromMCCEtoNAMD as vapbs
from pathlib import Path
import os


def arg_parsing_for_step_4():
    #--Argument setup
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "in_psfgen_script_paths", type=str,
        help="a list of paths to psfgen scripts generated in step 2."
    )
    parser.add_argument(
        "-s", "--salt_level", type=float, default=0.15,
        help="the salinity level for the system (unit: M)."
    )
    parser.add_argument(
        "-m", "--map_dimensions", type=str, default='100,100,100',
        help="the dimension of the EM map to be computed (a list of float numbers in Å for x, y, z dimension.)."
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
        "-o", "--out_dir", default=None,
        help="the directory to hold the output; if None, unique output directory will be generated automatically."
    )
    parser.add_argument(
        "-r", "--run_psfgen", default=False,
        help="to run psfgen in VMD to generate MD-ready structure."
    )
    parser.add_argument(
        "-p", "--run_apbs", default=False,
        help="to run APBS to get EM maps."
    )
    args = parser.parse_args()
    args.in_psfgen_script_paths = [elm.strip() for elm in args.in_psfgen_script_paths.split(",")]
    args.map_dimensions = [float(elm.strip()) for elm in args.map_dimensions.split(',')]
    args.run_psfgen = bool(args.run_psfgen)
    args.run_apbs = bool(args.run_apbs)
    return args



def main():

    args = arg_parsing_for_step_4()


    step_text = "step2_mcce_to_namd_Chun.py"
    print(f"Results from each MCCE run in {step_text} will be used to prepare VMD visualizations related to EM potentials.")
#    args.in_psfgen_script_paths = [elm.strip() for elm in args.in_psfgen_script_paths.split(',')]
    print(f"Input psfgen files: {args.in_psfgen_script_paths}")
#    args.map_dimensions = [float(elm.strip()) for elm in args.map_dimensions.split(',')]
    print(f"Map dimensions: {args.map_dimensions}")

    for ind, in_psfgen in enumerate(args.in_psfgen_script_paths):
        vEM_obj = vapbs(
            in_psfgen_path=in_psfgen,
            vmd_path=args.vmd_path,
            apbs_path=args.apbs_path,
            map_dimensions=args.map_dimensions,
            out_dir=args.out_dir,
            salinity=args.salt_level,
        )
        vEM_obj.get_APBS_EM_map_and_Vis()




if __name__ == "__main__":
    main()

