#!/usr/bin/env python3

"""
CLI entry point for the MCCE4 Topology Agent.

"""
import argparse
from datetime import datetime
import logging
import os
from pathlib import Path
import sys
from typing import Union

from mcce4 import CLI_EPILOG
from mcce4.pymol_ligand_cif2pdb import main as cif_ligand_converter
from mcce4.mcce_ftpl_agent import __version__
from mcce4.mcce_ftpl_agent.agent import run_agent
from mcce4.mcce_ftpl_agent.config import (DEFAULT_CHARGE_METHOD,
                                          DEFAULT_DIELECTRICS,
                                          DEFAULT_LLM_PROVIDER,
                                          DIMORPHITE_LABEL_STATES,
                                          DIMORPHITE_MAX_VARIANTS,
                                          DIMORPHITE_PH_MIN,
                                          DIMORPHITE_PH_MAX,
                                          DIMORPHITE_PRECISION,
                                          SUPPORTED_CHARGE_METHODS,
                                          SUPPORTED_LLM_PROVIDERS,
                                          USER_CONFIG_DICT,
                                          )
from mcce4.mcce_ftpl_agent.launch_gui import launch_gui
from mcce4.mcce_ftpl_agent.tools.mcce_tools import extract_lig_id_from_pdb


def setup_logging(log_file: str, verbose: bool = False):
    """Configure logging.
    """
    log_level = logging.DEBUG if verbose else logging.INFO
    file_fmt = logging.Formatter("%(asctime)s | %(levelname)-8s | %(message)s",
                                 datefmt="%Y-%m-%d %H:%M:%S")
    console_fmt = logging.Formatter("%(levelname)-8s | %(message)s")
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler(log_file, mode="a")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(file_fmt)
    logger.addHandler(fh)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(log_level)
    ch.setFormatter(console_fmt)
    logger.addHandler(ch)


USAGE ="""
EXAMPLES:
    %(prog)s -lig-code XYZ -lig-smiles C(C(C(=O)O)N        # Ligand code with smiles
    %(prog)s -lig-code EMH                                 # Ligand code only (query RCSB for smiles)
    %(prog)s EMH.pdb                                       # Full auto (PDB input)
    %(prog)s EMH.cif                                       # CIF input (auto-converts)
    %(prog)s EMH.pdb -lig-code EMH                         # PDB + ligand code
    %(prog)s -lig-code EMH --dry-run --no-llm              # Quick test
    %(prog)s EMH.pdb --gui                                 # Web GUI (requires Streamlit)
    %(prog)s EMH.pdb -state-pdbs EMH_01.pdb EMH_+1.pdb     # User state PDBs
    %(prog)s EMH.pdb --no-llm                              # No LLM reasoning
    %(prog)s EMH.pdb -charge-method am1bcc                 # Override charges
    %(prog)s EMH.pdb --dry-run                             # Skip calibration
    %(prog)s EMH.pdb -llm-provider claude -api-key sk-...  # Use Claude
    %(prog)s EMH.pdb -llm-provider chatgpt -api-key sk-... # Use ChatGPT
"""

EPILOG = USAGE + """
NOTE: Before using the GUI (--gui option), install its requirements:
      (mc4)> conda install -c conda-forge --file gui_requirements.txt
""" + CLI_EPILOG + "\n"


def cli_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mcce_ftpl",
        description="🤖 MCCE4 Topology File AI Agent",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=EPILOG
    )
    parser.add_argument("input_file",
                        nargs="?",
                        default=None,
                        help="""Ligand PDB or CIF file (e.g., EMH.pdb or EMH.cif).
Optional if lig_code is provided (default: %(default)s)"""
                        )
    parser.add_argument("-lig-code",
                        type=str,
                        default="",
                        help="""3-letter RCSB ligand code (e.g., EMH). If provided, SMILES and
metadata are fetched from RCSB. A PDB/CIF file is optional;
If omitted, the 3D structure is built from SMILES using RDKit (default: %(default)s)"""
                        )
    parser.add_argument("-lig-smiles",
                        type=str,
                        default="",
                        help="""The ligand SMILES string. If provided along with lig-code,
the 3D structure is built from SMILES using RDKit (default: %(default)s)"""
                        )
    parser.add_argument("-state-pdbs",
                        nargs="+",
                        default=None,
                        help="PDB files for specific states (e.g., EMH_01.pdb EMH_+1.pdb) (default: %(default)s)"
                        )
    parser.add_argument("--gui",
                        default=False,
                        action="store_true",
                        help="Launch the streamlit GUI for interactive review (default: %(default)s)"
                        )
    parser.add_argument("-ph",
                        type=float,
                        default=7.0,
                        help="Target pH (default: %(default)s)"
                        )
    parser.add_argument("-charge-method",
                        default=DEFAULT_CHARGE_METHOD,
                        choices=SUPPORTED_CHARGE_METHODS,
                        help="Charge method (default: %(default)s)"
                        )
    parser.add_argument("-d", "--dielectric",
                        nargs="+",
                        type=int,
                        default=DEFAULT_DIELECTRICS,
                        help="Dielectric constants (default: %(default)s)"
                        )
    parser.add_argument("--no-llm",
                        default=False,
                        action="store_true",
                        help="Disable LLM reasoning (default: %(default)s)"
                        )
    parser.add_argument("--dry-run",
                        default=False,
                        action="store_true",
                        help="Skip RXN calibration steps (default: %(default)s)"
                        )
    # Dimorphite-DL options
    dimorphite_group = parser.add_argument_group(
        "Dimorphite-DL options",
        "Control protonation state enumeration (Phase 2)"
    )
    dimorphite_group.add_argument("-ph-min",
                                  type=float,
                                  default=DIMORPHITE_PH_MIN,
                                  help="Minimum pH for protonation enumeration (default: %(default)s)"
                                  )
    dimorphite_group.add_argument("-ph-max",
                                  type=float,
                                  default=DIMORPHITE_PH_MAX,
                                  help="Maximum pH for protonation enumeration (default: %(default)s)"
                                  )
    dimorphite_group.add_argument("-precision",
                                  type=float,
                                  default=DIMORPHITE_PRECISION,
                                  help="""pKa precision factor: number of std deviations 
from mean pKa to consider (default: %(default)s)"""
                                  )
    dimorphite_group.add_argument("-max-variants",
                                  type=int,
                                  default=DIMORPHITE_MAX_VARIANTS,
                                  help=f"Max protonation variants per compound "
                                       f"(default: {DIMORPHITE_MAX_VARIANTS})"
                                  )
    dimorphite_group.add_argument("--label-states",
                                  default=DIMORPHITE_LABEL_STATES,
                                  action="store_true",
                                  help="Label output SMILES as (default: %(default)s)"
                                  )
    parser.add_argument("-work-dir",
                        default=".",
                        help="Working directory (default: %(default)s)"
                        )
    parser.add_argument("-o", "--output",
                        default=None,
                        help="Output .ftpl filename (default: %(default)s)"
                        )
    parser.add_argument("-v", "--verbose",
                        default=False,
                        action="store_true",
                        help="Use the DEBUG mode for logging (default: %(default)s)"
                        )
    # LLM provider options
    parser.add_argument("-llm-provider",
                        default=DEFAULT_LLM_PROVIDER,
                        choices=SUPPORTED_LLM_PROVIDERS,
                        help="LLM provider (default: %(default)s)"
                        )
    parser.add_argument("-api-key",
                        default=None,
                        help="API key for the chosen LLM provider (overrides env vars) (default: %(default)s)"
                        )

    return parser


def validate_args(args: argparse.Namespace):
    if args.input_file is None and args.lig_code is None:
        raise argparse.ArgumentTypeError("Either an input PDB/CIF file or -lig-code must be provided.")

    if args.input_file is not None:
        if not Path(args.input_file).exists():
            print(f"ERROR: Input file not found: {args.input_file}")
            sys.exit(1)
        if not args.input_file.lower().endswith((".cif", ".pdb")):
            print("ERROR: Input file must be either a .pdb or .cif file")
            sys.exit(1)


def _convert_cif_to_pdb(cif_path: str) -> str:
    """Convert a .cif file to .pdb using cif2pdb_PyMOL.

    Returns the path to the generated PDB file.
    """
    pdb_path = os.path.splitext(cif_path)[0] + ".pdb"
    print(f"🔄 Converting {cif_path} → {pdb_path} using cif2pdb_PyMOL ...")

    cif_ligand_converter({"input": cif_path,
                          "output": pdb_path,
                          "ligand_cache": None,
                          "offline": False})
    if not Path(pdb_path).exists():
        print(f"ERROR: Conversion completed but {pdb_path} was not created.")
        sys.exit(1)
    print(f"✅ Converted successfully: {pdb_path}")

    return pdb_path


def main(args: Union[argparse.Namespace, dict]):

    if isinstance(args, dict):
        args = argparse.Namespace(**args)

    validate_args(args)

    if not args.lig_code:
        lig_id = ""
    else:
        lig_code = args.lig_code.upper()
        lig_id = lig_code

    pdb_path = None

    if not args.lig_smiles:
        if args.input_file is not None:
            # ── CIF → PDB conversion ──
            pdb_path = args.input_file
            if args.input_file.lower().endswith(".cif"):
                pdb_path = _convert_cif_to_pdb(args.input_file)

            # Determine lig_id: -lig-code takes precedence, else extract from PDB
            if not lig_id:
                lig_id = extract_lig_id_from_pdb(pdb_path)

    # ── GUI mode: launch Streamlit ──
    if args.gui:
        if pdb_path is None:
            print("ERROR: GUI mode requires an input PDB/CIF file.")
            sys.exit(1)
        args.pdb = pdb_path  # GUI expects args.pdb
        launch_gui(args)
        return

    # ── CLI mode: run agent directly ──
    setup_logging(f"mcce_ftpl_{lig_id}.log", args.verbose)

    logging.info(f"{'='*60}")
    logging.info(f"  🤖 MCCE4 Topology Agent v{__version__} (LangGraph + per-state PDBs)")
    logging.info(f"{'='*60}")
    if pdb_path:
        logging.info(f"  Input:  {os.path.abspath(pdb_path)}")
        if args.input_file and args.input_file.lower().endswith(".cif"):
            logging.info(f"  (converted from {Path(args.input_file)!s})")
    else:
        logging.info(f"  Input:  -lig-code {lig_id} (no PDB file: will build from RCSB SMILES)")
        if args.lig_smiles:
            logging.info(f"  Input:  -lig-smiles (no PDB file: will build from given SMILES)")

    logging.info(f"  Ligand: {lig_id}   pH: {args.ph}   Method: {args.charge_method}")
    logging.info(f"  Dimorphite-DL: ph_min={args.ph_min}, ph_max={args.ph_max}, "
                 f"precision={args.precision}, max_variants={args.max_variants}"
                 f"{', label_states=True' if args.label_states else ''}")
    logging.info(f"  LLM:    {args.llm_provider}")
    logging.info(f"  Time:   {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logging.info(f"{'='*60}\n")

    # Set env to disable LLM if requested
    if args.no_llm:
        USER_CONFIG_DICT = None

    final_state = run_agent(
        pdb_path=pdb_path,
        use_gui=False,
        charge_method=args.charge_method,
        dielectrics=args.dielectric,
        ph=args.ph,
        work_dir=args.work_dir,
        dry_run=args.dry_run,
        user_state_pdbs=args.state_pdbs,
        output=args.output,
        llm_provider=args.llm_provider,
        api_key=args.api_key,
        lig_code=lig_code,
        lig_smiles=args.lig_smiles,
        ph_min=args.ph_min,
        ph_max=args.ph_max,
        precision=args.precision,
        max_variants=args.max_variants,
        label_states=args.label_states,
    )

    # Exit code based on errors
    if final_state.get("errors"):
        sys.exit(1)


def cli():
    p = cli_parser()
    args = p.parse_args()
    main(args)


if __name__ == "__main__":
    cli()
