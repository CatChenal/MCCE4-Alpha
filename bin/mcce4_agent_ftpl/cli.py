"""
CLI entry point for the MCCE4 Topology Agent.
"""

import argparse
import os
import sys
import textwrap
import logging
import subprocess
from datetime import datetime

from .config import SUPPORTED_CHARGE_METHODS, DEFAULT_CHARGE_METHOD, DEFAULT_DIELECTRICS


def setup_logging(log_file: str, verbose: bool = False):
    """Configure logging."""
    log_level = logging.DEBUG if verbose else logging.INFO
    file_fmt = logging.Formatter("%(asctime)s | %(levelname)-8s | %(message)s",
                                  datefmt="%Y-%m-%d %H:%M:%S")
    console_fmt = logging.Formatter("%(levelname)-8s | %(message)s")
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler(log_file, mode="w")
    fh.setLevel(logging.DEBUG); fh.setFormatter(file_fmt)
    logger.addHandler(fh)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(log_level); ch.setFormatter(console_fmt)
    logger.addHandler(ch)


def main():
    parser = argparse.ArgumentParser(
        description="🤖 MCCE4 Topology File AI Agent",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
        Examples:
          %(prog)s EMH.pdb                                         # Full auto
          %(prog)s EMH.pdb --gui                                   # Web GUI (Streamlit)
          %(prog)s EMH.pdb --state-pdbs EMH_01.pdb EMH_+1.pdb     # User state PDBs
          %(prog)s EMH.pdb --no-llm                                # No LLM reasoning
          %(prog)s EMH.pdb --charge-method am1bcc                  # Override charges
          %(prog)s EMH.pdb --dry-run                               # Skip calibration

        Install:
          pip install google-genai langgraph dimorphite-dl rdkit streamlit
          export GEMINI_API_KEY="your_free_key"  # from https://ai.google.dev
        """))

    parser.add_argument("pdb", help="Ligand PDB file (e.g., EMH.pdb)")
    parser.add_argument("--state-pdbs", nargs="+", default=None,
                        help="PDB files for specific states (e.g., EMH_01.pdb EMH_+1.pdb)")
    parser.add_argument("--gui", action="store_true",
                        help="Launch web GUI (Streamlit) for interactive review")
    parser.add_argument("--ph", type=float, default=7.4, help="Target pH (default: 7.4)")
    parser.add_argument("--charge-method", default=DEFAULT_CHARGE_METHOD,
                        choices=SUPPORTED_CHARGE_METHODS, help="Charge method")
    parser.add_argument("-d", "--dielectric", nargs="+", type=int, default=DEFAULT_DIELECTRICS,
                        help="Dielectric constants (default: 2 4 8)")
    parser.add_argument("--no-llm", action="store_true", help="Disable LLM reasoning")
    parser.add_argument("--dry-run", action="store_true", help="Skip RXN calibration")
    parser.add_argument("--work-dir", default=".", help="Working directory")
    parser.add_argument("-o", "--output", default=None, help="Output .ftpl filename")
    parser.add_argument("-v", "--verbose", action="store_true")

    args = parser.parse_args()

    # Validate input
    if not os.path.exists(args.pdb):
        print(f"ERROR: PDB file not found: {args.pdb}")
        sys.exit(1)

    # ── GUI mode: launch Streamlit ──
    if args.gui:
        launch_gui(args)
        return

    # ── CLI mode: run agent directly ──
    from .tools.mcce_tools import extract_lig_id_from_pdb
    lig_id = extract_lig_id_from_pdb(args.pdb)

    setup_logging(f"mcce4_agent_ftpl_{lig_id}.log", args.verbose)

    logging.info(f"{'='*60}")
    logging.info(f"  🤖 MCCE4 Topology Agent v3.0 (LangGraph + per-state PDBs)")
    logging.info(f"{'='*60}")
    logging.info(f"  Input:  {os.path.abspath(args.pdb)}")
    logging.info(f"  Ligand: {lig_id}   pH: {args.ph}   Method: {args.charge_method}")
    logging.info(f"  Time:   {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logging.info(f"{'='*60}\n")

    # Set env to disable LLM if requested
    if args.no_llm:
        os.environ.pop("GEMINI_API_KEY", None)
        os.environ.pop("GOOGLE_API_KEY", None)

    from .agent import run_agent

    final_state = run_agent(
        pdb_path=args.pdb,
        use_gui=False,
        charge_method=args.charge_method,
        dielectrics=args.dielectric,
        ph=args.ph,
        work_dir=args.work_dir,
        dry_run=args.dry_run,
        user_state_pdbs=args.state_pdbs,
        output=args.output,
    )

    # Exit code based on errors
    if final_state.get("errors"):
        sys.exit(1)


def launch_gui(args):
    """Launch the Streamlit web GUI."""
    gui_path = os.path.join(os.path.dirname(__file__), "gui", "app.py")

    if not os.path.exists(gui_path):
        print(f"ERROR: GUI app not found at {gui_path}")
        sys.exit(1)

    # Pass PDB path via environment so Streamlit can access it
    os.environ["MCCE_AGENT_PDB"] = os.path.abspath(args.pdb)
    os.environ["MCCE_AGENT_PH"] = str(args.ph)
    os.environ["MCCE_AGENT_CHARGE_METHOD"] = args.charge_method

    print(f"🌐 Launching MCCE4 Topology Agent GUI...")
    print(f"   PDB: {args.pdb}")
    print(f"   Open your browser to: http://localhost:8501")
    print(f"   (If remote: ssh -L 8501:localhost:8501 user@server)")
    print()

    try:
        subprocess.run([
            sys.executable, "-m", "streamlit", "run", gui_path,
            "--server.headless", "true",
            "--server.port", "8501",
            "--browser.gatherUsageStats", "false",
        ])
    except KeyboardInterrupt:
        print("\n  GUI stopped.")
    except FileNotFoundError:
        print("ERROR: Streamlit not installed — run: pip install streamlit")
        sys.exit(1)


if __name__ == "__main__":
    main()
