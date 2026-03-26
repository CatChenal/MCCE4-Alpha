#!/usr/bin/env python3
"""
mcce4_agent_ftpl.py — MCCE4 Topology File AI Agent
====================================================

Launcher script. Place this in MCCE4-Alpha/bin/ alongside the
mcce4_agent_ftpl/ package directory.

Usage:
  mcce4_agent_ftpl.py EMH.pdb                                     # Full auto
  mcce4_agent_ftpl.py EMH.pdb --gui                               # Web GUI
  mcce4_agent_ftpl.py EMH.pdb --state-pdbs EMH_01.pdb EMH_+1.pdb  # User states
  mcce4_agent_ftpl.py EMH.pdb --no-llm --dry-run                  # Minimal run

Install dependencies:
  pip install google-genai langgraph dimorphite-dl rdkit streamlit
  export GEMINI_API_KEY="your_free_key"   # from https://ai.google.dev
"""

import sys
import os

# Ensure the package is findable from bin/
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from mcce4_agent_ftpl.cli import main

if __name__ == "__main__":
    main()
