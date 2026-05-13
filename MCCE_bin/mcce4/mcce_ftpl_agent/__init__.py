#!/usr/bin/env python3

"""
mcce_ftpl_agent — MCCE4 Topology File AI Agent Package
=========================================================

An AI agent for creating MCCE4 topology files (.ftpl) from ligand PDB files.
Uses Dimorphite-DL for protonation state enumeration, RDKit for chemistry
validation, Google Gemini (free tier) for reasoning, and LangGraph for
agentic orchestration.

For GUI:
    conda install -c conda-forge --file gui_requirements.txt

Usage:
    # CLI
    mcce_ftpl EMH.pdb --gui

    # Python
    from mcce_ftpl_agent.agent import run_agent
    run_agent("EMH.pdb")
"""

__version__ = "1.0.5"
__author__ = "Gehan / MCCE4 Team (GunnerLab)"
