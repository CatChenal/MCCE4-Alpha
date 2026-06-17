#!/usr/bin/env python3

# definitely beta
__version__ = "0.1.6"
__authors__ = "Gehan / MCCE4 Team (GunnerLab)"

from pathlib import Path


GUI_REQS = "gui_requirements.txt"


def gui_reqs_cmd() -> str:
    gui_reqs_fp = Path(__file__).parent.joinpath(GUI_REQS)
    cmd_to_run = f"""
To install the GUI-specific packages run these commands:
 conda activate mc4
 conda install -c conda-forge -f {gui_reqs_fp!s}
"""
    return cmd_to_run
