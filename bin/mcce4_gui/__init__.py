"""
mcce4_gui — Web-based GUI for the MCCE4 protein electrostatics pipeline.

A subpackage of MCCE4-Alpha that provides a modern browser interface for
running MCCE4 steps, visualizing protein structures, analyzing titration
results, and submitting SLURM jobs.

Usage:
    python -m mcce4_gui [--port 8080] [--host 0.0.0.0] [--no-browser]

Then open http://localhost:8080 in your browser.
For remote servers, SSH tunnel: ssh -L 8080:localhost:8080 user@server
"""

__version__ = "1.0.0"
__author__ = "Gehan Ranepura"
__license__ = "MIT"
