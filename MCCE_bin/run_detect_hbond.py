#!/usr/bin/env python
"""
Created on Mar 12 06:00:00 2025

@author: Gehan Ranepura
"""

import os
import subprocess
import argparse
import shutil

# Set up argument parser
parser = argparse.ArgumentParser(description="Run detect_hbond.py for every PDB file in the specified directory.")
parser.add_argument("-pdb_dir",    type=str, default="pdb_output_mc",        help="Directory containing PDB files (default: pdb_output_mc).")
parser.add_argument("-output_dir", type=str, default="pdb_output_mc_hbonds", help="Directory to save the output files (default: output_files).")
args = parser.parse_args()

# Get the directory containing PDB files
pdb_directory    = args.pdb_dir
output_directory = args.output_dir

# Ensure the provided and output directory exists
if not os.path.isdir(pdb_directory):
    print(f"Error: The directory '{pdb_directory}' does not exist.")
    exit(1)

# Ensure the output directory exists, if not create it
if os.path.exists(output_directory):
    print(f"Output directory '{output_directory}' already exists. Deleting it and recreating.")
    shutil.rmtree(output_directory)  # Use shutil.rmtree() if there are files inside
    os.makedirs(output_directory)
else:
    os.makedirs(output_directory)
    print(f"Created output directory: {output_directory}")

# Print the absolute output directory path
abs_pdb_directory = f"{os.path.abspath(pdb_directory)}"
abs_output_directory = f"{os.path.abspath(output_directory)}"
print(f"Output files will be saved to: {abs_output_directory}")
print()

# Script to run
script_name = "detect_hbond.py"

# List all PDB files in the directory
pdb_files = [f for f in os.listdir(abs_pdb_directory) if f.endswith(".pdb")]

if not pdb_files:
    print(f"No PDB files found in the directory '{pdb_directory}'.")
    exit(1)

# Run detect_hbond.py for each PDB file
for pdb_file in pdb_files:
    pdb_path = os.path.abspath(os.path.join(abs_pdb_directory, pdb_file))
    print(f"Processing {pdb_file}...")
    
    # Change the current working directory to the output directory
    os.chdir(abs_output_directory)
    print(f"Output directory: {abs_output_directory}")
 
    # Run detect_hbond.py with the PDB file as an argument
    subprocess.run([script_name, "-inpdb", pdb_path], check=True)
    os.chdir(abs_pdb_directory)

print("Processing complete!")

