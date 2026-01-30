#!/bin/bash

# Parameter/Options for SLURM (Simple Linux Utility for Resource Management)
#SBATCH --job-name=mcce4_run
#SBATCH --output=submit_mcce4.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=12G                 # Adjust memory if needed
#SBATCH --export=ALL

#=============================================================================
#-----------------------------------------------------------------------------
# Input and Output:
input_pdb="prot.pdb"    # (INPDB) keep if you have soft-linked your pdb as prot.pdb, e.g.: ln -s 4lzt.pdb prot.pdb

# >>> Automated Parameters (best not change)
MCBIN=$(dirname $(which mcce))
MCCE_HOME=$(dirname "$MCBIN")

# the python executable will come from your virtual env if activated:
PYEX=$(python3 -c "import sys; print(sys.executable)")
PYENV="${PYEX%python3}"

APTNR=$(which apptainer) > /dev/null 2>&1
if [[ -n "$APTNR" ]]; then
    APPTAINER_BIN=$(dirname "$APTNR")
else
    echo "Error: apptainer program not found."
    exit 1
fi
# <<<

# Set MCCE4 Parameters: default USER_PARAM: "/param"->MCCE_HOME/param; default EXTRA is 'extra.tpl'->MCCE_HOME/extra.tpl
USER_PARAM="/param"     # could be "/path/to/new/topologies"
EXTRA="extra.tpl"       # could be "/path/to/test/extra.tpl" or "./extra.tpl" (local file), etc.
TMP="/tmp"
CPUS=1

# Step control flags
step1="t"               # STEP1: pre-run, pdb-> mcce pdb  (DO_PREMCCE)
step2="t"               # STEP2: make rotamers            (DO_ROTAMERS)
step3="t"               # STEP3: Energy calculations      (DO_ENERGY)
step4="t"               # STEP4: Monte Carlo Sampling     (DO_MONTE)
step_clean="t"          # Clean PBE data                  (BACKUP CLEAN) : Set to f if step3 --debug option is used

# Optional step controls
center="t"              # Center protein structure before MCCE run    : Set to f to skip centering and use input PDB as-is
# Optional, user-defined scripts. User MUST satisfy condidtions of their custom script
stepM="f"               # Generate Partial Membranes                  : Sample found in $MCCE_HOME/schedulers/stepM.sh
stepA="f"               # Run a custom script between step1 and step2
stepB="f"               # Run a custom script between step2 and step3
stepC="f"               # Run a custom script between step3 and step4

# MCCE Simulation
EPS=8                   # Protein dielectric constant

STEP1="$PYEX $MCBIN/step1.py -d $EPS --dry"
STEP2="$PYEX $MCBIN/step2.py -d $EPS -l 1"
STEP3="$PYEX $MCBIN/step3.py -d $EPS -s ngpb -p $CPUS -t $TMP"
STEP4="$PYEX $MCBIN/step4.py --xts -i 7 -n 1 --ms"

# Optional scripts locations
STEPM="/path/to/stepM.sh"         # Optional StepM: Bash script; Example in $MCCE_HOME/schedulers/stepM.sh
STEPA="/path/to/stepA_script.py"  # Optional StepA: Python script to run between step1 and step2.
STEPB="/path/to/stepB_script.py"  # Optional StepB: Python script to run between step2 and step3.
STEPC="/path/to/stepC_script.py"  # Optional StepC: Python script to run between step3 and step4.

# NO USER INPUT NECCESARY BELOW THIS LINE
#==============================================================================

# Check MCCE_HOME exists before PATH export
if [[ ! -d "$MCBIN" ]]; then
    echo -e "\033[0;31m[ERROR]\033[0m MCCE4 executable was not found or $MCCE_HOME/bin does not exist."
    echo "Please check your MCCE_HOME path in submit_mcce4.sh."
    exit 1
fi

if [[ ! -d "$MCCE_HOME" ]]; then
    echo "Error: MCCE_HOME is not set."
    exit 1
fi

# Initialize Apptainer to ensure job uses user-installed Apptainer and avoid systemd/cgroups (DBus)
# issues on compute nodes
export PATH="$APPTAINER_BIN:$PATH"
export APPTAINER_CONFIG_FILE="$HOME/.apptainer/apptainer.conf"
mkdir -p "$HOME/.apptainer"
cat > "$APPTAINER_CONFIG_FILE" <<'EOF'
systemd cgroups = no
EOF

# Remove any existing instance of mc_bin from PATH and prepend mc_bin to PATH
PATH=$(echo "$PATH" | tr ':' '\n' | grep -vx "$MCCE_HOME/MCCE_bin" | paste -sd ':' -)
export PATH="$MCCE_HOME/MCCE_bin:$PATH"

# Remove any existing instance of mc_bin from PATH and prepend mc_bin to PATH
PATH=$(echo "$PATH" | tr ':' '\n' | grep -vx "$MCBIN" | paste -sd ':' -)
export PATH="$MCBIN:$PATH"

# Remove any existing instance of PYENV from PATH and prepend PYENV to PATH
PATH=$(echo "$PATH" | tr ':' '\n' | grep -vx "$PYENV" | paste -sd ':' -)
export PATH="$PYENV:$PATH"

echo "============================================================"
echo "MCCE4 SUBMIT SHELL JOB ENVIRONMENT (startup diagnostics)"
echo "------------------------------------------------------------"
echo "Date:             $(date)"
echo "Host:             $(hostname)"
echo "User:             $(whoami)"
echo
echo -e "Apptainer:        $(which apptainer)"
echo -e "Config File:      $APPTAINER_CONFIG_FILE"
echo -e "MCCE_HOME:        $MCCE_HOME"
echo -e "MCBIN:            $MCBIN"
echo -e "Driver:           $MCBIN/driver_mcce4.sh"
echo -e "PYEX:             $PYEX"
echo -e "PATH:             $PATH"
echo "============================================================"
echo

# Export environment for downstream script
export PYEX
export input_pdb MCCE_HOME MCBIN USER_PARAM EXTRA TMP
export step1 step2 step3 step4 step_clean
export center stepM stepA stepB stepC
export STEP1 STEP2 STEP3 STEP4
export STEPM STEPA STEPB STEPC

# Inititiate MCCE_HOME PATH and call driver_mcce4.sh
"$MCBIN"/driver_mcce4.sh

# ==============================================================================
# Script Name   : submit_mcce4.sh
# Purpose       : Automate and control the full MCCE4 simulation pipeline including optional custom
#               : preprocessing steps.
#
# Description   :
#   This script manages the sequential execution of MCCE4 simulation steps (1 to 4), with optional
#   hooks (stepM, stepA, stepB, stepC) that allow the user to insert custom membrane generation and
#   intermediate processing scripts.
#   It records the timing and success/failure of each step in a detailed log file (`mcce_timing.log`).
#   The script supports flexible control through flags to enable/disable specific MCCE steps or custom steps.
#
# Main Features :
#   - Step 1: Convert standard PDB to MCCE-compatible format
#   - Step 2: Generate rotamers
#   - Step 3: Perform energy calculations
#   - Step 4: Run Monte Carlo sampling
#   - Step clean: Clean temporary pbe_data
#   - Optional StepM: Add membrane-specific conformers (e.g., using IPECE)
#   - Optional StepA/B/C: Insert custom preprocessing scripts between core steps
#   - Intelligent skip logic and output checking to prevent redundant work
#   - Runtime logging with timestamps for each phase
#
# Input Requirements :
#   - input_pdb       : A protein PDB file named `prot.pdb`
#   - MCCE_HOME       : Path to MCCE4 installation directory
#   - USER_PARAM      : Directory with user-defined MCCE parameters (optional)
#   - EXTRA           : Custom extra.tpl file (optional)
#   - TMP             : Path to store temporary pbe_data (default: /tmp)
#   - CPUS            : Multiprocessing for step3
#   - Optional scripts for stepM/A/B/C must exist and be executable if enabled
#
# Output Files :
#   - step1_out.pdb, step2_out.pdb, head3.lst, pK.out
#   - Timing report: mcce_timing.log
#   - Logs for each step: step1.log, step2.log, etc.
#
# Usage :
#   Set control flags (`step1`, `step2`, etc.) to "t" or "f" to enable/disable each step.
#   Set paths to optional scripts as needed.
#   Submit this script to a SLURM cluster or run locally if sbatch is not used.
#
# Author        : Gehan A. Ranepura
# Date Created  : July 2025
# ==============================================================================
