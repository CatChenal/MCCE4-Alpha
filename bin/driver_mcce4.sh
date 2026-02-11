#!/bin/bash

# Script Name   : driver_mcce4.sh
# Purpose       : Automate and control the full MCCE4 simulation pipeline
#
# Description   : This driver script orchestrates the complete MCCE4 workflow
# including:
#                 - Optional protein structure centering
#                 - Optional membrane generation (STEP M)
#                 - Core MCCE steps 1-4 (protonation, rotamer, MC sampling, analysis)
#                 - Custom user scripts at key pipeline stages (STEP A, B, C)
#                 - Comprehensive timing and logging
#                 - Automatic cleanup of temporary files
#
# Execution Flow:
#   1. Initialize timing log and verify MCCE environment
#   2. Setup parameters and merge -u flags if needed
#   3. [Optional] Center protein structure
#   4. [Optional] Generate membrane (STEP M)
#   5. Run STEP 1: Generate protonation states
#   6. [Optional] Run custom script after STEP 1 (STEP A)
#   7. Run STEP 2: Generate rotamers and conformers
#   8. [Optional] Append membrane to step2_out.pdb
#   9. [Optional] Run custom script after STEP 2 (STEP B)
#   10. Run STEP 3: Monte Carlo sampling
#   11. [Optional] Run custom script after STEP 3 (STEP C)
#   12. Run STEP 4: Analyze titration results
#   13. Clean up temporary PBE data directories
#
# Requirements  : MCCE_HOME, input_pdb, and step flags must be defined before 
#                 execution
# Output        : mcce_timing.log with detailed timing information for all steps
# ...............................................................................

#=================================================================================
# SECTION 1: INITIALIZATION AND LOGGING SETUP
#
# Initialize timing log file and configure error handling for critical operations
#echo ">> driver_mcce4.sh - Python path from submit_mcce4.sh: $PYEX"

TIMING_FILE="mcce_timing.log"
echo "MCCE Timing Report" > $TIMING_FILE
echo "====================================" >> $TIMING_FILE
set -e  # Exit immediately if any command fails

#=================================================================================
# SECTION 2: VERIFY AND LOG MCCE ENVIRONMENT VARIABLES
#
# Verify critical MCCE paths exist and log them to timing file

# Log MCCE_HOME (required)
echo "MCCE_HOME: $MCCE_HOME" >> $TIMING_FILE

# Verify EXTRA template file exists, use default if not provided
if [ -f "$EXTRA" ]; then
    EXTRA="$EXTRA"
else
    EXTRA="$MCCE_HOME/extra.tpl"  # Default fallback location
fi
echo "EXTRA: $EXTRA" >> $TIMING_FILE

# Check for custom user parameter directory
if [ -d "$USER_PARAM" ]; then
    echo "USER_PARAM: $USER_PARAM" >> $TIMING_FILE
else
    echo "USER_PARAM: N/A" >> $TIMING_FILE
fi

#=================================================================================
# SECTION 3: SETUP TEMPORARY DIRECTORY FOR PBE DATA
#
# Configure temp directory for Poisson-Boltzmann electrostatics calculations
# If TMP is not /tmp, create directory or clean existing contents

if [[ "$TMP" != "/tmp" ]]; then
    echo "TMP is not set to /tmp. Using custom path: $TMP"
    if [[ ! -d "$TMP" ]]; then
        echo "Creating directory: $TMP"
        mkdir -p "$TMP"
        echo "TMP: $TMP" >> $TIMING_FILE
    else
        echo "Cleaning existing directory: $TMP"
        rm -rf "$TMP"/*  # Remove all contents but keep directory
        echo "TMP: $TMP" >> $TIMING_FILE
    fi
fi

#=================================================================================
# SECTION 4: LOG OPTIONAL FEATURES AND START TIMING
#
# Document which optional features are enabled in the timing log

echo -e "====================================\n" >> $TIMING_FILE
echo "Optional Features Enabled:" >> $TIMING_FILE
[[ "$center" == "t" ]] && echo " - Center Protein Structure (center)" >> $TIMING_FILE
[[ "$stepM" == "t" ]]  && echo " - Partial Membrane Generation (stepM)" >> $TIMING_FILE
[[ "$stepA" == "t" ]]  && echo " - Custom Script After Step1 (stepA)" >> $TIMING_FILE
[[ "$stepB" == "t" ]]  && echo " - Custom Script After Step2 (stepB)" >> $TIMING_FILE
[[ "$stepC" == "t" ]]  && echo " - Custom Script After Step3 (stepC)" >> $TIMING_FILE
echo -e "====================================\n" >> $TIMING_FILE

# Record the start time for total runtime calculation
script_start_time=$(date +%s)
echo "Run started at: $(date)" >> $TIMING_FILE

#=================================================================================
# SECTION 5: CONFIGURE STEP FLAGS
#
# Append --norun flag to steps that are disabled (when step flag is "f")
# This allows the script to document skipped steps without executing them

[[ "$step1" == "f" ]] && STEP1="$STEP1 --norun"
[[ "$step2" == "f" ]] && STEP2="$STEP2 --norun"
[[ "$step3" == "f" ]] && STEP3="$STEP3 --norun"
[[ "$step4" == "f" ]] && STEP4="$STEP4 --norun"

#=================================================================================
# SECTION 6: BUILD PARAMETER STRING FOR -u FLAG
#
# Create the PARAM string containing MCCE environment variables that will be passed
# to all MCCE steps via the -u flag. This includes MCCE_HOME, EXTRA, and optionally
# USER_PARAM

if [ -d "$USER_PARAM" ]; then
    PARAM="-u MCCE_HOME=$MCCE_HOME,EXTRA=$EXTRA,USER_PARAM=$USER_PARAM"
else
    PARAM="-u MCCE_HOME=$MCCE_HOME,EXTRA=$EXTRA"
fi

#=================================================================================
# SECTION 7: DEFINE UTILITY FUNCTIONS
#
#---------------------------------------------------------------------------------
# FUNCTION: merge_u_params
# Purpose : Intelligently merge -u parameters when STEP commands already contain them
# Usage   : merge_u_params "$STEP_CMD" "$PARAM"
# Returns : Modified command with merged -u parameters
#
# Example : Input:  STEP="step1.py -u PARAM1=val1" PARAM="-u PARAM2=val2"
#           Output: "step1.py -u PARAM1=val1,PARAM2=val2"
#
# Notes   : - Handles both quoted (-u "value") and unquoted (-u value) formats
#           - If STEP has no -u flag, simply appends PARAM
#           - Merges duplicate -u flags into a single comma-separated list
#---------------------------------------------------------------------------------
merge_u_params() {
    local step_cmd="$1"
    local param="$2"

    # Check if step_cmd already contains a -u flag
    if [[ "$step_cmd" =~ (^|[[:space:]])-u[[:space:]] ]]; then
        # Split command at the -u flag to isolate components
        local before_u="${step_cmd%% -u *}"
        local after_u="${step_cmd#* -u }"

        # Extract the -u value from STEP (handles both quoted and unquoted)
        local step_u_value=""
        local remaining=""

        # Try to match quoted value first: -u "VALUE"
        if [[ "$after_u" =~ ^\"([^\"]*)\"(.*)$ ]]; then
            step_u_value="${BASH_REMATCH[1]}"
            remaining="${BASH_REMATCH[2]}"
        # Try unquoted value: -u VALUE
        elif [[ "$after_u" =~ ^([^[:space:]]+)(.*)$ ]]; then
            step_u_value="${BASH_REMATCH[1]}"
            remaining="${BASH_REMATCH[2]}"
        fi

        # Extract the parameter's -u value (remove "-u " prefix from PARAM)
        local param_u_value="${param#-u }"

        # Merge existing and new parameters with comma separator
        local merged_u="${step_u_value},${param_u_value}"

        # Return reconstructed command with merged -u parameters
        echo "${before_u} -u ${merged_u}${remaining}"
    else
        # No -u in step_cmd, simply append PARAM as-is
        echo "$step_cmd $param"
    fi
}

#---------------------------------------------------------------------------------
# FUNCTION: file_just_made
# Purpose : Check if a file was created within the last 5 minutes
# Usage   : file_just_made "$step_flag" "$filename"
# Returns : 0 (true) if file was recently created or step was disabled
#           1 (false) otherwise
#
# Notes   : - Used to verify that MCCE steps successfully produced output files
#           - If step_flag is "f" (disabled), always returns true (no check needed)
#           - Uses 'find -mmin -5' to check for files modified in last 5 minutes
#---------------------------------------------------------------------------------
function file_just_made {
    step_flag="$1"   # "t" (enabled) or "f" (disabled)
    file="$2"        # File path to check

    # If step was disabled, treat file as OK to proceed
    if [[ "$step_flag" == "f" ]]; then
        return 0
    fi

    # Check if file exists and was modified in the last 5 minutes
    if [[ -f "$file" ]] && [[ $(find "$file" -mmin -5) ]]; then
        return 0
    else
        return 1
    fi
}

#---------------------------------------------------------------------------------
# FUNCTION: format_time
# Purpose : Convert elapsed seconds into human-readable HH:MM:SS format
# Usage   : format_time $elapsed_seconds
# Returns : Formatted string like "01h:23m:45s"
#
# Example : format_time 3665  →  "01h:01m:05s"
#---------------------------------------------------------------------------------
format_time() {
    local elapsed=$1
    local hours=$((elapsed / 3600))
    local minutes=$(( (elapsed % 3600) / 60 ))
    local seconds=$((elapsed % 60))
    printf "%02dh:%02dm:%02ds" "$hours" "$minutes" "$seconds"
}

#---------------------------------------------------------------------------------
# FUNCTION: time_step
# Purpose : Execute and time an MCCE step, logging results to timing file
# Usage   : time_step "STEP_NAME" "$STEP_CMD" "output_file" "success_msg" "$step_flag"
#
# Parameters:
#   $1 - step_name       : Display name (e.g., "STEP1", "STEP2")
#   $2 - step_cmd        : Full command to execute
#   $3 - success_output  : Expected output file/directory to verify success
#   $4 - success_msg     : Message to log on successful completion
#   $5 - step_flag       : "t" (enabled) or "f" (disabled with --norun)
#
# Behavior:
#   - Times the step execution
#   - Logs "Running..." message to timing file
#   - Executes the command using eval
#   - Checks if expected output was created
#   - Logs success/warning/failure to timing file and console
#   - Formats elapsed time in HH:MM:SS
#
# Exit Codes:
#   - Success: Expected output exists and is recent (or --norun was used)
#   - Warning: Command succeeded but output not fresh
#   - Failure: Command failed or output not created
#---------------------------------------------------------------------------------
function time_step {
    step_name="$1"
    step_cmd="$2"
    success_output="$3"
    success_msg="$4"
    step_flag="$5"    # "t" or "f"

    echo "Running $step_name ..."
    printf "%-6s: Running $step_name ...\n" "$step_name" >> "$TIMING_FILE"
    start_time=$(date +%s)

    # Execute the step command (works even if --norun is present)
    if eval "$step_cmd"; then
        end_time=$(date +%s)
        elapsed=$((end_time - start_time))
        formatted_time=$(format_time "$elapsed")

        # Remove the "running..." placeholder line from timing file
        sed -i "/^${step_name}.*Running ${step_name}.*/d" "$TIMING_FILE"

        # Verify expected output was created
        if [[ -f "$success_output" ]] || [[ -d "$success_output" ]]; then
            if [[ "$step_flag" == "f" ]]; then
                # Step was run with --norun, don't check file freshness
                echo "$step_name completed in $formatted_time (Note: --norun was used)."
                printf "%-6s: %s   - %s\n" "$step_name" "$formatted_time" "(--norun): $success_msg" >> "$TIMING_FILE"
            elif file_just_made "$step_flag" "$success_output"; then
                # Output file was recently created - success
                echo "$step_name completed SUCCESSFULLY in $formatted_time."
                printf "%-6s: %s   - Success: %s\n" "$step_name" "$formatted_time" "$success_msg" >> "$TIMING_FILE"
            else
                # Output exists but wasn't recently created - warning
                echo "$step_name completed with warnings in $formatted_time (Output file not fresh or just created)."
                printf "%-6s: %s   - Warning: %s (not recently created)\n" "$step_name" "$formatted_time" "$success_msg" >> "$TIMING_FILE"
            fi
        else
            # Expected output file was not created - failure
            echo "$step_name FAILED in $formatted_time!"
            printf "%-6s: %s   - FAILED: Expected output not created.\n" "$step_name" "$formatted_time" >> "$TIMING_FILE"
        fi
    else
        # Command returned non-zero exit status - failure
        end_time=$(date +%s)
        elapsed=$((end_time - start_time))
        formatted_time=$(format_time "$elapsed")

        # Remove the "running..." placeholder line from timing file
        sed -i "/^${step_name}.*Running ${step_name}.*/d" "$TIMING_FILE"

        echo "$step_name FAILED in $formatted_time!"
        printf "%-6s: %s   - FAILED: Command returned non-zero exit status.\n" "$step_name" "$formatted_time" >> "$TIMING_FILE"
    fi
}

#=================================================================================
# SECTION 8: OPTIONAL PROTEIN CENTERING
#
# If enabled (center="t"), center the protein structure at the origin using orientation.py
# This ensures consistent coordinate system for simulations, especially important for
# membrane proteins and systems requiring specific geometric orientation

if [[ "$center" == "t" ]]; then
    echo "Centering protein structure..."
    start_time=$(date +%s)

    # Run orientation.py and capture output
    output=$("$PYEX" $MCBIN/orientation.py "$input_pdb")
    # Extract the new centered filename from orientation.py output
    centered_pdb=$(echo "$output" | awk '/Structure moved to origin/ {print $NF}')

    end_time=$(date +%s)
    elapsed=$((end_time - start_time))
    formatted_time=$(format_time "$elapsed")

    # Verify the centered file was successfully created
    if [[ ! -f "$centered_pdb" ]]; then
        echo -e "\033[0;31m[ERROR]\033[0m Centered file '$centered_pdb' was not created. Exiting."
        printf "%-6s: %s   - FAILED: Centered file not created.\n" "CENTER" "$formatted_time" >> "$TIMING_FILE"
        exit 1
    else
        # Update input_pdb to use the centered file for all subsequent steps
        input_pdb="$centered_pdb"
        echo "Structure centered SUCCESSFULLY in $formatted_time. Using file: $input_pdb"
        printf "%-6s: %s   - Success: Structure centered to '%s'.\n" "CENTER" "$formatted_time" "$centered_pdb" >> "$TIMING_FILE"
    fi
else
    # Centering disabled, use original input file
    echo "Skipping protein centering. Using original file: $input_pdb"
fi

echo "Final input PDB file: $input_pdb"

#=================================================================================
# SECTION 9: FINALIZE MCCE STEP COMMANDS
#
# Merge -u parameters from STEP variables with the global PARAM string
# This handles cases where STEP commands already contain -u flags with custom parameters

STEP1_MERGED="$(merge_u_params "$STEP1" "$PARAM")"
STEP2_MERGED="$(merge_u_params "$STEP2" "$PARAM")"
STEP3_MERGED="$(merge_u_params "$STEP3" "$PARAM")"
STEP4_MERGED="$(merge_u_params "$STEP4" "$PARAM")"

# Build final command strings with redirected output
STEP1_CMD="$STEP1_MERGED \$input_pdb > step1.log"
STEP2_CMD="$STEP2_MERGED > step2.log"
STEP3_CMD="$STEP3_MERGED > step3.log 2>&1"  # Redirect both stdout and stderr for step3
STEP4_CMD="$STEP4_MERGED > step4.log"

#=================================================================================
# SECTION 10: CONFIGURE OPTIONAL SCRIPT COMMANDS
#
# Setup commands for optional custom scripts (STEPM, STEPA, STEPB, STEPC)
# STEPM uses sbatch if available (for HPC clusters), otherwise runs with bash

if command -v sbatch &> /dev/null; then
    STEPM_CMD="sbatch --wait $STEPM"  # --wait ensures sbatch waits for job completion
else
    STEPM_CMD="bash $STEPM > stepM.log"
fi
STEPA_CMD="$PYEX $STEPA > stepA.log"
STEPB_CMD="$PYEX $STEPB > stepB.log"
STEPC_CMD="$PYEX $STEPC > stepC.log"

#=================================================================================
# SECTION 11: SETUP CLEANUP HANDLER FOR INTERRUPTS
#
# Define cleanup function to handle script interruptions (Ctrl+C, scancel, etc.)
# Ensures temporary PBE data directories are cleaned up even on premature exit

cleanup_on_exit() {
    echo ""
    echo ">>> Caught termination signal. Running STEP_CLEAN before exiting..."

    WDIR=$(pwd -L)
    WDIR_FLAT=$(echo "$WDIR" | cut -d'/' -f3- | sed 's|/|.|g')

    if [ "$step_clean" = "t" ] && [ -d "$TMP" ]; then
        # Find and remove temporary directories matching this working directory
        MATCHING_DIRS=$(find "$TMP" -mindepth 1 -maxdepth 1 -type d -user "$USER" -name "*${WDIR_FLAT}*")
        if [ -n "$MATCHING_DIRS" ]; then
            echo "$MATCHING_DIRS" | xargs -r rm -rf
            printf "%-6s: done.     - STEP_CLEAN after cancellation.\n" "STEP_CLEAN" >> "$TIMING_FILE"
        else
            printf "%-6s: skipped.  - No matching tmp dirs at cancel.\n" "STEP_CLEAN" >> "$TIMING_FILE"
        fi
    fi
}

# Register the cleanup function to run on interrupt signals
trap cleanup_on_exit SIGINT SIGTERM EXIT

#=================================================================================
#                         MAIN EXECUTION: MCCE4 SIMULATION PIPELINE
#
# This section executes the MCCE4 workflow in the following order:
#   1. Optional membrane generation (STEP M)
#   2. STEP 1: Generate protonation states (step1.py)
#   3. Optional custom script A
#   4. STEP 2: Generate rotamers and conformers (step2.py) + membrane append
#   5. Optional custom script B
#   6. STEP 3: Monte Carlo sampling for titration (step3.py)
#   7. Optional custom script C
#   8. STEP 4: Analyze results and generate output (step4.py)
#   9. Cleanup temporary PBE data directories
#
# Each step includes:
#   - Conditional execution based on flags and dependencies
#   - Timing measurements and logging
#   - Output file verification
#   - Dependency checking (ensures prerequisite files exist)
#=================================================================================

#---------------------------------------------------------------------------------
# OPTIONAL STEP M: Membrane Generation
#
# Generate partial membrane structure if requested (membrane proteins only)
# Runs before STEP 1 to prepare membrane environment
# Execution: Controlled by stepM flag

if [[ "$stepM" == "t" ]]; then
    script_name=$(basename "$STEPM")

    if [[ -f "$STEPM" ]]; then
        echo "Running custom stepM script \"$script_name\" (before STEP1)..."
        start_time=$(date +%s)

        if eval "$STEPM_CMD"; then
            end_time=$(date +%s)
            elapsed=$((end_time - start_time))
            formatted_time=$(format_time "$elapsed")

            echo "$script_name completed SUCCESSFULLY in $formatted_time."
            printf "%-6s: %s   - Success: \"%s\" executed.\n" "STEPM" "$formatted_time" "$script_name" >> "$TIMING_FILE"
        else
            end_time=$(date +%s)
            elapsed=$((end_time - start_time))
            formatted_time=$(format_time "$elapsed")

            echo "$script_name FAILED after $formatted_time!"
            printf "%-6s: %s   - FAILED: \"%s\" encountered an error.\n" "STEPM" "$formatted_time" "$script_name" >> "$TIMING_FILE"
        fi
    else
        echo "Warning: stepM script \"$STEPM\" not found. Skipping stepM."
        printf "%-6s: skipped.           - Skipped: Script not found \"%s\"\n" "STEPM" "$script_name" >> "$TIMING_FILE"
    fi
fi

#---------------------------------------------------------------------------------
# STEP 1: Generate Protonation States
#
# Creates initial protonation states and conformers from input PDB
# Output: step1_out.pdb (all atoms with all possible protonation states)
# Execution: Always runs if step1 flag is "t"

if [[ "$step1" == "t" ]]; then
    time_step "STEP1" "$STEP1_CMD" "step1_out.pdb" "step1_out.pdb created." "$step1"
else
    echo "STEP1 SKIPPED (flag was f)."
    printf "%-6s: skipped.      - Skipped (flag f): step not run.\n" "STEP1" >> "$TIMING_FILE"
fi

#---------------------------------------------------------------------------------
# OPTIONAL STEP A: Custom Script After STEP 1
#
# User-defined processing between STEP 1 and STEP 2
# Common uses: Modify step1_out.pdb, add custom residues, filter conformers
# Execution: Controlled by stepA flag

if [[ "$stepA" == "t" ]]; then
    script_name=$(basename "$STEPA")

    if [[ -f "$STEPA" ]]; then
        echo "Running custom stepA script \"$script_name\" (between STEP1 and STEP2)..."
        start_time=$(date +%s)

        if eval "$STEPA_CMD"; then
            end_time=$(date +%s)
            elapsed=$((end_time - start_time))
            formatted_time=$(format_time "$elapsed")

            echo "$script_name completed SUCCESSFULLY in $formatted_time."
            printf "%-6s: %s   - Success: \"%s\" executed.\n" "STEPA" "$formatted_time" "$script_name" >> "$TIMING_FILE"
        else
            end_time=$(date +%s)
            elapsed=$((end_time - start_time))
            formatted_time=$(format_time "$elapsed")

            echo "$script_name FAILED after $formatted_time!"
            printf "%-6s: %s   - FAILED: \"%s\" encountered an error.\n" "STEPA" "$formatted_time" "$script_name" >> "$TIMING_FILE"
        fi
    else
        echo "Warning: stepA script \"$STEPA\" not found. Skipping stepA."
        printf "%-6s: skipped.           - Skipped: Script not found \"%s\"\n" "STEPA" "$script_name" >> "$TIMING_FILE"
    fi
fi

#---------------------------------------------------------------------------------
# STEP 2: Generate Rotamers and Conformers
#
# Creates rotamers and energy calculations for all conformers
# Output: step2_out.pdb (optimized conformers with energy data)
# Dependency: Requires step1_out.pdb from STEP 1
# Execution: Runs if step2 flag is "t" AND step1_out.pdb exists/was recently created
if [[ "$step1" == "t" && "$step2" == "t" ]]; then
    # Both step1 and step2 enabled: check if step1_out.pdb was recently created
    if file_just_made "$step1" "step1_out.pdb"; then
        time_step "STEP2" "$STEP2_CMD" "step2_out.pdb" "step2_out.pdb updated." "$step2"
    else
        echo "STEP2 SKIPPED — step1_out.pdb was not created recently."
        printf "%-6s: skipped.       - Reason: step1_out.pdb was not created recently.\n" "STEP2" >> "$TIMING_FILE"
    fi
elif [[ "$step1" == "f" && "$step2" == "t" ]]; then
    # Step1 disabled but step2 enabled: check if step1_out.pdb exists from previous run
    if [[ -f "step1_out.pdb" ]]; then
        time_step "STEP2" "$STEP2_CMD" "step2_out.pdb" "step2_out.pdb updated." "$step2"
    else
        echo "STEP2 SKIPPED — step1_out.pdb not found."
        printf "%-6s: skipped.       - Reason: step1_out.pdb not found.\n" "STEP2" >> "$TIMING_FILE"
    fi
else
    echo "STEP2 SKIPPED (flag was f)."
    printf "%-6s: skipped.      - Skipped (flag f): step not run.\n" "STEP2" >> "$TIMING_FILE"
fi

#---------------------------------------------------------------------------------
# Membrane Append Operation
#
# If membrane generation (STEP M) was run and STEP 2 completed successfully,
# append the membrane structure to step2_out.pdb
# This integrates the membrane environment with the protein conformers

if [[ "$step2" == "t" && "$stepM" == "t" ]]; then
    if [[ -f "PROT_MEM/MEM_step2_out.pdb" ]]; then
        # Backup original step2_out.pdb and append membrane
        mv "step2_out.pdb" "BK_step2_out.pdb"
        cat "BK_step2_out.pdb" "PROT_MEM/MEM_step2_out.pdb" > "step2_out.pdb"
        echo "MEM successfully appended to step2_out.pdb."
        printf "    : MEM appendment. - Success: PROT_MEM/MEM_step2_out.pdb appended onto step2_out.pdb.\n" >> "$TIMING_FILE"
    else
        echo "Warning: MEM_step2_out.pdb not found. Skipping MEM append."
        printf "    : MEM appendment. -  FAILED: File missing: PROT_MEM/MEM_step2_out.pdb\n" >> "$TIMING_FILE"
    fi
fi

#---------------------------------------------------------------------------------
# OPTIONAL STEP B: Custom Script After STEP 2
#
# User-defined processing between STEP 2 and STEP 3
# Common uses: Modify step2_out.pdb, adjust conformers, prepare for MC sampling
# Execution: Controlled by stepB flag
if [[ "$stepB" == "t" ]]; then
    script_name=$(basename "$STEPB")

    if [[ -f "$STEPB" ]]; then
        echo "Running custom stepB script \"$script_name\" (between STEP2 and STEP3)..."
        start_time=$(date +%s)

        if eval "$STEPB_CMD"; then
            end_time=$(date +%s)
            elapsed=$((end_time - start_time))
            formatted_time=$(format_time "$elapsed")

            echo "$script_name completed SUCCESSFULLY in $formatted_time."
            printf "%-6s: %s   - Success: \"%s\" executed.\n" "STEPB" "$formatted_time" "$script_name" >> "$TIMING_FILE"
        else
            end_time=$(date +%s)
            elapsed=$((end_time - start_time))
            formatted_time=$(format_time "$elapsed")

            echo "$script_name FAILED after $formatted_time!"
            printf "%-6s: %s   - FAILED: \"%s\" encountered an error.\n" "STEPB" "$formatted_time" "$script_name" >> "$TIMING_FILE"
        fi
    else
        echo "Warning: stepB script \"$STEPB\" not found. Skipping stepB."
        printf "%-6s: skipped.           - Skipped: Script not found \"%s\"\n" "STEPB" "$script_name" >> "$TIMING_FILE"
    fi
fi

#---------------------------------------------------------------------------------
# STEP 3: Monte Carlo Sampling
#
# Performs Monte Carlo sampling to determine protonation states at different pH values
# Output: head3.lst (MC results), energies/ directory (detailed energy files)
# Dependency: Requires step2_out.pdb from STEP 2
# Execution: Runs if step3 flag is "t" AND step2_out.pdb exists/was recently created
# Note: This is typically the longest-running step (hours to days depending on system size)

if [[ "$step2" == "t" && "$step3" == "t" ]]; then
    # Both step2 and step3 enabled: check if step2_out.pdb was recently created
    if file_just_made "$step2" "step2_out.pdb"; then
        time_step "STEP3" "$STEP3_CMD" "head3.lst" "head3.lst and energies directory created." "$step3"
    else
        echo "STEP3 SKIPPED — step2_out.pdb was not created recently."
        printf "%-6s: skipped.      - Reason: step2_out.pdb was not created recently.\n" "STEP3" >> "$TIMING_FILE"
    fi
elif [[ "$step2" == "f" && "$step3" == "t" ]]; then
    # Step2 disabled but step3 enabled: check if step2_out.pdb exists from previous run
    if [[ -f step2_out.pdb ]]; then
        time_step "STEP3" "$STEP3_CMD" "head3.lst" "head3.lst and energies directory created." "$step3"
    else
        echo "STEP3 SKIPPED — step2_out.pdb not found."
        printf "%-6s: skipped.      - Reason: step2_out.pdb not found.\n" "STEP3" >> "$TIMING_FILE"
    fi
else
    echo "STEP3 SKIPPED (flag was f)."
    printf "%-6s: skipped.      - Skipped (flag f): step not run.\n" "STEP3" >> "$TIMING_FILE"
fi

#---------------------------------------------------------------------------------
# OPTIONAL STEP C: Custom Script After STEP 3
#
# User-defined processing between STEP 3 and STEP 4
# Common uses: Parse head3.lst, extract specific pH results, custom analysis
# Execution: Controlled by stepC flag

if [[ "$stepC" == "t" ]]; then
    script_name=$(basename "$STEPC")

    if [[ -f "$STEPC" ]]; then
        echo "Running custom stepC script \"$script_name\" (between STEP3 and STEP4)..."
        start_time=$(date +%s)

        if eval "$STEPC_CMD"; then
            end_time=$(date +%s)
            elapsed=$((end_time - start_time))
            formatted_time=$(format_time "$elapsed")

            echo "$script_name completed SUCCESSFULLY in $formatted_time."
            printf "%-6s: %s   - Success: \"%s\" executed.\n" "STEPC" "$formatted_time" "$script_name" >> "$TIMING_FILE"
        else
            end_time=$(date +%s)
            elapsed=$((end_time - start_time))
            formatted_time=$(format_time "$elapsed")

            echo "$script_name FAILED after $formatted_time!"
            printf "%-6s: %s   - FAILED: \"%s\" encountered an error.\n" "STEPC" "$formatted_time" "$script_name" >> "$TIMING_FILE"
        fi
    else
        echo "Warning: stepC script \"$STEPC\" not found. Skipping stepC."
        printf "%-6s: skipped.           - Skipped: Script not found \"%s\"\n" "STEPC" "$script_name" >> "$TIMING_FILE"
    fi
fi

#---------------------------------------------------------------------------------
# STEP 4: Analyze and Generate Output
#
# Analyzes MC sampling results and generates final titration curves and pKa values
# Output: sum_crg.out (charge vs pH), pK.out (calculated pKa values), and other analysis files
# Dependency: Requires head3.lst and energies/ directory from STEP 3
# Execution: Runs if step4 flag is "t" AND head3.lst/energies exist/were recently created

if [[ "$step3" == "t" && "$step4" == "t" ]]; then
    # Both step3 and step4 enabled: check if head3.lst and energies were recently created
    if [[ -f "head3.lst" && -d "energies" ]] && file_just_made "$step3" "head3.lst"; then
        time_step "STEP4" "$STEP4_CMD" "sum_crg.out" "sum_crg.out created." "$step4"
    else
        echo "STEP4 SKIPPED — head3.lst was not created recently."
        printf "%-6s: skipped.      - Reason: head3.lst was not created recently.\n" "STEP4" >> "$TIMING_FILE"
    fi
elif [[ "$step3" == "f" && "$step4" == "t" ]]; then
    # Step3 disabled but step4 enabled: check if head3.lst and energies exist from previous run
    if [[ -f "head3.lst" && -d "energies" ]]; then
        time_step "STEP4" "$STEP4_CMD" "sum_crg.out" "sum_crg.out created." "$step4"
    else
        echo "STEP4 SKIPPED — head3.lst or energies directory not found."
        printf "%-6s: skipped.      - Reason: head3.lst or energies directory not found.\n" "STEP4" >> "$TIMING_FILE"
    fi
else
    echo "STEP4 SKIPPED (flag was f)."
    printf "%-6s: skipped.      - Skipped (flag f): step not run.\n" "STEP4" >> "$TIMING_FILE"
fi

#---------------------------------------------------------------------------------
# STEP_CLEAN: Remove Temporary PBE Data Directories
#
# Clean up temporary Poisson-Boltzmann electrostatics data from $TMP directory
# These files can be several GB in size and are safe to remove after successful completion
# Only removes directories matching the current working directory pattern
# Execution: Controlled by step_clean flag ("t" to clean, "f" to keep for debugging)
# Note: To keep PBE data for debugging, set step_clean="f" and add --debug to STEP3
# Generate flattened working directory path pattern for cleanup matching
WDIR=$(pwd -L)  # Get full working directory path (follow symlinks)
WDIR_FLAT=$(echo "$WDIR" | cut -d'/' -f3- | sed 's|/|.|g')  # Convert path to flat pattern

# Display cleanup configuration for user awareness
printf "%-65s %s\n" "Temporary PBE data will be stored in:" "$TMP"
printf "%-65s %s\n" "Full path of the working directory is:" "$WDIR"
printf "%-65s %s\n" "Subdirectories matching this pattern will be removed from temp dir:" "$WDIR_FLAT"

if [[ "$step_clean" == "t" ]]; then
    if [[ -d "$TMP" ]]; then
        echo "STEP_CLEAN: cleaning matching pbe_data folders in $TMP..."

        # Find all directories in $TMP belonging to current user matching working directory pattern
        MATCHING_DIRS=$(find "$TMP" -mindepth 1 -maxdepth 1 -type d -user "$USER" -name "*${WDIR_FLAT}*")

        if [[ -n "$MATCHING_DIRS" ]]; then
            # Remove matching directories
            echo "$MATCHING_DIRS" | xargs -r rm -rf
            printf "%-6s: done.     - Matching /tmp dirs removed:\n" "STEP_CLEAN" >> "$TIMING_FILE"
            echo "$MATCHING_DIRS" | sed 's/^/          - /' >> "$TIMING_FILE"
        else
            echo "STEP_CLEAN: no matching /tmp pbe_data folders found."
            printf "%-6s: skipped.  - No pbe_data dirs matching this run. ($WDIR_FLAT)\n" "STEP_CLEAN" >> "$TIMING_FILE"
        fi
    else
        echo "STEP_CLEAN: /tmp directory not found."
        printf "%-6s: skipped.  - /tmp directory not found.\n" "STEP_CLEAN" >> "$TIMING_FILE"
    fi
else
    echo "STEP_CLEAN: skipped (flag was f)."
    printf "%-6s: skipped.  - Skipped (flag f): /tmp cleanup not run.\n" "STEP_CLEAN" >> "$TIMING_FILE"
fi

#=================================================================================
# FINALIZATION: Calculate Total Runtime and Close Timing Log
#=================================================================================
# Record completion time, calculate total runtime, and finalize timing report

script_end_time=$(date +%s)
total_elapsed=$((script_end_time - script_start_time))
formatted_total_elapsed=$(format_time "$total_elapsed")

# Write completion information to timing log
echo -e "\nRun ended at: $(date)" >> "$TIMING_FILE"
printf "Total script runtime: %s\n" "$formatted_total_elapsed" >> "$TIMING_FILE"

# Brief pause before final messages
sleep 5

# Display completion status to user
echo "Script complete."
echo "Total runtime: $formatted_total_elapsed"
echo "Timing report written to $TIMING_FILE"

# END OF SCRIPT
