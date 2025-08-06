#!/bin/bash

# =========================================================================================
# Script Name   : driver_mcce4.sh
# Purpose       : Automate and control the full MCCE4 simulation pipeline including optional custom preprocessing steps.
#
# MAIN EXECUTION CALL — Run the core MCCE4 workflow
#
# This script defines all necessary environment variables and configuration parameters 
# (input files, flags, paths, simulation commands, and optional step scripts), and then 
# hands off execution to a secondary driver script.
#
# The driver script (`mcce4_driver.sh`) performs the actual step-by-step execution of 
# the MCCE4 simulation pipeline using the parameters defined above.
#
# NOTE:
# - All variables are exported to the environment so they can be accessed by the driver.
# - The driver script assumes all configurations are pre-validated by this script.
# - This modular structure improves maintainability and separates setup from execution.
#
# =========================================================================================

# Inititiate timing log and set to exit on errors for critical parts
TIMING_FILE="mcce_timing.log"
echo "MCCE Timing Report" > $TIMING_FILE
echo "====================================" >> $TIMING_FILE
set -e

# Finalize MCCE Parameters Options:
# Check if EXTRA exists; if not, use fallback
echo "MCCE_HOME: $MCCE_HOME" >> $TIMING_FILE

# Check if EXTRA exists; if not, use fallback
if [ -f "$EXTRA" ]; then
    EXTRA="$EXTRA"
else
    EXTRA="$MCCE_HOME/extra.tpl"
fi
echo "EXTRA: $EXTRA" >> $TIMING_FILE

# Check if USER_PARAM exists; if not, print N/A
if [ -d "$USER_PARAM" ]; then
    echo "USER_PARAM: $USER_PARAM" >> $TIMING_FILE
else
    echo "USER_PARAM: N/A" >> $TIMING_FILE
fi

# If TMP is not /tmp, create the directory if it doesn't exist; otherwise, empty its contents
if [[ "$TMP" != "/tmp" ]]; then
    echo "TMP is not set to /tmp. Using custom path: $TMP"
    if [[ ! -d "$TMP" ]]; then
        echo "Creating directory: $TMP"
        mkdir -p "$TMP"
        echo "TMP: $TMP" >> $TIMING_FILE
    else
        echo "Cleaning existing directory: $TMP"
        rm -rf "$TMP"/*
        echo "TMP: $TMP" >> $TIMING_FILE
    fi
fi

echo -e "====================================\n" >> $TIMING_FILE
echo "Optional Features Enabled:" >> $TIMING_FILE
[[ "$stepM" == "t" ]] && echo " - Partial Membrane Generation (stepM)" >> $TIMING_FILE
[[ "$stepA" == "t" ]] && echo " - Custom Script After Step1 (stepA)" >> $TIMING_FILE
[[ "$stepB" == "t" ]] && echo " - Custom Script After Step2 (stepB)" >> $TIMING_FILE
[[ "$stepC" == "t" ]] && echo " - Custom Script After Step3 (stepC)" >> $TIMING_FILE
echo -e "====================================\n" >> $TIMING_FILE
script_start_time=$(date +%s)
echo "Run started at: $(date)" >> $TIMING_FILE

# Append --norun if corresponding step is disabled (flag == "f")
[[ "$step1" == "f" ]] && STEP1="$STEP1 --norun"
[[ "$step2" == "f" ]] && STEP2="$STEP2 --norun"
[[ "$step3" == "f" ]] && STEP3="$STEP3 --norun"
[[ "$step4" == "f" ]] && STEP4="$STEP4 --norun"

#-----------------------------------------------------------------------------------------
#=========================================================================================
# Finalize MCCE step commands to run
if [ -d "$USER_PARAM" ]; then
    PARAM="-u MCCE_HOME=$MCCE_HOME,EXTRA=$EXTRA,USER_PARAM=$USER_PARAM"
else
    PARAM="-u MCCE_HOME=$MCCE_HOME,EXTRA=$EXTRA"
fi


# Run orientation.py to center the structure and extract the new filename from its output.
# Update input_pdb, print the final filename, and exit if the file doesn't exist.
output=$(orientation.py "$input_pdb")
input_pdb=$(echo "$output" | grep "Structure moved to origin" | awk '{print $NF}')
if [[ ! -f "$input_pdb" ]]; then
    echo "Error: Centered file '$input_pdb' was not created. Exiting."
    exit 1
else
    echo "Final input PDB file: $input_pdb"
fi


# Finalize MCCE simulation steps
STEP1_CMD="$STEP1 \$input_pdb $PARAM > step1.log"
STEP2_CMD="$STEP2 $PARAM > step2.log"
STEP3_CMD="$STEP3 $PARAM > step3.log 2>&1"  # Redirect stderr to stdout for step3
STEP4_CMD="$STEP4 $PARAM > step4.log"

# Finalize Optional Scripts
if command -v sbatch &> /dev/null; then
    STEPM_CMD="sbatch --wait $STEPM"  # --wait makes sbatch wait for job completion
else
    STEPM_CMD="bash $STEPM > stepM.log"
fi
STEPA_CMD="python $STEPA > stepA.log"
STEPB_CMD="python $STEPB > stepB.log"
STEPC_CMD="python $STEPC > stepC.log"

# Define cleanup function to remove temporary pbe_data directories
cleanup_on_exit() {
    echo ""
    echo ">>> Caught termination signal. Running STEP_CLEAN before exiting..."

    WDIR=$(pwd -L)
    WDIR_FLAT=$(echo "$WDIR" | cut -d'/' -f3- | sed 's|/|.|g')

    if [ "$step_clean" = "t" ] && [ -d "$TMP" ]; then
        MATCHING_DIRS=$(find "$TMP" -mindepth 1 -maxdepth 1 -type d -user "$USER" -name "*${WDIR_FLAT}*")
        if [ -n "$MATCHING_DIRS" ]; then
            echo "$MATCHING_DIRS" | xargs -r rm -rf
            printf "%-6s: done.     - STEP_CLEAN after cancellation.\n" "STEP_CLEAN" >> "$TIMING_FILE"
        else
            printf "%-6s: skipped.  - No matching tmp dirs at cancel.\n" "STEP_CLEAN" >> "$TIMING_FILE"
        fi
    fi
}

# Trap script termination signals (e.g., Ctrl+C or scancel)
trap cleanup_on_exit SIGINT SIGTERM EXIT
#-----------------------------------------------------------------------------------------
#=========================================================================================

# Function to check if file was just made
function file_just_made {
    step_flag="$1"   # "t" or "f"
    file="$2"

    if [[ "$step_flag" == "f" ]]; then
        # If step was disabled, treat file as OK to proceed.
        return 0
    fi

    if [[ -f "$file" ]] && [[ $(find "$file" -mmin -5) ]]; then
        return 0
    else
        return 1
    fi
}

# Function to format elapsed time
format_time() {
    local elapsed=$1
    local hours=$((elapsed / 3600))
    local minutes=$(( (elapsed % 3600) / 60 ))
    local seconds=$((elapsed % 60))
    printf "%02dh:%02dm:%02ds" "$hours" "$minutes" "$seconds"
}

# Helper to time and record a step
function time_step {
    step_name="$1"
    step_cmd="$2"
    success_output="$3"
    success_msg="$4"
    step_flag="$5"    # "t" or "f"

    echo "Running $step_name ..."
    printf "%-6s: Running $step_name ...\n" "$step_name" >> "$TIMING_FILE"
    start_time=$(date +%s)

    # Run the step (even if --norun is present)
    if eval "$step_cmd"; then
        end_time=$(date +%s)
        elapsed=$((end_time - start_time))
        formatted_time=$(format_time "$elapsed")

        # Remove the "running..." line
        grep -v "^$step_name *: Running $step_name ..." "$TIMING_FILE" > "${TIMING_FILE}.tmp" && mv "${TIMING_FILE}.tmp" "$TIMING_FILE"

        # Log final status to mcce_timing.log:
        # - If step was skipped, mark as "Skipped"
        # - If output file was correctly updated, mark as "Success"
        # - If step ran but expected output was not updated, mark as "Failed"
        if [[ "$step_flag" == "f" ]]; then
            # If the step was intentionally skipped, mark as Skipped.
            echo "$step_name SKIPPED (flag was f)."
            printf "%-6s: %s   - Skipped (flag f): step not run.\n" "$step_name" "$formatted_time" >> "$TIMING_FILE"
        elif file_just_made "$step_flag" "$success_output"; then
            # Success case
            echo "$step_name completed SUCCESSFULLY in $formatted_time."
            printf "%-6s: %s   - Success: %s\n" "$step_name" "$formatted_time" "$success_msg" >> "$TIMING_FILE"
        else
            # Failed — but only if step_flag == "t"
            echo "$step_name completed, but expected output $success_output was NOT updated!"
            printf "%-6s: %s   - Failed: expected output %s not updated.\n" "$step_name" "$formatted_time" "$success_output" >> "$TIMING_FILE"
        fi
    else
        # eval failed
        end_time=$(date +%s)
        elapsed=$((end_time - start_time))
        formatted_time=$(format_time "$elapsed")

        if [[ "$step_flag" == "f" ]]; then
            echo "$step_name SKIPPED (flag was f)."
            printf "%-6s: %s   - Skipped (flag f): step not run.\n" "$step_name" "$formatted_time" >> "$TIMING_FILE"
        else
            echo "$step_name FAILED after $formatted_time!"
            printf "%-6s: %s   - Failed: see console output.\n" "$step_name" "$formatted_time" >> "$TIMING_FILE"
        fi
    fi
}


#--------------------------------------------------------------------------------------------------------------
# START: RUN MCCE4 Simulations
#==============================================================================================================

# Optional: STEP M — Run if stepM is enabled (flag = "t") and STEPM exists and is executable
if [[ "$stepM" == "t" ]]; then
    if [[ -x "$STEPM" ]]; then
        echo "Running optional stepM (partial membrane generation)..."
        time_step "STEPM" "$STEPM" "PROT_MEM/MEM_step2_out.pdb" "MEM_step2_out.pdb updated." "$stepM"
        input_pdb="prot.pdb"
        echo "New Centered Input PDB to have MEM appended to is: \"$input_pdb\""
        [[ -f "prot_center.pdb" ]] && rm "prot_center.pdb"
    else
        echo "STEP M SKIPPED — $STEPM not found or not executable."
        printf "%-6s: skipped.      - Reason: STEPM not found or not executable (%s)\n" "STEPM" "$STEPM" >> "$TIMING_FILE"
    fi
fi


# STEP 1 — Run if step1 is enabled (flag = "t")
if [[ "$step1" == "t" ]]; then
    time_step "STEP1" "$STEP1_CMD" "step1_out.pdb" "step1_out.pdb updated." "$step1"
else
    echo "STEP1 SKIPPED (flag was f)."
    printf "%-6s: skipped.      - Skipped (flag f): step not run.\n" "STEP1" >> $TIMING_FILE
fi


# Optional: STEP A — Run if stepA is enabled (flag = "t") and STEPA exists and is executable
if [[ "$stepA" == "t" ]]; then
    script_name=$(basename "$STEPA")

    if [[ -f "$STEPA" ]]; then
        echo "Running custom stepC script \"$script_name\" (between STEP1 and STEP2)..."
        start_time=$(date +%s)

        if python "$STEPA_CMD"; then
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


# STEP 2 — Run if step2 is enabled (flag = "t") AND:
#           - If step1 was enabled (flag = "t"), step1_out.pdb was recently created
#           - If step1 was skipped (flag = "f"), step1_out.pdb exists
if [[ "$step1" == "t" && "$step2" == "t" ]]; then
    if file_just_made "$step1" "step1_out.pdb"; then
        time_step "STEP2" "$STEP2_CMD" "step2_out.pdb" "step2_out.pdb updated." "$step2"
    else
        echo "STEP2 SKIPPED — step1_out.pdb was not created recently."
        printf "%-6s: skipped.       - Reason: step1_out.pdb was not created recently.\n" "STEP2" >> "$TIMING_FILE"
    fi
elif [[ "$step1" == "f" && "$step2" == "t" ]]; then
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
# MEM appendment — only if stepM = t and step2 ran successfully
if [[ "$step2" == "t" && "$stepM" == "t" ]]; then
    if [[ -f "PROT_MEM/MEM_step2_out.pdb" ]]; then
        mv "step2_out.pdb" "BK_step2_out.pdb"
        cat "BK_step2_out.pdb" "PROT_MEM/MEM_step2_out.pdb" > "step2_out.pdb"
        echo "MEM successfully appended to step2_out.pdb."
        printf "    : MEM appendment. - Success: PROT_MEM/MEM_step2_out.pdb appended onto step2_out.pdb.\n" >> "$TIMING_FILE"
    else
        echo "Warning: MEM_step2_out.pdb not found. Skipping MEM append."
        printf "    : MEM appendment. -  FAILED: File missing: PROT_MEM/MEM_step2_out.pdb\n" >> "$TIMING_FILE"
    fi
fi


# Optional: STEP B — Run if stepB is enabled (flag = "t") and STEPB exists and is executable
if [[ "$stepB" == "t" ]]; then
    script_name=$(basename "$STEPB")

    if [[ -f "$STEPB" ]]; then
        echo "Running custom stepB script \"$script_name\" (between STEP2 and STEP3)..."
        start_time=$(date +%s)

        if python "$STEPB_CMD"; then
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


# STEP 3 — Run if step3 is enabled (flag = "t") AND:
#           - If step2 was enabled (flag = "t"), step2_out.pdb was recently created
#           - If step2 was skipped (flag = "f"), step2_out.pdb exists
if [[ "$step2" == "t" && "$step3" == "t" ]]; then
    if file_just_made "$step2" "step2_out.pdb"; then
        time_step "STEP3" "$STEP3_CMD" "head3.lst" "head3.lst and energies directory created." "$step3"
    else
        echo "STEP3 SKIPPED — step2_out.pdb was not created recently."
        printf "%-6s: skipped.      - Reason: step2_out.pdb was not created recently.\n" "STEP3" >> "$TIMING_FILE"
    fi
elif [[ "$step2" == "f" && "$step3" == "t" ]]; then
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


# Optional: STEP C — Run if stepC is enabled (flag = "t") and STEPC exists and is executable
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


# STEP 4 — Run if step4 is enabled (flag = "t") AND:
#           - If step3 was enabled (flag = "t"), head3.lst and energies directory was recently created
#           - If step3 was skipped (flag = "f"), head3.lst and energies directory exists
if [[ "$step3" == "t" && "$step4" == "t" ]]; then
    if [[ -f "head3.lst" && -d "energies" ]] && file_just_made "$step3" "head3.lst"; then
        time_step "STEP4" "$STEP4_CMD" "sum_crg.out" "sum_crg.out created." "$step4"
    else
        echo "STEP4 SKIPPED — head3.lst was not created recently."
        printf "%-6s: skipped.      - Reason: head3.lst was not created recently.\n" "STEP4" >> "$TIMING_FILE"
    fi
elif [[ "$step3" == "f" && "$step4" == "t" ]]; then
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


# STEP_CLEAN — Remove temporary pbe_data folders from /tmp related to this specific working directory.
#              Ensures no leftover pbe_data from this run remains.
#              To keep these files for debugging, set step_clean="f" and add --debug option to STEP3.
WDIR=$(pwd -L)
WDIR_FLAT=$(echo "$WDIR" | cut -d'/' -f3- | sed 's|/|.|g')
printf "%-65s %s\n" "Temporary PBE data will be stored in:" "$TMP"
printf "%-65s %s\n" "Full path of the working directory is:" "$WDIR"
printf "%-65s %s\n" "Subdirectories matching this pattern will be removed from temp dir:" "$WDIR_FLAT"

if [[ "$step_clean" == "t" ]]; then
    if [[ -d "$TMP" ]]; then
        echo "STEP_CLEAN: cleaning matching pbe_data folders in $TMP..."

        MATCHING_DIRS=$(find "$TMP" -mindepth 1 -maxdepth 1 -type d -user "$USER" -name "*${WDIR_FLAT}*")

        if [[ -n "$MATCHING_DIRS" ]]; then
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


#--------------------------------------------------------------------------------------------------------------
# END: Run MCCE4 Simulations
#==============================================================================================================

# Finish Run Prompt and Wait 10 seconds
script_end_time=$(date +%s)
total_elapsed=$((script_end_time - script_start_time))
formatted_total_elapsed=$(format_time "$total_elapsed")

echo -e "\nRun ended at: $(date)" >> "$TIMING_FILE"
printf "Total script runtime: %s\n" "$formatted_total_elapsed" >> "$TIMING_FILE"

sleep 5
echo "Script complete."
echo "Total runtime: $formatted_total_elapsed"
echo "Timing report written to $TIMING_FILE"

