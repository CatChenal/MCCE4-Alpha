#!/bin/bash
#SBATCH --job-name=mcce_run
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4         # Adjust number of cores if needed
#SBATCH --mem=24G                 # Adjust memory if needed

# MCCE Simulation
STEP1="step1.py prot.pdb --dry -d 4"
STEP2="step2.py -l 1 -d 4"
STEP3="step3.py -d 4 -s zap -salt 0.05 --new_vdw"
STEP4="step4.py --xts"

# Set MCCE4 Parameters
MCCE_HOME="/home/granepura/MCCE4"
USER_PARAM="./user_param"
EXTRA="./user_param/extra.tpl"

# Check if EXTRA exists; if not, use fallback
if [ -f "$EXTRA" ]; then
    EXTRA="$EXTRA"
else
    EXTRA="$MCCE_HOME/extra.tpl"
fi

# Inititiate timing log and set to exit on errors for critical parts
set -e
TIMING_FILE="mcce_timing.log"
echo "MCCE Timing Report" > $TIMING_FILE
echo "==================" >> $TIMING_FILE

# Print MCCE Parameters used
echo "MCCE_HOME: $MCCE_HOME" >> $TIMING_FILE
echo "EXTRA: $EXTRA" >> $TIMING_FILE 
if [ -d "$USER_PARAM" ]; then
    echo "USER_PARAM: $USER_PARAM" >> $TIMING_FILE
else
    echo "USER_PARAM: N/A" >> $TIMING_FILE
fi

echo -e "==================\n" >> $TIMING_FILE

# Define commands
PARAM="-u MCCE_HOME=$MCCE_HOME,EXTRA=$EXTRA"
STEP1_CMD="$STEP1 $PARAM > step1.log"
STEP2_CMD="$STEP2 $PARAM > step2.log"
STEP3_CMD="$STEP3 $PARAM > step3.log"
STEP4_CMD="$STEP4 $PARAM > step4.log"


# Function to check if file was just made
function file_just_made {
    file="$1"
    if [[ -f "$file" ]] && [[ $(find "$file" -mmin -1) ]]; then
        return 0
    else
        return 1
    fi
}

# Helper to time and record a step
function time_step {
    step_name="$1"
    step_cmd="$2"
    success_output="$3"    # what output file(s) we expect
    success_msg="$4"       # what to write in the success message

    echo "Running $step_name ..."
    start_time=$(date +%s)

    # Run the step; catch any errors but don't stop the script
    if eval $step_cmd; then
        end_time=$(date +%s)
        elapsed=$((end_time - start_time))

        # Check if expected output was updated
        if file_just_made "$success_output"; then
            echo "$step_name completed SUCCESSFULLY in $elapsed seconds."
            printf "%-6s: %2d seconds.   - Success: %s\n" "$step_name" "$elapsed" "$success_msg" >> $TIMING_FILE
        else
            echo "$step_name completed, but expected output $success_output was NOT updated!"
            printf "%-6s: %2d seconds.   - Failed: expected output $success_output not updated.\n" "$step_name" "$elapsed" >> $TIMING_FILE
        fi
    else
        end_time=$(date +%s)
        elapsed=$((end_time - start_time))
        echo "$step_name FAILED after $elapsed seconds!"
        printf "%-6s: %2d seconds.   - Failed: see console output.\n" "$step_name" "$elapsed" >> $TIMING_FILE
    fi
}

# STEP 1
time_step "STEP1" "$STEP1_CMD" "step1_out.pdb" "step1_out.pdb updated."

# STEP 2 — only if step1_out.pdb was just made
if file_just_made "step1_out.pdb"; then
    time_step "STEP2" "$STEP2_CMD" "step2_out.pdb" "step2_out.pdb updated."
else
    echo "Skipping step2.py — step1_out.pdb not updated."
    printf "%-6s: skipped.      - Reason: step1_out.pdb not updated.\n" "STEP2" >> $TIMING_FILE
fi

# STEP 3 — only if step2_out.pdb was just made
if file_just_made "step2_out.pdb"; then
    time_step "STEP3" "$STEP3_CMD" "head3.lst" "head3.lst and energies updated."
else
    echo "Skipping step3.py — step2_out.pdb not updated."
    printf "%-6s: skipped.      - Reason: step2_out.pdb not updated.\n" "STEP3" >> $TIMING_FILE
fi

# STEP 4 — always run if head3.lst and energies/ just made; ensure pK.out updated
if [[ -f head3.lst ]] && [[ -d energies ]] && file_just_made "head3.lst"; then
    time_step "STEP4" "$STEP4_CMD" "pK.out" "pK.out updated."
else
    echo "Skipping step4.py — head3.lst or energies directory not newly created."
    printf "%-6s: skipped.      - Reason: head3.lst or energies not updated.\n" "STEP4" >> $TIMING_FILE
fi

# Wait 10 seconds
sleep 10

echo "Script complete."
echo "Timing report written to $TIMING_FILE"

