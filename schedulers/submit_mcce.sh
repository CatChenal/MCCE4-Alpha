#!/bin/bash
#SBATCH --job-name=mcce_run
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4         # Adjust number of cores if needed
#SBATCH --mem=24G                 # Adjust memory if needed

#=============================================================================
# Input and Output:
input_pdb="prot.pdb"    # (INPDB)

# Step control flags
step1="t"               # STEP1: pre-run, pdb-> mcce pdb  (DO_PREMCCE)
step2="t"               # STEP2: make rotamers            (DO_ROTAMERS)
step3="f"               # STEP3: Energy calculations      (DO_ENERGY)
step4="f"               # STEP4: Monte Carlo Sampling     (DO_MONTE)

# Set MCCE4 Parameters
MCCE_HOME="/home/granepura/MCCE4"
USER_PARAM="./user_param"
EXTRA="./user_param/extra.tpl"

# MCCE Simulation
STEP1="step1.py $input_pdb -d 4 --noter --dry"
STEP2="step2.py -l 1 -d 4"
STEP3="step3.py -d 4 -s zap -salt 0.05"
STEP4="step4.py --xts"
#==============================================================================

# Inititiate timing log and set to exit on errors for critical parts
set -e
TIMING_FILE="mcce_timing.log"
echo "MCCE Timing Report" > $TIMING_FILE
echo "====================================" >> $TIMING_FILE

# Print MCCE Parameters used
# Check if EXTRA exists; if not, use fallback
# Check if USER_PARAM exists; if not, print N/A
echo "MCCE_HOME: $MCCE_HOME" >> $TIMING_FILE

if [ -f "$EXTRA" ]; then
    EXTRA="$EXTRA"
else
    EXTRA="$MCCE_HOME/extra.tpl"
fi
echo "EXTRA: $EXTRA" >> $TIMING_FILE 

if [ -d "$USER_PARAM" ]; then
    echo "USER_PARAM: $USER_PARAM" >> $TIMING_FILE
else
    echo "USER_PARAM: N/A" >> $TIMING_FILE
fi

echo -e "====================================\n" >> $TIMING_FILE
script_start_time=$(date +%s)
echo "Run started at: $(date)" >> $TIMING_FILE

# Append --norun if corresponding step is disabled (flag == "f")
[[ "$step1" == "f" ]] && STEP1="$STEP1 --norun"
[[ "$step2" == "f" ]] && STEP2="$STEP2 --norun"
[[ "$step3" == "f" ]] && STEP3="$STEP3 --norun"
[[ "$step4" == "f" ]] && STEP4="$STEP4 --norun"

# Finalize MCCE step commands to run
PARAM="-u MCCE_HOME=$MCCE_HOME,EXTRA=$EXTRA"
STEP1_CMD="$STEP1 $PARAM > step1.log"
STEP2_CMD="$STEP2 $PARAM > step2.log"
STEP3_CMD="$STEP3 $PARAM > step3.log"
STEP4_CMD="$STEP4 $PARAM > step4.log"


# Function to check if file was just made
function file_just_made {
    step_flag="$1"   # "t" or "f"
    file="$2"

    if [[ "$step_flag" == "f" ]]; then
        # If step was disabled, treat file as OK to proceed.
        return 0
    fi

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
    success_output="$3"
    success_msg="$4"
    step_flag="$5"    # NEW argument — "t" or "f"

    echo "Running $step_name ..."
    start_time=$(date +%s)

    # Run the step (even if --norun is present)
    if eval "$step_cmd"; then
        end_time=$(date +%s)
        elapsed=$((end_time - start_time))

        if [[ "$step_flag" == "f" ]]; then
            # If the step was intentionally skipped, mark as Skipped.
            echo "$step_name SKIPPED (flag was f)."
            printf "%-6s: %2d seconds.   - Skipped (flag f): step not run.\n" "$step_name" "$elapsed" >> $TIMING_FILE
        elif file_just_made "$step_flag" "$success_output"; then
            # Success case
            echo "$step_name completed SUCCESSFULLY in $elapsed seconds."
            printf "%-6s: %2d seconds.   - Success: %s\n" "$step_name" "$elapsed" "$success_msg" >> $TIMING_FILE
        else
            # Failed — but only if step_flag == "t"
            echo "$step_name completed, but expected output $success_output was NOT updated!"
            printf "%-6s: %2d seconds.   - Failed: expected output $success_output not updated.\n" "$step_name" "$elapsed" >> $TIMING_FILE
        fi
    else
        # eval failed
        end_time=$(date +%s)
        elapsed=$((end_time - start_time))

        if [[ "$step_flag" == "f" ]]; then
            echo "$step_name SKIPPED (flag was f)."
            printf "%-6s: %2d seconds.   - Skipped (flag f): step not run.\n" "$step_name" "$elapsed" >> $TIMING_FILE
        else
            echo "$step_name FAILED after $elapsed seconds!"
            printf "%-6s: %2d seconds.   - Failed: see console output.\n" "$step_name" "$elapsed" >> $TIMING_FILE
        fi
    fi
}


# STEP 1 — run only if step1="t", else log Skipped
if [[ "$step1" == "t" ]]; then
    time_step "STEP1" "$STEP1_CMD" "step1_out.pdb" "step1_out.pdb updated." "$step1"
else
    echo "STEP1 SKIPPED (flag was f)."
    printf "%-6s: skipped.      - Skipped (flag f): step not run.\n" "STEP1" >> $TIMING_FILE
fi

# STEP 2 — run only if STEP1 was "t" and step1_out.pdb updated, OR if STEP1 was "f"
if file_just_made "$step1" "step1_out.pdb"; then
    if [[ "$step2" == "t" ]]; then
        time_step "STEP2" "$STEP2_CMD" "step2_out.pdb" "step2_out.pdb updated." "$step2"
    else
        echo "STEP2 SKIPPED (flag was f)."
        printf "%-6s: skipped.      - Skipped (flag f): step not run.\n" "STEP2" >> $TIMING_FILE
    fi
else
    if [[ "$step2" == "f" ]]; then
        echo "STEP2 SKIPPED (flag was f)."
        printf "%-6s: skipped.      - Skipped (flag f): step not run.\n" "STEP2" >> $TIMING_FILE
    else
        echo "Skipping STEP2 — step1_out.pdb not updated."
        printf "%-6s: skipped.      - Reason: step1_out.pdb not updated.\n" "STEP2" >> $TIMING_FILE
    fi
fi

# STEP 3 — run only if STEP2 was "t" and step2_out.pdb updated, OR if STEP2 was "f"
if file_just_made "$step2" "step2_out.pdb"; then
    if [[ "$step3" == "t" ]]; then
        time_step "STEP3" "$STEP3_CMD" "head3.lst" "head3.lst and energies updated." "$step3"
    else
        echo "STEP3 SKIPPED (flag was f)."
        printf "%-6s: skipped.      - Skipped (flag f): step not run.\n" "STEP3" >> $TIMING_FILE
    fi
else
    if [[ "$step3" == "f" ]]; then
        echo "STEP3 SKIPPED (flag was f)."
        printf "%-6s: skipped.      - Skipped (flag f): step not run.\n" "STEP3" >> $TIMING_FILE
    else
        echo "Skipping STEP3 — step2_out.pdb not updated."
        printf "%-6s: skipped.      - Reason: step2_out.pdb not updated.\n" "STEP3" >> $TIMING_FILE
    fi
fi

# STEP 4 — run only if STEP3 was "t" and head3.lst updated, OR if STEP3 was "f"
if [[ -f head3.lst ]] && [[ -d energies ]] && file_just_made "$step3" "head3.lst"; then
    if [[ "$step4" == "t" ]]; then
        time_step "STEP4" "$STEP4_CMD" "pK.out" "pK.out updated." "$step4"
    else
        echo "STEP4 SKIPPED (flag was f)."
        printf "%-6s: skipped.      - Skipped (flag f): step not run.\n" "STEP4" >> $TIMING_FILE
    fi
else
    if [[ "$step4" == "f" ]]; then
        echo "STEP4 SKIPPED (flag was f)."
        printf "%-6s: skipped.      - Skipped (flag f): step not run.\n" "STEP4" >> $TIMING_FILE
    else
        echo "Skipping STEP4 — head3.lst or energies directory not newly created."
        printf "%-6s: skipped.      - Reason: head3.lst or energies not updated.\n" "STEP4" >> $TIMING_FILE
    fi
fi

# Finsish Run Prompt and Wait 10 seconds
script_end_time=$(date +%s)
total_elapsed=$((script_end_time - script_start_time))
echo -e "\nRun ended at: $(date)" >> $TIMING_FILE
printf "Total script runtime: %d seconds\n" "$total_elapsed" >> $TIMING_FILE

sleep 10
echo "Script complete."
echo "Timing report written to $TIMING_FILE"

