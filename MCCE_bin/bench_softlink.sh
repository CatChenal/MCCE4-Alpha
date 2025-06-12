#!/bin/bash


# MCCE4 Benchmarking app utility script.


read -r -d '' desc << EndOfText
  Description: MCCE4 Benchmarking app utility script.
  ...................................................
  Purpose:
    Soft-link customized run.prm, extra.tpl or name.txt files in the current bench_dir
    as 'run.prm.user', 'extra.tpl' or 'name.txt' in each of the proteins' subfolders of
    bench_dir/runs/.

  Context:
    Each of the MCCE4 stepX.py programs has a -load_runprm command-line option.
    The value of this option is a file name. As MCCE4 retains the hard-coding of
    resource files of the previous versions, this file is expected to be in the
    protein directory were the program is launched.

  Use case for the benchmarking app:
    An existing bench_dir needs recalculated & the run script in the bench_dir/runs
    subfolder has been modified to accept the -load_runprm option on the stepX.py
    command(s) defined.
    In this case, the user run.prm file along with any of these resource files
    customized therein, i.e: extra.tpl, name.txt must be linked into the proteins
    folders.

  Warnings:
    1. The file paths in a customized run.prm file in a bench_dir folder must refer
       to the names - not paths! - of the linked files.
    2. After running this script, ensure that the runs/ script to be used passes the value
       of 'run.prm.user' to the -load_runprm option, e.g.:

       step1.py --dry -load_runprm run.prm.user

  Note:
    A future version of the benchmarking app will do this linking automatically
    if the -load_runprm option is used, or the -u option includes EXTRA or RENAME_RULES.
    Status as of 11-09-2024: WIP.

  Script version: 11-11-24.
EndOfText

help_msg() {
  echo "Usage:"
  echo " $(basename "$0") -h"
  echo " $(basename "$0") target_file target_type [ p | e | n ]"
  echo ""
  echo "  -h              Display this help message only"
  echo "  target_file     Name of the customized resource file in the current folder"
  echo "  target_type     Type of customized file: 'p': run.prm, 'e': extra.tpl, 'n': name.txt"
  echo "                  Used to determines the name of the linked files: 'run.prm.user', 'extra.tpl', 'name.txt'"
  echo ""
  echo "Examples:"
  echo " $(basename "$0") customized_extra.tpl e"
  echo " $(basename "$0") customized_name.txt n"
  echo ""
  echo "$desc"
  exit 0
}

# Process options using getopts
while getopts ":h" opt; do
  case $opt in
    h) help_msg ;;
    # Exit on invalid options:
    \?) printf 'Invalid option: -%s\n' "$OPTARG" >&2
      exit 1
      ;;
  esac
done

# Shift arguments to access target_file and link_name after options processing
shift $((OPTIND-1))

# Get arguments from command line
target_file="$1"
target_type="$2"

valid_types=("p" "e" "n")

if [[ ! " ${valid_types[@]} " =~ " $target_type " ]]; then
  printf 'Invalid target_type: %s; must be one of: p, e, n\n' "$target_type" >&2
  exit 1
fi

# Check if target exists
if [[ ! -f "$target_file" ]]; then
  printf 'Error: Customized file '%s' not found in %s\n' "$target_file" "$(pwd)"
  exit 1
 fi

# Check if standard bench_dir/runs dir exists:
runs_dir="runs"
# Check if runs directory exists
if [[ ! -d "$runs_dir" ]]; then
  printf 'Error: Subdirectory '%s' not found in %s\n' "$runs_dir" "$(pwd)"
  exit 1
fi

# Set the link_name according to target_type:
case $target_type in
  p) link_name="run.prm.user" ;;
  e) link_name="extra.tpl" ;;
  n) link_name="name.txt" ;;
  \?) printf 'Unrecognized type: %s\n' "$target_type" >&2
   exit 1
   ;;
esac

# Loop through subdirectories in runs
for subdir in "$runs_dir"/*; do
  # Skip non-directories
  if [[ ! -d "$subdir" ]]; then
    continue
  fi
  # Check if file already exists:
  linked="$subdir/$link_name"
  if [[ -f $linked ]]; then
    if [[ -L $linked ]]; then
      printf 'File %s is already linked.\n' $linked
    else
      printf 'WARNING: File %s exists but is not soft-linked!\n' $linked
    fi
    continue
  else
    # Create symlink using ln -s with relative path
    ln -s "../../$target_file" "$linked"
    printf 'Soft-linked %s as %s\n' $target_file $linked
  fi
done

printf '\a\n%s' "Processing over."
printf '\n%s' "If the runs/ script is using -load_runprm, make sure it passes 'run.prm.user'"
printf '\n%s' "as a value, e.g.: step1.py --dry -load_runprm run.prm.user"
printf '\n\n%s' "Alternatively*, if passing just a few customized file paths, the -u option should"
printf '\n%s' "used, e.g.: step2.py -d8 -u EXTRA=extra.tpl,RENAME_RULES=name.txt"
printf '\n*  %s\n\n' "As of 11-11-24, this alternative MUST be used until Issue #524 is resolved."
 