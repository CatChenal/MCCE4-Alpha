#!/bin/bash


# Job management utility script.

# set desc variable with text:
read -r -d '' desc << EndOfText
  Description: Job management utility script.
  ......................................
  Purpose:
    Given a folder path, list the python jobs running in it.
    If a second argument is is given, it is assumed to be a specifc script name.

  Context:
    Debugging failed jobs becomes an intractable task when the user keeps resubmitting
    jobs in the problematic run folder, which increases the queue insteasd of clearing it.

  Additional considerations if using the cron scheduler:
    If the problematic run folder appears in the output of this command:
    > crontab -l
    then the related lines must be deleted as well. Use the edit command:
    > crontab -e

  Script version: 12-03-24."
EndOfText

# ..............................................................
# functions
# ..............................................................

help_msg() {
  echo ""
  echo "Usage:"
  echo " Note: Options come first!"
  echo ""
  echo " $(basename "$0") [-h]"
  echo " $(basename "$0") [-K] folder_path [command or script name]"
  echo ""
  echo "  -h              Display this help message & exit."
  echo "  -K              List, then kill the jobs running in the folder_path."
  echo "  folder_path     Path to the target folder."
  echo ""
  echo "Examples:"
  echo " running_jobs_in_folder.sh .             # List all python jobs running in the current folder."
  echo " running_jobs_in_folder.sh . new2.sh     # List python AND given script jobs running in the current folder"
  echo " running_jobs_in_folder.sh -K .          # Kill python jobs running in current folder"
  echo " running_jobs_in_folder.sh -K . new2.sh  # Kill python AND given script jobs running in current folder"
  echo ""
  echo "$desc"
  exit 0
}

aborted() {
  echo "Caught user interrupt: exiting."
  exit
}

trap cleanup INT

process_cmd() {
  local folder="$1"
  local cmd="$2"
  local kill_switch="$3"

  if [[ $kill_switch == 1 ]]; then
    which="Processing"
  else
    which="Listing"
  fi

  printf "\n%s user's '%s' jobs, if any:\n" "$which" "$cmd"
  # list only the pid (no -a option):
  command="pgrep -u $USER $cmd"
  exec 3< <(eval "$command")
  i=0;
  while read -r pid <&3;
  do
    d=$(pwdx $pid)
    if [[ "$d" == "$folder" ]]; then
      printf 'pid: %d ; dir: %s' "$pid" "$d";
      if [[ $kill_switch == 1 ]]; then
        kill $pid;
      fi
      i=$(($i+1));
    fi
  done
  # Close file descriptor 3
  exec 3<&-

  if [[ $i == 0 ]]; then
    printf "\nNo '%s' jobs found running in target folder\n" "$cmd"
  else
    # show total:
    i=$(($i+1));
    if [[ $kill_switch == 1 ]]; then
        printf "\nKilled %d '%s' jobs running in target folder\n" "$i" "$cmd"
    else
        printf "\nFound %d '%s' jobs running in target folder\n" "$i" "$cmd"
    fi
  fi
}

# ..............................................................
# main
# ..............................................................
# preset boolean to indicate whether to also kill the jobs:
kill_switch=0

# Process options
while getopts "hK" opt; do
  case $opt in
    h) help_msg ;;
    K) kill_switch=1 ;;
    # Exit on invalid options:
    \?) printf "Invalid option: -%s\n" "$OPTARG" >&2
      exit 1
  esac
done
# Shift arguments to access folder_path after options processing
shift $((OPTIND-1))

# Get arguments from command line
if [[ -z "$1" ]]; then
  echo "Error: Missing folder path."
  exit 1
fi

if [[ "$1" == "." ]]; then
  folder_path="$(pwd)"
else
  folder_path="$1"
fi

printf "Target folder: %s\n" "$folder_path"
# always list python jobs:
cmd_name="python"
process_cmd $folder_path $cmd_name $kill_switch

if [[ "$#" == 2 ]]; then
  cmd_name="$2"
  if [[ $cmd_name != "python " ]]; then
    process_cmd $folder_path $cmd_name $kill_switch
  fi
fi

if [[ $kill_switch == 1 ]]; then
  printf "\nProcessing over.\n\n"
else
  printf "\nListing over.\n\n"
fi
