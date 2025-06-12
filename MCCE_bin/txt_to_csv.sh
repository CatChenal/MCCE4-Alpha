#!/bin/bash

# Check for -h flag
for arg in "$@"; do
    if [[ $arg == "-h" ]]; then
        echo "Usage: $(basename "$0") my_file.lst"
        echo "Note: A file extension is not required."
        echo "Options:"
        echo "  -h              Show this help message"
        exit 0
    fi
done

# Check if a file was provided
if [[ -z $1 ]]; then
    echo "Error: Missing input file."
    exit 1
fi

if [[ -x $1 ]]; then
    echo "Not handling executable files."
    exit 1
fi

if ! [[ -f $1 && -r $1 ]]; then
    printf "Input does not exist or is not readable\nUsage: $0 <input_file.xyz>"
    exit 1
fi

# check if input file is empty:
if ! [[ -s $1 ]]; then
    echo "Error: Input file is empty."
    exit 1
fi

in_file="$1"
# Replace last extension to '.csv'
out_file="${in_file%.*}.csv"
# Convert spaces or tabs to commas and write to output file
awk '{gsub(/[[:space:]]+/, ","); print}' "$in_file" > "$out_file"

printf "Converted $in_file to $out_file\n\n"
