#!/usr/bin/env python
"""
Created on Apr 24 05:15:00 2025

@author: Gehan Ranepura
"""

import sys
import os

def modify_hohdm_by_fixed_char(input_file):
    if not input_file.endswith('.lst'):
        print("Error: Input file must be a .lst file (e.g., head3.lst)")
        sys.exit(1)

    output_file = os.path.splitext(input_file)[0] + '_HOHDMoff.lst'

    with open(input_file, 'r') as f:
        lines = f.readlines()

    modified_lines = []
    for line in lines:
        if 'HOHDM' in line:
            if len(line) > 21 and line[21] == 'f':
                line = line[:21] + 't' + line[22:]
        elif 'WATDM' in line:  # Added this condition
            if len(line) > 21 and line[21] == 'f':
                line = line[:21] + 't' + line[22:]
        modified_lines.append(line)

    with open(output_file, 'w') as f:
        f.writelines(modified_lines)

    print(f"Modified file saved as: {output_file}")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python script.py head3.lst")
        sys.exit(1)

    input_file = sys.argv[1]
    modify_hohdm_by_fixed_char(input_file)

