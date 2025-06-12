#!/usr/bin/env python

"""
Module: msout_file.py

Contains functions to process the msout files in a run folder
for testing purposes.

Functions:

* `msout_with_mc_lines`:
   Reduce an 'msout file' to smaller files with a known number of 
   accepted state lines: 10, 20, 50, 100, 200, 500, 1000.
 
"""

from subprocess import run as subp_run, CalledProcessError, CompletedProcess
from pathlib import Path

from mcce4.io_utils import get_mcce_filepaths


def msout_with_mc_lines(mcce_dir: str, ph: str, eh: str):
    """
    Create smaller msout files with these number of accepted state lines:
    10, 20, 50, 100, 200, 500, 1000.
    The smaller files have the count in their extension: ".txt.sm<target_mc>",
    i.e. pH7eH0ms.txt.sm10, and are save in the parent file location.
    """
    mcce_dir = Path(mcce_dir)

    # h3_fp & s2_fp output paths not used here:
    _, _, msout_fp = get_mcce_filepaths(mcce_dir, ph, eh)

    N = 23  # N :: number of msout file lines to get 10 mc accepted state lines

    # how many mc_lines to output; will be used in smaller file extension:
    target_mc = [10, 20, 50, 100, 200, 500, 1000]
    # head_n :: integers to pass to `head` command -n option:
    # sorted in descending order so that the largest small file wil be reused as
    # input for the others:
    head_n = sorted([N + x for x in target_mc], reverse=True)

    largest_fp = str(msout_fp.with_suffix(f".txt.sm{target_mc[-1]}"))

    # create smaller msout files
    for i, hn in enumerate(head_n, start=1):
        msout_in = largest_fp
        if i == 1:
            msout_in = str(msout_fp)
 
        outname = str(msout_fp.with_suffix(f".txt.sm{target_mc[-i]}"))
        # create sed command:
        cmd_str = "sed '/^MC\:1/q'" + f" {msout_in} | head -n{hn} > {outname}"
        #print(cmd_str)
        try:
            subp_run(cmd_str, cwd=mcce_dir, shell=True, check=True)
        except CalledProcessError as e:
            print(f"Command failed with exit code {e.returncode}:\n{e.stderr}")
