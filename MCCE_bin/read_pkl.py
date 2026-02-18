#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


def main():
    # Setup argument parsing
    parser = argparse.ArgumentParser(
        description="Load and display a spatial features pickle file.",
        usage="read_pkl.py path/to/file.pkl"
        )

    # Adding arguments
    parser.add_argument("file_path",
                        type=str,
                        help="Path to the .pkl file"
                        )
    parser.add_argument("-v", "--verbose",
                        default=False,
                        action="store_true",
                        help="Display the entire DataFrame"
                        )
    parser.add_argument("-i", "--info",
                        default=False,
                        action="store_true",
                        help="Display DataFrame structure and memory usage"
                        )

    args = parser.parse_args()

    pkl_fp = Path(args.file_path)
    if not pkl_fp.exists():
        print(f"Error: File '{pkl_fp!s}' not found.")
        return

    try:
        # Load the pickle file into pandas
        df = pd.read_pickle(pkl_fp)
    except Exception as e:
        print(f"Could not load the pickle file into pandas: {e}")
        return
    
    # Global display settings for better terminal viewing
    pd.set_option('display.width', 1000)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_colwidth', 50) 

    if args.info:
        print(f"--- DataFrame Information ({pkl_fp!s}) ---")
        print(df.info())
        return

    n_rows, n_cols = df.shape
    if args.verbose:
        print(f"--- Full Dataset ({n_rows:,} rows) ---")
        print(df.to_string())
    else:
        # Show first 15 and last 15 with a manual separator
        # This avoids the Dtype concatenation error you encountered
        print(f"\n--- Summary: First 15 & Last 15 of {n_rows:,} rows ---")
        print(f"File: {pkl_fp.name}\n")

        if n_rows <= 30:  # small: print it all
            print(df)
        else:
            print(df.head(15))
            # Print a visual break that won't conflict with DataTypes
            # We use a mix of symbols to make it very obvious in the scroll
            if n_rows-30 > 15:
                # indicate gap between head & tail rows:
                print(f"\n{'='*20} [ DATA GAP: {n_rows-30} ROWS ] {'='*20}\n")
            print(df.tail(15))

    print(f"\nShape: {n_rows:,} rows, {n_cols:,} columns")


if __name__ == "__main__":
    main()
