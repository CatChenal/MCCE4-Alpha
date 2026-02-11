#!/usr/bin/env python

import pandas as pd
import argparse
import os

def main():
    # Setup argument parsing
    parser = argparse.ArgumentParser(description="Load and display a spatial features pickle file.")

    # Adding arguments
    parser.add_argument("file_path", nargs="?", help="Path to the .pkl file")
    parser.add_argument("-v", "--verbose", action="store_true", help="Display the entire DataFrame")
    parser.add_argument("-i", "--info", action="store_true", help="Display DataFrame structure and memory usage")

    args = parser.parse_args()

    # 1. Check if the argument was provided
    if not args.file_path:
        print("Error: Please specify a .pkl file as an argument.")
        print("Usage: python read_pkl.py path/to/file.pkl")
        return

    try:
        # 2. Check if the file actually exists
        if not os.path.exists(args.file_path):
            print(f"Error: File '{args.file_path}' not found.")
            return

        # Load the pickle file
        df = pd.read_pickle(args.file_path)

        # Global display settings for better terminal viewing
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', 1000)
        pd.set_option('display.max_colwidth', 50) 

        # Option: Just show info
        if args.info:
            print(f"--- DataFrame Information ({args.file_path}) ---")
            print(df.info())
            return

        if args.verbose:
            # Show everything
            print(f"--- Full Dataset ({len(df)} rows) ---")
            print(df.to_string())
        else:
            # Show first 15 and last 15 with a manual separator
            # This avoids the Dtype concatenation error you encountered
            print(f"\n--- Summary: First 15 & Last 15 of {len(df)} rows ---")
            print(f"File: {os.path.basename(args.file_path)}\n")

            if len(df) <= 30:
                # If the file is small, just print it all
                print(df)
            else:
                # Print the top slice
                print(df.head(15))
                
                # Print a visual break that won't conflict with DataTypes
                # We use a mix of symbols to make it very obvious in the scroll
                print(f"\n{'='*20} [ DATA GAP: {len(df)-30} ROWS ] {'='*20}\n")
                
                # Print the bottom slice
                print(df.tail(15))

        print(f"\nShape: {df.shape[0]} rows, {df.shape[1]} columns")

    except Exception as e:
        print(f"An error occurred while reading the pickle: {e}")

if __name__ == "__main__":
    main()
