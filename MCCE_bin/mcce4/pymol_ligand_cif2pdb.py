#!/usr/bin/env python3

"""
Module: pymol_ligand_cif2pdb.py

Created on Aug 19 19:20:00 2025
@author: Gehan Ranepura and Raihan Uddin

Convert ligands in .cif files to PDB format using PyMOL.

Usage:
Convert with explicit output filename:
    ligand_cif2pdb 4lzt.cif 4lzt_from_cif.pdb

Convert and let the script choose output (input name with .pdb extension):
    ligand_cif2pdb 4lzt.cif

Optional flags:
    --offline
        Prevents PyMOL from fetching missing ligand definitions from the internet.
        Use this if you see messages like:
            "Downloading https://files.rcsb.org/ligands/download/UNK.cif -> ./UNK.cif"

    --ligand-cache /path/to/dir
        Sets a local cache directory for ligand CIF files. If PyMOL needs a ligand
        definition (e.g., UNK) and it exists in this directory, it will use the local
        copy instead of downloading it.

Background:
    When loading some CIF files, PyMOL may encounter non-standard residues/ligands
    (e.g., "UNK"). To interpret them correctly, PyMOL sometimes fetches their
    definitions from the RCSB ligand dictionary. That's why you might see:
        Downloading https://files.rcsb.org/ligands/download/UNK.cif
        -> ./UNK.cif
    This is normal and does not affect the conversion; it just helps PyMOL parse the
    structure more accurately. Use --offline or a --ligand-cache to control this.
"""
import argparse
from pathlib import Path
import shutil
import sys
import tempfile
from typing import Union

try:
    from pymol import cmd, finish_launching
except ImportError:
    sys.exit("Please, visit https://pymol.org/conda/ to install PyMOL on your system.")


def list_nonpolymer_ligands(obj_name: str="molecule") -> list:
    """Return a sorted list of unique ligand 3-letter codes present (HETATM, excluding HOH).
    """
    resnames = set()
    model = cmd.get_model(f"(hetatm and not resn HOH) and model {obj_name}")
    for a in model.atom:
        if a.resn:
            resnames.add(a.resn.upper())
    return sorted(resnames)


def snapshot_cif_files(directory: str) -> set:
    """Return a set of 2-tuples: (uppercase file stem, filepath) 
    for all *.cif files in directory.
    """
    dir_fp = Path(directory)
    if directory and dir_fp.is_dir():
       out = set([(fp.stem.upper(), str(fp)) 
                  for fp in dir_fp.iterdir()
                  if fp.suffix == ".cif"])
    else:
        out = set()
    return out


def cif_only(filename: str):
    if not filename.endswith((".cif", ".CIF")):
        raise argparse.ArgumentTypeError(f"File '{filename}' must be a .cif file")
    return filename


def cli_parser():
    p = argparse.ArgumentParser(
        prog="ligand_cif2pdb",
        description="Convert ligands in .cif file to PDB using PyMOL",
        usage=("\n  %(prog)s 4lzt.cif 4lzt.pdb\n"
            "  %(prog)s 4lzt.cif\n"
            "  %(prog)s 8esw.cif --offline\n"
            "  %(prog)s 8esw.cif --ligand-cache ./ligands\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "input",
        type=cif_only,
        help="Input CIF file"
    )
    p.add_argument(
        "output",
        nargs="?",
        help="Output PDB file (optional)"
    )
    p.add_argument(
        "-ligand-cache",
        default=None,
        metavar="DIR",
        help="Directory to use as local cache for ligand CIF files"
    )
    p.add_argument(
        "--offline",
        default=False,
        action="store_true",
        help="Disable network fetching of missing ligand CIFs"
    )

    return p


def main(args: Union[argparse.Namespace, dict]):

    if isinstance(args, dict):
        args = argparse.Namespace(**args)
        cif_only(args.input)

    # Default output name if not provided
    output_file = args.output if args.output else Path(args.input).stem + ".pdb"

    # Start PyMOL without GUI
    finish_launching(['pymol', '-cq'])

    # Decide fetch_path target (for detecting downloads)
    created_temp_cache = False
    fetch_dir = None
    if args.ligand_cache:
        fetch_dir = Path(args.ligand_cache)
        fetch_dir.mkdir(exist_ok=True)
        fetch_dir = str(fetch_dir)
    elif not args.offline:
        fetch_dir = tempfile.mkdtemp(prefix="pymol_ligand_cache_")
        created_temp_cache = True

    if fetch_dir:
        cmd.set("fetch_path", fetch_dir)

    if args.offline:
        cmd.set("fetch_host", "")
        if not fetch_dir:
            cmd.set("fetch_path", "/dev/null")

    fileset_before = snapshot_cif_files(fetch_dir) if fetch_dir else set()

    # This is where PyMOL itself may print "Downloading ..." lines
    cmd.load(args.input, "molecule")

    ligands_present = list_nonpolymer_ligands("molecule")
    if ligands_present:
        print("\n🔎 Non-polymer ligands present in structure (excluding water):")
        print("   " + ", ".join(ligands_present))
    else:
        print("\n🔎 No non-polymer ligands (excluding water) detected.")

    fileset_after = snapshot_cif_files(fetch_dir) if fetch_dir else set()
    new_files = {code: path for (code, path) in fileset_after - fileset_before}

    if args.offline:
        print("\n🌐 Offline mode: network fetching disabled.")
        if new_files:
            print("⚠️  Unexpected: ligand files appeared despite --offline:")
            for code, path in sorted(new_files.items()):
                print(f"   - {code} -> {path}")
    else:
        if new_files:
            print("\n⬇️  Detected ligand definitions fetched during load:")
            for code in sorted(new_files):
                print(f"   - {code} -> {new_files[code]}")
        else:
            print("✅ No new ligand definitions were fetched during this conversion.")

    cmd.save(output_file, "molecule")
    print(f"\n✅ Converted {args.input} -> {output_file}")
    cmd.quit()
    # cleanup
    if created_temp_cache and fetch_dir and Path(fetch_dir).is_dir():
        try:
            shutil.rmtree(fetch_dir)
        except Exception:
            pass

    return


def cli():
    p = cli_parser()
    args = p.parse_args()
    main(args)


if __name__ == "__main__":
    cli()
