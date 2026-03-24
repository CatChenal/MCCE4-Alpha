"""
Analysis module for MCCE4 output parsing.

Parses pK.out, sum_crg.out, and other Step 4 outputs into structured
data suitable for Plotly visualization in the frontend.
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Optional


def parse_pk_out(path: str = "pK.out") -> dict:
    """
    Parse pK.out to extract pKa values and titration curves.

    Returns:
        {
            "residues": [{"name": "GLUA0035", "pka": 4.2, "n_crg": -1.0, "hill": 1.0}, ...],
            "ph_values": [0.0, 1.0, 2.0, ...],
            "titration_curves": {
                "GLUA0035": [occ_at_ph0, occ_at_ph1, ...],
                ...
            }
        }
    """
    if not os.path.isfile(path):
        return {"error": f"File not found: {path}"}

    with open(path, "r") as f:
        lines = f.readlines()

    if not lines:
        return {"error": "pK.out is empty"}

    residues = []
    ph_values = []
    titration_curves = {}

    # First line is typically a header with pH values
    header = lines[0].strip()
    header_parts = header.split()

    # Try to identify pH columns from header
    # Format varies: sometimes "Residue pKa n 1000*Hill  0.0  1.0  2.0 ..."
    ph_start_col = None
    for i, part in enumerate(header_parts):
        try:
            val = float(part)
            if ph_start_col is None:
                ph_start_col = i
            ph_values.append(val)
        except ValueError:
            if ph_start_col is not None:
                break  # Stop after continuous float block

    # Parse residue lines
    for line in lines[1:]:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        parts = line.split()
        if len(parts) < 2:
            continue

        res_name = parts[0]

        # Try to parse pKa (second column)
        try:
            pka = float(parts[1])
        except ValueError:
            pka = None

        # Parse n_crg (third column) if present
        n_crg = None
        if len(parts) > 2:
            try:
                n_crg = float(parts[2])
            except ValueError:
                pass

        # Parse Hill coefficient (sometimes column 3 or 4, as 1000*Hill)
        hill = None
        if len(parts) > 3:
            try:
                hill_raw = float(parts[3])
                # If stored as 1000*Hill
                if abs(hill_raw) > 10:
                    hill = hill_raw / 1000.0
                else:
                    hill = hill_raw
            except ValueError:
                pass

        residues.append({
            "name": res_name,
            "pka": pka,
            "n_crg": n_crg,
            "hill": hill,
        })

        # Extract titration curve data (occupancy at each pH)
        if ph_start_col is not None and len(parts) > ph_start_col:
            curve = []
            for j in range(ph_start_col, min(len(parts), ph_start_col + len(ph_values))):
                try:
                    curve.append(float(parts[j]))
                except ValueError:
                    curve.append(None)
            if curve:
                titration_curves[res_name] = curve

    return {
        "residues": residues,
        "ph_values": ph_values,
        "titration_curves": titration_curves,
    }


def parse_sum_crg(path: str = "sum_crg.out") -> dict:
    """
    Parse sum_crg.out for net charge vs pH/Eh data.

    Returns:
        {
            "ph_values": [...],
            "residues": {"GLUA0035": [crg_at_ph0, ...], ...},
            "total_charge": [total_at_ph0, ...],
        }
    """
    if not os.path.isfile(path):
        return {"error": f"File not found: {path}"}

    with open(path, "r") as f:
        lines = f.readlines()

    if not lines:
        return {"error": "sum_crg.out is empty"}

    ph_values = []
    residue_charges = {}
    total_charge = []

    # First line is header with pH values
    header = lines[0].strip()
    header_parts = header.split()

    # Find where pH values start
    for part in header_parts:
        try:
            ph_values.append(float(part))
        except ValueError:
            continue

    # Parse each residue line
    for line in lines[1:]:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        parts = line.split()
        if len(parts) < 2:
            continue

        res_name = parts[0]
        charges = []
        for p in parts[1:]:
            try:
                charges.append(float(p))
            except ValueError:
                charges.append(0.0)

        if res_name.lower() in ("total", "net", "sum"):
            total_charge = charges
        else:
            residue_charges[res_name] = charges

    return {
        "ph_values": ph_values,
        "residues": residue_charges,
        "total_charge": total_charge,
    }


def parse_head3_lst(path: str = "head3.lst") -> dict:
    """
    Parse head3.lst for conformer information.

    Returns:
        {"conformers": [{"name": ..., "occ": ..., "crg": ..., ...}, ...]}
    """
    if not os.path.isfile(path):
        return {"error": f"File not found: {path}"}

    conformers = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("iConf"):
                continue

            parts = line.split()
            if len(parts) < 5:
                continue

            try:
                conformers.append({
                    "index": int(parts[0]) if parts[0].isdigit() else 0,
                    "name": parts[1],
                    "occ": float(parts[2]) if len(parts) > 2 else 0.0,
                    "crg": float(parts[3]) if len(parts) > 3 else 0.0,
                    "em": float(parts[4]) if len(parts) > 4 else 0.0,
                })
            except (ValueError, IndexError):
                continue

    return {"conformers": conformers}


def get_available_outputs() -> dict:
    """Check which MCCE4 output files exist in the current directory."""
    files = {
        "pK.out": os.path.isfile("pK.out"),
        "sum_crg.out": os.path.isfile("sum_crg.out"),
        "head3.lst": os.path.isfile("head3.lst"),
        "step1_out.pdb": os.path.isfile("step1_out.pdb"),
        "step2_out.pdb": os.path.isfile("step2_out.pdb"),
        "energies": os.path.isdir("energies"),
        "run.prm": os.path.isfile("run.prm"),
        "run.log": os.path.isfile("run.log"),
    }
    return {"files": files, "workdir": os.getcwd()}
