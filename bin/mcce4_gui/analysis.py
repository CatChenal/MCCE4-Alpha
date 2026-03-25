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
    Parse pK.out which has TWO sections:

    Section 1 — pKa summary:
        pH             pKa/Em  n(slope) 1000*chi2
        NTR+A0001_        7.948     1.012     0.001
        ARG+A0001_        >14.0
        ASP-A0003_        3.008     1.004     0.005

    Section 2 — titration curves:
         ph              0.0   1.0   2.0   3.0  ...  14.0
        NTR+A0001_      1.00  1.00  1.00  ...   0.00
        ASP-A0003_      0.00 -0.01 -0.09  ...  -1.00

    Returns:
        {
            "residues": [{"name": ..., "pka": ..., "n_slope": ..., "chi2": ...}, ...],
            "ph_values": [0.0, 1.0, 2.0, ...],
            "titration_curves": {"NTR+A0001_": [1.00, 1.00, ...], ...}
        }
    """
    if not os.path.isfile(path):
        return {"error": f"File not found: {path}"}

    with open(path, "r") as f:
        content = f.read()

    if not content.strip():
        return {"error": "pK.out is empty"}

    lines = content.strip().split("\n")

    residues = []
    ph_values = []
    titration_curves = {}

    # ── Find where section 2 starts ──
    # Section 2 header: a line starting with "ph" followed by float values
    # e.g. " ph              0.0   1.0   2.0   3.0 ..."
    # Section 1 header starts with "pH" followed by text labels like "pKa/Em"
    section2_start = None
    for i, line in enumerate(lines):
        stripped = line.strip()
        # Match: starts with "ph" (case-insensitive), then a float number
        m = re.match(r'^ph\s+([\d.-]+)', stripped, re.IGNORECASE)
        if m:
            section2_start = i
            break

    # ── Parse Section 1: pKa summary ──
    section1_end = section2_start if section2_start is not None else len(lines)
    for line in lines[:section1_end]:
        stripped = line.strip()
        if not stripped:
            continue
        # Skip header line ("pH   pKa/Em   n(slope)  1000*chi2")
        if stripped.lower().startswith("ph"):
            continue

        parts = stripped.split()
        if len(parts) < 2:
            continue

        res_name = parts[0]

        # Parse pKa — handle ">14.0", "<0.0", "nan", etc.
        pka_str = parts[1]
        pka = None
        try:
            cleaned = pka_str.lstrip("><")
            pka = float(cleaned)
        except ValueError:
            pka = None

        # Parse n(slope)
        n_slope = None
        if len(parts) > 2:
            try:
                n_slope = float(parts[2])
            except ValueError:
                pass

        # Parse 1000*chi2
        chi2 = None
        if len(parts) > 3:
            try:
                chi2 = float(parts[3])
            except ValueError:
                pass

        residues.append({
            "name": res_name,
            "pka": pka,
            "pka_display": pka_str,
            "n_slope": n_slope,
            "chi2": chi2,
        })

    # ── Parse Section 2: titration curves ──
    if section2_start is not None:
        # Parse pH values from header
        header_line = lines[section2_start].strip()
        header_parts = header_line.split()
        for part in header_parts[1:]:  # Skip "ph" label
            try:
                ph_values.append(float(part))
            except ValueError:
                break

        # Parse residue curves
        for line in lines[section2_start + 1:]:
            stripped = line.strip()
            if not stripped:
                continue

            parts = stripped.split()
            if len(parts) < 2:
                continue

            res_name = parts[0]
            curve = []
            for val_str in parts[1:]:
                try:
                    curve.append(float(val_str))
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

    Format:
         ph         0.0    1.0    2.0  ...  14.0
        NTR+A0001_  1.00   1.00   1.00 ...   0.00
        ...
        Total       6.00   6.00   5.91 ...  -9.05

    Returns:
        {
            "ph_values": [...],
            "residues": {"NTR+A0001_": [crg_at_ph0, ...], ...},
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

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue

        parts = stripped.split()

        # Detect header line: starts with "ph" followed by numbers
        if parts[0].lower() == "ph" and len(parts) > 1:
            for p in parts[1:]:
                try:
                    ph_values.append(float(p))
                except ValueError:
                    break
            continue

        if len(parts) < 2:
            continue

        res_name = parts[0]

        # Skip separator lines (e.g. "------...")
        if res_name.startswith("---"):
            continue

        charges = []
        for p in parts[1:]:
            try:
                charges.append(float(p))
            except ValueError:
                charges.append(0.0)

        # Detect summary rows
        name_lower = res_name.lower().replace("_", "")
        if name_lower in ("total", "net", "netcharge", "sum", "total:"):
            total_charge = charges
        elif name_lower in ("protons", "electrons"):
            # Store but don't add to residues
            pass
        else:
            residue_charges[res_name] = charges

    return {
        "ph_values": ph_values,
        "residues": residue_charges,
        "total_charge": total_charge,
    }


def parse_head3_lst(path: str = "head3.lst") -> dict:
    """Parse head3.lst for conformer information."""
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
