"""
parsers.py - Read MCCE4 output files for dipole analysis.

MCCE4 step2_out.pdb format (0-indexed columns):
    [0:6]    record type ("ATOM  " / "HETATM")
    [6:11]   atom serial number
    [12:16]  atom name (" CA ", " OD1", etc.)
    [17:20]  residue name (3 chars: "ARG", "NTR", "GLU")
    [20]     space
    [21]     chain ID ("A", "B", etc.)
    [22:26]  residue sequence number (4 digits: "0001")
    [26]     "_"
    [27:30]  conformer number (3 digits: "001")
    [30:38]  x coordinate (8.3f)
    [38:46]  y coordinate (8.3f)
    [46:54]  z coordinate (8.3f)
    [54:62]  radius
    [62:80]  charge (embedded partial charge from topology)
    [80:90]  history field — first 2 chars = conformer type ("01","BK","+1","-1")

    Full conformer ID = RES[17:20] + CONFTYPE[80:82] + CHAIN[21] + SEQ[22:26] + "_" + NUM[27:30]
    Example: "NTR01A0001_001", "ARGBKA0001_000", "ASP-1A0003_005"

fort.38 format (Monte Carlo occupancies):
    Line 1:  " ph     0.0   1.0   2.0  ..."   (pH values as column headers)
    Line 2+: "CONF_ID  occ_pH0  occ_pH1  ..."  (one row per conformer)

mcce.tpl format (concatenated topology):
    CHARGE lines: "CHARGE   ARGBK  CA   0.100"
    conf_key = RES(3) + CONFTYPE(2) = "ARGBK", "ARG01", "NTR+1"
"""

import os
import numpy as np
from pathlib import Path
from collections import defaultdict


# ---------------------------------------------------------------------------
# step2_out.pdb parser
# ---------------------------------------------------------------------------
class Atom:
    """Single atom parsed from step2_out.pdb."""
    __slots__ = ("serial", "name", "res_name", "conf_type", "chain",
                 "res_seq", "conf_id", "coords", "charge", "radius")

    def __init__(self, serial, name, res_name, conf_type, chain,
                 res_seq, conf_id, coords, charge, radius):
        self.serial = serial
        self.name = name          # e.g. "CA", "OD1"
        self.res_name = res_name  # e.g. "GLU", "NTR"
        self.conf_type = conf_type  # e.g. "BK", "01", "-1", "+1"
        self.chain = chain
        self.res_seq = res_seq    # int
        self.conf_id = conf_id    # full 14-char conformer ID, e.g. "NTR01A0001_001"
        self.coords = coords      # np.array([x, y, z])
        self.charge = charge      # partial charge (float)
        self.radius = radius

    @property
    def conf_key(self):
        """5-char conformer type key: RES(3)+CONFTYPE(2), e.g. 'ARG01', 'ARGBK'."""
        return self.res_name + self.conf_type

    @property
    def is_backbone(self):
        return self.conf_type == "BK"

    def __repr__(self):
        return f"Atom({self.serial}, {self.name}, {self.conf_id}, crg={self.charge:.3f})"


def parse_step2_pdb(path):
    """
    Parse step2_out.pdb into conformer groups.

    Reads coordinates, charges, and conformer identity for every atom.
    Conformer type is extracted from the history field (columns 80-81).
    Charges are read directly from the embedded charge field (columns 62-80).

    Returns:
        conformers: dict {conf_id: [Atom, ...]}
        all_atoms:  flat list of all Atom objects
    """
    conformers = defaultdict(list)
    all_atoms = []

    with open(path) as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            if len(line) < 54:
                continue

            try:
                serial = int(line[6:11].strip())
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain = line[21] if len(line) > 21 else " "
                res_seq_str = line[22:26].strip()
                res_seq = int(res_seq_str) if res_seq_str else 0
                conf_num = line[27:30].strip()

                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                radius = float(line[54:62]) if len(line) >= 62 else 0.0

                # Charge: embedded in columns 62-80
                charge = 0.0
                if len(line) >= 74:
                    charge_str = line[62:80].strip()
                    if charge_str:
                        try:
                            charge = float(charge_str)
                        except ValueError:
                            # If charge field has non-numeric data, try first token
                            parts = charge_str.split()
                            if parts:
                                try:
                                    charge = float(parts[0])
                                except ValueError:
                                    charge = 0.0

                # Conformer type: first 2 chars of history field (col 80+)
                conf_type = "??"
                if len(line) >= 82:
                    history = line[80:].strip()
                    if len(history) >= 2:
                        conf_type = history[:2]
                        # Clean up: "BK" stays "BK", "+1" stays "+1", "01" stays "01"
                        # Handle edge case where history might start with underscore
                        if conf_type[0] == "_":
                            conf_type = "BK"

                # Build full conformer ID to match head3.lst / fort.38 format:
                # RES(3) + CONFTYPE(2) + CHAIN(1) + SEQ(4) + "_" + NUM(3)
                conf_id = f"{res_name}{conf_type}{chain}{res_seq_str}_{conf_num}"

            except (ValueError, IndexError) as e:
                continue

            atom = Atom(serial, atom_name, res_name, conf_type, chain,
                        res_seq, conf_id, np.array([x, y, z]), charge, radius)
            conformers[conf_id].append(atom)
            all_atoms.append(atom)

    return dict(conformers), all_atoms


# ---------------------------------------------------------------------------
# head3.lst parser
# ---------------------------------------------------------------------------
def parse_head3lst(path):
    """
    Parse head3.lst for conformer metadata.

    Returns:
        dict {conf_id: {"iconf": int, "flag": str, "occ": float, "crg": float, ...}}
    """
    result = {}
    with open(path) as fh:
        for line in fh:
            line_stripped = line.strip()
            if not line_stripped or line_stripped.startswith("iConf") or line_stripped.startswith("-"):
                continue
            parts = line_stripped.split()
            if len(parts) < 5:
                continue
            try:
                iconf = int(parts[0])
                conf_id = parts[1]
                result[conf_id] = {
                    "iconf": iconf,
                    "flag": parts[2],
                    "occ": float(parts[3]),
                    "crg": float(parts[4]),
                    "em0": float(parts[5]) if len(parts) > 5 else 0.0,
                    "pka0": float(parts[6]) if len(parts) > 6 else 0.0,
                    "ne": int(parts[7]) if len(parts) > 7 else 0,
                    "nh": int(parts[8]) if len(parts) > 8 else 0,
                }
            except (ValueError, IndexError):
                continue
    return result


# ---------------------------------------------------------------------------
# fort.38 parser
# ---------------------------------------------------------------------------
def parse_fort38(path):
    """
    Parse fort.38 (Monte Carlo occupancy output).

    Format:
        Line 1:  " ph      0.0   1.0   2.0  ..."   (pH column headers)
        Line 2+: "CONF_ID  occ0  occ1  occ2  ..."   (one row per conformer)

    Returns:
        ph_values:   np.array of pH values
        conf_ids:    list of conformer ID strings
        occupancies: np.array of shape (n_ph, n_conf)
    """
    with open(path) as fh:
        lines = fh.readlines()

    if not lines:
        raise ValueError(f"fort.38 is empty: {path}")

    # --- Parse header: extract pH values ---
    header = lines[0].strip()
    header_parts = header.split()
    # First token should be "ph", rest are pH values
    ph_start = 0
    for i, token in enumerate(header_parts):
        try:
            float(token)
            ph_start = i
            break
        except ValueError:
            continue

    ph_values = np.array([float(x) for x in header_parts[ph_start:]])
    n_ph = len(ph_values)

    if n_ph == 0:
        raise ValueError(f"No pH values found in fort.38 header: {header}")

    # --- Parse data rows: each row is one conformer ---
    conf_ids = []
    occ_rows = []  # will be (n_conf, n_ph)

    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) < 2:
            continue

        conf_id = parts[0]
        try:
            occs = [float(x) for x in parts[1:]]
        except ValueError:
            continue

        # Pad or truncate to match number of pH values
        if len(occs) < n_ph:
            occs.extend([0.0] * (n_ph - len(occs)))
        occs = occs[:n_ph]

        conf_ids.append(conf_id)
        occ_rows.append(occs)

    # occ_rows is (n_conf, n_ph) — transpose to (n_ph, n_conf)
    occ_matrix = np.array(occ_rows)  # shape: (n_conf, n_ph)
    occupancies = occ_matrix.T        # shape: (n_ph, n_conf)

    return ph_values, conf_ids, occupancies


# ---------------------------------------------------------------------------
# Charge parser — supports both mcce.tpl and *.ftpl
# ---------------------------------------------------------------------------
def parse_tpl_charges(param_path):
    """
    Read CHARGE entries from MCCE4 topology files.

    Supports:
      - Single concatenated file: param/mcce.tpl (MCCE4 default)
      - Directory of .ftpl files: param/*.ftpl (MCCE3 style)
      - Direct path to a .tpl file

    CHARGE line format:
        CHARGE   ARGBK  CA   0.100
        conf_key = RES(3) + CONFTYPE(2) = "ARGBK"

    Returns:
        charges: dict {(conf_key, atom_name): float}
    """
    charges = {}
    param = Path(param_path)

    # Determine what we're looking at
    files_to_parse = []
    if param.is_file():
        # Direct path to a .tpl file
        files_to_parse = [param]
    elif param.is_dir():
        # Look for mcce.tpl first (MCCE4 default), then *.ftpl
        mcce_tpl = param / "mcce.tpl"
        if mcce_tpl.exists():
            files_to_parse = [mcce_tpl]
        else:
            files_to_parse = list(param.glob("*.ftpl"))
            if not files_to_parse:
                files_to_parse = list(param.glob("*.tpl"))
    else:
        raise FileNotFoundError(f"Parameter path not found: {param_path}")

    if not files_to_parse:
        raise FileNotFoundError(
            f"No topology files found in {param_path}. "
            f"Expected mcce.tpl or *.ftpl files."
        )

    for tpl_file in files_to_parse:
        with open(tpl_file) as fh:
            for line in fh:
                if not line.startswith("CHARGE"):
                    continue
                parts = line.split()
                if len(parts) < 4:
                    continue
                try:
                    conf_key = parts[1].strip()
                    atom_name = parts[2].strip()
                    charge = float(parts[3])
                    charges[(conf_key, atom_name)] = charge
                except (ValueError, IndexError):
                    continue

    return charges


# ---------------------------------------------------------------------------
# PQR parser (microstate / single state)
# ---------------------------------------------------------------------------
def parse_pqr(path):
    """
    Parse a PQR file for single-state dipole calculation.

    PQR format variations:
        ATOM serial name resname [chain] resseq x y z charge radius

    Returns:
        list of dicts with keys: name, res_name, chain, res_seq, coords, charge, radius
    """
    atoms = []
    with open(path) as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            parts = line.split()
            if len(parts) < 9:
                continue
            try:
                # PQR can have 10 fields (with chain) or 9 fields (without)
                if len(parts) >= 11:
                    # serial name resname chain resseq x y z charge radius
                    atom = {
                        "serial": int(parts[1]),
                        "name": parts[2],
                        "res_name": parts[3],
                        "chain": parts[4],
                        "res_seq": int(parts[5]),
                        "coords": np.array([float(parts[6]), float(parts[7]), float(parts[8])]),
                        "charge": float(parts[9]),
                        "radius": float(parts[10]),
                    }
                elif len(parts) >= 10:
                    # Try: might have chain or not
                    # If parts[4] is numeric → no chain column
                    try:
                        int(parts[4])
                        # No chain column
                        atom = {
                            "serial": int(parts[1]),
                            "name": parts[2],
                            "res_name": parts[3],
                            "chain": "",
                            "res_seq": int(parts[4]),
                            "coords": np.array([float(parts[5]), float(parts[6]), float(parts[7])]),
                            "charge": float(parts[8]),
                            "radius": float(parts[9]),
                        }
                    except ValueError:
                        # Has chain column
                        atom = {
                            "serial": int(parts[1]),
                            "name": parts[2],
                            "res_name": parts[3],
                            "chain": parts[4],
                            "res_seq": int(parts[5]),
                            "coords": np.array([float(parts[6]), float(parts[7]), float(parts[8])]),
                            "charge": float(parts[9]),
                            "radius": float(parts[-1]),
                        }
                else:
                    # Minimal: serial name resname resseq x y z charge radius
                    atom = {
                        "serial": int(parts[1]),
                        "name": parts[2],
                        "res_name": parts[3],
                        "chain": "",
                        "res_seq": int(parts[4]),
                        "coords": np.array([float(parts[5]), float(parts[6]), float(parts[7])]),
                        "charge": float(parts[8]),
                        "radius": float(parts[-1]) if len(parts) > 9 else 0.0,
                    }
                atoms.append(atom)
            except (ValueError, IndexError):
                continue
    return atoms
