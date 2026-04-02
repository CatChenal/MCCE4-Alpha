"""
MCCE4 pipeline tools — template generation, calibration, parsing.
"""

import os
import re
import shutil
import logging
import subprocess
import json
from pathlib import Path
from typing import Optional
from ..config import DIELECTRIC_MAP, RCSB_SMILES_URL


def run_cmd(cmd, description="", capture=False, cwd=None):
    """Run a shell command with logging."""
    if description:
        logging.info(f"▶ {description}")
    logging.debug(f"  CMD: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=capture, text=True, cwd=cwd)
    if result.returncode != 0:
        logging.error(f"  Command failed (exit code {result.returncode})")
        if capture and result.stderr:
            logging.error(f"  STDERR: {result.stderr.strip()}")
        return None
    return result


def extract_lig_id_from_pdb(pdb_path: str) -> str:
    """Extract 3-letter ligand code from PDB HETATM/ATOM records."""
    res_counts = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("HETATM", "ATOM  ")):
                resname = line[17:20].strip()
                if resname and resname not in ("HOH", "WAT", "TIP"):
                    res_counts[resname] = res_counts.get(resname, 0) + 1
    if not res_counts:
        return Path(pdb_path).stem.upper()[:3]
    return max(res_counts, key=res_counts.get)


def fetch_rcsb_info(lig_id: str) -> dict:
    """Fetch SMILES and metadata from RCSB chemical component dictionary."""
    logging.info(f"🔍 Fetching info for '{lig_id}' from RCSB...")
    url = RCSB_SMILES_URL.format(lig_id=lig_id)
    result = run_cmd(f'curl -s "{url}"', description="Querying RCSB", capture=True)
    if result is None or not result.stdout.strip():
        return {}

    try:
        data = json.loads(result.stdout)
        chem = data.get("chem_comp", {})
        info = {
            "name": chem.get("name", "Unknown"),
            "formula": chem.get("formula", "Unknown"),
            "formal_charge": chem.get("pdbx_formal_charge", 0) or 0,
            "smiles": None,
        }
        # Get SMILES
        descriptors = data.get("rcsb_chem_comp_descriptor", [])
        if isinstance(descriptors, list):
            for desc in descriptors:
                if isinstance(desc, dict) and "SMILES" in desc.get("type", ""):
                    info["smiles"] = desc.get("descriptor")
                    if "CANONICAL" in desc.get("type", ""):
                        break
        for k, v in info.items():
            if v is not None:
                logging.info(f"  ✓ {k}: {v}")
        return info
    except (json.JSONDecodeError, KeyError) as e:
        logging.warning(f"  RCSB parse error: {e}")
        return {}


def generate_ftpl_template(lig_id: str, pdb_path: str,
                            conformer_labels: list, ftpl_path: str) -> int:
    """Generate .ftpl template via pdb2ftpl.py. Returns count of unfilled entries."""
    conf_str = " ".join(conformer_labels)
    logging.info(f"📝 Generating .ftpl: {lig_id} with conformers = {conf_str}")
    result = run_cmd(f"pdb2ftpl.py -p {pdb_path} -c {conf_str}",
                     description="pdb2ftpl.py", capture=True)
    if result is None:
        logging.error("  pdb2ftpl.py failed!")
        return -1
    with open(ftpl_path, "w") as f:
        f.write(result.stdout)
    unfilled = result.stdout.count("to_be_filled")
    logging.info(f"  ✓ Template: {ftpl_path} ({unfilled} charges to fill)")
    return unfilled


def fill_ftpl_charges(ftpl_path: str, charges: dict, lig_id: str) -> int:
    """Replace 'to_be_filled' with charges. Returns count of unfilled entries."""
    logging.info(f"✏️  Filling charges in {ftpl_path}")
    with open(ftpl_path) as f:
        lines = f.readlines()

    pat = re.compile(
        r'^(CHARGE,\s*' + re.escape(lig_id) +
        r'(\S+),\s*"(.{4})":\s*)to_be_filled(.*)$'
    )
    filled = unfilled = 0
    new_lines = []
    for line in lines:
        m = pat.match(line)
        if m:
            prefix, conf_id, atom_name, suffix = m.groups()
            matched = False
            # Exact match
            if atom_name in charges:
                new_lines.append(f"{prefix}{charges[atom_name]:7.3f} # auto-filled{suffix}\n")
                filled += 1; matched = True
            # Stripped match
            if not matched:
                stripped = atom_name.strip()
                for cn, cv in charges.items():
                    if cn.strip() == stripped:
                        new_lines.append(f"{prefix}{cv:7.3f} # auto-filled{suffix}\n")
                        filled += 1; matched = True; break
            if not matched:
                new_lines.append(line); unfilled += 1
        else:
            new_lines.append(line)

    with open(ftpl_path, "w") as f:
        f.writelines(new_lines)
    logging.info(f"  ✓ Filled {filled}, unfilled {unfilled}")
    return unfilled


def update_conformer_params(ftpl_path: str, states: list, lig_id: str):
    """Update nH and pKa0 in CONFORMER lines."""
    logging.info("📋 Updating CONFORMER parameters...")
    with open(ftpl_path) as f:
        content = f.read()

    for s in states:
        ct = f"{lig_id}{s['label']}" if isinstance(s, dict) else f"{lig_id}{s.label}"
        nH = (s.get('nH', 0) if isinstance(s, dict) else s.nH) or 0
        try:
            nH = int(nH)
        except (ValueError, TypeError):
            nH = 0

        pka_raw = (s.get('pka', 0) if isinstance(s, dict) else s.pka)
        try:
            pka = float(str(pka_raw).lstrip("~≈><≥≤ ")) if pka_raw is not None else 0.0
        except (ValueError, TypeError):
            pka = 0.0

        p = re.compile(rf"(CONFORMER,\s*{re.escape(ct)}:.*?nH=)\s*(\d+)")
        if p.search(content):
            content = p.sub(rf"\g<1>{nH}", content)
        p2 = re.compile(rf"(CONFORMER,\s*{re.escape(ct)}:.*?pKa0=)\s*([-\d.]+)")
        if p2.search(content):
            content = p2.sub(rf"\g<1>{pka:.2f}", content)
        logging.info(f"  ✓ {ct}: nH={nH}, pKa0={pka:.2f}")

    with open(ftpl_path, "w") as f:
        f.write(content)


def run_rxn_calibration(lig_id: str, ftpl_path: str, pdb_path: str,
                         conformer_labels: list, dielectrics: list, work_dir: str):
    """Run MCCE4 steps 1-3, calibrate rxn values for all dielectrics."""
    user_param = os.path.join(work_dir, "user_param")
    os.makedirs(user_param, exist_ok=True)

    # Symlink ftpl
    link = os.path.join(user_param, os.path.basename(ftpl_path))
    if os.path.exists(link):
        os.remove(link)
    os.symlink(os.path.abspath(ftpl_path), link)
    logging.info(f"  🔗 Linked {ftpl_path} → user_param/")

    # Ensure PDB accessible
    pdb_base = os.path.basename(pdb_path)
    pdb_wd = os.path.join(work_dir, pdb_base)
    if not os.path.exists(pdb_wd):
        shutil.copy2(os.path.abspath(pdb_path), pdb_wd)

    # Steps 1-2 once
    if run_cmd(f"step1.py {pdb_base}", description="MCCE4 Step 1", cwd=work_dir) is None:
        return {"error": "Step 1 failed"}
    if run_cmd("step2.py", description="MCCE4 Step 2", cwd=work_dir) is None:
        return {"error": "Step 2 failed"}

    rxn_results = {}

    for eps in dielectrics:
        rxn_key = DIELECTRIC_MAP.get(eps)
        if not rxn_key:
            continue
        logging.info(f"\n{'='*60}")
        logging.info(f"  🔬 RXN Calibration ε={eps} ({rxn_key})")
        logging.info(f"{'='*60}")

        run_cmd(f"step3.py -d {eps}", description=f"Step 3 (ε={eps})", cwd=work_dir)
        head3 = os.path.join(work_dir, "head3.lst")
        dsolv = parse_head3_dsolv(head3, lig_id, conformer_labels)
        if not dsolv:
            logging.error("  Could not parse dsolv"); continue

        for ct, v in dsolv.items():
            logging.info(f"  📊 {ct}: dsolv = {v:.3f}")
        update_rxn(ftpl_path, dsolv, rxn_key)

        # Validate
        logging.info("  🔄 Validating...")
        run_cmd(f"step3.py -d {eps}", description=f"Validation (ε={eps})", cwd=work_dir)
        check = parse_head3_dsolv(head3, lig_id, conformer_labels)
        ok = all(abs(v) <= 0.01 for v in check.values())
        for ct, v in check.items():
            sym = "✓" if abs(v) <= 0.01 else "⚠"
            logging.info(f"  {sym} {ct}: dsolv = {v:.3f}")
        if ok:
            logging.info(f"  🎉 {rxn_key} calibration successful!")

        rxn_results[rxn_key] = dsolv

    return rxn_results


def parse_head3_dsolv(head3_path: str, lig_id: str, labels: list) -> dict:
    """Parse most negative dsolv per conformer type from head3.lst.

    head3.lst format:
    iConf CONFORMER     FL  occ  crg  Em0 pKa0 ne nH vdw0 vdw1 tors epol dsolv extra history
    00001 EMH01_0000_001 f 0.00 0.000  0  0.00  0  0 -15.9 0.0  0.0  0.0  -8.1  0.0 01O000 t

    For each conformer type (e.g., EMH01, EMH+1), take the most negative dsolv.
    """
    if not os.path.exists(head3_path):
        return {}
    with open(head3_path) as f:
        lines = f.readlines()

    types = [f"{lig_id}{l}" for l in labels]
    logging.info(f"  Looking for conformer types: {types}")

    # Find dsolv column from header
    header = next((l for l in lines if "dsolv" in l.lower() and "iConf" in l), None)
    if not header:
        logging.warning("  No header with 'dsolv' found in head3.lst")
        return {}
    try:
        col = header.split().index("dsolv")
    except ValueError:
        logging.warning("  'dsolv' column not found in header")
        return {}
    logging.info(f"  dsolv is column {col} in head3.lst")

    vals = {t: 0.0 for t in types}
    found_any = {t: False for t in types}

    for line in lines:
        if line.strip().startswith("iConf") or not line.strip():
            continue
        parts = line.split()
        if len(parts) <= col:
            continue
        # CONFORMER name is parts[1], e.g., "EMH01_0000_001" or "EMH+1_0000_003"
        if len(parts) < 2:
            continue
        conf_name = parts[1]

        for t in types:
            # Match conformer type: "EMH01" is prefix of "EMH01_0000_001"
            if conf_name.startswith(t):
                try:
                    v = float(parts[col])
                    found_any[t] = True
                    if v < vals[t]:
                        vals[t] = v
                        logging.info(f"    {t}: dsolv={v:.3f} (from {conf_name})")
                except (ValueError, IndexError):
                    pass

    for t in types:
        if not found_any[t]:
            logging.warning(f"  ⚠ No conformers found for {t} in head3.lst")
        else:
            logging.info(f"  📊 {t}: most negative dsolv = {vals[t]:.3f}")

    return vals


def update_rxn(ftpl_path: str, dsolv_values: dict, rxn_key: str):
    """Update rxn values in CONFORMER lines."""
    with open(ftpl_path) as f:
        content = f.read()
    for ct, dsolv in dsolv_values.items():
        p = re.compile(rf"(CONFORMER,\s*{re.escape(ct)}:.*?{rxn_key}=\s*)([-\d.]+)")
        if p.search(content):
            content = p.sub(rf"\g<1>{dsolv:8.3f}", content)
            logging.info(f"    {ct}: {rxn_key} = {dsolv:.3f}")
    with open(ftpl_path, "w") as f:
        f.write(content)


# ═════════════════════════════════════════════════════════════════════════════
# v3: Per-State Template Generation & Merge
# ═════════════════════════════════════════════════════════════════════════════

def generate_ftpl_for_state(lig_id, state_pdb, label, output_dir="."):
    """Call pdb2ftpl.py for a SINGLE protonation state PDB.

    Command: pdb2ftpl.py -p EMH_01.pdb -c 01
    Produces a .ftpl with CONNECT and CHARGE for that state only.

    Returns path to the generated ftpl, or None on failure.
    """
    try:
        pdb_base = os.path.basename(state_pdb)
        cmd = f"pdb2ftpl.py -p {pdb_base} -c {label}"
        logging.info(f"  Running: {cmd}")
        result = run_cmd(cmd, description=f"pdb2ftpl for {label}",
                         capture=True, cwd=output_dir)
        if result is None:
            logging.error(f"  pdb2ftpl.py failed for state {label}")
            return None

        # pdb2ftpl.py writes to stdout — save it
        ftpl_name = f"{lig_id}_{label.replace('+', 'p').replace('-', 'm')}.ftpl"
        ftpl_path = os.path.join(output_dir, ftpl_name)

        # If pdb2ftpl wrote to stdout, save it
        if result.stdout.strip():
            with open(ftpl_path, "w") as f:
                f.write(result.stdout)
        else:
            # Check if it wrote a file directly
            default_ftpl = os.path.join(output_dir, f"{lig_id}.ftpl")
            if os.path.exists(default_ftpl):
                shutil.move(default_ftpl, ftpl_path)
            else:
                logging.error(f"  No ftpl output for state {label}")
                return None

        logging.info(f"  → {ftpl_path}")
        return ftpl_path

    except Exception as e:
        logging.error(f"  pdb2ftpl.py failed: {e}")
        return None


def generate_all_state_ftpls(lig_id, state_pdbs, output_dir="."):
    """Call pdb2ftpl.py for each protonation state PDB.

    Args:
        lig_id: 3-letter ligand code
        state_pdbs: {label: pdb_path}
        output_dir: Working directory

    Returns:
        {label: ftpl_path} for successful states
    """
    ftpls = {}
    for label, pdb_path in state_pdbs.items():
        logging.info(f"\n  Generating ftpl for state {label}...")
        # Copy state PDB to output_dir so pdb2ftpl can find it
        pdb_dest = os.path.join(output_dir, os.path.basename(pdb_path))
        if not os.path.exists(pdb_dest):
            shutil.copy2(pdb_path, pdb_dest)

        ftpl = generate_ftpl_for_state(lig_id, pdb_path, label, output_dir)
        if ftpl:
            ftpls[label] = ftpl
        else:
            logging.error(f"  FAILED: No ftpl for state {label}")
    return ftpls


def merge_ftpl_files(per_state_ftpls, lig_id, output_path):
    """Merge per-state ftpl files into a single master .ftpl.

    Each state's ftpl has its own CONNECT, CHARGE, RADIUS, and CONFORMER
    blocks.  The merge:
      1. Rebuilds CONFLIST to list ALL conformer states (not just neutral).
      2. Combines CONNECT, CHARGE, RADIUS, and CONFORMER blocks from every
         state in sorted label order.

    Returns True on success.
    """
    if not per_state_ftpls:
        logging.error("No ftpl files to merge")
        return False

    neutral_label = "01" if "01" in per_state_ftpls else list(per_state_ftpls.keys())[0]
    all_labels = sorted(per_state_ftpls.keys())

    parsed = {}
    for label, ftpl_path in per_state_ftpls.items():
        parsed[label] = _parse_ftpl_sections(ftpl_path, lig_id)

    base = parsed[neutral_label]
    merged = []

    # ── Header: write section comment + merged CONFLIST ──
    conftypes = ", ".join(f"{lig_id}{l}" for l in all_labels)
    merged_conflist = f"CONFLIST, {lig_id}: {lig_id}BK, {conftypes}\n"
    # Preserve any file-level comments (before the first CONFLIST/section line)
    # and skip the per-section header comments that we re-emit below.
    section_comments = {
        "# conformer definition", "# atom name and bonds", "# atom charges",
        "# atom radius", "# conformer parameters",
    }
    for line in base["header"]:
        stripped = line.strip().lower()
        if stripped.startswith("conflist"):
            merged.append(merged_conflist)
        elif stripped in section_comments or stripped == "":
            pass  # drop — we emit our own section headers below
        else:
            merged.append(line)
    merged.append("\n")

    # ── CONNECT blocks from all states ──
    merged.append(f"# ATOM name and bonds\n")
    for label in all_labels:
        if parsed[label]["connect"]:
            merged.extend(parsed[label]["connect"])
            merged.append("\n")

    # ── CHARGE blocks from all states ──
    merged.append(f"# ATOM charges\n")
    for label in all_labels:
        if parsed[label]["charge"]:
            merged.extend(parsed[label]["charge"])
            merged.append("\n")

    # ── RADIUS blocks from all states ──
    merged.append(f"# Atom radius, dielectric boundary radius, VDW radius, and energy well depth\n")
    for label in all_labels:
        if parsed[label].get("radius"):
            merged.extend(parsed[label]["radius"])
            merged.append("\n")

    # ── CONFORMER lines from all states ──
    merged.append(f"# Conformer parameters that appear in head3.lst: ne, Em0, nH, pKa0, rxn\n")
    for label in all_labels:
        merged.extend(parsed[label]["conformer"])

    with open(output_path, "w") as f:
        f.writelines(merged)

    logging.info(f"  Merged ftpl: {output_path}")
    logging.info(f"    States: {', '.join(all_labels)}")
    return True


def _parse_ftpl_sections(ftpl_path, lig_id):
    """Parse a single-state ftpl into sections: header, connect, charge, radius, conformer."""
    header, connect, charge, radius, conformer = [], [], [], [], []
    past_header = False

    with open(ftpl_path) as f:
        for line in f:
            s = line.strip()
            if s.startswith("CONNECT"):
                connect.append(line)
                past_header = True
            elif s.startswith("CHARGE"):
                charge.append(line)
                past_header = True
            elif s.startswith("RADIUS"):
                radius.append(line)
                past_header = True
            elif s.startswith("CONFORMER"):
                conformer.append(line)
                past_header = True
            elif not past_header:
                header.append(line)
            # blank/comment lines between sections are intentionally dropped

    return {"header": header, "connect": connect, "charge": charge,
            "radius": radius, "conformer": conformer}


def fill_ftpl_charges_per_state(ftpl_path, per_state_charges, lig_id):
    """Fill CHARGE blocks in merged ftpl for each conformer state.

    Each conformer (e.g., EMH01, EMH+1) gets charges from its own PDB.
    'to_be_filled' values are replaced with actual charges.

    Returns count of unfilled entries.
    """
    from .charge_tools import to_mcce_atom_name

    with open(ftpl_path) as f:
        lines = f.readlines()

    pat = re.compile(
        r'^(CHARGE,\s*' + re.escape(lig_id) +
        r'(\S+),\s*"(.{4})":\s*)to_be_filled(.*)$'
    )
    filled = unfilled = 0
    new_lines = []

    for line in lines:
        m = pat.match(line)
        if m:
            prefix, conf_suffix, atom_name, suffix = m.groups()
            # conf_suffix is like "01", "+1", etc.
            # Find matching charges
            charges = per_state_charges.get(conf_suffix, {})
            matched = False

            # Try exact match
            if atom_name in charges:
                new_lines.append(f"{prefix}{charges[atom_name]:7.3f} # auto-filled{suffix}\n")
                filled += 1; matched = True

            # Try stripped match
            if not matched:
                stripped = atom_name.strip()
                mcce = to_mcce_atom_name(stripped)
                for cn, cv in charges.items():
                    if cn.strip() == stripped or to_mcce_atom_name(cn) == mcce:
                        new_lines.append(f"{prefix}{cv:7.3f} # auto-filled{suffix}\n")
                        filled += 1; matched = True; break

            if not matched:
                new_lines.append(line)
                unfilled += 1
        else:
            new_lines.append(line)

    with open(ftpl_path, "w") as f:
        f.writelines(new_lines)

    logging.info(f"  ✓ Filled {filled}, unfilled {unfilled} (per-state)")
    return unfilled


# ═════════════════════════════════════════════════════════════════════════════
# v3: PyMOL visualization script generator
# ═════════════════════════════════════════════════════════════════════════════

def generate_pymol_script(state_pdbs, h_diffs, lig_id, output_path):
    """Generate a PyMOL .pml script for side-by-side protonation state comparison.

    Shows all states side by side with:
      - H atoms visible (show sticks, h_add)
      - Added H highlighted in GREEN
      - Removed H sites highlighted in RED
      - Each state labeled with its conformer name

    Usage: pymol compare_states.pml

    Args:
        state_pdbs: {label: pdb_path}
        h_diffs: {label: {"added": [...], "removed": [...], ...}}
        lig_id: Ligand code
        output_path: Where to write the .pml file

    Returns:
        Path to the .pml file.
    """
    labels = sorted(state_pdbs.keys())
    spacing = 25.0  # Angstroms between states

    lines = [
        f"# PyMOL script: {lig_id} protonation state comparison",
        f"# Generated by mcce4_agent_ftpl v3",
        f"# Usage: pymol {os.path.basename(output_path)}",
        "",
        "# Settings for H atom visibility",
        "set valence, 1",
        "set h_bond_max_angle, 30",
        "set label_size, 14",
        "set label_color, white",
        "set stick_radius, 0.12",
        "set sphere_scale, 0.15",
        "bg_color white",
        "",
    ]

    for i, label in enumerate(labels):
        pdb_path = state_pdbs[label]
        obj_name = f"{lig_id}_{label.replace('+', 'p').replace('-', 'm')}"
        diff = h_diffs.get(label, {})
        h_added = diff.get("added", [])
        h_removed = diff.get("removed", [])

        lines.append(f"# ── State {label} ──")
        lines.append(f"load {os.path.abspath(pdb_path)}, {obj_name}")
        lines.append(f"show sticks, {obj_name}")
        lines.append(f"show spheres, {obj_name} and elem H")
        lines.append(f"set sphere_scale, 0.2, {obj_name} and elem H")
        lines.append(f"color gray80, {obj_name}")
        lines.append(f"color skyblue, {obj_name} and elem N")
        lines.append(f"color red, {obj_name} and elem O")
        lines.append(f"util.cbag {obj_name}")

        # Highlight added H in green
        for h_name in h_added:
            lines.append(f"color green, {obj_name} and name {h_name}")
            lines.append(f"set sphere_scale, 0.35, {obj_name} and name {h_name}")
            lines.append(f"label {obj_name} and name {h_name}, name")

        # Highlight parent atoms of removed H in red
        removed_from = diff.get("removed_from", {})
        for h_name, parent in removed_from.items():
            lines.append(f"color red, {obj_name} and name {parent}")

        # Show all H atom labels
        lines.append(f"label {obj_name} and elem H, name")

        # Position states side by side
        if i > 0:
            lines.append(f"translate [{spacing * i}, 0, 0], {obj_name}")

        # Add state label
        charge_sign = "+" if "p" in obj_name or "+" in label else ("-" if "m" in obj_name or "-" in label else "")
        lines.append(f'pseudoatom label_{obj_name}, label="{lig_id}{label}", '
                     f'pos=[{spacing * i}, -10, 0]')
        lines.append(f"set label_size, 20, label_{obj_name}")
        lines.append(f"set label_color, black, label_{obj_name}")
        lines.append("")

    lines.extend([
        "# Show all H atoms",
        "set h_bond_max_angle, 30",
        "show sticks, elem H",
        "",
        "# Zoom to fit",
        "orient",
        "zoom",
        "",
        f"# States: {', '.join(f'{lig_id}{l}' for l in labels)}",
        f"# Green H = added vs neutral, Red site = H removed vs neutral",
    ])

    with open(output_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    logging.info(f"  📄 PyMOL script: {output_path}")
    return output_path
