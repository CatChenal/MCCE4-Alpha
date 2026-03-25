#!/usr/bin/env python3
"""
mcce4_agent_ftpl.py — MCCE4 Topology File AI Agent
====================================================

An AI agent for creating MCCE4 topology files (.ftpl) from ligand PDB files.

Input:  A PDB file of the ligand molecule.
Output: A complete .ftpl topology file with charges and calibrated rxn values.

The agent can:
  - Automatically enumerate protonation states (via Dimorphite-DL + LLM)
  - Accept user-supplied PDB files for specific protonation states
  - Display a GUI for reviewing/editing proposed states before proceeding
  - Validate hybridizations (via RDKit)
  - Generate and fill atomic charges
  - Calibrate rxn02, rxn04, rxn08 via MCCE4 steps 1-3

Usage:
  mcce4_agent_ftpl.py EMH.pdb                                   # Auto — agent decides states
  mcce4_agent_ftpl.py EMH.pdb --state-pdbs EMH_01.pdb EMH_+1.pdb  # User provides state PDBs
  mcce4_agent_ftpl.py EMH.pdb --gui                              # GUI to review states
  mcce4_agent_ftpl.py EMH.pdb --no-llm                           # No LLM, Dimorphite-DL only
  mcce4_agent_ftpl.py EMH.pdb --interactive                      # Terminal confirmation mode
  mcce4_agent_ftpl.py EMH.pdb --dry-run                          # Generate ftpl, skip calibration

Requirements:
  pip install dimorphite-dl rdkit google-genai PyQt5
  export GEMINI_API_KEY="your_free_key"    # from https://ai.google.dev (optional)

Author:  Gehan / MCCE4 Team (GunnerLab)
Version: 2.1.0
"""

import argparse
import os
import sys
import re
import subprocess
import shutil
import logging
import json
import textwrap
from pathlib import Path
from datetime import datetime
from typing import Optional

# ──────────────────────────────────────────────────────────────────────────────
# Constants
# ──────────────────────────────────────────────────────────────────────────────
RCSB_SMILES_URL = "https://data.rcsb.org/rest/v1/core/chemcomp/{lig_id}"
LIGAND_EXPO_URL = "http://ligand-expo.rcsb.org/reports/{first_char}/{lig_id}/{lig_id}_ideal.pdb"

SUPPORTED_CHARGE_METHODS = [
    "mmff94", "am1bcc", "am1bccelf10", "am1bccnosymspt",
    "amber", "amberff94", "antechamber", "file",
]
DEFAULT_CHARGE_METHOD = "mmff94"
DIELECTRIC_MAP = {2: "rxn02", 4: "rxn04", 8: "rxn08"}

CHARGE_TO_CONF = {0: "01", 1: "+1", -1: "-1", 2: "+2", -2: "-2"}


# ──────────────────────────────────────────────────────────────────────────────
# Logging
# ──────────────────────────────────────────────────────────────────────────────
def setup_logging(log_file: str = "mcce4_agent_ftpl.log", verbose: bool = False):
    """Configure logging to file and console."""
    log_level = logging.DEBUG if verbose else logging.INFO
    file_fmt = logging.Formatter("%(asctime)s | %(levelname)-8s | %(message)s",
                                  datefmt="%Y-%m-%d %H:%M:%S")
    console_fmt = logging.Formatter("%(levelname)-8s | %(message)s")
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler(log_file, mode="w")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(file_fmt)
    logger.addHandler(fh)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(log_level)
    ch.setFormatter(console_fmt)
    logger.addHandler(ch)
    return logger


# ──────────────────────────────────────────────────────────────────────────────
# Utility
# ──────────────────────────────────────────────────────────────────────────────
def run_cmd(cmd, description="", capture=False, cwd=None):
    """Run a shell command with logging and error handling."""
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


def to_mcce_atom_name(name: str) -> str:
    """Convert an atom name to MCCE's 4-character padded format."""
    name = name.strip()
    if len(name) >= 4:
        return name[:4]
    elif len(name) == 3:
        return f" {name}"
    elif len(name) == 2:
        return f" {name} "
    elif len(name) == 1:
        return f" {name}  "
    return name.ljust(4)


def extract_lig_id_from_pdb(pdb_path: str) -> str:
    """Extract the 3-letter ligand/residue code from a PDB file.

    Looks at HETATM records first, then ATOM records.
    Returns the most common residue name found.
    """
    res_counts = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("HETATM", "ATOM  ")):
                resname = line[17:20].strip()
                if resname and resname not in ("HOH", "WAT", "TIP"):
                    res_counts[resname] = res_counts.get(resname, 0) + 1

    if not res_counts:
        # Fallback: use filename
        return Path(pdb_path).stem.upper()[:3]

    return max(res_counts, key=res_counts.get)


def parse_state_pdb_label(filename: str, lig_id: str) -> Optional[str]:
    """Extract conformer state label from a filename like EMH_01.pdb or EMH_+1.pdb.

    Supported patterns:
      EMH_01.pdb  -> "01"
      EMH_+1.pdb  -> "+1"
      EMH_-1.pdb  -> "-1"
      EMH01.pdb   -> "01"
    """
    stem = Path(filename).stem.upper()
    lig_upper = lig_id.upper()

    # Try EMH_01 or EMH_+1 pattern
    pattern = re.compile(
        re.escape(lig_upper) + r'[_\-]?([0+-]\d*|[+-]\d+)$', re.IGNORECASE
    )
    m = pattern.match(stem)
    if m:
        return m.group(1)

    # Try just the suffix after removing lig_id
    if stem.startswith(lig_upper):
        suffix = stem[len(lig_upper):].lstrip("_-")
        if suffix:
            return suffix

    return None


# ──────────────────────────────────────────────────────────────────────────────
# Conformer State — data class
# ──────────────────────────────────────────────────────────────────────────────
class ConformerState:
    """Represents a single protonation/redox conformer state."""

    def __init__(self, label: str, charge: int = 0, nH: int = 0,
                 pka: Optional[float] = None, smiles: str = "",
                 pdb_path: Optional[str] = None, source: str = "auto",
                 site_atom: Optional[str] = None,
                 rationale: str = ""):
        self.label = label          # e.g., "01", "+1", "-1"
        self.charge = charge        # formal charge
        self.nH = nH                # protons relative to neutral
        self.pka = pka              # solution pKa (if known)
        self.smiles = smiles        # SMILES for this state
        self.pdb_path = pdb_path    # PDB file for this state (if user-supplied)
        self.source = source        # "user", "dimorphite", "llm", "auto"
        self.site_atom = site_atom  # atom where protonation occurs
        self.rationale = rationale  # why this state was included

    def to_dict(self):
        return {
            "label": self.label, "charge": self.charge, "nH": self.nH,
            "pka": self.pka, "smiles": self.smiles, "pdb_path": self.pdb_path,
            "source": self.source, "site_atom": self.site_atom,
            "rationale": self.rationale,
        }

    def __repr__(self):
        return (f"ConformerState({self.label}, charge={self.charge:+d}, "
                f"nH={self.nH:+d}, source={self.source})")


# ──────────────────────────────────────────────────────────────────────────────
# PHASE 1: Molecule Intelligence
# ──────────────────────────────────────────────────────────────────────────────
class MoleculeIntelligence:
    """Gathers chemical information about the target ligand."""

    def __init__(self, lig_id: str, pdb_path: str):
        self.lig_id = lig_id.upper()
        self.pdb_path = pdb_path
        self.smiles = None
        self.name = None
        self.formula = None
        self.formal_charge = 0

    def fetch_from_rcsb(self) -> bool:
        """Fetch SMILES and metadata from RCSB chemical component dictionary."""
        logging.info(f"🔍 Fetching chemical info for '{self.lig_id}' from RCSB...")
        url = RCSB_SMILES_URL.format(lig_id=self.lig_id)
        result = run_cmd(f'curl -s "{url}"', description="Querying RCSB API", capture=True)

        if result is None or not result.stdout.strip():
            logging.warning(f"  Could not fetch from RCSB — will try RDKit")
            return False
        try:
            data = json.loads(result.stdout)
            chem = data.get("chem_comp", {})
            self.name = chem.get("name", "Unknown")
            self.formula = chem.get("formula", "Unknown")
            self.formal_charge = chem.get("pdbx_formal_charge", 0) or 0

            # Get SMILES from descriptors
            descriptors = data.get("rcsb_chem_comp_descriptor", [])
            if isinstance(descriptors, list):
                for desc in descriptors:
                    if isinstance(desc, dict) and "SMILES" in desc.get("type", ""):
                        self.smiles = desc.get("descriptor")
                        if "CANONICAL" in desc.get("type", ""):
                            break

            logging.info(f"  ✓ Name:    {self.name}")
            logging.info(f"  ✓ Formula: {self.formula}")
            logging.info(f"  ✓ Charge:  {self.formal_charge}")
            if self.smiles:
                logging.info(f"  ✓ SMILES:  {self.smiles}")
            return True
        except (json.JSONDecodeError, KeyError) as e:
            logging.warning(f"  RCSB parse error: {e}")
            return False

    def get_smiles_from_pdb(self) -> Optional[str]:
        """Extract SMILES from PDB via RDKit."""
        try:
            from rdkit import Chem
            mol = Chem.MolFromPDBFile(self.pdb_path, removeHs=False)
            if mol:
                self.smiles = Chem.MolToSmiles(mol)
                logging.info(f"  ✓ SMILES from PDB: {self.smiles}")
                return self.smiles
        except ImportError:
            logging.warning("  RDKit not available for SMILES extraction")
        except Exception as e:
            logging.warning(f"  RDKit error: {e}")
        return None


# ──────────────────────────────────────────────────────────────────────────────
# PHASE 2: Protonation State Enumeration
# ──────────────────────────────────────────────────────────────────────────────
class ProtonationAnalyzer:
    """Enumerate protonation states using Dimorphite-DL."""

    def __init__(self, smiles: str, ph: float = 7.4, ph_tolerance: float = 1.0):
        self.smiles = smiles
        self.ph = ph
        self.ph_tolerance = ph_tolerance

    def enumerate_states(self) -> list:
        """Return list of ConformerState from Dimorphite-DL enumeration."""
        logging.info(f"🧪 Enumerating protonation states at pH {self.ph} "
                     f"(± {self.ph_tolerance})...")

        state_smiles = [self.smiles]  # fallback

        try:
            from dimorphite_dl.protonate import protonate_smiles
            results = protonate_smiles(
                self.smiles,
                min_ph=self.ph - self.ph_tolerance,
                max_ph=self.ph + self.ph_tolerance,
            )
            if results:
                state_smiles = list(set(results))
            logging.info(f"  ✓ Dimorphite-DL found {len(state_smiles)} state(s)")
        except ImportError:
            logging.warning("  Dimorphite-DL not installed — using input SMILES only")
        except Exception as e:
            logging.warning(f"  Dimorphite-DL error: {e}")

        # Compute formal charges via RDKit
        states = []
        seen_charges = set()
        try:
            from rdkit import Chem
            for smi in state_smiles:
                mol = Chem.MolFromSmiles(smi)
                charge = Chem.GetFormalCharge(mol) if mol else 0
                if charge in seen_charges:
                    continue
                seen_charges.add(charge)
                label = CHARGE_TO_CONF.get(charge, f"{charge:+d}"[-2:])
                nH = charge  # relative to neutral
                states.append(ConformerState(
                    label=label, charge=charge, nH=nH if charge != 0 else 0,
                    smiles=smi, source="dimorphite",
                    rationale=f"Dimorphite-DL at pH {self.ph}"
                ))
        except ImportError:
            # No RDKit — basic fallback
            states.append(ConformerState(
                label="01", charge=0, nH=0, smiles=self.smiles,
                source="fallback", rationale="No RDKit, assumed neutral"
            ))

        for s in states:
            logging.info(f"  📋 {s.label}: charge={s.charge:+d}, nH={s.nH:+d}")

        return states


# ──────────────────────────────────────────────────────────────────────────────
# PHASE 3: LLM Reasoning (Gemini Free API)
# ──────────────────────────────────────────────────────────────────────────────
class AgentBrain:
    """LLM-powered chemistry reasoning using Google Gemini free tier."""

    def __init__(self, api_key: Optional[str] = None, enabled: bool = True):
        self.available = False
        self.client = None
        self.model_id = "gemini-2.5-flash"

        if not enabled:
            logging.info("🧠 Agent brain: Rule-based mode (--no-llm)")
            return

        api_key = api_key or os.environ.get("GEMINI_API_KEY") or os.environ.get("GOOGLE_API_KEY")
        if not api_key:
            logging.info("🧠 Agent brain: Rule-based mode (no GEMINI_API_KEY)")
            return

        try:
            from google import genai
            self.client = genai.Client(api_key=api_key)
            self.available = True
            logging.info(f"🧠 Agent brain: {self.model_id} (free tier)")
        except ImportError:
            logging.warning("  google-genai not installed — run: pip install google-genai")
        except Exception as e:
            logging.warning(f"  Gemini init failed: {e}")

    def _ask(self, prompt: str) -> Optional[str]:
        if not self.available or not self.client:
            return None
        try:
            response = self.client.models.generate_content(
                model=self.model_id,
                contents=prompt,
            )
            return response.text
        except Exception as e:
            logging.warning(f"  Gemini API error: {e}")
            return None

    def refine_states(self, lig_name, smiles, formula,
                       states: list) -> list:
        """Ask LLM to review and refine conformer states.

        May adjust pKa values, add warnings, suggest additional states.
        Returns refined list of ConformerState.
        """
        states_desc = "\n".join(
            f"  {s.label}: charge={s.charge:+d}, nH={s.nH:+d}, SMILES={s.smiles}"
            for s in states
        )

        prompt = f"""You are a computational chemistry expert for MCCE simulations.

Analyze this ligand for MCCE topology file creation:
  Name: {lig_name}
  Formula: {formula}
  SMILES: {smiles}

Proposed protonation states:
{states_desc}

For each state, provide:
1. Is it chemically reasonable at pH 7.4?
2. Estimated solution pKa for the protonation/deprotonation
3. Which atom is the protonation site?
4. Any missing states?

Respond ONLY in JSON (no markdown fences):
{{
  "states": [
    {{"label": "01", "charge": 0, "nH": 0, "pka": null, "site_atom": null, "rationale": "neutral"}}
  ],
  "warnings": [],
  "summary": "brief analysis"
}}"""

        response = self._ask(prompt)
        if response:
            try:
                clean = re.sub(r'^```json\s*', '', response.strip())
                clean = re.sub(r'\s*```$', '', clean)
                data = json.loads(clean)

                refined = []
                for s in data.get("states", []):
                    refined.append(ConformerState(
                        label=s["label"], charge=s.get("charge", 0),
                        nH=s.get("nH", 0), pka=s.get("pka"),
                        site_atom=s.get("site_atom"),
                        source="llm", rationale=s.get("rationale", "")
                    ))

                if data.get("warnings"):
                    for w in data["warnings"]:
                        logging.warning(f"  ⚠ LLM: {w}")

                logging.info(f"  🧠 {data.get('summary', 'Analysis complete')}")
                return refined if refined else states

            except (json.JSONDecodeError, KeyError) as e:
                logging.warning(f"  Could not parse LLM response: {e}")

        return states  # Return original if LLM fails

    def validate_hybridizations(self, ftpl_content, smiles, lig_id) -> list:
        """Ask LLM to review CONNECT hybridization assignments."""
        connect_lines = [l for l in ftpl_content.splitlines()
                         if l.strip().startswith("CONNECT")][:30]
        if not connect_lines:
            return []

        prompt = f"""Review these MCCE CONNECT records for {lig_id} (SMILES: {smiles}).
Orbital types: s, sp, sp2, sp3, d2sp3

{chr(10).join(connect_lines)}

Are hybridizations correct? Respond with JSON list of corrections only:
[{{"atom": "C7", "current": "sp3", "correct": "sp2", "reason": "aromatic carbon"}}]
If all correct: []"""

        response = self._ask(prompt)
        if response:
            try:
                clean = re.sub(r'^```json\s*', '', response.strip())
                clean = re.sub(r'\s*```$', '', clean)
                corrections = json.loads(clean)
                if corrections:
                    logging.info(f"  🧠 {len(corrections)} hybridization correction(s)")
                return corrections
            except (json.JSONDecodeError, KeyError):
                pass
        return []


# ──────────────────────────────────────────────────────────────────────────────
# PHASE 4: GUI Review Interface
# ──────────────────────────────────────────────────────────────────────────────
def launch_gui_review(lig_id: str, pdb_path: str,
                       states: list) -> list:
    """Launch PyQt5 GUI for reviewing and editing conformer states.

    Shows:
      - 2D molecule depiction (via RDKit if available)
      - Table of proposed conformer states with label, charge, nH, pKa, source
      - Buttons: Approve, Edit, Remove, Add State
      - File picker to assign/change PDB files per state

    Returns:
        Updated list of ConformerState (user-approved).
    """
    try:
        from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout,
                                      QHBoxLayout, QTableWidget, QTableWidgetItem,
                                      QPushButton, QLabel, QFileDialog, QHeaderView,
                                      QMessageBox, QGroupBox, QInputDialog, QComboBox,
                                      QAbstractItemView, QFrame)
        from PyQt5.QtGui import QPixmap, QImage, QFont, QColor
        from PyQt5.QtCore import Qt
    except ImportError:
        logging.warning("  PyQt5 not available — falling back to terminal mode")
        return _terminal_review(lig_id, states)

    approved_states = list(states)  # Will be modified by GUI
    app = QApplication.instance() or QApplication(sys.argv)

    # Try to generate 2D depiction via RDKit (QApplication must exist first)
    mol_pixmap = None
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw
        import io

        mol = Chem.MolFromPDBFile(pdb_path, removeHs=True)
        if mol is None:
            mol = Chem.MolFromPDBFile(pdb_path, removeHs=True, sanitize=False)

        if mol:
            img = Draw.MolToImage(mol, size=(400, 300))
            buf = io.BytesIO()
            img.save(buf, format="PNG")
            buf.seek(0)
            qimg = QImage.fromData(buf.read())
            mol_pixmap = QPixmap.fromImage(qimg)
    except Exception as e:
        logging.debug(f"  RDKit 2D depiction failed: {e}")

    class ReviewWindow(QMainWindow):
        def __init__(self):
            super().__init__()
            self.setWindowTitle(f"🤖 MCCE4 Agent — Review Conformer States for {lig_id}")
            self.setMinimumSize(800, 600)
            self.approved = False

            central = QWidget()
            self.setCentralWidget(central)
            layout = QVBoxLayout(central)

            # ── Header ──
            header = QLabel(f"<h2>🧬 {lig_id} — Conformer State Review</h2>"
                           f"<p>Input PDB: <b>{pdb_path}</b></p>")
            header.setAlignment(Qt.AlignCenter)
            layout.addWidget(header)

            # ── Molecule + Table side by side ──
            content = QHBoxLayout()

            # Left: molecule depiction
            mol_group = QGroupBox("Molecule Structure")
            mol_layout = QVBoxLayout(mol_group)
            mol_label = QLabel()
            if mol_pixmap:
                mol_label.setPixmap(mol_pixmap.scaled(380, 280, Qt.KeepAspectRatio,
                                                       Qt.SmoothTransformation))
            else:
                mol_label.setText("(RDKit not available for 2D depiction)\n\n"
                                  "Install rdkit for molecule visualization")
                mol_label.setAlignment(Qt.AlignCenter)
            mol_layout.addWidget(mol_label)
            content.addWidget(mol_group)

            # Right: states table
            table_group = QGroupBox("Proposed Conformer States")
            table_layout = QVBoxLayout(table_group)

            self.table = QTableWidget()
            self.table.setColumnCount(7)
            self.table.setHorizontalHeaderLabels(
                ["Label", "Charge", "nH", "pKa", "Source", "PDB File", "Rationale"]
            )
            self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
            self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
            self._populate_table()
            table_layout.addWidget(self.table)

            # Table buttons
            btn_row = QHBoxLayout()
            btn_add = QPushButton("➕ Add State")
            btn_add.clicked.connect(self._add_state)
            btn_row.addWidget(btn_add)

            btn_remove = QPushButton("❌ Remove Selected")
            btn_remove.clicked.connect(self._remove_state)
            btn_row.addWidget(btn_remove)

            btn_pdb = QPushButton("📂 Assign PDB")
            btn_pdb.clicked.connect(self._assign_pdb)
            btn_row.addWidget(btn_pdb)

            table_layout.addLayout(btn_row)
            content.addWidget(table_group)
            layout.addLayout(content)

            # ── Bottom: Approve / Cancel ──
            sep = QFrame()
            sep.setFrameShape(QFrame.HLine)
            layout.addWidget(sep)

            bottom = QHBoxLayout()
            btn_cancel = QPushButton("✖ Cancel")
            btn_cancel.setStyleSheet("background-color: #d9534f; color: white; "
                                      "padding: 8px 20px; font-size: 14px;")
            btn_cancel.clicked.connect(self._cancel)
            bottom.addWidget(btn_cancel)

            bottom.addStretch()

            btn_approve = QPushButton("✔ Approve && Proceed")
            btn_approve.setStyleSheet("background-color: #5cb85c; color: white; "
                                       "padding: 8px 20px; font-size: 14px; font-weight: bold;")
            btn_approve.clicked.connect(self._approve)
            bottom.addWidget(btn_approve)

            layout.addLayout(bottom)

        def _populate_table(self):
            self.table.setRowCount(len(approved_states))
            for i, s in enumerate(approved_states):
                self.table.setItem(i, 0, QTableWidgetItem(s.label))
                self.table.setItem(i, 1, QTableWidgetItem(f"{s.charge:+d}"))
                self.table.setItem(i, 2, QTableWidgetItem(str(s.nH)))
                self.table.setItem(i, 3, QTableWidgetItem(
                    f"{s.pka:.1f}" if s.pka is not None else "—"))
                self.table.setItem(i, 4, QTableWidgetItem(s.source))
                self.table.setItem(i, 5, QTableWidgetItem(
                    os.path.basename(s.pdb_path) if s.pdb_path else "—"))
                self.table.setItem(i, 6, QTableWidgetItem(s.rationale))

                # Color-code by source
                color = {"user": QColor(200, 230, 255),
                         "dimorphite": QColor(220, 255, 220),
                         "llm": QColor(255, 240, 200)}.get(s.source, QColor(245, 245, 245))
                for col in range(7):
                    item = self.table.item(i, col)
                    if item:
                        item.setBackground(color)

        def _add_state(self):
            label, ok = QInputDialog.getText(self, "Add State",
                                              "Conformer label (e.g., 01, +1, -1):")
            if ok and label:
                charge, ok2 = QInputDialog.getInt(self, "Charge",
                                                    "Formal charge:", 0, -5, 5)
                if ok2:
                    nH = charge  # default assumption
                    new_state = ConformerState(
                        label=label, charge=charge, nH=nH,
                        source="user", rationale="Manually added by user"
                    )
                    # Optionally assign a PDB right away
                    reply = QMessageBox.question(
                        self, "Assign PDB",
                        f"Assign a PDB file for state '{label}'?",
                        QMessageBox.Yes | QMessageBox.No
                    )
                    if reply == QMessageBox.Yes:
                        fpath, _ = QFileDialog.getOpenFileName(
                            self, f"PDB for {label}", "", "PDB files (*.pdb);;All (*)"
                        )
                        if fpath:
                            new_state.pdb_path = fpath

                    approved_states.append(new_state)
                    self._populate_table()

        def _remove_state(self):
            rows = set(idx.row() for idx in self.table.selectedIndexes())
            if not rows:
                QMessageBox.information(self, "Remove", "Select a row first.")
                return
            for row in sorted(rows, reverse=True):
                if row < len(approved_states):
                    removed = approved_states.pop(row)
                    logging.info(f"  ✖ Removed state: {removed.label}")
            self._populate_table()

        def _assign_pdb(self):
            rows = set(idx.row() for idx in self.table.selectedIndexes())
            if not rows:
                QMessageBox.information(self, "Assign PDB", "Select a row first.")
                return
            row = min(rows)
            fpath, _ = QFileDialog.getOpenFileName(
                self, f"PDB for {approved_states[row].label}", "",
                "PDB files (*.pdb);;All (*)"
            )
            if fpath:
                approved_states[row].pdb_path = fpath
                logging.info(f"  📂 Assigned {fpath} to state {approved_states[row].label}")
                self._populate_table()

        def _approve(self):
            # Read back any edited cells
            for i in range(self.table.rowCount()):
                if i < len(approved_states):
                    label_item = self.table.item(i, 0)
                    if label_item:
                        approved_states[i].label = label_item.text()
                    pka_item = self.table.item(i, 3)
                    if pka_item and pka_item.text() != "—":
                        try:
                            approved_states[i].pka = float(pka_item.text())
                        except ValueError:
                            pass
            self.approved = True
            self.close()

        def _cancel(self):
            self.approved = False
            self.close()

    window = ReviewWindow()
    window.show()
    app.exec_()

    if not window.approved:
        logging.info("  ✖ User cancelled — aborting.")
        sys.exit(0)

    logging.info(f"  ✅ User approved {len(approved_states)} state(s):")
    for s in approved_states:
        pdb_info = f", PDB={os.path.basename(s.pdb_path)}" if s.pdb_path else ""
        logging.info(f"     {s.label}: charge={s.charge:+d}{pdb_info}")

    return approved_states


def _terminal_review(lig_id: str, states: list) -> list:
    """Fallback terminal-based review when GUI is unavailable."""
    print(f"\n  === Conformer States for {lig_id} ===")
    for i, s in enumerate(states):
        pdb_info = f" (PDB: {os.path.basename(s.pdb_path)})" if s.pdb_path else ""
        print(f"  [{i+1}] {s.label}: charge={s.charge:+d}, nH={s.nH:+d}, "
              f"source={s.source}{pdb_info}")

    print(f"\n  Options: [A]pprove / [E]dit labels / [C]ancel")
    choice = input("  > ").strip().lower()

    if choice == "c":
        print("  Aborted.")
        sys.exit(0)
    elif choice == "e":
        new_labels = input("  Enter conformer labels (space-separated, e.g., 01 +1 -1): ")
        labels = new_labels.split()
        states = [ConformerState(label=l, source="user", rationale="User-specified")
                  for l in labels]

    return states


# ──────────────────────────────────────────────────────────────────────────────
# PHASE 5: Charge Generation
# ──────────────────────────────────────────────────────────────────────────────
def generate_charges(pdb_path: str, method: str, lig_id: str) -> dict:
    """Generate charges using specified method. Returns MCCE atom name -> charge."""
    if method == "antechamber":
        return _charges_antechamber(pdb_path, lig_id)
    return _charges_openeye(pdb_path, method, lig_id)


def _charges_openeye(pdb_path, method, lig_id):
    logging.info(f"⚡ Generating charges via OpenEye (method: {method})")
    result = run_cmd(
        f"oe_assigncharges_QuacpakTK.py -method {method} -in {pdb_path}",
        description=f"OpenEye {method}", capture=True
    )
    if result is None:
        return {}
    charges = {}
    pattern = re.compile(
        r"Atom Name:\s+(\S+)\s+\|\s+Symbol:\s+(\S+)\s+\|\s+Charge:\s+([-\d.]+)")
    atom_counts = {}
    raw = []
    for line in result.stdout.splitlines():
        m = pattern.search(line)
        if m:
            name, sym, charge = m.group(1), m.group(2), float(m.group(3))
            raw.append((name, charge))
            atom_counts[name] = atom_counts.get(name, 0) + 1
    water_atoms = {n for n, c in atom_counts.items() if c > 5}
    seen = set()
    for name, charge in raw:
        if name in water_atoms:
            continue
        mcce_name = to_mcce_atom_name(name)
        if mcce_name not in seen:
            charges[mcce_name] = charge
            seen.add(mcce_name)
    non_zero = sum(1 for v in charges.values() if abs(v) > 1e-6)
    if non_zero == 0 and method != "mmff94":
        logging.warning(f"  ⚠ All zeros — falling back to mmff94")
        return _charges_openeye(pdb_path, "mmff94", lig_id)
    logging.info(f"  ✓ {len(charges)} charges (total: {sum(charges.values()):.3f})")
    return charges


def _charges_antechamber(pdb_path, lig_id, nc=0):
    logging.info(f"⚡ Generating charges via antechamber (nc={nc})")
    mol2 = f"{lig_id}_antechamber.mol2"
    result = run_cmd(
        f"antechamber -i {pdb_path} -fi pdb -o {mol2} -fo mol2 -c bcc -s 2 -nc {nc} -at gaff",
        description="antechamber AM1-BCC", capture=True)
    if result is None or not os.path.exists(mol2):
        return {}
    charges = {}
    in_atoms = False
    with open(mol2) as f:
        for line in f:
            l = line.strip()
            if l.startswith("@<TRIPOS>ATOM"):
                in_atoms = True; continue
            elif l.startswith("@<TRIPOS>"):
                in_atoms = False; continue
            if in_atoms and l:
                p = l.split()
                if len(p) >= 9:
                    charges[to_mcce_atom_name(p[1])] = float(p[8])
    logging.info(f"  ✓ {len(charges)} charges (total: {sum(charges.values()):.3f})")
    return charges


# ──────────────────────────────────────────────────────────────────────────────
# PHASE 6: FTPL Assembly & Charge Filling
# ──────────────────────────────────────────────────────────────────────────────
def generate_ftpl_template(lig_id, pdb_path, conformer_labels, ftpl_path):
    conf_str = " ".join(conformer_labels)
    logging.info(f"📝 Generating .ftpl: conformers = {conf_str}")
    result = run_cmd(f"pdb2ftpl.py -p {pdb_path} -c {conf_str}",
                     description="pdb2ftpl.py", capture=True)
    if result is None:
        logging.error("  pdb2ftpl.py failed!"); sys.exit(1)
    with open(ftpl_path, "w") as f:
        f.write(result.stdout)
    unfilled = result.stdout.count("to_be_filled")
    logging.info(f"  ✓ Template: {ftpl_path} ({unfilled} charges to fill)")
    return unfilled


def fill_ftpl_charges(ftpl_path, charges, lig_id):
    logging.info(f"✏️  Filling charges in {ftpl_path}")
    with open(ftpl_path) as f:
        lines = f.readlines()
    pat = re.compile(
        r'^(CHARGE,\s*' + re.escape(lig_id) + r'(\S+),\s*"(.{4})":\s*)to_be_filled(.*)$')
    filled = unfilled = 0
    new_lines = []
    for line in lines:
        m = pat.match(line)
        if m:
            prefix, conf_id, atom_name, suffix = m.groups()
            matched = False
            if atom_name in charges:
                new_lines.append(f"{prefix}{charges[atom_name]:7.3f} # auto-filled{suffix}\n")
                filled += 1; matched = True
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


def update_conformer_params(ftpl_path, states, lig_id):
    logging.info("📋 Updating CONFORMER parameters...")
    with open(ftpl_path) as f:
        content = f.read()
    for s in states:
        ct = f"{lig_id}{s.label}"
        nH = s.nH if s.nH else 0
        pka = s.pka if s.pka else 0.0
        # Update nH
        p = re.compile(rf"(CONFORMER,\s*{re.escape(ct)}:.*?nH=)\s*(\d+)")
        if p.search(content):
            content = p.sub(rf"\g<1>{nH}", content)
        # Update pKa0
        p2 = re.compile(rf"(CONFORMER,\s*{re.escape(ct)}:.*?pKa0=)\s*([-\d.]+)")
        if p2.search(content):
            content = p2.sub(rf"\g<1>{pka:.2f}", content)
        logging.info(f"  ✓ {ct}: nH={nH}, pKa0={pka:.2f}")
    with open(ftpl_path, "w") as f:
        f.write(content)


# ──────────────────────────────────────────────────────────────────────────────
# PHASE 7: RXN Calibration
# ──────────────────────────────────────────────────────────────────────────────
def run_rxn_calibration(lig_id, ftpl_path, pdb_path, conformer_labels,
                         dielectrics, work_dir):
    user_param = os.path.join(work_dir, "user_param")
    os.makedirs(user_param, exist_ok=True)
    link = os.path.join(user_param, os.path.basename(ftpl_path))
    if os.path.exists(link):
        os.remove(link)
    os.symlink(os.path.abspath(ftpl_path), link)

    pdb_base = os.path.basename(pdb_path)
    pdb_in_wd = os.path.join(work_dir, pdb_base)
    if not os.path.exists(pdb_in_wd):
        shutil.copy2(os.path.abspath(pdb_path), pdb_in_wd)

    # Steps 1-2 once
    if run_cmd(f"step1.py {pdb_base}", description="MCCE4 Step 1", cwd=work_dir) is None:
        logging.error("  Step 1 failed"); return
    if run_cmd("step2.py", description="MCCE4 Step 2", cwd=work_dir) is None:
        logging.error("  Step 2 failed"); return

    for eps in dielectrics:
        rxn_key = DIELECTRIC_MAP.get(eps)
        if not rxn_key:
            continue
        logging.info(f"\n{'='*60}")
        logging.info(f"  🔬 RXN Calibration ε={eps} ({rxn_key})")
        logging.info(f"{'='*60}")

        run_cmd(f"step3.py -d {eps}", description=f"Step 3 (ε={eps})", cwd=work_dir)
        head3 = os.path.join(work_dir, "head3.lst")
        dsolv = _parse_dsolv(head3, lig_id, conformer_labels)
        if not dsolv:
            logging.error("  Could not parse dsolv"); continue
        for ct, v in dsolv.items():
            logging.info(f"  📊 {ct}: dsolv = {v:.3f}")
        _update_rxn(ftpl_path, dsolv, rxn_key)

        logging.info("  🔄 Validating...")
        run_cmd(f"step3.py -d {eps}", description=f"Validation (ε={eps})", cwd=work_dir)
        check = _parse_dsolv(head3, lig_id, conformer_labels)
        ok = all(abs(v) <= 0.01 for v in check.values())
        for ct, v in check.items():
            sym = "✓" if abs(v) <= 0.01 else "⚠"
            logging.info(f"  {sym} {ct}: dsolv = {v:.3f}")
        if ok:
            logging.info(f"  🎉 {rxn_key} calibration successful!")


def _parse_dsolv(head3_path, lig_id, labels):
    if not os.path.exists(head3_path):
        return {}
    with open(head3_path) as f:
        lines = f.readlines()
    types = [f"{lig_id}{l}" for l in labels]
    header = next((l for l in lines if "dsolv" in l.lower() and "iConf" in l), None)
    if not header:
        return {}
    try:
        col = header.split().index("dsolv")
    except ValueError:
        return {}
    vals = {t: 0.0 for t in types}
    for line in lines:
        for t in types:
            if t in line and not line.strip().startswith("iConf"):
                parts = line.split()
                if len(parts) > col:
                    try:
                        v = float(parts[col])
                        if v < vals[t]:
                            vals[t] = v
                    except (ValueError, IndexError):
                        pass
    return vals


def _update_rxn(ftpl_path, dsolv_values, rxn_key):
    with open(ftpl_path) as f:
        content = f.read()
    for ct, dsolv in dsolv_values.items():
        p = re.compile(rf"(CONFORMER,\s*{re.escape(ct)}:.*?{rxn_key}=\s*)([-\d.]+)")
        if p.search(content):
            content = p.sub(rf"\g<1>{dsolv:8.3f}", content)
            logging.info(f"    {ct}: {rxn_key} = {dsolv:.3f}")
    with open(ftpl_path, "w") as f:
        f.write(content)


# ──────────────────────────────────────────────────────────────────────────────
# MAIN AGENT LOOP
# ──────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="🤖 MCCE4 Topology File AI Agent",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
        Examples:
          %(prog)s EMH.pdb                                          # Full auto
          %(prog)s EMH.pdb --state-pdbs EMH_01.pdb EMH_+1.pdb      # User state PDBs
          %(prog)s EMH.pdb --gui                                    # GUI review
          %(prog)s EMH.pdb --interactive                            # Terminal review
          %(prog)s EMH.pdb --no-llm                                 # No LLM
          %(prog)s EMH.pdb --charge-method am1bcc                   # Override charges
          %(prog)s EMH.pdb --dry-run                                # Skip calibration
        """))

    parser.add_argument("pdb", help="Ligand PDB file (e.g., EMH.pdb)")
    parser.add_argument("--state-pdbs", nargs="+", default=None,
                        help="PDB files for specific states (e.g., EMH_01.pdb EMH_+1.pdb)")
    parser.add_argument("--gui", action="store_true",
                        help="Launch GUI to review/edit conformer states")
    parser.add_argument("--interactive", action="store_true",
                        help="Terminal-based review before proceeding")
    parser.add_argument("--ph", type=float, default=7.4, help="Target pH (default: 7.4)")
    parser.add_argument("--ph-tolerance", type=float, default=1.0,
                        help="pH tolerance (default: ±1.0)")
    parser.add_argument("--charge-method", default=DEFAULT_CHARGE_METHOD,
                        choices=SUPPORTED_CHARGE_METHODS, help="Charge method")
    parser.add_argument("-d", "--dielectric", nargs="+", type=int, default=[2, 4, 8],
                        help="Dielectric constants (default: 2 4 8)")
    parser.add_argument("--no-llm", action="store_true", help="Disable LLM reasoning")
    parser.add_argument("--dry-run", action="store_true", help="Skip RXN calibration")
    parser.add_argument("--work-dir", default=".", help="Working directory")
    parser.add_argument("-o", "--output", default=None, help="Output .ftpl filename")
    parser.add_argument("-v", "--verbose", action="store_true")

    args = parser.parse_args()

    # Validate input PDB
    if not os.path.exists(args.pdb):
        print(f"ERROR: PDB file not found: {args.pdb}")
        sys.exit(1)

    pdb_path = os.path.abspath(args.pdb)
    lig_id = extract_lig_id_from_pdb(pdb_path)
    ftpl_path = args.output or f"{lig_id}.ftpl"
    work_dir = os.path.abspath(args.work_dir)

    setup_logging(f"mcce4_agent_ftpl_{lig_id}.log", args.verbose)

    logging.info(f"{'='*60}")
    logging.info(f"  🤖 MCCE4 Topology Agent — mcce4_agent_ftpl.py v2.1")
    logging.info(f"{'='*60}")
    logging.info(f"  Input PDB: {pdb_path}")
    logging.info(f"  Ligand ID: {lig_id}")
    logging.info(f"  pH: {args.ph}   Charge method: {args.charge_method}")
    logging.info(f"  Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logging.info(f"{'='*60}\n")

    # ── PHASE 1: Molecule Intelligence ──
    logging.info(f"{'─'*60}")
    logging.info(f"  PHASE 1: Molecule Intelligence")
    logging.info(f"{'─'*60}")
    mol = MoleculeIntelligence(lig_id, pdb_path)
    mol.fetch_from_rcsb()
    if not mol.smiles:
        mol.get_smiles_from_pdb()

    # ── PHASE 2: Determine conformer states ──
    logging.info(f"\n{'─'*60}")
    logging.info(f"  PHASE 2: Conformer State Determination")
    logging.info(f"{'─'*60}")

    states = []

    # Option A: User provided state-specific PDBs
    if args.state_pdbs:
        logging.info("  📂 Using user-supplied state PDB files:")
        for sp in args.state_pdbs:
            if not os.path.exists(sp):
                logging.error(f"  State PDB not found: {sp}")
                sys.exit(1)
            label = parse_state_pdb_label(sp, lig_id)
            if not label:
                logging.error(f"  Could not extract state label from filename: {sp}")
                logging.error(f"  Expected format: {lig_id}_01.pdb, {lig_id}_+1.pdb, etc.")
                sys.exit(1)
            charge = int(label.replace("+", "").replace("-", "-") or "0") if label != "01" else 0
            # Better charge parsing
            if label.startswith("+"):
                charge = int(label)
            elif label.startswith("-"):
                charge = int(label)
            elif label == "01":
                charge = 0
            else:
                charge = 0

            states.append(ConformerState(
                label=label, charge=charge,
                nH=charge,  # default: nH = charge relative to neutral
                pdb_path=os.path.abspath(sp),
                source="user",
                rationale=f"User-supplied PDB: {os.path.basename(sp)}"
            ))
            logging.info(f"     {label}: {sp} (charge={charge:+d})")

    # Option B: Auto-enumerate via Dimorphite-DL
    if not states and mol.smiles:
        prot = ProtonationAnalyzer(mol.smiles, ph=args.ph,
                                    ph_tolerance=args.ph_tolerance)
        states = prot.enumerate_states()

    # Option C: Fallback — neutral only
    if not states:
        logging.warning("  No SMILES and no state PDBs — defaulting to neutral")
        states = [ConformerState(label="01", charge=0, nH=0,
                                 source="fallback", rationale="Default neutral")]

    # ── PHASE 3: LLM Reasoning ──
    logging.info(f"\n{'─'*60}")
    logging.info(f"  PHASE 3: Agent Reasoning")
    logging.info(f"{'─'*60}")
    brain = AgentBrain(enabled=not args.no_llm)
    if brain.available and mol.smiles:
        states = brain.refine_states(
            mol.name or lig_id, mol.smiles, mol.formula or "",
            states
        )

    # ── PHASE 4: Review (GUI or terminal) ──
    if args.gui or args.interactive:
        logging.info(f"\n{'─'*60}")
        logging.info(f"  PHASE 4: User Review")
        logging.info(f"{'─'*60}")
        if args.gui:
            states = launch_gui_review(lig_id, pdb_path, states)
        else:
            states = _terminal_review(lig_id, states)

    conformer_labels = [s.label for s in states]
    logging.info(f"\n  ✅ Final conformers: {conformer_labels}")

    # ── PHASE 5: Generate template ──
    logging.info(f"\n{'─'*60}")
    logging.info(f"  PHASE 5: Template Generation")
    logging.info(f"{'─'*60}")
    generate_ftpl_template(lig_id, pdb_path, conformer_labels, ftpl_path)

    # Hybridization validation
    if brain.available and mol.smiles:
        with open(ftpl_path) as f:
            brain.validate_hybridizations(f.read(), mol.smiles, lig_id)

    # ── PHASE 6: Charges ──
    logging.info(f"\n{'─'*60}")
    logging.info(f"  PHASE 6: Charge Assignment")
    logging.info(f"{'─'*60}")
    charges = generate_charges(pdb_path, args.charge_method, lig_id)
    if not charges:
        logging.error("  No charges generated!"); sys.exit(1)
    unfilled = fill_ftpl_charges(ftpl_path, charges, lig_id)

    # ── PHASE 6b: Update CONFORMER params ──
    update_conformer_params(ftpl_path, states, lig_id)

    # ── PHASE 7: RXN Calibration ──
    if args.dry_run:
        logging.info(f"\n  ⏩ Dry run — skipping RXN calibration")
    elif unfilled > 0:
        logging.warning(f"\n  ⏩ Skipping calibration ({unfilled} unfilled charges)")
    else:
        logging.info(f"\n{'─'*60}")
        logging.info(f"  PHASE 7: RXN Calibration")
        logging.info(f"{'─'*60}")
        run_rxn_calibration(lig_id, ftpl_path, pdb_path, conformer_labels,
                             args.dielectric, work_dir)

    # ── Summary ──
    logging.info(f"\n{'='*60}")
    logging.info(f"  🤖 MCCE4 Topology Agent — COMPLETE")
    logging.info(f"{'='*60}")
    logging.info(f"  ✅ Topology file: {ftpl_path}")
    logging.info(f"  📋 Log:           mcce4_agent_ftpl_{lig_id}.log")
    logging.info(f"  🧪 Conformers:    {conformer_labels}")
    logging.info(f"  ⚡ Charges:       {args.charge_method}")
    logging.info(f"  🌊 Dielectrics:   {args.dielectric}")
    logging.info(f"{'='*60}")


if __name__ == "__main__":
    main()
