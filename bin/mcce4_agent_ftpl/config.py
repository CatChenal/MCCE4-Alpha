"""
Configuration constants and settings for the MCCE4 Topology Agent.
"""

import os

# ──────────────────────────────────────────────────────────────────────────────
# URLs
# ──────────────────────────────────────────────────────────────────────────────
RCSB_SMILES_URL = "https://data.rcsb.org/rest/v1/core/chemcomp/{lig_id}"
LIGAND_EXPO_URL = "http://ligand-expo.rcsb.org/reports/{first_char}/{lig_id}/{lig_id}_ideal.pdb"

# ──────────────────────────────────────────────────────────────────────────────
# Charge methods
# ──────────────────────────────────────────────────────────────────────────────
SUPPORTED_CHARGE_METHODS = [
    "mmff94", "am1bcc", "am1bccelf10", "am1bccnosymspt",
    "amber", "amberff94", "antechamber",
]
DEFAULT_CHARGE_METHOD = "am1bccelf10"

# ──────────────────────────────────────────────────────────────────────────────
# Dielectric constants for RXN calibration
# ──────────────────────────────────────────────────────────────────────────────
DIELECTRIC_MAP = {2: "rxn02", 4: "rxn04", 8: "rxn08"}
DEFAULT_DIELECTRICS = [2, 4, 8]

# ──────────────────────────────────────────────────────────────────────────────
# MCCE conformer naming
# ──────────────────────────────────────────────────────────────────────────────
CHARGE_TO_CONF = {
    0: "01",    # neutral
    1: "+1",    # protonated (+1 charge)
    -1: "-1",   # deprotonated (-1 charge)
    2: "+2",    # doubly protonated
    -2: "-2",   # doubly deprotonated
}

# ──────────────────────────────────────────────────────────────────────────────
# LLM
# ──────────────────────────────────────────────────────────────────────────────
GEMINI_MODEL = "gemini-2.5-flash"

def get_gemini_api_key() -> str:
    """Get API key from environment."""
    return os.environ.get("GEMINI_API_KEY") or os.environ.get("GOOGLE_API_KEY", "")

# ──────────────────────────────────────────────────────────────────────────────
# GUI
# ──────────────────────────────────────────────────────────────────────────────
STREAMLIT_PORT = 8501
GUI_TITLE = "🤖 MCCE4 Topology Agent"

# 7-phase pipeline definitions for the GUI progress tracker
AGENT_PHASES = [
    {"id": "phase1", "num": 1, "name": "Molecule Intelligence",
     "desc": "Fetch RCSB info, detect ionizable sites"},
    {"id": "phase2", "num": 2, "name": "Protonation State Enumeration",
     "desc": "Dimorphite-DL: enumerate states, apply naming (01, +1, -1, +a, +b)"},
    {"id": "phase3", "num": 3, "name": "Agent Reasoning (LLM)",
     "desc": "LLM validates states, estimates pKa, checks chemistry"},
    {"id": "phase4", "num": 4, "name": "Per-State PDB Generation",
     "desc": "Generate PDB with correct H atoms for each protonation state"},
    {"id": "phase5", "num": 5, "name": "Template Generation & Validation",
     "desc": "pdb2ftpl.py per state, merge into single .ftpl"},
    {"id": "phase6", "num": 6, "name": "Charge Assignment",
     "desc": "OpenEye QuacPac TK partial charges per state"},
    {"id": "phase7", "num": 7, "name": "RXN Calibration",
     "desc": "MCCE step1/step2/step3 for each dielectric (2, 4, 8)"},
]

# Map agent internal phase names to phase IDs
PHASE_NAME_MAP = {
    "init": None,
    "molecule_intel": "phase1",
    "enumerate_states": "phase2",
    "llm_reasoning": "phase3",
    "generate_state_pdbs": "phase4",
    "user_review": "phase4",  # review happens after phase 4
    "generate_template": "phase5",
    "assign_charges": "phase6",
    "rxn_calibration": "phase7",
}

# ──────────────────────────────────────────────────────────────────────────────
# v3: Per-state PDB / RXN calibration
# ──────────────────────────────────────────────────────────────────────────────
DSOLV_TOLERANCE = 0.5     # kcal/mol — consider converged if |dsolv| < this
MAX_RXN_ITERATIONS = 3    # Max re-runs of step3 per dielectric
PDB2FTPL_CMD = "pdb2ftpl.py"
STEP1_CMD = "step1.py"
STEP3_CMD = "step3.py"
