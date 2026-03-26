"""
Data models for the MCCE4 Topology Agent v3.

v3 additions to ConformerState:
  - h_added / h_removed: track which H atoms differ from neutral
  - per_state_charges: charges computed on this state's PDB
  - rxn_values / dsolv_values: per-state calibration results
"""

from dataclasses import dataclass, field, asdict
from typing import Optional, TypedDict


@dataclass
class ConformerState:
    """Represents a single protonation/redox conformer state.

    v3: Each state has its own PDB with correct H atoms.
    """
    label: str              # e.g., "01", "+1", "-1"
    charge: int = 0         # formal charge
    nH: int = 0             # protons relative to neutral
    pka: Optional[float] = None  # solution pKa
    smiles: str = ""        # SMILES for this state
    pdb_path: Optional[str] = None  # PDB file for this state (v3: per-state)
    source: str = "auto"    # "user", "dimorphite", "llm", "auto"
    site_atom: Optional[str] = None  # atom where protonation occurs
    rationale: str = ""     # why this state was included

    # ── v3: per-state PDB fields ──
    h_added: list = field(default_factory=list)     # H atom names added vs neutral
    h_removed: list = field(default_factory=list)   # H atom names removed vs neutral
    per_state_charges: dict = field(default_factory=dict)  # atom_name → charge
    rxn_values: dict = field(default_factory=dict)  # {rxn02: v, rxn04: v, rxn08: v}
    dsolv_values: dict = field(default_factory=dict)  # {2: dsolv, 4: dsolv, 8: dsolv}

    # ── v3: provenance & proton exchange info ──
    proton_exchange: str = ""   # Which protons differ vs neutral (e.g., "+H on N1 piperidine")
    llm_model: str = ""         # LLM that produced this state (e.g., "gemini-2.5-flash")
    references: list = field(default_factory=list)  # Literature refs (DOIs, URLs, titles)

    def to_dict(self) -> dict:
        return asdict(self)

    def __repr__(self):
        return (f"ConformerState({self.label}, charge={self.charge:+d}, "
                f"nH={self.nH:+d}, source={self.source})")


class AgentState(TypedDict, total=False):
    """State object for the LangGraph agent.

    This is passed between nodes in the agent graph.
    Each node reads/writes fields as needed.
    """
    # ── Inputs ──
    pdb_path: str                   # Input ligand PDB file
    lig_id: str                     # 3-letter ligand code
    ph: float                       # Target pH
    charge_method: str              # Charge calculation method
    dielectrics: list               # Dielectric constants for RXN
    work_dir: str                   # Working directory

    # ── Molecule info ──
    smiles: str                     # Canonical SMILES
    name: str                       # Molecule name
    formula: str                    # Molecular formula
    formal_charge: int              # Net formal charge from RCSB

    # ── Conformer states ──
    states: list                    # List of ConformerState dicts
    conformer_labels: list          # List of label strings (e.g., ["01", "+1"])
    user_state_pdbs: list           # User-supplied state PDB paths

    # ── v3: Per-state PDBs ──
    state_pdbs: dict                # label → pdb_path (each with correct H atoms)
    h_diffs: dict                   # label → {"added": [...], "removed": [...]}
    per_state_connects: dict        # label → ftpl_path from pdb2ftpl

    # ── Generated files ──
    ftpl_path: str                  # Output .ftpl file path
    charges: dict                   # atom_name -> charge mapping
    per_state_charges: dict         # label → {atom_name: charge}

    # ── Validation ──
    hybridization_issues: list      # List of hybridization mismatches
    charge_issues: list             # List of charge validation issues
    rxn_values: dict                # conformer_type -> {rxn02: v, rxn04: v, rxn08: v}

    # ── Agent control ──
    phase: str                      # Current phase name
    errors: list                    # Accumulated errors
    warnings: list                  # Accumulated warnings
    messages: list                  # LLM conversation messages
    needs_user_review: bool         # Whether to pause for GUI review
    user_approved: bool             # Whether user approved states
    complete: bool                  # Whether agent is done
    dry_run: bool                   # Skip RXN calibration
