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


def make_labels_unique(states: list) -> list:
    """Ensure all conformer state labels are distinct and ≤ 2 characters.

    MCCE conformer type names are lig_id (3 chars) + label (≤ 2 chars),
    so labels must not exceed 2 characters.

    When multiple states share the same base label the disambiguation scheme
    replaces the numeric part with a letter to stay within 2 chars:

        '+1', '+1'  →  '+a', '+b'   (sign + letter)
        '-1', '-1'  →  '-a', '-b'
        '01', '01'  →  '0a', '0b'
        '+2', '+2'  →  '2a', '2b'   (digit + letter, drop sign to fit)
        '-2', '-2'  →  'Ma', 'Mb'   (M = minus, + letter)

    A lone state keeps its original label.

    Works with both ConformerState objects (mutated in place) and dicts
    (replaced with a shallow copy).  Returns the modified list.
    """
    import logging
    from collections import Counter

    def _disambig(base: str, idx: int) -> str:
        """Return a unique 2-char label for the idx-th duplicate of base."""
        letter = chr(ord('a') + idx)
        if base in ('+1', '-1', '01'):
            return base[0] + letter          # '+a'/'-a'/'0a'
        sign = base[0] if base[0] in ('+', '-') else ''
        digit = base[1] if len(base) > 1 and base[1].isdigit() else base[0]
        if sign == '+':
            return digit + letter            # '2a', '3a', …
        if sign == '-':
            return 'M' + letter             # 'Ma', 'Mb', … (Minus)
        return base[0] + letter             # fallback

    label_counts = Counter(
        s['label'] if isinstance(s, dict) else s.label for s in states
    )
    label_index: dict = {}
    result = []
    for s in states:
        label = s['label'] if isinstance(s, dict) else s.label
        if label_counts[label] > 1:
            idx = label_index.get(label, 0)
            new_label = _disambig(label, idx)
            label_index[label] = idx + 1
            logging.info(f"  Label disambiguated: '{label}' → '{new_label}'")
        else:
            new_label = label
        if isinstance(s, dict):
            s = {**s, 'label': new_label}
        else:
            s.label = new_label
        result.append(s)
    return result


def sort_conformer_labels(labels):
    """Sort conformer labels: neutral first, then positive, then negative.

    Label categories:
      Neutral:  starts with '0'  (01, 0a, 0b)
      Positive: starts with '+' or digit 1-9 (e.g., +1, +a, 2a, 3b)
      Negative: starts with '-' or 'M'  (e.g., -1, -a, Ma, Mb)

    Within each category, labels are sorted alphabetically.
    Example: ['01', '+a', '+b', '2a', '-a', 'Ma']
    """
    def _sort_key(label):
        if label.startswith('0'):
            return (0, label)  # neutral
        elif label.startswith('+') or (label[0].isdigit() and label[0] != '0'):
            return (1, label)  # positive
        elif label.startswith('-') or label.startswith('M'):
            return (2, label)  # negative
        else:
            return (3, label)  # unknown — put last
    return sorted(labels, key=_sort_key)


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
