"""
RDKit-based tools for molecule intelligence and validation.
"""

import logging
import os
import re
import io
from typing import Optional
from pathlib import Path


def get_smiles_from_pdb(pdb_path: str) -> Optional[str]:
    """Extract SMILES from a PDB file using RDKit."""
    try:
        from rdkit import Chem
        mol = Chem.MolFromPDBFile(pdb_path, removeHs=False)
        if mol:
            smiles = Chem.MolToSmiles(mol)
            logging.info(f"  ✓ SMILES from PDB: {smiles}")
            return smiles
    except ImportError:
        logging.warning("  RDKit not available")
    except Exception as e:
        logging.warning(f"  RDKit SMILES extraction failed: {e}")
    return None


def get_mol_from_pdb(pdb_path: str, remove_hs: bool = True):
    """Load an RDKit Mol object from PDB, with fallback for sanitization."""
    try:
        from rdkit import Chem
        mol = Chem.MolFromPDBFile(pdb_path, removeHs=remove_hs, sanitize=True)
        if mol is None:
            mol = Chem.MolFromPDBFile(pdb_path, removeHs=remove_hs, sanitize=False)
        return mol
    except ImportError:
        return None
    except Exception:
        return None


def get_mol_from_smiles(smiles: str):
    """Load an RDKit Mol from SMILES."""
    try:
        from rdkit import Chem
        return Chem.MolFromSmiles(smiles)
    except ImportError:
        return None


def get_formal_charge(smiles: str) -> int:
    """Compute formal charge from SMILES."""
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Chem.GetFormalCharge(mol)
    except ImportError:
        pass
    return 0


def mol_to_png_bytes(mol_or_pdb, size=(450, 350), highlight_atoms=None) -> Optional[bytes]:
    """Render a molecule to PNG bytes for GUI display.

    Args:
        mol_or_pdb: RDKit Mol object or path to PDB file.
        size: Image dimensions (width, height).
        highlight_atoms: List of atom indices to highlight (for protonation sites).

    Returns:
        PNG image bytes, or None if rendering fails.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw, AllChem

        if isinstance(mol_or_pdb, str):
            mol = get_mol_from_pdb(mol_or_pdb, remove_hs=True)
        else:
            mol = mol_or_pdb

        if mol is None:
            return None

        # Generate 2D coordinates for clean depiction
        AllChem.Compute2DCoords(mol)

        # Draw with optional highlighting
        drawer = Draw.MolDraw2DPNG(size[0], size[1])
        if highlight_atoms:
            colors = {idx: (1.0, 0.4, 0.4) for idx in highlight_atoms}  # red
            drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms,
                                highlightAtomColors=colors)
        else:
            drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()

    except ImportError:
        logging.debug("  RDKit not available for molecule rendering")
        return None
    except Exception as e:
        logging.debug(f"  Molecule rendering failed: {e}")
        return None


def mol_to_svg(mol_or_pdb, size=(450, 350), highlight_atoms=None) -> Optional[str]:
    """Render a molecule to SVG string for web GUI.

    Args:
        mol_or_pdb: RDKit Mol object or path to PDB file.
        size: Image dimensions.
        highlight_atoms: List of atom indices to highlight.

    Returns:
        SVG string, or None.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw, AllChem

        if isinstance(mol_or_pdb, str):
            mol = get_mol_from_pdb(mol_or_pdb, remove_hs=True)
        else:
            mol = mol_or_pdb

        if mol is None:
            return None

        AllChem.Compute2DCoords(mol)

        drawer = Draw.MolDraw2DSVG(size[0], size[1])
        if highlight_atoms:
            colors = {idx: (1.0, 0.3, 0.3) for idx in highlight_atoms}
            drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms,
                                highlightAtomColors=colors)
        else:
            drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()

    except Exception as e:
        logging.debug(f"  SVG rendering failed: {e}")
        return None


def find_protonation_site_atoms(smiles_neutral: str, smiles_protonated: str) -> list:
    """Find atom indices that differ between two protonation states.

    Compares hydrogen counts per heavy atom to identify where protons were added/removed.

    Returns:
        List of atom indices (in the neutral molecule) where protonation changed.
    """
    try:
        from rdkit import Chem

        mol_n = Chem.MolFromSmiles(smiles_neutral)
        mol_p = Chem.MolFromSmiles(smiles_protonated)
        if mol_n is None or mol_p is None:
            return []

        mol_n = Chem.AddHs(mol_n)
        mol_p = Chem.AddHs(mol_p)

        # Count H neighbors per heavy atom
        def h_counts(mol):
            counts = {}
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() != 1:
                    n_h = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 1)
                    counts[atom.GetIdx()] = n_h
            return counts

        hc_n = h_counts(mol_n)
        hc_p = h_counts(mol_p)

        # Find atoms where H count changed
        changed = []
        for idx in hc_n:
            if idx in hc_p and hc_n[idx] != hc_p[idx]:
                changed.append(idx)

        return changed

    except Exception:
        return []


RDKIT_TO_MCCE_HYB = {
    "SP": "sp", "SP2": "sp2", "SP3": "sp3",
    "SP3D2": "d2sp3", "S": "s",
}


# ═════════════════════════════════════════════════════════════════════════════
# v3: Per-State PDB Generation
# ═════════════════════════════════════════════════════════════════════════════

def get_atom_info(mol) -> list:
    """Get atom names, symbols, and indices from PDB molecule.

    Returns list of dicts: [{"idx": 0, "name": "C1", "symbol": "C", "is_h": False}, ...]
    """
    atoms = []
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        name = info.GetName().strip() if info else atom.GetSymbol()
        atoms.append({
            "idx": atom.GetIdx(),
            "name": name,
            "symbol": atom.GetSymbol(),
            "is_h": atom.GetAtomicNum() == 1,
        })
    return atoms


def get_ionizable_sites(mol) -> list:
    """Identify potential ionizable sites (N, O, S with lone pairs).

    Returns list of dicts with atom index, name, symbol, and site type.
    """
    try:
        from rdkit import Chem
    except ImportError:
        return []

    sites = []
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol not in ("N", "O", "S"):
            continue
        info = atom.GetPDBResidueInfo()
        name = info.GetName().strip() if info else f"{symbol}{atom.GetIdx()}"
        n_hs = atom.GetTotalNumHs()
        fc = atom.GetFormalCharge()

        if symbol == "N" and (n_hs > 0 or fc == 0):
            sites.append({"idx": atom.GetIdx(), "name": name, "symbol": symbol,
                          "type": "protonatable_N", "current_hs": n_hs, "formal_charge": fc})
        elif symbol == "O" and n_hs > 0:
            sites.append({"idx": atom.GetIdx(), "name": name, "symbol": symbol,
                          "type": "deprotonatable_OH", "current_hs": n_hs, "formal_charge": fc})
        elif symbol == "S" and n_hs > 0:
            sites.append({"idx": atom.GetIdx(), "name": name, "symbol": symbol,
                          "type": "deprotonatable_SH", "current_hs": n_hs, "formal_charge": fc})
    return sites


def generate_state_pdb(neutral_pdb, state_smiles, lig_id, label, output_dir="."):
    """Generate a PDB file for a specific protonation state.

    Strategy:
      1. Load neutral PDB as RDKit mol (keep all H)
      2. Compare H counts per heavy atom between neutral and target SMILES
      3. Add/remove H at specific sites
      4. Write new PDB preserving atom naming conventions

    Args:
        neutral_pdb: Path to neutral PDB file
        state_smiles: SMILES for this protonation state
        lig_id: 3-letter ligand code
        label: Conformer label (e.g., "+1", "-1")
        output_dir: Where to write the per-state PDB

    Returns:
        (pdb_path, h_added_names, h_removed_names)
        pdb_path is None on failure.
    """
    import shutil
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, rdFMCS
    except ImportError:
        logging.error("RDKit required for per-state PDB generation")
        return None, [], []

    # If neutral state, just copy
    if label in ("01", "00"):
        out_path = os.path.join(output_dir, f"{lig_id}_{label}.pdb")
        shutil.copy2(neutral_pdb, out_path)
        return out_path, [], []

    try:
        # 1. Load neutral mol from PDB (with H)
        neutral_mol = Chem.MolFromPDBFile(neutral_pdb, removeHs=False, sanitize=True)
        if neutral_mol is None:
            neutral_mol = Chem.MolFromPDBFile(neutral_pdb, removeHs=False, sanitize=False)
        if neutral_mol is None:
            logging.error(f"Cannot load neutral PDB: {neutral_pdb}")
            return None, [], []

        # 2. Build target mol from SMILES
        target_mol = Chem.MolFromSmiles(state_smiles)
        if target_mol is None:
            logging.error(f"Cannot parse state SMILES: {state_smiles}")
            return None, [], []
        target_mol = Chem.AddHs(target_mol)

        # 3. Compare H counts per heavy atom
        neutral_hc = _count_h_per_heavy(neutral_mol)
        target_hc = _map_target_h_counts(neutral_mol, target_mol)

        # 4. Determine what to add/remove
        h_to_add = []
        h_to_remove = []
        for heavy_idx, n_count in neutral_hc.items():
            t_count = target_hc.get(heavy_idx, n_count)
            diff = t_count - n_count
            if diff > 0:
                h_to_add.append((heavy_idx, diff))
            elif diff < 0:
                h_to_remove.append((heavy_idx, abs(diff)))

        # 5. Apply modifications
        state_mol, added_names, removed_names = _modify_mol_hydrogens(
            neutral_mol, neutral_pdb, h_to_add, h_to_remove, lig_id, label
        )
        if state_mol is None:
            return None, [], []

        # 6. Write PDB
        safe_label = label.replace('+', 'p').replace('-', 'm')
        out_path = os.path.join(output_dir, f"{lig_id}_{safe_label}.pdb")
        _write_state_pdb(state_mol, neutral_pdb, out_path, lig_id,
                         added_names, removed_names)

        logging.info(f"  State {label}: wrote {out_path} "
                     f"(+{len(added_names)}H, -{len(removed_names)}H)")
        return out_path, added_names, removed_names

    except Exception as e:
        logging.error(f"Per-state PDB generation failed for {label}: {e}")
        import traceback
        traceback.print_exc()
        return None, [], []


def compute_h_diff(neutral_pdb, state_pdb):
    """Compare H atoms between neutral and state PDBs.

    Returns:
        {"added": ["HN1"], "removed": ["HO3"],
         "added_on": {"HN1": "N1"}, "removed_from": {"HO3": "O3"}}
    """
    neutral_hs = _extract_h_atoms_from_pdb(neutral_pdb)
    state_hs = _extract_h_atoms_from_pdb(state_pdb)

    neutral_names = set(neutral_hs.keys())
    state_names = set(state_hs.keys())

    added = sorted(state_names - neutral_names)
    removed = sorted(neutral_names - state_names)

    return {
        "added": added,
        "removed": removed,
        "added_on": {h: state_hs[h].get("bonded_to", "?") for h in added},
        "removed_from": {h: neutral_hs[h].get("bonded_to", "?") for h in removed},
    }


def mol_to_svg_with_h_diff(mol_or_pdb, h_added_names, h_removed_names,
                           size=(450, 350)):
    """Render molecule SVG highlighting added (green) and removed (red) H sites.

    Shows PDB atom name labels on ALL H atoms near the changed sites
    so the user can see exactly which hydrogens differ.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw, AllChem
        from rdkit.Chem.Draw import rdMolDraw2D

        if isinstance(mol_or_pdb, str):
            mol = get_mol_from_pdb(mol_or_pdb, remove_hs=False)
        else:
            mol = mol_or_pdb
        if mol is None:
            return None

        AllChem.Compute2DCoords(mol)

        # Map PDB atom names to indices
        added_indices = []
        all_h_indices = []
        name_map = {}  # idx → PDB name
        for atom in mol.GetAtoms():
            info = atom.GetPDBResidueInfo()
            if info:
                name = info.GetName().strip()
                name_map[atom.GetIdx()] = name
                if name in h_added_names:
                    added_indices.append(atom.GetIdx())
                if atom.GetAtomicNum() == 1:
                    all_h_indices.append(atom.GetIdx())

        # Find parent heavy atoms of removed H (they won't be in this mol)
        removed_parent_indices = []
        # h_removed_names are H atoms that existed in neutral but not here
        # We can't show them, but we can mark their parent sites

        # Build highlight + color maps
        highlight_atoms = list(added_indices)
        colors = {}
        for idx in added_indices:
            colors[idx] = (0.3, 0.9, 0.3)  # green

        # Set atom labels: show PDB names on H atoms and changed sites
        drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
        opts = drawer.drawOptions()
        opts.addAtomIndices = False
        opts.addStereoAnnotation = False

        # Label ALL H atoms with their PDB name so user can see them
        for h_idx in all_h_indices:
            opts.atomLabels[h_idx] = name_map.get(h_idx, "H")

        # Also label the parent heavy atoms of added H
        for h_idx in added_indices:
            atom = mol.GetAtomWithIdx(h_idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() != 1:
                    parent_name = name_map.get(nbr.GetIdx(), "")
                    if parent_name:
                        opts.atomLabels[nbr.GetIdx()] = parent_name

        if highlight_atoms:
            drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms,
                                highlightAtomColors=colors)
        else:
            drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()

    except Exception as e:
        logging.debug(f"  SVG H-diff rendering failed: {e}")
        return None


def render_protonation_site_view(pdb_path, site_atom_name, h_atom_names,
                                  lig_id, label, size=(350, 300),
                                  highlight_color=(0.3, 0.9, 0.3)):
    """Render a FOCUSED view of a protonation site and its neighborhood.

    Shows only the ionizable site atom + atoms within 2 bonds, WITH
    all H atoms explicitly drawn and labeled with PDB names. This gives
    a clear, readable view of exactly which H atoms are present.

    Args:
        pdb_path: Path to per-state PDB
        site_atom_name: PDB name of the ionizable site (e.g., "N1")
        h_atom_names: List of H atom names to highlight (added H)
        lig_id: Ligand code
        label: Conformer label
        size: SVG dimensions
        highlight_color: RGB for highlighted H atoms

    Returns:
        SVG string or None.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit.Chem.Draw import rdMolDraw2D
    except ImportError:
        return None

    mol = Chem.MolFromPDBFile(pdb_path, removeHs=False, sanitize=False)
    if mol is None:
        return None

    # Find the site atom index
    site_idx = None
    name_map = {}
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info:
            name = info.GetName().strip()
            name_map[atom.GetIdx()] = name
            if name == site_atom_name:
                site_idx = atom.GetIdx()

    if site_idx is None:
        # Can't find the site — fall back to full molecule
        return mol_to_svg_with_h_diff(mol, h_atom_names, [], size)

    # Collect atoms within 2 bonds of the site
    neighborhood = {site_idx}
    frontier = {site_idx}
    for depth in range(2):
        next_frontier = set()
        for idx in frontier:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() not in neighborhood:
                    neighborhood.add(nbr.GetIdx())
                    next_frontier.add(nbr.GetIdx())
        frontier = next_frontier

    # Also include H neighbors of all atoms in the neighborhood
    h_to_add = set()
    for idx in list(neighborhood):
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 1:
                h_to_add.add(nbr.GetIdx())
    neighborhood.update(h_to_add)

    AllChem.Compute2DCoords(mol)

    # Highlight the H atoms that were added
    highlight = []
    colors = {}
    for idx in neighborhood:
        name = name_map.get(idx, "")
        if name in h_atom_names:
            highlight.append(idx)
            colors[idx] = highlight_color

    # Also highlight the site atom
    highlight.append(site_idx)
    colors[site_idx] = (0.7, 0.85, 1.0)  # light blue for the site

    # Draw with atom labels for the neighborhood
    drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
    opts = drawer.drawOptions()
    opts.addAtomIndices = False

    # Label atoms in the neighborhood with their PDB names
    for idx in neighborhood:
        name = name_map.get(idx, "")
        if name:
            opts.atomLabels[idx] = name

    # Draw only the neighborhood atoms highlighted, but show whole mol faded
    # (RDKit doesn't support substructure-only drawing easily,
    #  so we highlight the neighborhood and let the rest fade)
    drawer.DrawMolecule(mol, highlightAtoms=highlight,
                        highlightAtomColors=colors)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def generate_state_comparison_table(states, h_diffs, lig_id):
    """Generate a text-based comparison of H atoms across all states.

    Returns a list of dicts suitable for a Streamlit dataframe:
    [{"State": "EMH01", "Charge": "+0", "H atoms at N1": "HN1",
      "H atoms at O3": "HO3", "Δ vs neutral": "none"}, ...]
    """
    if not states:
        return []

    # Collect all ionizable site names across all states
    all_sites = set()
    for label, diff in h_diffs.items():
        for h_name in diff.get("added", []):
            parent = diff.get("added_on", {}).get(h_name, "")
            if parent:
                all_sites.add(parent)
        for h_name in diff.get("removed", []):
            parent = diff.get("removed_from", {}).get(h_name, "")
            if parent:
                all_sites.add(parent)

    rows = []
    for s in states:
        sd = s if isinstance(s, dict) else s.to_dict()
        label = sd.get("label", "?")
        charge = sd.get("charge", 0)
        diff = h_diffs.get(label, {})

        row = {
            "State": f"{lig_id}{label}",
            "Charge": f"{charge:+d}",
        }

        # Show H atoms at each ionizable site
        added = diff.get("added", [])
        added_on = diff.get("added_on", {})
        removed = diff.get("removed", [])
        removed_from = diff.get("removed_from", {})

        # Build delta description
        deltas = []
        for h in added:
            parent = added_on.get(h, "?")
            deltas.append(f"+{h} on {parent}")
        for h in removed:
            parent = removed_from.get(h, "?")
            deltas.append(f"-{h} from {parent}")

        row["Δ vs neutral"] = "; ".join(deltas) if deltas else "— (reference)"
        row["H added"] = ", ".join(added) if added else "—"
        row["H removed"] = ", ".join(removed) if removed else "—"

        rows.append(row)

    return rows


# ── v3 internal helpers ──────────────────────────────────────────────────────

def _count_h_per_heavy(mol) -> dict:
    """Count H atoms bonded to each heavy atom. Returns {heavy_idx: count}."""
    counts = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        n_h = sum(1 for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 1)
        counts[atom.GetIdx()] = n_h
    return counts


def _map_target_h_counts(neutral_mol, target_mol) -> dict:
    """Map target SMILES H counts back to neutral PDB atom indices.

    Uses substructure matching on heavy atoms to align the two molecules.
    Returns {neutral_heavy_idx: target_h_count}.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import rdFMCS

        neutral_noH = Chem.RemoveHs(Chem.RWMol(neutral_mol))
        target_noH = Chem.RemoveHs(Chem.RWMol(target_mol))

        match = neutral_noH.GetSubstructMatch(target_noH)
        if not match and target_noH.GetNumAtoms() <= neutral_noH.GetNumAtoms():
            match_rev = target_noH.GetSubstructMatch(neutral_noH)
            if match_rev:
                # Build forward mapping
                mapping = {}
                for neut_i, tgt_i in enumerate(match_rev):
                    mapping[neut_i] = tgt_i
                result = {}
                for neut_noH_idx, tgt_noH_idx in mapping.items():
                    tgt_atom = target_noH.GetAtomWithIdx(tgt_noH_idx)
                    n_h = tgt_atom.GetTotalNumHs()
                    neut_full_idx = _noH_idx_to_full(neutral_mol, neut_noH_idx)
                    result[neut_full_idx] = n_h
                return result

        if not match:
            # MCS fallback
            mcs = rdFMCS.FindMCS([neutral_noH, target_noH], timeout=10)
            mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
            if mcs_mol is None:
                return _count_h_per_heavy(neutral_mol)  # give up
            m1 = neutral_noH.GetSubstructMatch(mcs_mol)
            m2 = target_noH.GetSubstructMatch(mcs_mol)
            result = {}
            for i, j in zip(m1, m2):
                n_h = target_noH.GetAtomWithIdx(j).GetTotalNumHs()
                neut_full_idx = _noH_idx_to_full(neutral_mol, i)
                result[neut_full_idx] = n_h
            return result

        result = {}
        for neut_noH_idx in range(neutral_noH.GetNumAtoms()):
            if neut_noH_idx >= len(match):
                continue
            tgt_noH_idx = match[neut_noH_idx]
            n_h = target_noH.GetAtomWithIdx(tgt_noH_idx).GetTotalNumHs()
            neut_full_idx = _noH_idx_to_full(neutral_mol, neut_noH_idx)
            result[neut_full_idx] = n_h
        return result

    except Exception as e:
        logging.warning(f"  H-count mapping failed: {e}")
        return _count_h_per_heavy(neutral_mol)


def _noH_idx_to_full(full_mol, noH_idx):
    """Map no-H atom index to full mol (with H) atom index."""
    heavy_count = 0
    for atom in full_mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        if heavy_count == noH_idx:
            return atom.GetIdx()
        heavy_count += 1
    return noH_idx


def _modify_mol_hydrogens(mol, pdb_path, h_to_add, h_to_remove, lig_id, label):
    """Add/remove H atoms on a copy of the neutral mol.

    Returns (modified_mol, added_atom_names, removed_atom_names).
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        return None, [], []

    added_names = []
    removed_names = []
    rw_mol = Chem.RWMol(mol)

    # Remove H first (high→low index to preserve indices)
    h_indices_to_remove = []
    for heavy_idx, count in h_to_remove:
        heavy_atom = rw_mol.GetAtomWithIdx(heavy_idx)
        heavy_info = heavy_atom.GetPDBResidueInfo()
        heavy_name = heavy_info.GetName().strip() if heavy_info else "?"
        h_neighbors = [n for n in heavy_atom.GetNeighbors() if n.GetAtomicNum() == 1]
        for h_atom in h_neighbors[:count]:
            h_info = h_atom.GetPDBResidueInfo()
            h_name = h_info.GetName().strip() if h_info else f"H_{heavy_name}"
            h_indices_to_remove.append(h_atom.GetIdx())
            removed_names.append(h_name)
            logging.info(f"    Remove H: {h_name} (from {heavy_name})")

    for idx in sorted(h_indices_to_remove, reverse=True):
        rw_mol.RemoveAtom(idx)

    # Add H atoms
    for heavy_idx_orig, count in h_to_add:
        heavy_info_orig = mol.GetAtomWithIdx(heavy_idx_orig).GetPDBResidueInfo()
        heavy_name = heavy_info_orig.GetName().strip() if heavy_info_orig else None

        # Find in modified mol by PDB name
        heavy_idx_new = None
        for atom in rw_mol.GetAtoms():
            info = atom.GetPDBResidueInfo()
            if info and info.GetName().strip() == heavy_name:
                heavy_idx_new = atom.GetIdx()
                break
        if heavy_idx_new is None:
            logging.warning(f"Cannot find {heavy_name} after removal step")
            continue

        for i in range(count):
            h_idx = rw_mol.AddAtom(Chem.Atom(1))
            rw_mol.AddBond(heavy_idx_new, h_idx, Chem.BondType.SINGLE)
            h_name = _make_h_name(heavy_name, rw_mol)
            added_names.append(h_name)

            # Set PDB info
            ri = Chem.AtomPDBResidueInfo()
            ri.SetName(f" {h_name:<3s}")
            ri.SetResidueName(lig_id)
            ri.SetResidueNumber(1)
            ri.SetChainId("L")
            ri.SetIsHeteroAtom(True)
            rw_mol.GetAtomWithIdx(h_idx).SetPDBResidueInfo(ri)
            logging.info(f"    Add H: {h_name} (on {heavy_name})")

    try:
        Chem.SanitizeMol(rw_mol)
    except Exception as e:
        logging.warning(f"Sanitization warning for {label}: {e}")

    try:
        AllChem.EmbedMolecule(rw_mol, useRandomCoords=True, randomSeed=42)
        if mol.GetNumConformers() > 0:
            AllChem.AlignMol(rw_mol, mol)
    except Exception:
        pass

    return rw_mol.GetMol(), added_names, removed_names


def _make_h_name(heavy_name, mol):
    """Generate MCCE-compatible H atom name: H on N1 → HN1, then HN1A, etc."""
    base = f"H{heavy_name}"[:4]
    existing = set()
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info:
            existing.add(info.GetName().strip())
    if base not in existing:
        return base
    for suffix in "ABCDEFGH":
        candidate = f"{base[:3]}{suffix}"[:4]
        if candidate not in existing:
            return candidate
    return base


def _extract_h_atoms_from_pdb(pdb_path):
    """Extract H atom info from PDB. Returns {name: {"serial": N, "bonded_to": "X"}}."""
    atoms = {}
    all_atoms = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                name = line[12:16].strip()
                serial = int(line[6:11].strip())
                element = line[76:78].strip() if len(line) > 76 else ""
                all_atoms[serial] = name
                if element == "H" or name.startswith("H"):
                    atoms[name] = {"serial": serial, "bonded_to": "?"}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("CONECT"):
                parts = line[6:].split()
                if not parts:
                    continue
                serial = int(parts[0])
                name = all_atoms.get(serial, "")
                if name in atoms:
                    for nbr_str in parts[1:]:
                        try:
                            nbr_name = all_atoms.get(int(nbr_str), "")
                            if nbr_name and not nbr_name.startswith("H"):
                                atoms[name]["bonded_to"] = nbr_name
                                break
                        except ValueError:
                            continue
    return atoms


def _write_state_pdb(state_mol, orig_pdb, out_path, lig_id,
                     added_names, removed_names):
    """Write per-state PDB preserving original formatting for unchanged atoms."""
    try:
        from rdkit import Chem
    except ImportError:
        return

    with open(orig_pdb) as f:
        orig_lines = f.readlines()

    removed_set = set(removed_names)
    kept_lines = []
    for line in orig_lines:
        if line.startswith(("ATOM", "HETATM")):
            atom_name = line[12:16].strip()
            if atom_name in removed_set:
                removed_set.discard(atom_name)
                continue
            kept_lines.append(line)
        elif line.startswith(("CONECT", "END")):
            continue
        else:
            kept_lines.append(line)

    # Add new H atoms
    next_serial = max(
        (int(l[6:11].strip()) for l in kept_lines if l.startswith(("ATOM", "HETATM"))),
        default=0
    ) + 1

    conf = state_mol.GetConformer(0) if state_mol.GetNumConformers() > 0 else None
    for h_name in added_names:
        pos = (0.0, 0.0, 0.0)
        for atom in state_mol.GetAtoms():
            info = atom.GetPDBResidueInfo()
            if info and info.GetName().strip() == h_name and conf:
                p = conf.GetAtomPosition(atom.GetIdx())
                pos = (p.x, p.y, p.z)
                break
        hetatm = (f"HETATM{next_serial:5d}  {h_name:<3s} {lig_id:>3s} L   1    "
                  f"{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}"
                  f"  1.00  0.00           H  \n")
        kept_lines.append(hetatm)
        next_serial += 1

    # Regenerate CONECT from RDKit bonds
    name_to_serial = {}
    for line in kept_lines:
        if line.startswith(("ATOM", "HETATM")):
            name_to_serial[line[12:16].strip()] = int(line[6:11].strip())

    idx_to_serial = {}
    for atom in state_mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info:
            name = info.GetName().strip()
            if name in name_to_serial:
                idx_to_serial[atom.GetIdx()] = name_to_serial[name]

    conect = {}
    for bond in state_mol.GetBonds():
        s1 = idx_to_serial.get(bond.GetBeginAtomIdx())
        s2 = idx_to_serial.get(bond.GetEndAtomIdx())
        if s1 and s2:
            conect.setdefault(s1, []).append(s2)
            conect.setdefault(s2, []).append(s1)

    for serial in sorted(conect.keys()):
        nbrs = sorted(conect[serial])[:4]
        line = f"CONECT{serial:5d}" + "".join(f"{n:5d}" for n in nbrs) + "\n"
        kept_lines.append(line)

    kept_lines.append("END\n")
    with open(out_path, "w") as f:
        f.writelines(kept_lines)


def validate_hybridizations(pdb_path: str, ftpl_content: str,
                              lig_id: str) -> list:
    """Compare RDKit-perceived hybridizations against CONNECT records.

    Returns:
        List of (atom_name, ftpl_hyb, rdkit_hyb) for mismatches.
    """
    try:
        from rdkit import Chem
    except ImportError:
        return []

    logging.info("🔬 Validating hybridizations with RDKit...")
    mol = get_mol_from_pdb(pdb_path, remove_hs=False)
    if mol is None:
        logging.warning("  RDKit could not parse PDB")
        return []

    # Build RDKit hybridization map
    rdkit_hybs = {}
    for atom in mol.GetAtoms():
        pdb_info = atom.GetPDBResidueInfo()
        if pdb_info:
            name = pdb_info.GetName().strip()
            hyb = str(atom.GetHybridization()).split(".")[-1]
            rdkit_hybs[name] = RDKIT_TO_MCCE_HYB.get(hyb, hyb.lower())

    # Parse CONNECT lines
    pat = re.compile(
        r'CONNECT,\s*"(.{4})",\s*' + re.escape(lig_id) + r'\S+:\s*(sp3|sp2|sp|s|d2sp3)'
    )

    mismatches = []
    for line in ftpl_content.splitlines():
        m = pat.search(line)
        if m:
            atom_name = m.group(1).strip()
            ftpl_hyb = m.group(2)
            rdkit_hyb = rdkit_hybs.get(atom_name)
            if rdkit_hyb and rdkit_hyb != ftpl_hyb and not atom_name.startswith("H"):
                mismatches.append((atom_name, ftpl_hyb, rdkit_hyb))

    if mismatches:
        logging.info(f"  ⚠ {len(mismatches)} hybridization mismatch(es)")
    else:
        logging.info(f"  ✓ All hybridizations match")

    return mismatches


# ═════════════════════════════════════════════════════════════════════════════
# v3: Side-by-side H-diff grid visualization
# ═════════════════════════════════════════════════════════════════════════════

def generate_h_diff_grid_svg(state_pdbs, h_diffs, lig_id,
                             cell_size=(400, 350), max_cols=3):
    """Generate a single SVG showing all protonation states side by side.

    Each cell shows the molecule for that state with:
      - Added H atoms highlighted in GREEN with atom name labels
      - Removed H sites highlighted in RED on the parent heavy atom
      - State label and charge as caption

    Args:
        state_pdbs: {label: pdb_path}
        h_diffs: {label: {"added": [...], "removed": [...], ...}}
        lig_id: 3-letter ligand code
        cell_size: (width, height) per cell
        max_cols: Max columns in the grid

    Returns:
        SVG string, or None on failure.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Draw
        from rdkit.Chem.Draw import rdMolDraw2D
    except ImportError:
        return None

    labels = sorted(state_pdbs.keys())
    if not labels:
        return None

    n = len(labels)
    n_cols = min(n, max_cols)
    n_rows = (n + n_cols - 1) // n_cols
    cw, ch = cell_size
    total_w = n_cols * cw
    total_h = n_rows * (ch + 30)  # +30 for label text

    svg_parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{total_w}" height="{total_h}" '
        f'viewBox="0 0 {total_w} {total_h}">'
    ]

    for i, label in enumerate(labels):
        col = i % n_cols
        row = i // n_cols
        x_off = col * cw
        y_off = row * (ch + 30)

        pdb_path = state_pdbs[label]
        diff = h_diffs.get(label, {})
        h_added = diff.get("added", [])
        h_removed = diff.get("removed", [])

        # Load mol from per-state PDB
        mol = None
        if pdb_path and os.path.exists(str(pdb_path)):
            mol = Chem.MolFromPDBFile(str(pdb_path), removeHs=False, sanitize=False)

        if mol is None:
            # Placeholder
            svg_parts.append(
                f'<text x="{x_off + cw//2}" y="{y_off + ch//2}" '
                f'text-anchor="middle" font-size="14" fill="#999">'
                f'No PDB for {label}</text>'
            )
        else:
            try:
                AllChem.Compute2DCoords(mol)
            except Exception:
                pass

            # Find atom indices to highlight
            added_indices = []
            removed_parent_indices = []

            for atom in mol.GetAtoms():
                info = atom.GetPDBResidueInfo()
                if info:
                    name = info.GetName().strip()
                    if name in h_added:
                        added_indices.append(atom.GetIdx())

            # For removed H, find parent heavy atoms in the NEUTRAL mol
            if h_removed and "01" in state_pdbs:
                neutral_mol = Chem.MolFromPDBFile(
                    str(state_pdbs["01"]), removeHs=False, sanitize=False
                )
                if neutral_mol:
                    removed_from = diff.get("removed_from", {})
                    for h_name, parent_name in removed_from.items():
                        # Find parent in this state's mol
                        for atom in mol.GetAtoms():
                            info = atom.GetPDBResidueInfo()
                            if info and info.GetName().strip() == parent_name:
                                removed_parent_indices.append(atom.GetIdx())
                                break

            # Build highlight colors
            highlight_atoms = added_indices + removed_parent_indices
            colors = {}
            for idx in added_indices:
                colors[idx] = (0.3, 0.85, 0.3)  # green
            for idx in removed_parent_indices:
                colors[idx] = (0.9, 0.3, 0.3)  # red

            # Draw this cell
            drawer = rdMolDraw2D.MolDraw2DSVG(cw, ch)
            opts = drawer.drawOptions()
            opts.addAtomIndices = False
            opts.addStereoAnnotation = False

            if highlight_atoms:
                drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms,
                                    highlightAtomColors=colors)
            else:
                drawer.DrawMolecule(mol)
            drawer.FinishDrawing()

            cell_svg = drawer.GetDrawingText()
            # Strip the outer <svg> tags and embed as a group
            import re as _re
            inner = _re.sub(r'<\?xml[^>]*\?>', '', cell_svg)
            inner = _re.sub(r'<svg[^>]*>', '', inner)
            inner = inner.replace('</svg>', '')

            svg_parts.append(
                f'<g transform="translate({x_off},{y_off})">'
                f'{inner}</g>'
            )

        # Label below each cell
        charge_str = ""
        # Try to find charge from h_diffs structure
        n_added = len(diff.get("added", []))
        n_removed = len(diff.get("removed", []))
        if label in ("01", "00"):
            charge_str = "0"
        elif n_added > n_removed:
            charge_str = f"+{n_added - n_removed}"
        elif n_removed > n_added:
            charge_str = f"-{n_removed - n_added}"

        # H-diff annotation
        h_text = ""
        if h_added:
            h_text += f" +H: {', '.join(h_added[:3])}"
            if len(h_added) > 3:
                h_text += "..."
        if h_removed:
            h_text += f" -H: {', '.join(h_removed[:3])}"

        svg_parts.append(
            f'<text x="{x_off + cw//2}" y="{y_off + ch + 15}" '
            f'text-anchor="middle" font-size="13" font-weight="bold" '
            f'font-family="sans-serif">'
            f'{lig_id}{label} (charge: {charge_str}){h_text}</text>'
        )

    svg_parts.append('</svg>')
    return '\n'.join(svg_parts)
