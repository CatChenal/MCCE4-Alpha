"""
LangGraph Agent for MCCE4 topology file creation.

v3: Per-state PDB pipeline — each protonation state gets its own PDB,
CONNECT records, charges, and rxn calibration.

Uses a state graph where each node is an agent phase.
The LLM reasons at decision points and the graph routes accordingly.
"""

import os
import sys
import logging
from typing import Literal
from datetime import datetime

from .models import AgentState, ConformerState, make_labels_unique, sort_conformer_labels
from .config import DEFAULT_CHARGE_METHOD, DEFAULT_DIELECTRICS, CHARGE_TO_CONF
from .llm import create_llm
from .tools.mcce_tools import (
    extract_lig_id_from_pdb, fetch_rcsb_info, generate_ftpl_template,
    fill_ftpl_charges, update_conformer_params, run_rxn_calibration,
    # v3:
    generate_all_state_ftpls, merge_ftpl_files, fill_ftpl_charges_per_state,
    ensure_connect_hybridizations,
)
from .tools.rdkit_tools import (
    get_smiles_from_pdb, validate_hybridizations, mol_to_svg,
    # v3:
    get_mol_from_pdb, get_atom_info, get_ionizable_sites,
    generate_state_pdb, compute_h_diff, mol_to_svg_with_h_diff,
)
from .tools.dimorphite_tool import enumerate_protonation_states
from .tools.charge_tools import generate_charges, generate_per_state_charges, get_ionizable_sites_oe


# ──────────────────────────────────────────────────────────────────────────────
# Agent Node Functions
# ──────────────────────────────────────────────────────────────────────────────

def node_molecule_intel(state: AgentState) -> AgentState:
    """PHASE 1: Gather molecule intelligence."""
    logging.info(f"\n{'─'*60}")
    logging.info(f"  PHASE 1: Molecule Intelligence")
    logging.info(f"{'─'*60}")

    pdb_path = state["pdb_path"]
    lig_id = state.get("lig_id") or extract_lig_id_from_pdb(pdb_path)

    # Fetch RCSB info
    info = fetch_rcsb_info(lig_id)
    smiles = info.get("smiles") or get_smiles_from_pdb(pdb_path)

    state["lig_id"] = lig_id
    state["smiles"] = smiles or ""
    state["name"] = info.get("name", lig_id)
    state["formula"] = info.get("formula", "")
    state["formal_charge"] = info.get("formal_charge", 0)

    # v3: detect ionizable sites for GUI protonation editor
    # Try RDKit first; fall back to OEChem when RDKit cannot parse the molecule.
    try:
        mol = get_mol_from_pdb(pdb_path, remove_hs=False)
        ionizable = get_ionizable_sites(mol) if mol else []
        atom_info = get_atom_info(mol) if mol else []
    except Exception as e:
        logging.warning(f"  RDKit ionizable site detection failed: {e}")
        ionizable = []
        atom_info = []

    if not ionizable:
        logging.info("  RDKit found no ionizable sites — trying OEChem fallback")
        ionizable = get_ionizable_sites_oe(pdb_path)

    state["_ionizable_sites"] = ionizable
    state["_atom_info"] = atom_info
    logging.info(f"  🔍 {len(ionizable)} ionizable site(s) detected:")
    for site in ionizable:
        logging.info(f"     {site['name']} ({site['type']})")

    state["phase"] = "molecule_intel"

    return state


def node_enumerate_states(state: AgentState) -> AgentState:
    """PHASE 2: Enumerate protonation states."""
    logging.info(f"\n{'─'*60}")
    logging.info(f"  PHASE 2: Protonation State Enumeration")
    logging.info(f"{'─'*60}")

    # If user provided state PDBs, use those
    user_pdbs = state.get("user_state_pdbs", [])
    if user_pdbs:
        logging.info("  📂 Using user-supplied state PDB files")
        states = _parse_user_state_pdbs(user_pdbs, state["lig_id"])
    elif state.get("smiles"):
        states = enumerate_protonation_states(
            state["smiles"], ph=state.get("ph", 7.4),
            pdb_path=state["pdb_path"],
        )
    else:
        logging.warning("  No SMILES — defaulting to neutral")
        states = [ConformerState(label="01", charge=0, nH=0,
                                 source="fallback", rationale="Default neutral")]

    states = make_labels_unique(states)
    state["states"] = [s.to_dict() for s in states]
    state["conformer_labels"] = sort_conformer_labels([s.label for s in states])
    state["phase"] = "enumerate_states"
    return state


def node_llm_reasoning(state: AgentState) -> AgentState:
    """PHASE 3: LLM reasons about states, estimates pKa, validates chemistry."""
    logging.info(f"\n{'─'*60}")
    logging.info(f"  PHASE 3: Agent Reasoning (LLM)")
    logging.info(f"{'─'*60}")

    llm_provider = state.get("llm_provider", "gemini")
    llm_api_key = state.get("llm_api_key")
    llm = create_llm(provider=llm_provider, api_key=llm_api_key)
    if not llm.available:
        logging.info("  Skipping LLM (not configured)")
        state["phase"] = "llm_reasoning"
        return state

    llm_model_name = llm.model_name

    states = state.get("states", [])
    states_desc = "\n".join(
        f"  {s['label']}: charge={s['charge']:+d}, nH={s['nH']:+d}, smiles={s.get('smiles','')}"
        for s in states
    )

    # Include ionizable sites so LLM can reference PDB atom names
    ionizable = state.get("_ionizable_sites", [])
    sites_desc = ""
    if ionizable:
        sites_desc = "\nIonizable sites detected (PDB atom names):\n"
        for site in ionizable:
            sites_desc += (
                f"  - {site['name']} ({site['symbol']}): {site['type']}, "
                f"current H count={site['current_hs']}\n"
            )

    result = llm.ask_json(f"""You are a computational chemistry expert for MCCE simulations.

Analyze this ligand for MCCE topology file creation:
  Name: {state.get('name', '')}
  Formula: {state.get('formula', '')}
  SMILES: {state.get('smiles', '')}
{sites_desc}
Proposed protonation states from Dimorphite-DL:
{states_desc}

IMPORTANT RULES:
  - Use MCCE label conventions: "01" for neutral, "+1" for protonated,
    "-1" for deprotonated. Do NOT use "02", "03", "04" etc.
  - For the site_atom field, use the EXACT PDB atom name from the
    ionizable sites list above (e.g., "N19" not "morpholine nitrogen").
  - For proton_exchange, specify which PDB atom name gains/loses a proton
    (e.g., "+H on N19" or "-H from O20").

For each state, provide:
  1. Estimated pKa and the PDB atom name of the protonation site
  2. EXACTLY which proton(s) are added or removed vs the neutral form
     using the PDB atom name (e.g., "+H on N19 (morpholine N)")
  3. Literature references: cite academic papers, textbooks, or trusted
     databases (PubChem, DrugBank, IUPHAR). Provide DOI or URL.

Do NOT invent new protonation states beyond what Dimorphite-DL proposed.
You may remove states that are chemically unreasonable, but do not add new ones.
Flag any warnings (metals, unusual chemistry, etc.).

Respond ONLY in JSON (no markdown fences):
{{
  "states": [
    {{
      "label": "01", "charge": 0, "nH": 0, "pka": null,
      "site_atom": null,
      "rationale": "neutral form — reference state",
      "proton_exchange": "none (neutral reference)",
      "references": ["PubChem CID 12345"]
    }},
    {{
      "label": "+1", "charge": 1, "nH": 1, "pka": 8.5,
      "site_atom": "N19",
      "rationale": "protonation at morpholine nitrogen",
      "proton_exchange": "+H on N19 (morpholine nitrogen)",
      "references": ["doi:10.xxxx/yyyy"]
    }}
  ],
  "warnings": [],
  "summary": "brief analysis"
}}""")

    if result and "states" in result:
        # Build lookup of original Dimorphite-DL data by label and by SMILES
        orig_states = state.get("states", [])
        smiles_by_label = {}
        smiles_by_charge = {}  # charge → list of SMILES
        for os_dict in orig_states:
            smi = os_dict.get("smiles", "")
            lbl = os_dict.get("label", "")
            chg = os_dict.get("charge", 0)
            if smi:
                smiles_by_label[lbl] = smi
                smiles_by_charge.setdefault(chg, []).append(smi)

        refined = []
        for s in result["states"]:
            charge = int(s.get("charge", 0) or 0)

            # Assign standard MCCE label based on charge
            # make_labels_unique will handle duplicates later (+1 → +a/+b)
            from .config import CHARGE_TO_CONF
            label = CHARGE_TO_CONF.get(charge, f"{charge:+d}"[-2:])

            # Preserve SMILES from Dimorphite-DL — LLM rarely returns them
            llm_smiles = s.get("smiles", "")
            preserved_smiles = llm_smiles
            if not preserved_smiles:
                # Try to find matching SMILES from Dimorphite by charge
                charge_smiles = smiles_by_charge.get(charge, [])
                if charge_smiles:
                    preserved_smiles = charge_smiles.pop(0)  # consume one
                else:
                    # Try label-based lookup
                    for lbl_key, smi_val in smiles_by_label.items():
                        if smi_val:
                            preserved_smiles = smi_val
                            break

            refined.append({
                "label": label, "charge": charge,
                "nH": int(s.get("nH", 0) or 0), "pka": s.get("pka"),
                "site_atom": s.get("site_atom"),
                "smiles": preserved_smiles,
                "source": f"llm ({llm_model_name})", "pdb_path": None,
                "rationale": s.get("rationale", ""),
                "proton_exchange": s.get("proton_exchange", ""),
                "llm_model": llm_model_name,
                "references": s.get("references", []),
                "h_added": [], "h_removed": [],
            })
            if preserved_smiles:
                logging.info(f"  State {label} (charge={charge:+d}): SMILES={preserved_smiles[:60]}...")
            else:
                logging.warning(f"  State {label} (charge={charge:+d}): NO SMILES — per-state PDB may use site-based generation")

        # Apply naming convention: single per charge → 01/+1/-1, multiple → +a/+b
        refined = make_labels_unique(refined)

        # Ensure neutral (01) state is always present
        has_neutral = any(
            int(s.get("charge", 0) or 0) == 0 for s in refined
        )
        if not has_neutral:
            logging.info("  Adding neutral (01) state — LLM did not include one")
            neutral_smiles = smiles_by_charge.get(0, [""])[0] if smiles_by_charge.get(0) else ""
            refined.insert(0, {
                "label": "01", "charge": 0, "nH": 0, "pka": None,
                "site_atom": None, "smiles": neutral_smiles,
                "source": "auto", "pdb_path": None,
                "rationale": "Neutral reference state (auto-added)",
                "proton_exchange": "none (neutral reference)",
                "llm_model": "", "references": [],
                "h_added": [], "h_removed": [],
            })

        state["states"] = refined
        state["conformer_labels"] = sort_conformer_labels([s["label"] for s in refined])

        if result.get("warnings"):
            warnings = state.get("warnings", [])
            warnings.extend(result["warnings"])
            state["warnings"] = warnings
            for w in result["warnings"]:
                logging.warning(f"  ⚠ LLM: {w}")

        logging.info(f"  🧠 {result.get('summary', 'Analysis complete')}")

    state["phase"] = "llm_reasoning"
    return state


def _read_pdb_atom_names(pdb_path: str) -> set:
    """Return the set of atom names (stripped) present in a PDB file."""
    names = set()
    try:
        with open(pdb_path) as f:
            for line in f:
                if line.startswith(("ATOM  ", "HETATM")):
                    names.add(line[12:16].strip())
    except Exception:
        pass
    return names


def node_generate_state_pdbs(state: AgentState) -> AgentState:
    """PHASE 4: Generate per-state PDBs with correct H atoms.

    For each protonation state:
      - Start from the PREVIOUS state's PDB (incremental protonation)
      - Validate LLM site_atom against actual PDB atom names
      - Fall back to detected ionizable sites when LLM gives a wrong name
      - Write per-state PDB (e.g., EMH_01.pdb, EMH_p1.pdb)
      - Compute H-diff vs neutral

    Incremental approach: +2 is built from +1, not from neutral.
    This ensures each state correctly carries forward all previous protonations.
    """
    logging.info(f"\n{'─'*60}")
    logging.info(f"  PHASE 4: Per-State PDB Generation")
    logging.info(f"{'─'*60}")

    pdb_path = state["pdb_path"]
    lig_id = state.get("lig_id", "LIG")
    work_dir = state.get("work_dir", ".")

    state_pdbs_dir = os.path.join(work_dir, f"{lig_id}_state_pdbs")
    os.makedirs(state_pdbs_dir, exist_ok=True)

    state_pdbs = {}
    h_diffs = {}
    states = state.get("states", [])

    # Trust Dimorphite-DL + LLM reasoning for protonation states.
    # Do NOT auto-add a -1 state — use only what was enumerated and validated.

    # Sort states by charge so we build each from the previous one
    sorted_states = sorted(states, key=lambda s: int(s.get('charge', 0) if isinstance(s, dict) else s.charge))

    # sites_used_per_base: {base_label → [site_atoms used]} for same-base disambiguation
    sites_used_per_base: dict = {}

    for s in sorted_states:
        sd = s if isinstance(s, dict) else s.to_dict()
        label = sd["label"]
        smiles = sd.get("smiles", "")
        charge = int(sd.get('charge', 0) or 0)
        llm_site_atom = (sd.get("site_atom") or "").strip()

        # Determine action from charge
        if charge > 0:
            action = "protonate"
        elif charge < 0:
            action = "deprotonate"
        else:
            action = None

        # Neutral state (charge 0): just copy the neutral PDB
        if charge == 0:
            logging.info(f"\n  Generating PDB for state {label} (neutral)...")
            pdb_out, added, removed = generate_state_pdb(
                neutral_pdb=pdb_path, state_smiles="",
                lig_id=lig_id, label=label, output_dir=state_pdbs_dir,
            )
            if pdb_out:
                state_pdbs[label] = pdb_out
                if isinstance(s, dict):
                    s["pdb_path"] = pdb_out
                    s["h_added"] = []; s["h_removed"] = []; s["nH"] = 0
                    s["proton_exchange"] = "none (neutral reference)"
                h_diffs[label] = {"added": [], "removed": [], "added_on": {}, "removed_from": {}}
            continue

        # ── Non-neutral state: determine base PDB and site ──

        # Base PDB: use charge±1 state if available, else neutral
        prev_charge = charge - 1 if charge > 0 else charge + 1
        base_label = next((
            l for l, p in state_pdbs.items()
            if int(next((x.get('charge', 0) for x in sorted_states
                         if (x.get('label') if isinstance(x, dict) else x.label) == l), 0)) == prev_charge
        ), None)
        base_pdb = state_pdbs.get(base_label, pdb_path)
        base_key  = base_label or "neutral"

        # Detect available ionizable sites from the SPECIFIC base PDB using OEChem.
        # This correctly handles multi-step protonation: a base that is already
        # singly protonated (e.g., +a) will show fewer available tertiary amines
        # than the original neutral, so the next protonation goes to the right atom.
        base_oe_sites = get_ionizable_sites_oe(base_pdb)
        if action == "protonate":
            candidates = [
                s2 for s2 in base_oe_sites
                if s2.get("type") == "protonatable_N" and s2.get("current_hs", 0) == 0
            ]
        else:  # deprotonate
            candidates = [
                s2 for s2 in base_oe_sites
                if s2.get("current_hs", 0) > 0
            ]

        base_atom_names = _read_pdb_atom_names(base_pdb)
        site_atom = None

        # 1) Validate LLM suggestion against BOTH base PDB and OEChem candidates.
        #    The LLM often confuses SMILES numbering with PDB atom names
        #    (e.g. saying "N19" for nitrile N when it means N25 the piperidine).
        #    We require the suggested atom to actually be an ionizable candidate.
        candidate_names = {c["name"] for c in candidates}
        if llm_site_atom and llm_site_atom in base_atom_names and llm_site_atom in candidate_names:
            site_atom = llm_site_atom
            logging.info(f"  ✓ LLM site_atom '{site_atom}' validated against base PDB "
                         f"and OEChem ionizable sites")
        else:
            if llm_site_atom:
                logging.warning(f"  ⚠ LLM site_atom '{llm_site_atom}' not in OEChem "
                                 f"ionizable candidates {candidate_names} — "
                                 "using OEChem-detected sites")
            # 2) From OEChem candidates, pick the first not yet used for this base
            used_for_base = sites_used_per_base.get(base_key, [])
            for cand in candidates:
                if cand["name"] not in used_for_base and cand["name"] in base_atom_names:
                    site_atom = cand["name"]
                    logging.info(f"  Using OEChem-detected site: {site_atom} "
                                 f"(from base {base_key})")
                    break

        if not site_atom:
            logging.warning(f"  No valid site_atom for state {label} (charge={charge:+d}) "
                             f"— skipping (will attempt charge assignment from existing PDB)")
            warns = state.get("warnings", [])
            warns.append(f"No valid protonation site for state {label} — state PDB skipped")
            state["warnings"] = warns
            continue

        logging.info(f"\n  Generating PDB for state {label} (charge={charge:+d}, "
                     f"site={site_atom}, base={os.path.basename(base_pdb)})...")

        pdb_out, added, removed = generate_state_pdb(
            neutral_pdb=base_pdb,
            state_smiles=smiles,
            lig_id=lig_id,
            label=label,
            output_dir=state_pdbs_dir,
            site_atom=site_atom,
            action=action,
        )

        if pdb_out:
            sites_used_per_base.setdefault(base_key, []).append(site_atom)
            state_pdbs[label] = pdb_out
            # nH is vs neutral (use h_diff count, not just added-removed for this step)
            if isinstance(s, dict):
                s["pdb_path"] = pdb_out
                s["h_added"] = added
                s["h_removed"] = removed
                s["nH"] = charge  # nH = charge for simple protonation
            else:
                s.pdb_path = pdb_out
                s.h_added = added
                s.h_removed = removed
                s.nH = charge

            # H-diff vs neutral
            if "01" in state_pdbs:
                diff = compute_h_diff(state_pdbs["01"], pdb_out)
                h_diffs[label] = diff
                logging.info(f"    H-diff vs neutral: +{diff['added']}  -{diff['removed']}")

                exchange_parts = []
                for h_name in diff.get("added", []):
                    parent = diff.get("added_on", {}).get(h_name, "?")
                    exchange_parts.append(f"+{h_name} on {parent}")
                for h_name in diff.get("removed", []):
                    parent = diff.get("removed_from", {}).get(h_name, "?")
                    exchange_parts.append(f"-{h_name} from {parent}")
                proton_ex = "; ".join(exchange_parts) if exchange_parts else ""
                if isinstance(s, dict):
                    s["proton_exchange"] = proton_ex
                else:
                    s.proton_exchange = proton_ex
            else:
                h_diffs[label] = {"added": [], "removed": [], "added_on": {}, "removed_from": {}}
        else:
            warns = state.get("warnings", [])
            warns.append(f"Failed to generate PDB for state {label} — will use existing PDB if available")
            state["warnings"] = warns

    # Re-sort state_pdbs back to original label order
    state["state_pdbs"] = state_pdbs
    state["h_diffs"] = h_diffs

    logging.info(f"\n  Generated {len(state_pdbs)}/{len(states)} state PDBs")

    # Generate PyMOL comparison script
    if len(state_pdbs) > 1:
        try:
            from .tools.mcce_tools import generate_pymol_script
            pymol_script = os.path.join(
                state.get("work_dir", "."),
                f"{lig_id}_compare_states.pml"
            )
            generate_pymol_script(state_pdbs, h_diffs, lig_id, pymol_script)
            state["pymol_script"] = pymol_script
        except Exception as e:
            logging.warning(f"  PyMOL script generation failed: {e}")

    state["phase"] = "generate_state_pdbs"
    return state


def node_user_review(state: AgentState) -> AgentState:
    """Pause for user review (GUI or terminal) between Phase 4 and Phase 5.

    This node sets needs_user_review=True.
    The GUI/CLI layer handles the actual review and sets user_approved.
    """
    logging.info(f"\n{'─'*60}")
    logging.info(f"  Awaiting User Review (between Phase 4 & 5)")
    logging.info(f"{'─'*60}")

    state["needs_user_review"] = True
    state["phase"] = "user_review"
    return state


def node_generate_template(state: AgentState) -> AgentState:
    """PHASE 5: Generate .ftpl template and validate hybridizations.

    v3: If state_pdbs are available, call pdb2ftpl.py once per state PDB
    then merge. Falls back to v2 single-call if no state_pdbs.
    """
    logging.info(f"\n{'─'*60}")
    logging.info(f"  PHASE 5: Template Generation & Validation")
    logging.info(f"{'─'*60}")

    lig_id = state["lig_id"]
    pdb_path = state["pdb_path"]
    labels = state["conformer_labels"]
    ftpl_path = state.get("ftpl_path", f"{lig_id}.ftpl")
    state_pdbs = state.get("state_pdbs", {})

    if state_pdbs:
        # ── v3 path: per-state ftpl generation & merge ──
        logging.info("  Using per-state PDB pipeline (v3)")
        work_dir = state.get("work_dir", ".")
        ftpl_dir = os.path.join(work_dir, f"{lig_id}_ftpls")
        os.makedirs(ftpl_dir, exist_ok=True)

        per_state_ftpls = generate_all_state_ftpls(lig_id, state_pdbs, ftpl_dir)

        if not per_state_ftpls:
            state["errors"] = state.get("errors", []) + [
                "No ftpl files generated — pdb2ftpl.py may not be available"
            ]
            state["complete"] = True
            return state

        merged_path = os.path.join(work_dir, ftpl_path)
        success = merge_ftpl_files(per_state_ftpls, lig_id, merged_path)
        if not success:
            state["errors"] = state.get("errors", []) + ["ftpl merge failed"]
            state["complete"] = True
            return state

        state["ftpl_path"] = merged_path
        state["per_state_connects"] = per_state_ftpls
    else:
        # ── v2 fallback: single pdb2ftpl call with all labels ──
        logging.info("  Using single-PDB pipeline (v2 fallback)")
        unfilled = generate_ftpl_template(lig_id, pdb_path, labels, ftpl_path)
        if unfilled < 0:
            state["errors"] = state.get("errors", []) + ["pdb2ftpl.py failed"]
            state["complete"] = True
            return state
        state["ftpl_path"] = ftpl_path

    # Ensure every CONNECT line has a valid hybridization (s, sp, sp2, sp3)
    try:
        n_fixed = ensure_connect_hybridizations(
            state["ftpl_path"], pdb_path=pdb_path, lig_id=lig_id
        )
        if n_fixed:
            logging.info(f"  Fixed {n_fixed} CONNECT line(s) with missing/invalid hybridization")
    except Exception as e:
        logging.warning(f"  Hybridization enforcement skipped: {e}")

    # Hybridization validation
    if state.get("smiles"):
        try:
            with open(state["ftpl_path"]) as f:
                issues = validate_hybridizations(pdb_path, f.read(), lig_id)
            state["hybridization_issues"] = [
                {"atom": a, "ftpl": ft, "rdkit": r} for a, ft, r in issues
            ]
        except Exception as e:
            logging.warning(f"  Hybridization validation skipped: {e}")

    state["phase"] = "generate_template"
    return state


def node_assign_charges(state: AgentState) -> AgentState:
    """PHASE 6: Generate and fill atomic charges.

    v3: If state_pdbs are available, compute charges per state PDB
    and fill per-conformer. Falls back to v2 single-charge if no state_pdbs.
    """
    logging.info(f"\n{'─'*60}")
    logging.info(f"  PHASE 6: Charge Assignment")
    logging.info(f"{'─'*60}")

    lig_id = state["lig_id"]
    method = state.get("charge_method", DEFAULT_CHARGE_METHOD)
    state_pdbs = dict(state.get("state_pdbs", {}))  # copy — we may augment it

    # Augment state_pdbs with any existing per-state PDB files not yet registered.
    # State PDB filenames follow: {lig_id}_{label.replace('+','p').replace('-','m')}.pdb
    work_dir = state.get("work_dir", ".")
    state_pdbs_dir = os.path.join(work_dir, f"{lig_id}_state_pdbs")
    for label in state.get("conformer_labels", []):
        if label not in state_pdbs:
            safe = label.replace('+', 'p').replace('-', 'm')
            candidate = os.path.join(state_pdbs_dir, f"{lig_id}_{safe}.pdb")
            if os.path.exists(candidate):
                logging.info(f"  Found existing state PDB for {label}: {candidate}")
                state_pdbs[label] = candidate

    if state_pdbs:
        # ── v3 path: per-state charges ──
        logging.info("  Using per-state charge pipeline (v3)")
        per_state_charges = generate_per_state_charges(state_pdbs, method, lig_id)
        state["per_state_charges"] = per_state_charges

        if state.get("ftpl_path"):
            unfilled = fill_ftpl_charges_per_state(
                state["ftpl_path"], per_state_charges, lig_id
            )
            if unfilled > 0:
                state["warnings"] = state.get("warnings", []) + [
                    f"{unfilled} charge entries still unfilled"
                ]
    else:
        # ── v2 fallback: single charge calculation ──
        charges = generate_charges(state["pdb_path"], method, lig_id)
        if not charges:
            state["errors"] = state.get("errors", []) + ["No charges generated"]
            state["complete"] = True
            return state
        state["charges"] = charges
        unfilled = fill_ftpl_charges(state["ftpl_path"], charges, lig_id)
        if unfilled > 0:
            state["warnings"] = state.get("warnings", []) + [
                f"{unfilled} charge entries unfilled"
            ]

    # Update conformer params (nH, pKa0)
    update_conformer_params(state["ftpl_path"], state["states"], lig_id)

    state["phase"] = "assign_charges"
    return state


def node_rxn_calibration(state: AgentState) -> AgentState:
    """PHASE 7: RXN calibration via MCCE4 steps 1-3."""
    logging.info(f"\n{'─'*60}")
    logging.info(f"  PHASE 7: RXN Calibration")
    logging.info(f"{'─'*60}")

    results = run_rxn_calibration(
        state["lig_id"], state["ftpl_path"], state["pdb_path"],
        state["conformer_labels"],
        state.get("dielectrics", DEFAULT_DIELECTRICS),
        state.get("work_dir", ".")
    )

    state["rxn_values"] = results
    state["phase"] = "rxn_calibration"
    state["complete"] = True
    return state


def node_done(state: AgentState) -> AgentState:
    """Final node — print summary."""
    logging.info(f"\n{'='*60}")
    logging.info(f"  🤖 MCCE4 Topology Agent — COMPLETE")
    logging.info(f"{'='*60}")
    logging.info(f"  ✅ Topology file: {state.get('ftpl_path', '?')}")
    logging.info(f"  🧪 Conformers:    {state.get('conformer_labels', [])}")
    logging.info(f"  ⚡ Charges:       {state.get('charge_method', '?')}")
    if state.get("errors"):
        for e in state["errors"]:
            logging.error(f"  ❌ {e}")
    if state.get("warnings"):
        for w in state["warnings"]:
            logging.warning(f"  ⚠ {w}")
    logging.info(f"{'='*60}")
    state["complete"] = True
    return state


# ──────────────────────────────────────────────────────────────────────────────
# Graph Construction
# ──────────────────────────────────────────────────────────────────────────────

def should_review(state: AgentState) -> Literal["user_review", "generate_template"]:
    """Decide whether to pause for user review."""
    if state.get("needs_user_review"):
        return "user_review"
    return "generate_template"


def after_review(state: AgentState) -> Literal["generate_template", "done"]:
    """After user review, proceed or abort."""
    if state.get("user_approved", True):
        return "generate_template"
    return "done"


def should_calibrate(state: AgentState) -> Literal["rxn_calibration", "done"]:
    """Decide whether to run RXN calibration."""
    if state.get("errors"):
        return "done"
    dry_run = state.get("dry_run", False)
    if dry_run:
        logging.info("  ⏩ Dry run — skipping calibration")
        return "done"
    return "rxn_calibration"


def build_agent_graph(use_gui: bool = False):
    """Build the LangGraph StateGraph for the agent.

    v3: Includes node_generate_state_pdbs between llm_reasoning and
    generate_template.

    Args:
        use_gui: If True, the graph includes a user_review node.

    Returns:
        Compiled LangGraph graph.
    """
    try:
        from langgraph.graph import StateGraph, END
    except ImportError:
        logging.error("langgraph not installed — run: pip install langgraph")
        sys.exit(1)

    graph = StateGraph(AgentState)

    # Add nodes
    graph.add_node("molecule_intel", node_molecule_intel)
    graph.add_node("enumerate_states", node_enumerate_states)
    graph.add_node("llm_reasoning", node_llm_reasoning)
    graph.add_node("generate_state_pdbs", node_generate_state_pdbs)  # Phase 4
    graph.add_node("generate_template", node_generate_template)
    graph.add_node("assign_charges", node_assign_charges)
    graph.add_node("rxn_calibration", node_rxn_calibration)
    graph.add_node("done", node_done)

    if use_gui:
        graph.add_node("user_review", node_user_review)

    # Set entry point
    graph.set_entry_point("molecule_intel")

    # Add edges
    graph.add_edge("molecule_intel", "enumerate_states")
    graph.add_edge("enumerate_states", "llm_reasoning")
    graph.add_edge("llm_reasoning", "generate_state_pdbs")  # Phase 3 → Phase 4

    if use_gui:
        graph.add_conditional_edges("generate_state_pdbs", should_review,
                                     {"user_review": "user_review",
                                      "generate_template": "generate_template"})
        graph.add_conditional_edges("user_review", after_review,
                                     {"generate_template": "generate_template",
                                      "done": "done"})
    else:
        graph.add_edge("generate_state_pdbs", "generate_template")

    graph.add_edge("generate_template", "assign_charges")
    graph.add_conditional_edges("assign_charges", should_calibrate,
                                 {"rxn_calibration": "rxn_calibration",
                                  "done": "done"})
    graph.add_edge("rxn_calibration", "done")
    graph.add_edge("done", END)

    return graph.compile()


# ──────────────────────────────────────────────────────────────────────────────
# Runner
# ──────────────────────────────────────────────────────────────────────────────

def run_agent(pdb_path: str, use_gui: bool = False,
               charge_method: str = DEFAULT_CHARGE_METHOD,
               dielectrics: list = None, ph: float = 7.4,
               work_dir: str = ".", dry_run: bool = False,
               user_state_pdbs: list = None,
               output: str = None,
               llm_provider: str = "gemini",
               api_key: str = None) -> AgentState:
    """Run the full agent pipeline.

    Args:
        pdb_path: Path to the ligand PDB file.
        use_gui: Whether to launch Streamlit GUI for review.
        charge_method: Charge calculation method.
        dielectrics: List of dielectric constants.
        ph: Target pH for protonation.
        work_dir: MCCE4 working directory.
        dry_run: Skip RXN calibration.
        user_state_pdbs: User-supplied state PDB paths.
        output: Output .ftpl filename.
        llm_provider: LLM provider ("gemini", "claude", or "chatgpt").
        api_key: Optional API key for the LLM provider.

    Returns:
        Final AgentState dict.
    """
    lig_id = extract_lig_id_from_pdb(pdb_path)

    initial_state: AgentState = {
        "pdb_path": os.path.abspath(pdb_path),
        "lig_id": lig_id,
        "ph": ph,
        "charge_method": charge_method,
        "dielectrics": dielectrics or DEFAULT_DIELECTRICS,
        "work_dir": os.path.abspath(work_dir),
        "ftpl_path": output or f"{lig_id}.ftpl",
        "user_state_pdbs": user_state_pdbs or [],
        "states": [],
        "conformer_labels": [],
        "errors": [],
        "warnings": [],
        "messages": [],
        "needs_user_review": use_gui,
        "user_approved": not use_gui,  # auto-approve if no GUI
        "complete": False,
        "dry_run": dry_run,
        "smiles": "",
        "name": "",
        "formula": "",
        "formal_charge": 0,
        "charges": {},
        "hybridization_issues": [],
        "charge_issues": [],
        "rxn_values": {},
        "phase": "init",
        # v3 fields:
        "state_pdbs": {},
        "h_diffs": {},
        "per_state_charges": {},
        "per_state_connects": {},
        # LLM configuration
        "llm_provider": llm_provider,
        "llm_api_key": api_key,
    }

    # Build and run graph
    graph = build_agent_graph(use_gui=use_gui)

    if use_gui:
        # For GUI mode: run phases 1-3b directly (no LangGraph streaming issues)
        # Then pause for user review. GUI calls resume_agent() after approval.
        state = initial_state
        state = node_molecule_intel(state)
        state = node_enumerate_states(state)
        state = node_llm_reasoning(state)
        state = node_generate_state_pdbs(state)  # v3: generate per-state PDBs
        state = node_user_review(state)
        return state
    else:
        # CLI mode: run full pipeline via LangGraph
        final_state = graph.invoke(initial_state)
        return final_state


def resume_agent(state: AgentState) -> AgentState:
    """Resume agent after user review (GUI sets user_approved=True).

    Runs phases 5-7: template generation, charges, and RXN calibration.
    """
    state["needs_user_review"] = False
    state["user_approved"] = True

    state = node_generate_template(state)
    if state.get("errors"):
        state = node_done(state)
        return state

    state = node_assign_charges(state)
    if state.get("errors"):
        state = node_done(state)
        return state

    if not state.get("dry_run", False):
        state = node_rxn_calibration(state)

    state = node_done(state)
    return state


# ──────────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────────

def _parse_user_state_pdbs(pdb_paths: list, lig_id: str) -> list:
    """Parse user-supplied state PDB filenames into ConformerState objects."""
    import re
    from pathlib import Path

    states = []
    for sp in pdb_paths:
        stem = Path(sp).stem.upper()
        lig_upper = lig_id.upper()
        label = None

        # Try EMH_01, EMH_+1, EMH01 patterns
        pat = re.compile(re.escape(lig_upper) + r'[_\-]?([0+-]\d*|[+-]\d+)$')
        m = pat.match(stem)
        if m:
            label = m.group(1)
        elif stem.startswith(lig_upper):
            suffix = stem[len(lig_upper):].lstrip("_-")
            if suffix:
                label = suffix

        if not label:
            logging.warning(f"  Could not parse state from: {sp}")
            continue

        # Parse charge from label
        if label.startswith("+"):
            charge = int(label)
        elif label.startswith("-"):
            charge = int(label)
        elif label == "01":
            charge = 0
        else:
            charge = 0

        states.append(ConformerState(
            label=label, charge=charge, nH=charge,
            pdb_path=os.path.abspath(sp), source="user",
            rationale=f"User PDB: {os.path.basename(sp)}"
        ))
        logging.info(f"     {label}: {sp} (charge={charge:+d})")

    return make_labels_unique(states)
