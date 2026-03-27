"""
Streamlit Web GUI for the MCCE4 Topology Agent.

Features:
  - 2D molecule depiction with protonation site highlighting
  - Conformer state table with add/edit/remove
  - Molecule editor via Ketcher (if streamlit-ketcher installed)
  - Runs in browser — works over SSH tunnel, no X11 needed

Launch:
  streamlit run mcce4_agent_ftpl/gui/app.py -- EMH.pdb
  # Or via the CLI:
  mcce4_agent_ftpl.py EMH.pdb --gui
"""

import os
import sys
import json
import base64
import logging
import tempfile
from pathlib import Path

# Ensure package is importable
PACKAGE_DIR = str(Path(__file__).resolve().parent.parent.parent)
if PACKAGE_DIR not in sys.path:
    sys.path.insert(0, PACKAGE_DIR)

import streamlit as st

from mcce4_agent_ftpl.models import ConformerState, AgentState
from mcce4_agent_ftpl.config import (
    GUI_TITLE, SUPPORTED_CHARGE_METHODS, DEFAULT_CHARGE_METHOD, DEFAULT_DIELECTRICS
)
from mcce4_agent_ftpl.tools.rdkit_tools import mol_to_svg, get_mol_from_pdb, get_mol_from_smiles
from mcce4_agent_ftpl.tools.mcce_tools import extract_lig_id_from_pdb


# ──────────────────────────────────────────────────────────────────────────────
# Page config
# ──────────────────────────────────────────────────────────────────────────────
st.set_page_config(
    page_title=GUI_TITLE,
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)


# ──────────────────────────────────────────────────────────────────────────────
# Session state init
# ──────────────────────────────────────────────────────────────────────────────
def init_session():
    """Initialize session state variables."""
    defaults = {
        "agent_state": None,
        "pdb_path": None,
        "lig_id": None,
        "states": [],
        "phase": "upload",
        "approved": False,
        "running": False,
        "analysis_done": False,
        "log_output": "",
    }
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v


init_session()


# ──────────────────────────────────────────────────────────────────────────────
# Sidebar
# ──────────────────────────────────────────────────────────────────────────────
with st.sidebar:
    st.title("🤖 MCCE4 Agent")
    st.markdown("---")

    # Settings
    st.subheader("⚙️ Settings")
    ph = st.slider("Target pH", 0.0, 14.0, 7.4, 0.1)
    charge_method = st.selectbox("Charge Method", SUPPORTED_CHARGE_METHODS,
                                  index=SUPPORTED_CHARGE_METHODS.index(DEFAULT_CHARGE_METHOD))

    dielectric_opts = st.multiselect("Dielectric Constants", [2, 4, 8], default=[2, 4, 8])
    dry_run = st.checkbox("Dry Run (skip RXN calibration)", value=False)

    st.markdown("---")
    st.subheader("📋 Agent Status")
    phase = st.session_state.get("phase", "upload")
    phase_info = {
        "upload": ("🔵", "Waiting for input"),
        "analyzing": ("🟡", "Analyzing molecule..."),
        "review": ("🟠", "Review conformer states"),
        "running": ("🟢", "Generating topology file..."),
        "complete": ("✅", "Complete!"),
        "error": ("🔴", "Error occurred"),
    }
    icon, desc = phase_info.get(phase, ("⚪", phase))
    st.write(f"{icon} **{desc}**")

    # Show conformer count if available
    states = st.session_state.get("states", [])
    if states:
        st.caption(f"{len(states)} conformer state(s) proposed")
    lig_id = st.session_state.get("lig_id")
    if lig_id:
        st.caption(f"Ligand: {lig_id}")


# ──────────────────────────────────────────────────────────────────────────────
# Main content
# ──────────────────────────────────────────────────────────────────────────────
st.title("🧬 MCCE4 Topology File Agent")
st.markdown("Create MCCE4 topology files (.ftpl) with AI-powered protonation state analysis")

# ── Tab layout (v3: merged Editor into Conformer States) ──
tab_input, tab_states, tab_output = st.tabs([
    "📂 Input", "🧪 Conformer States & Editor", "📄 Output"
])

# ──────────────────────────────────────────────────────────────────────────────
# TAB 1: Input
# ──────────────────────────────────────────────────────────────────────────────
with tab_input:
    col1, col2 = st.columns([1, 1])

    with col1:
        st.subheader("Upload Ligand PDB")

        uploaded_pdb = st.file_uploader("Ligand PDB file", type=["pdb"],
                                         key="pdb_upload")
        if uploaded_pdb:
            # Save to temp file
            tmp_dir = tempfile.mkdtemp()
            pdb_path = os.path.join(tmp_dir, uploaded_pdb.name)
            with open(pdb_path, "wb") as f:
                f.write(uploaded_pdb.getvalue())
            st.session_state["pdb_path"] = pdb_path
            st.session_state["lig_id"] = extract_lig_id_from_pdb(pdb_path)
            st.success(f"✓ Loaded: {uploaded_pdb.name} (Ligand: {st.session_state['lig_id']})")

        # Or specify path directly
        pdb_path_input = st.text_input("Or enter PDB file path on server:",
                                        placeholder="/path/to/EMH.pdb")
        if pdb_path_input and os.path.exists(pdb_path_input):
            st.session_state["pdb_path"] = pdb_path_input
            st.session_state["lig_id"] = extract_lig_id_from_pdb(pdb_path_input)
            st.success(f"✓ Found: {pdb_path_input} (Ligand: {st.session_state['lig_id']})")

        st.markdown("---")
        st.subheader("Optional: State-specific PDBs")
        st.caption("Upload PDB files for specific protonation states (e.g., EMH_01.pdb, EMH_+1.pdb)")
        state_pdbs = st.file_uploader("State PDB files", type=["pdb"],
                                       accept_multiple_files=True, key="state_pdbs")

    with col2:
        st.subheader("Molecule Preview")
        if st.session_state.get("pdb_path"):
            svg = mol_to_svg(st.session_state["pdb_path"], size=(500, 400))
            if svg:
                st.image(svg, use_container_width=True)
            else:
                st.info("RDKit not available for 2D preview. Install: `pip install rdkit`")
        else:
            st.info("Upload or specify a PDB file to see the molecule.")

    # ── Analyze button ──
    if st.session_state.get("pdb_path"):
        if st.button("🧠 Analyze Protonation States", type="primary", use_container_width=True):
            st.session_state["phase"] = "analyzing"

            from mcce4_agent_ftpl.agent import run_agent

            # Show real-time progress below the button
            status = st.status("🤖 Agent is analyzing the molecule...", expanded=True)

            with status:
                st.write("**Phase 1:** Fetching molecule info from RCSB...")

                # Capture log output for display
                import io
                log_capture = io.StringIO()
                log_handler = logging.StreamHandler(log_capture)
                log_handler.setLevel(logging.INFO)
                log_handler.setFormatter(logging.Formatter("%(message)s"))
                logging.getLogger().addHandler(log_handler)

                agent_state = run_agent(
                    st.session_state["pdb_path"],
                    use_gui=True,
                    charge_method=charge_method,
                    dielectrics=dielectric_opts,
                    ph=ph,
                    dry_run=dry_run,
                )

                # Remove the capture handler
                logging.getLogger().removeHandler(log_handler)

                # Show what happened
                log_text = log_capture.getvalue()
                if agent_state.get("smiles"):
                    st.write(f"✅ **SMILES:** `{agent_state['smiles'][:80]}...`"
                             if len(agent_state.get('smiles', '')) > 80
                             else f"✅ **SMILES:** `{agent_state.get('smiles', '')}`")
                if agent_state.get("name"):
                    st.write(f"✅ **Name:** {agent_state['name']}")

                st.write("**Phase 2:** Enumerating protonation states...")
                states = agent_state.get("states", [])
                if states:
                    st.write(f"✅ Found **{len(states)}** conformer state(s):")
                    for s in states:
                        label = s.get('label', '?')
                        charge = int(s.get('charge', 0) or 0)
                        pka = s.get('pka')
                        pka_str = ""
                        if pka is not None:
                            try:
                                # LLM may return "~8.7", ">10", "≈4.5" etc.
                                pka_clean = str(pka).lstrip("~≈><≥≤ ")
                                pka_str = f", pKa≈{float(pka_clean):.1f}"
                            except (ValueError, TypeError):
                                pka_str = f", pKa={pka}"
                        source = s.get('source', '?')
                        st.write(f"   • **{label}** (charge={charge:+d}{pka_str}) — {source}")

                # v3: show per-state PDB generation results
                state_pdbs = agent_state.get("state_pdbs", {})
                if state_pdbs:
                    st.write(f"**Phase 3b:** Generated **{len(state_pdbs)}** per-state PDB(s)")
                    h_diffs = agent_state.get("h_diffs", {})
                    for label, pdb_path in state_pdbs.items():
                        diff = h_diffs.get(label, {})
                        added = diff.get("added", [])
                        removed = diff.get("removed", [])
                        diff_str = ""
                        if added:
                            diff_str += f" +H: {', '.join(added)}"
                        if removed:
                            diff_str += f" -H: {', '.join(removed)}"
                        st.write(f"   • **{label}**: `{os.path.basename(pdb_path)}`{diff_str}")

                if agent_state.get("warnings"):
                    st.write("**⚠️ Warnings:**")
                    for w in agent_state["warnings"]:
                        st.warning(w)

                # Show captured log in an expander
                if log_text.strip():
                    with st.expander("📋 Full agent log", expanded=False):
                        st.code(log_text, language="text")

            status.update(label="✅ Analysis complete — switch to **Conformer States** tab to review",
                          state="complete", expanded=False)

            st.session_state["agent_state"] = agent_state
            st.session_state["states"] = agent_state.get("states", [])
            st.session_state["phase"] = "review"
            st.session_state["analysis_done"] = True

            # Rerun so other tabs pick up the new state
            st.rerun()

    # Show success banner (persists after rerun)
    if st.session_state.get("analysis_done"):
        st.success("🎉 Analysis complete! Click the **Conformer States** tab above to review and approve.")


# ──────────────────────────────────────────────────────────────────────────────
# TAB 2: Conformer States & Editor (v3)
# ──────────────────────────────────────────────────────────────────────────────
with tab_states:
    states = st.session_state.get("states", [])
    agent_state = st.session_state.get("agent_state")

    if not states:
        st.info("Run analysis first to see proposed conformer states.")
    else:
        lig_id = st.session_state.get("lig_id", "?")
        h_diffs = agent_state.get("h_diffs", {}) if agent_state else {}
        state_pdbs = agent_state.get("state_pdbs", {}) if agent_state else {}

        st.subheader(f"Protonation States for {lig_id}")

        # ══════════════════════════════════════════════════════════════════
        # SECTION 1: Hydrogen comparison table + per-state focused views
        # ══════════════════════════════════════════════════════════════════

        # 1a. Hydrogen comparison table — the clearest way to show differences
        if h_diffs:
            st.markdown("#### Hydrogen Differences vs Neutral")
            try:
                from mcce4_agent_ftpl.tools.rdkit_tools import generate_state_comparison_table
                comp_rows = generate_state_comparison_table(states, h_diffs, lig_id)
                if comp_rows:
                    import pandas as pd
                    comp_df = pd.DataFrame(comp_rows)
                    st.dataframe(comp_df, use_container_width=True, hide_index=True)
            except Exception as e:
                logging.warning(f"Comparison table failed: {e}")

        # 1b. Per-state views side by side — H atoms labeled with PDB names
        if state_pdbs and len(state_pdbs) > 0:
            st.markdown("---")
            st.markdown(
                "#### Per-State Structures\n"
                "H atoms labeled with PDB names. "
                "🟢 **Green** = H added vs neutral. "
                "🔵 **Blue** = ionizable site."
            )

            cols_per_row = min(3, len(states))
            for row_start in range(0, len(states), cols_per_row):
                row_cols = st.columns(cols_per_row)
                for i, col in enumerate(row_cols):
                    idx = row_start + i
                    if idx >= len(states):
                        break
                    s = states[idx]
                    sd = s if isinstance(s, dict) else s.to_dict()
                    label = sd.get("label", "?")
                    charge = int(sd.get("charge", 0) or 0)
                    diff = h_diffs.get(label, {})
                    h_added = diff.get("added", []) or sd.get("h_added", [])
                    h_removed = diff.get("removed", []) or sd.get("h_removed", [])
                    site_atom = sd.get("site_atom", "")

                    with col:
                        st.markdown(f"**{lig_id}{label}** (charge: {charge:+d})")

                        pdb_for_state = state_pdbs.get(label) or sd.get("pdb_path")
                        rendered = False

                        if pdb_for_state and os.path.exists(str(pdb_for_state)):
                            try:
                                from mcce4_agent_ftpl.tools.rdkit_tools import (
                                    render_protonation_site_view,
                                    mol_to_svg_with_h_diff,
                                )
                                # If we know the site atom, show focused view
                                if site_atom and (h_added or h_removed):
                                    svg = render_protonation_site_view(
                                        pdb_for_state, site_atom,
                                        h_added, lig_id, label,
                                        size=(380, 320),
                                    )
                                    if svg:
                                        st.image(svg, use_container_width=True)
                                        rendered = True

                                # Fallback: full molecule with H labels
                                if not rendered:
                                    svg = mol_to_svg_with_h_diff(
                                        pdb_for_state, h_added, h_removed,
                                        size=(380, 320),
                                    )
                                    if svg:
                                        st.image(svg, use_container_width=True)
                                        rendered = True
                            except Exception as e:
                                st.warning(f"Render: {e}")

                        # Fallback: SMILES rendering
                        if not rendered:
                            smi = sd.get("smiles", "")
                            if smi:
                                mol = get_mol_from_smiles(smi)
                                if mol:
                                    svg_s = mol_to_svg(mol, size=(350, 280))
                                    if svg_s:
                                        st.image(svg_s, use_container_width=True)

                        # H-diff annotations below each molecule
                        if h_added:
                            added_on = diff.get("added_on", {})
                            for h in h_added:
                                parent = added_on.get(h, "?")
                                st.markdown(f"🟢 **{h}** on {parent}")
                        if h_removed:
                            removed_from = diff.get("removed_from", {})
                            for h in h_removed:
                                parent = removed_from.get(h, "?")
                                st.markdown(f"🔴 **{h}** from {parent}")
                        if label in ("01", "00"):
                            st.caption("_(neutral reference)_")
                        proton_ex = sd.get("proton_exchange", "")
                        if proton_ex:
                            st.caption(f"Exchange: {proton_ex}")

        else:
            # v2 fallback: single molecule view
            pdb = st.session_state.get("pdb_path")
            if pdb:
                svg = mol_to_svg(pdb, size=(450, 350))
                if svg:
                    st.image(svg, use_container_width=True)
                if len(states) > 1:
                    state_cols = st.columns(min(len(states), 3))
                    for i, s in enumerate(states):
                        with state_cols[i % 3]:
                            smi = s.get("smiles", "") if isinstance(s, dict) else s.smiles
                            if smi:
                                mol = get_mol_from_smiles(smi)
                                if mol:
                                    svg_s = mol_to_svg(mol, size=(250, 200))
                                    if svg_s:
                                        lbl = s.get("label", "?") if isinstance(s, dict) else s.label
                                        chg = s.get("charge", 0) if isinstance(s, dict) else s.charge
                                        st.image(svg_s, caption=f"{lbl} ({chg:+d})",
                                                 use_container_width=True)

        # ══════════════════════════════════════════════════════════════════
        # SECTION 1b: PyMOL visualization (full 3D with all H atoms)
        # ══════════════════════════════════════════════════════════════════
        if state_pdbs and len(state_pdbs) > 1:
            st.markdown("---")
            st.subheader("🔬 3D Visualization (PyMOL)")
            st.markdown(
                "For **full control** over hydrogen visualization, use PyMOL "
                "with the generated per-state PDB files. Each PDB has the "
                "correct H atoms for that protonation state."
            )

            pymol_script = agent_state.get("pymol_script") if agent_state else None

            col_pymol, col_files = st.columns([1, 1])

            with col_pymol:
                if pymol_script and os.path.exists(pymol_script):
                    with open(pymol_script) as f:
                        pml_content = f.read()
                    st.download_button(
                        "📥 Download PyMOL script (.pml)",
                        data=pml_content,
                        file_name=os.path.basename(pymol_script),
                        mime="text/plain",
                    )
                    st.caption("Run: `pymol " + os.path.basename(pymol_script) + "`")

                    with st.expander("Preview PyMOL script"):
                        st.code(pml_content[:2000], language="python")
                else:
                    st.info("PyMOL script will be generated after analysis completes.")

            with col_files:
                st.markdown("**Per-state PDB files:**")
                for label, pdb_path in sorted(state_pdbs.items()):
                    if os.path.exists(str(pdb_path)):
                        with open(pdb_path) as f:
                            pdb_content = f.read()
                        n_atoms = pdb_content.count("\nHETATM") + pdb_content.count("\nATOM")
                        n_h = sum(1 for line in pdb_content.splitlines()
                                  if line.startswith(("ATOM", "HETATM")) and
                                  line[76:78].strip() == "H")
                        st.download_button(
                            f"📥 {os.path.basename(pdb_path)} ({n_atoms} atoms, {n_h} H)",
                            data=pdb_content,
                            file_name=os.path.basename(pdb_path),
                            mime="chemical/x-pdb",
                            key=f"dl_{label}",
                        )

        # ══════════════════════════════════════════════════════════════════
        # SECTION 2: Enriched conformer state table
        # ══════════════════════════════════════════════════════════════════
        st.markdown("---")
        st.subheader("Conformer State Table")

        import pandas as pd
        table_data = []
        for s in states:
            sd = s if isinstance(s, dict) else s.to_dict()
            row = {
                "Label": sd.get("label", "?"),
                "Charge": sd.get("charge", 0),
                "nH": sd.get("nH", 0),
                "pKa": sd.get("pka"),
                "Proton Exchange": sd.get("proton_exchange", ""),
                "H Added": ", ".join(sd.get("h_added", [])) or "—",
                "H Removed": ", ".join(sd.get("h_removed", [])) or "—",
                "Source": sd.get("source", ""),
                "LLM Model": sd.get("llm_model", ""),
                "Rationale": sd.get("rationale", ""),
                "References": "; ".join(sd.get("references", [])) if sd.get("references") else "—",
            }
            table_data.append(row)

        df = pd.DataFrame(table_data)
        edited_df = st.data_editor(
            df,
            num_rows="dynamic",
            use_container_width=True,
            column_config={
                "Label": st.column_config.TextColumn("Label", help="e.g., 01, +1, -1"),
                "Charge": st.column_config.NumberColumn("Charge", format="%+d"),
                "nH": st.column_config.NumberColumn("nH", help="Protons relative to neutral"),
                "pKa": st.column_config.NumberColumn("pKa", format="%.1f"),
                "Proton Exchange": st.column_config.TextColumn(
                    "Proton Exchange", help="Which protons differ from neutral",
                    width="medium"),
                "H Added": st.column_config.TextColumn("H Added", disabled=True),
                "H Removed": st.column_config.TextColumn("H Removed", disabled=True),
                "Source": st.column_config.TextColumn("Source"),
                "LLM Model": st.column_config.TextColumn("LLM Model", disabled=True),
                "Rationale": st.column_config.TextColumn("Rationale", width="medium"),
                "References": st.column_config.TextColumn(
                    "References", help="Literature sources",
                    width="large"),
            },
            key="state_editor"
        )

        # Merge edits back into state dicts — INCLUDING label changes
        if edited_df is not None:
            old_labels = [
                (s.get("label", "?") if isinstance(s, dict) else s.label)
                for s in states
            ]
            label_remap = {}  # old_label → new_label

            for i, row in edited_df.iterrows():
                if i < len(states):
                    orig = states[i] if isinstance(states[i], dict) else states[i].to_dict()
                    new_label = str(row.get("Label", orig.get("label", "?")))
                    old_label = orig.get("label", "?")

                    if new_label != old_label:
                        label_remap[old_label] = new_label
                        logging.info(f"  Label change: {old_label} → {new_label}")

                    orig["label"] = new_label
                    try:
                        orig["charge"] = int(row.get("Charge", orig.get("charge", 0)) or 0)
                    except (ValueError, TypeError):
                        pass
                    try:
                        orig["nH"] = int(row.get("nH", orig.get("nH", 0)) or 0)
                    except (ValueError, TypeError):
                        pass
                    orig["pka"] = row.get("pKa")
                    orig["rationale"] = row.get("Rationale", "")
                    orig["proton_exchange"] = row.get("Proton Exchange", "")
                    states[i] = orig

            st.session_state["states"] = states

            # Propagate label renames to agent_state keys
            if label_remap and agent_state:
                # Update state_pdbs keys
                if agent_state.get("state_pdbs"):
                    new_pdbs = {}
                    for old_lbl, pdb_path in agent_state["state_pdbs"].items():
                        new_lbl = label_remap.get(old_lbl, old_lbl)
                        new_pdbs[new_lbl] = pdb_path
                    agent_state["state_pdbs"] = new_pdbs

                # Update h_diffs keys
                if agent_state.get("h_diffs"):
                    new_diffs = {}
                    for old_lbl, diff in agent_state["h_diffs"].items():
                        new_lbl = label_remap.get(old_lbl, old_lbl)
                        new_diffs[new_lbl] = diff
                    agent_state["h_diffs"] = new_diffs

                # Update conformer_labels
                agent_state["conformer_labels"] = [
                    (s.get("label") if isinstance(s, dict) else s.label)
                    for s in states
                ]

                st.session_state["agent_state"] = agent_state

        # Apply changes button — refreshes the page with updated labels
        if st.button("🔄 Apply Table Changes", help="Refresh page after editing labels or parameters"):
            st.rerun()

        # ══════════════════════════════════════════════════════════════════
        # SECTION 3: Add custom state via ionizable site selector (not Ketcher)
        # ══════════════════════════════════════════════════════════════════
        st.markdown("---")
        st.subheader("🧪 Add / Remove Protonation State")
        st.caption(
            "Select an ionizable site on the molecule to create a new "
            "protonation state. RDKit handles the chemistry automatically."
        )

        # Get ionizable sites from agent state
        ionizable = []
        if agent_state:
            ionizable = agent_state.get("_ionizable_sites", [])
            if not ionizable:
                # Try to compute from PDB
                try:
                    from mcce4_agent_ftpl.tools.rdkit_tools import (
                        get_ionizable_sites, get_mol_from_pdb as _gm
                    )
                    pdb = st.session_state.get("pdb_path")
                    if pdb:
                        mol = _gm(pdb, remove_hs=False)
                        if mol:
                            ionizable = get_ionizable_sites(mol)
                except Exception:
                    pass

        col_proto, col_depro = st.columns(2)

        with col_proto:
            st.markdown("**➕ Protonate a site** (add H)")
            proto_sites = [s for s in ionizable if "protonatable" in s.get("type", "")]
            if proto_sites:
                options = [
                    f"{s['name']} ({s['symbol']}) — {s['type'].replace('_', ' ')}"
                    for s in proto_sites
                ]
                selected = st.selectbox("Protonate at:", options, key="proto_site")
                new_label = st.text_input("Label:", value="+1", key="proto_label")
                if st.button("➕ Add protonated state", key="btn_proto"):
                    site_idx = options.index(selected)
                    site = proto_sites[site_idx]
                    new_state = {
                        "label": new_label, "charge": 1,
                        "nH": 1, "pka": None,
                        "smiles": "", "source": "user",
                        "site_atom": site["name"],
                        "rationale": f"User-added: protonate {site['name']}",
                        "proton_exchange": f"+H on {site['name']} ({site['type'].replace('_', ' ')})",
                        "pdb_path": None, "h_added": [], "h_removed": [],
                        "llm_model": "", "references": [],
                    }
                    states.append(new_state)
                    st.session_state["states"] = states
                    st.success(f"Added state '{new_label}': +H on {site['name']}")
                    st.rerun()
            else:
                st.info("No protonatable sites detected (run analysis first)")

        with col_depro:
            st.markdown("**➖ Deprotonate a site** (remove H)")
            depro_sites = [s for s in ionizable if "deprotonatable" in s.get("type", "")]
            if depro_sites:
                options = [
                    f"{s['name']} ({s['symbol']}) — {s['type'].replace('_', ' ')}"
                    for s in depro_sites
                ]
                selected = st.selectbox("Deprotonate at:", options, key="depro_site")
                new_label = st.text_input("Label:", value="-1", key="depro_label")
                if st.button("➖ Add deprotonated state", key="btn_depro"):
                    site_idx = options.index(selected)
                    site = depro_sites[site_idx]
                    new_state = {
                        "label": new_label, "charge": -1,
                        "nH": -1, "pka": None,
                        "smiles": "", "source": "user",
                        "site_atom": site["name"],
                        "rationale": f"User-added: deprotonate {site['name']}",
                        "proton_exchange": f"-H from {site['name']} ({site['type'].replace('_', ' ')})",
                        "pdb_path": None, "h_added": [], "h_removed": [],
                        "llm_model": "", "references": [],
                    }
                    states.append(new_state)
                    st.session_state["states"] = states
                    st.success(f"Added state '{new_label}': -H from {site['name']}")
                    st.rerun()
            else:
                st.info("No deprotonatable sites detected (run analysis first)")

        # Optional SMILES fallback for advanced users
        with st.expander("Advanced: add state from SMILES", expanded=False):
            manual_smiles = st.text_input(
                "SMILES for custom state:",
                placeholder="e.g., [NH2+]1CCCCC1",
                key="manual_smiles"
            )
            if manual_smiles:
                try:
                    from rdkit import Chem
                    mol = Chem.MolFromSmiles(manual_smiles)
                    if mol:
                        charge = Chem.GetFormalCharge(mol)
                        svg_custom = mol_to_svg(mol, size=(300, 200))
                        if svg_custom:
                            st.image(svg_custom, use_container_width=False)
                        st.info(f"Formal charge: {charge:+d}")
                        smiles_label = st.text_input("Label:",
                                                      value=f"{charge:+d}" if charge else "01",
                                                      key="smiles_label")
                        if st.button("➕ Add", key="add_smiles"):
                            new_state = {
                                "label": smiles_label, "charge": charge,
                                "nH": charge, "pka": None,
                                "smiles": manual_smiles,
                                "source": "user-smiles",
                                "rationale": "Custom SMILES entry",
                                "proton_exchange": "",
                                "pdb_path": None, "site_atom": None,
                                "h_added": [], "h_removed": [],
                                "llm_model": "", "references": [],
                            }
                            states.append(new_state)
                            st.session_state["states"] = states
                            st.success(f"Added '{smiles_label}' (charge={charge:+d})")
                            st.rerun()
                    else:
                        st.error("Invalid SMILES")
                except ImportError:
                    st.error("RDKit required — pip install rdkit")

        # ══════════════════════════════════════════════════════════════════
        # Approve / Run
        # ══════════════════════════════════════════════════════════════════
        st.markdown("---")
        col_cancel, col_spacer, col_approve = st.columns([1, 2, 1])

        with col_cancel:
            if st.button("✖ Cancel", type="secondary", use_container_width=True):
                st.session_state["phase"] = "upload"
                st.session_state["states"] = []
                st.rerun()

        with col_approve:
            if st.button("✔ Approve & Generate .ftpl", type="primary", use_container_width=True):
                st.session_state["phase"] = "running"
                st.session_state["approved"] = True

                from mcce4_agent_ftpl.agent import resume_agent

                agent_state = st.session_state["agent_state"]
                agent_state["states"] = st.session_state["states"]
                agent_state["conformer_labels"] = [
                    (s["label"] if isinstance(s, dict) else s.label)
                    for s in st.session_state["states"]
                ]
                agent_state["user_approved"] = True

                with st.spinner("Agent is generating topology file..."):
                    final = resume_agent(agent_state)

                st.session_state["agent_state"] = final
                st.session_state["phase"] = "complete"
                st.rerun()


# ──────────────────────────────────────────────────────────────────────────────
# TAB 3: Output
# ──────────────────────────────────────────────────────────────────────────────
with tab_output:
    st.subheader("📄 Generated Topology File")

    agent_state = st.session_state.get("agent_state")
    if agent_state and agent_state.get("ftpl_path"):
        ftpl_path = agent_state["ftpl_path"]

        if os.path.exists(ftpl_path):
            with open(ftpl_path) as f:
                ftpl_content = f.read()

            # Download button
            st.download_button(
                "⬇️ Download .ftpl",
                data=ftpl_content,
                file_name=os.path.basename(ftpl_path),
                mime="text/plain",
                type="primary"
            )

            # Show content
            st.code(ftpl_content, language="text", line_numbers=True)

            # Show summary
            if agent_state.get("rxn_values"):
                st.subheader("RXN Calibration Results")
                st.json(agent_state["rxn_values"])

            if agent_state.get("warnings"):
                st.subheader("⚠️ Warnings")
                for w in agent_state["warnings"]:
                    st.warning(w)
        else:
            st.warning(f"File not found: {ftpl_path}")
    else:
        st.info("Complete the analysis and approve states to generate the topology file.")


# ──────────────────────────────────────────────────────────────────────────────
# Footer
# ──────────────────────────────────────────────────────────────────────────────
st.markdown("---")
st.caption("MCCE4 Topology Agent v3.0 | GunnerLab | "
           "[Tutorial](https://gunnerlab.github.io/mcce4_tutorial/docs/topology/)")
