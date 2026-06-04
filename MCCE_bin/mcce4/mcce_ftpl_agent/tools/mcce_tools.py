#!/usr/bin/env python3

"""
MCCE4 pipeline tools — template generation, calibration, parsing.
"""
from collections import Counter
import json
import logging
from pathlib import Path
import re
import shutil
import subprocess
import sys

from mcce4.pymol_ligand_cif2pdb import convert as cif_ligand_converter
from mcce4.mcce_ftpl_agent.config import DIELECTRIC_MAP, RCSB_GRAPHQL_URL, RCSB_GRAPHQL_QUERY


def convert_cif_to_pdb(cif_path: str) -> str:
    """Convert a .cif file to .pdb using pymol_ligand_cif2pdb.

    Returns the path to the generated PDB file.
    """
    pdb_path = Path(cif_path).with_suffix(".pdb")
    print(f"🔄 Converting {cif_path} → {pdb_path!s} using pymol_ligand_cif2pdb...")

    cif_ligand_converter({"input": cif_path,
                          "output": pdb_path,
                          "ligand_cache": None,
                          "offline": False})
    if not pdb_path.exists():
        print(f"ERROR: Conversion completed but {pdb_path.name} was not created.")
        sys.exit(1)
    print(f"✅ Converted successfully: {pdb_path.name}")

    return pdb_path


def _copy_pdb_with_chain_id(src_pdb: str,
                            dst_pdb: str,
                            default_chain: str = "A"):
    """Copy a PDB file, assigning default_chain to any atom line with a blank chain ID.
    """
    with open(src_pdb) as f:
        lines = f.readlines()
    with open(dst_pdb, "w") as f:
        for line in lines:
            if line.startswith(("ATOM", "HETATM")) and len(line) > 21 and line[21] == " ":
                line = line[:21] + default_chain + line[22:]
            f.write(line)


def run_cmd(cmd, description="", capture=False, cwd=None):
    """Run a shell command with logging.
    """
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
    """Extract 3-letter ligand code from PDB HETATM records."""
    ligs = []
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith("HETATM"):
                lig = line[17:20].strip()
                if lig in ("HOH", "WAT", "TIP", "UNK", "UNL"):
                    continue
                ligs.append(lig)
            if line.startswith(("CONECT")):
                break
    if not ligs:
        return ""

    return Counter(ligs).most_common(n=1)[0][0]


def fetch_rcsb_info(lig_id: str, work_dir: str = ".") -> dict:
    """Fetch SMILES and comprehensive metadata from RCSB via GraphQL.

    Uses the RCSB GraphQL API to retrieve all available ligand metadata
    including SMILES, InChI, DrugBank info, synonyms, and related resources.
    Logs everything and saves the SMILES string to {lig_id}.smi.

    Args:
        lig_id: 3-letter RCSB ligand code (e.g., "EMH").
        work_dir: Directory to save the .smi file.

    Returns:
        Dict with keys: name, formula, formal_charge, molecular_weight,
        smiles, smiles_stereo, type, drugbank, etc.
    """
    logging.info(f"🔍 Fetching info for '{lig_id}' from RCSB GraphQL...")
    logging.info(f"   Ligand page: https://www.rcsb.org/ligand/{lig_id}")

    # Build GraphQL request
    # FIX: Use requests; run_cmd with logging is listing fields that are not
    # requested, possibly bc of curl cmd.
    payload = json.dumps({"query": RCSB_GRAPHQL_QUERY, "variables": {"id": lig_id}})
    cmd = ("curl -s -X POST -H 'Content-Type: application/json' -d "
           f"'{payload}' '{RCSB_GRAPHQL_URL}'")
    result = run_cmd(cmd,
                     description="Querying RCSB GraphQL",
                     capture=True
                     )
    if result is None or not result.stdout.strip():
        logging.warning("  RCSB GraphQL query returned empty response")
        return {}

    try:
        response = json.loads(result.stdout)
        data = response.get("data", {}).get("chem_comp", {})
        if not data:
            logging.warning(f"  No data returned for ligand '{lig_id}'")
            return {}

        # Core chemical component info
        chem = data.get("chem_comp") or {}
        info = {
            "name": chem.get("name", "Unknown"),
            "formula": chem.get("formula", "Unknown"),
            "formal_charge": chem.get("pdbx_formal_charge", 0),
            "molecular_weight": chem.get("formula_weight"),
            "type": chem.get("type"),
            "smiles": None,
            "smiles_stereo": None,
        }

        # Descriptors (SMILES)
        desc = data.get("rcsb_chem_comp_descriptor") or {}
        if desc.get("SMILES"):
            info["smiles"] = desc["SMILES"]
        if desc.get("SMILES_stereo"):
            info["smiles_stereo"] = desc["SMILES_stereo"]

        # Additional descriptors
        pdbx_desc = data.get("pdbx_chem_comp_descriptor") or []
        for d in pdbx_desc:
            if not isinstance(d, dict):
                continue
            dtype = d.get("type", "")
            descriptor = d.get("descriptor", "")
            program = d.get("program", "")
            version = d.get("program_version", "")
            if descriptor:
                logging.info(f"  ✓ {dtype} ({program} {version}): {descriptor[:120]}")
            # Fallback: if SMILES wasn't in rcsb_chem_comp_descriptor
            # if "SMILES" in dtype and not info["smiles"]:
            #     info["smiles"] = descriptor

        # Compound info (atom/bond counts, dates)
        comp_info = data.get("rcsb_chem_comp_info") or {}
        if comp_info:
            info["atom_count"] = comp_info.get("atom_count")
            info["bond_count"] = comp_info.get("bond_count")
            info["bond_count_aromatic"] = comp_info.get("bond_count_aromatic")
            info["atom_count_chiral"] = comp_info.get("atom_count_chiral")
            info["initial_deposition_date"] = comp_info.get("initial_deposition_date")
            info["revision_date"] = comp_info.get("revision_date")

        # Identifiers (IUPAC name, etc.)
        identifiers = data.get("pdbx_chem_comp_identifier") or []
        for ident in identifiers:
            if isinstance(ident, dict) and ident.get("identifier"):
                logging.info(f"  ✓ Identifier ({ident.get('program', '?')}): "
                             f"{ident['identifier'][:120]}")

        # Synonyms
        synonyms = data.get("rcsb_chem_comp_synonyms") or []
        synonym_names = [s.get("name") for s in synonyms
                         if isinstance(s, dict) and s.get("name")]
        if synonym_names:
            info["synonyms"] = synonym_names

        # DrugBank info
        drugbank = data.get("drugbank") or {}
        db_info = drugbank.get("drugbank_info") or {}
        if db_info:
            info["drugbank_id"] = db_info.get("drugbank_id")
            info["drugbank_name"] = db_info.get("name")
            info["cas_number"] = db_info.get("cas_number")
            info["drug_categories"] = db_info.get("drug_categories")
            info["mechanism_of_action"] = db_info.get("mechanism_of_action")
            info["drug_groups"] = db_info.get("drug_groups")
            info["indication"] = db_info.get("indication")

        # Reference molecule (PRD info)
        ref_mol = data.get("pdbx_reference_molecule") or {}
        if ref_mol and ref_mol.get("prd_id"):
            info["prd_id"] = ref_mol.get("prd_id")
            info["prd_class"] = ref_mol.get("class")
            info["representative_pdb"] = ref_mol.get("representative_PDB_id_code")

        # ── Log all retrieved metadata ──
        logging.info(f"\n  {'─'*55}")
        logging.info(f"  RCSB Ligand Summary: {lig_id}")
        logging.info(f"  {'─'*55}")
        for k, v in info.items():
            if v is not None and k not in ("external_ids", "synonyms",
                                           "drug_categories", "mechanism_of_action",
                                           "indication"):
                display = str(v)[:100]
                logging.info(f"    {k:30s}: {display}")
        if info.get("synonyms"):
            logging.info(f"    {'synonyms':30s}: {', '.join(info['synonyms'][:5])}")
        if info.get("external_ids"):
            for res, acc in info["external_ids"].items():
                logging.info(f"    {'external: ' + res:30s}: {acc}")
        if info.get("drugbank_id"):
            logging.info(f"    {'drugbank_id':30s}: {info['drugbank_id']}")
            if info.get("mechanism_of_action"):
                logging.info(f"    {'mechanism_of_action':30s}: "
                             f"{str(info['mechanism_of_action'])[:100]}")
        logging.info(f"  {'─'*55}")

        # ── Save SMILES to .smi file ──
        smiles = info.get("smiles") or info.get("smiles_stereo")
        if smiles:
            smi_path = Path(work_dir).joinpath(f"{lig_id}.smi")
            with open(smi_path, "w") as f:
                f.write(f"{smiles} {lig_id}\n")
            logging.info(f"  💾 SMILES saved to: {smi_path!s}")
            info["smi_path"] = str(smi_path)

        return info

    except (json.JSONDecodeError, KeyError) as e:
        logging.warning(f"  RCSB GraphQL parse error: {e}")
        return {}


def generate_ftpl_template(lig_id: str,
                           pdb_path: str,
                           conformer_labels: list,
                           ftpl_path: str) -> int:
    """Generate .ftpl template via pdb2ftpl.py. Returns count of unfilled entries.
    """
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
                filled += 1
                matched = True
            # Stripped match
            if not matched:
                stripped = atom_name.strip()
                for cn, cv in charges.items():
                    if cn.strip() == stripped:
                        new_lines.append(f"{prefix}{cv:7.3f} # auto-filled{suffix}\n")
                        filled += 1
                        matched = True
                        break
            if not matched:
                new_lines.append(line)
                unfilled += 1
        else:
            new_lines.append(line)

    with open(ftpl_path, "w") as f:
        f.writelines(new_lines)
    logging.info(f"  ✓ Filled {filled}, unfilled {unfilled}")
    return unfilled


def update_conformer_params(ftpl_path: str, states: list, lig_id: str):
    """Update nH and pKa0 in CONFORMER lines.
    """
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
        charge = (s.get('charge', 0) if isinstance(s, dict) else getattr(s, 'charge', 0)) or 0
        logging.info(f"  📝 {ct}: pka_raw={pka_raw!r} (type={type(pka_raw).__name__}), "
                     f"from state dict pka key={s.get('pka') if isinstance(s, dict) else s.pka!r}")
        try:
            pka = float(str(pka_raw).lstrip("~≈><≥≤ ")) if pka_raw is not None else 0.0
        except (ValueError, TypeError):
            pka = 0.0

        # Neutral reference states should have pKa0=0.0; charged states must not
        is_neutral = (int(charge) == 0)
        if not is_neutral and (pka == 0.0 or pka_raw is None):
            logging.warning(f"  ⚠ {ct}: charged conformer (charge={charge}) has no pKa estimate "
                            f"(pka_raw={pka_raw!r}). A pKa0 value is required for charged states.")
        logging.info(f"  📝 {ct}: pka_raw={pka_raw!r} → final pKa0={pka:.2f}")

        p = re.compile(rf"(CONFORMER,\s*{re.escape(ct)}:.*?nH=)\s*(\d+)")
        p2 = re.compile(rf"(CONFORMER,\s*{re.escape(ct)}:.*?pKa0=)\s*([-\d.]+)")
        if p.search(content) and p2.search(content):
            content = p.sub(rf"\g<1>{nH}", content)
            content = p2.sub(rf"\g<1>{pka:.2f}", content)
            logging.info(f"  ✓ {ct}: nH={nH}, pKa0={pka:.2f}")
        else:
            logging.warning(f"  ⚠ {ct}: CONFORMER line not found in ftpl — skipping nH/pKa0 update")

    with open(ftpl_path, "w") as f:
        f.write(content)


def run_rxn_calibration(lig_id: str,
                        ftpl_path: str,
                        pdb_path: str,
                        conformer_labels: list,
                        dielectrics: list,
                        work_dir: str):
    """Run MCCE4 steps 1-3, calibrate rxn values for all dielectrics.

    Workflow:
      1. Link ftpl to user_param/ so MCCE uses the generated topology
      2. Run step1.py and step2.py once
      3. Verify all expected conformers are in step2_out.pdb
      4. For each dielectric (2, 4, 8):
         a. Run step3.py -d <eps>
         b. Parse head3.lst — pick lowest dsolv per conformer type
         c. Update rxn<eps> in the ftpl CONFORMER lines
         d. Re-run step3.py -d <eps> to validate (dsolv should be ~0)
    """
    work_dir = Path(work_dir)
    user_param = work_dir.joinpath("user_param")
    user_param.mkdir(exist_ok=True)

    # ── Link ftpl to user_param BEFORE step1 ──
    ftpl_abs = Path(ftpl_path).resolve()
    link = user_param.joinpath(ftpl_abs.name)
    if link.exists():
        link.unlink()
    # if os.path.islink(link) or os.path.exists(link):
    #     os.remove(link)
    link.symlink_to(ftpl_abs)
    # os.symlink(ftpl_abs, link)
    logging.info(f"  🔗 Linked {ftpl_abs!s} → {link!s}")

    # Ensure PDB is accessible in work_dir, with chain ID defaulting to 'A'
    pdb_path = Path(pdb_path).resolve()
    pdb_base = pdb_path.name
    _copy_pdb_with_chain_id(str(pdb_path), str(work_dir.joinpath(pdb_base)))

    # ── Pre-flight: verify ftpl CONFLIST matches expected states ──
    expected_types = [f"{lig_id}{cl}" for cl in conformer_labels]
    ftpl_content = ftpl_abs.read_text()
    conflist_match = re.search(r'CONFLIST,\s*' + re.escape(lig_id) + r':\s*(.*)', ftpl_content)
    if conflist_match:
        conflist_str = conflist_match.group(1)
        for ct in expected_types:
            if ct not in conflist_str:
                logging.error(f"  ✗ {ct} missing from ftpl CONFLIST — step2 will not create it")
        # Also verify each conftype has CONNECT records
        for ct in expected_types:
            connect_count = len(re.findall(
                r'^CONNECT,.*' + re.escape(ct) + r':',
                ftpl_content, re.MULTILINE
            ))
            if connect_count == 0:
                logging.error(f"  ✗ {ct} has NO CONNECT records in ftpl")
            else:
                logging.info(f"  ✓ {ct}: {connect_count} CONNECT records in ftpl")
    else:
        logging.error(f"  No CONFLIST found for {lig_id} in ftpl!")

    # ── Step 1: Initialize MCCE ──
    if run_cmd(f"step1.py {pdb_base}", description="MCCE4 Step 1", cwd=work_dir) is None:
        return {"error": "Step 1 failed"}

    # ── Step 2: Generate conformers ──
    if run_cmd("step2.py", description="MCCE4 Step 2", cwd=work_dir) is None:
        return {"error": "Step 2 failed"}

    # ── Verify all conformers are in step2_out.pdb ──
    step2_out = work_dir.joinpath("step2_out.pdb")
    if step2_out.exists():
        found_types = _check_conformers_in_step2(str(step2_out), lig_id, conformer_labels)
        expected_types = [f"{lig_id}{cl}" for cl in conformer_labels]
        missing = [t for t in expected_types if t not in found_types]
        if missing:
            logging.warning(f"  ⚠ Missing conformers in step2_out.pdb: {missing}")
            logging.warning(f"  Found: {sorted(found_types)}")
            # Diagnose: check ftpl for CONFLIST and CONNECT records
            _diagnose_missing_conformers(str(ftpl_abs), lig_id, conformer_labels, missing)
        else:
            logging.info(f"  ✅ All {len(expected_types)} conformer types found in step2_out.pdb")
            for ct, count in sorted(found_types.items()):
                logging.info(f"    {ct}: {count} conformer(s)")
    else:
        logging.warning("  step2_out.pdb not found — cannot verify conformers")

    # ── Step 3: RXN calibration for each dielectric ──
    rxn_results = {}

    for eps in dielectrics:
        rxn_key = DIELECTRIC_MAP.get(eps)
        if not rxn_key:
            continue
        logging.info(f"\n{'='*60}")
        logging.info(f"  🔬 RXN Calibration ε={eps} ({rxn_key})")
        logging.info(f"{'='*60}")

        # Run step3.py -d <eps>
        result = run_cmd(f"step3.py -d {eps}", description=f"Step 3 (ε={eps})", cwd=work_dir)
        if result is None:
            logging.error(f"  step3.py -d {eps} failed")
            continue

        head3 = str(work_dir.joinpath("head3.lst"))
        dsolv = parse_head3_dsolv(head3, lig_id, conformer_labels)
        if not dsolv:
            logging.error("  Could not parse dsolv from head3.lst")
            continue

        for ct, v in dsolv.items():
            logging.info(f"  📊 {ct}: dsolv = {v:.3f} (lowest across all conformers)")

        # Update rxn value in ftpl
        update_rxn(str(ftpl_abs), dsolv, rxn_key)

        # Re-link ftpl after update (symlink target content changed in-place, so OK)
        logging.info(f"  Updated {rxn_key} in {ftpl_abs.name}")

        # Validate: re-run step3 to check dsolv is now ~0
        logging.info("  🔄 Validating calibration...")
        run_cmd(f"step3.py -d {eps}", description=f"Validation (ε={eps})", cwd=work_dir)
        check = parse_head3_dsolv(head3, lig_id, conformer_labels)
        all_ok = True
        for ct, v in check.items():
            ok = abs(v) <= 0.05
            sym = "✓" if ok else "⚠"
            if not ok:
                all_ok = False
            logging.info(f"  {sym} {ct}: dsolv = {v:.3f}")
        if all_ok:
            logging.info(f"  🎉 {rxn_key} calibration successful!")
        else:
            logging.warning(f"  ⚠ {rxn_key} calibration has residual dsolv > 0.05")

        rxn_results[rxn_key] = dsolv

    return rxn_results


def _check_conformers_in_step2(step2_path: str, lig_id: str, labels: list) -> dict:
    """Check which conformer types are present in step2_out.pdb.

    MCCE PDB format:
      - Columns 17:20 = resName (3-char, e.g., "EMH")
      - Columns 80:82 = 2-char conformer type suffix (e.g., "01", "+1", "-1")
      - confType = resName + suffix (e.g., "EMH01", "EMH+1")
      - Full confID = confType + chainID + resSeq + iCode + confNum

    Returns {conformer_type: count} where count = number of unique conformers.
    """
    type_counts = {}
    # Track unique conformers by their full identity (confType + chain + seq + confNum)
    seen_conformers = set()

    with open(step2_path) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if len(line) < 82:
                continue
            resname = line[17:20].strip()
            if resname != lig_id:
                continue

            # MCCE PDB: confType = resName(3) + cols 80:82(2)
            conf_suffix = line[80:82]
            conf_type = f"{resname}{conf_suffix}".rstrip()

            # Build a unique conformer key from chain + resSeq + iCode + history
            chain_id = line[21]
            res_seq = line[22:26].strip()
            icode = line[26]
            history = line[80:].strip()
            conf_key = f"{conf_type}_{chain_id}_{res_seq}_{icode}_{history}"

            if conf_key not in seen_conformers:
                seen_conformers.add(conf_key)
                type_counts[conf_type] = type_counts.get(conf_type, 0) + 1

    return type_counts


def _diagnose_missing_conformers(ftpl_path: str, lig_id: str, labels: list, missing: list):
    """Diagnose why conformers are missing from step2_out.pdb.

    Checks the ftpl file for:
      1. CONFLIST line — does it include all expected conformer types?
      2. CONNECT records — does each conformer type have atom definitions?
    """
    if not Path(ftpl_path).exists():
        logging.error(f"  ftpl file not found: {ftpl_path}")
        return

    content = Path(ftpl_path).read_text()
    # Check CONFLIST
    conflist_match = re.search(r'CONFLIST,\s*' + re.escape(lig_id) + r':\s*(.*)', content)
    if conflist_match:
        conflist_str = conflist_match.group(1)
        logging.info(f"  CONFLIST in ftpl: {conflist_str.strip()}")
        for m in missing:
            if m not in conflist_str:
                logging.error(f"    ✗ {m} is NOT in CONFLIST — step2 cannot create it")
            else:
                logging.info(f"    ✓ {m} is in CONFLIST")
    else:
        logging.error(f"  No CONFLIST line found for {lig_id} in ftpl!")

    # Check CONNECT records per conformer type
    for m in missing:
        connect_count = len(re.findall(
            r'^CONNECT,.*' + re.escape(m) + r':',
            content, re.MULTILINE
        ))
        if connect_count == 0:
            logging.error(f"    ✗ {m} has NO CONNECT records in ftpl — "
                          f"step2 rot_ionization skips types without CONNECT")
        else:
            logging.info(f"    ✓ {m} has {connect_count} CONNECT record(s)")

    # Check CHARGE records per conformer type
    for m in missing:
        charge_count = len(re.findall(
            r'^CHARGE,\s*' + re.escape(m) + r',',
            content, re.MULTILINE
        ))
        if charge_count == 0:
            logging.warning(f"    ⚠ {m} has NO CHARGE records in ftpl")
        else:
            logging.info(f"    ✓ {m} has {charge_count} CHARGE record(s)")


def parse_head3_dsolv(head3_path: str, lig_id: str, labels: list) -> dict:
    """Parse most negative dsolv per conformer type from head3.lst.

    head3.lst format:
    iConf CONFORMER     FL  occ  crg  Em0 pKa0 ne nH vdw0 vdw1 tors epol dsolv extra history
    00001 EMH01_0000_001 f 0.00 0.000  0  0.00  0  0 -15.9 0.0  0.0  0.0  -8.1  0.0 01O000 t

    For each conformer type (e.g., EMH01, EMH+1), take the most negative dsolv.
    """
    if not Path(head3_path).exists():
        return {}
    with open(head3_path) as fh:
        lines = fh.readlines()

    types = [f"{lig_id}{lb}" for lb in labels]
    logging.info(f"  Looking for conformer types: {types}")

    # Find dsolv column from header
    header = next((line for line in lines if "dsolv" in line.lower() and "iConf" in line), None)
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
        pdb_base = Path(state_pdb).name
        cmd = f"pdb2ftpl.py -p {pdb_base} -c {label}"
        logging.info(f"  Running: {cmd}")
        result = run_cmd(cmd, description=f"pdb2ftpl for {label}",
                         capture=True, cwd=output_dir)
        if result is None:
            logging.error(f"  pdb2ftpl.py failed for state {label}")
            return None

        output_dir = Path(output_dir)
        # pdb2ftpl.py writes to stdout — save it
        ftpl_name = f"{lig_id}_{label.replace('+', 'p').replace('-', 'm')}.ftpl"
        ftpl_path = output_dir.joinpath(ftpl_name)

        # If pdb2ftpl wrote to stdout, save it
        if result.stdout.strip():
            with open(ftpl_path, "w") as fh:
                fh.write(result.stdout)
        else:
            # Check if it wrote a file directly
            default_ftpl = output_dir.joinpath(f"{lig_id}.ftpl")
            if default_ftpl.exists():
                shutil.move(default_ftpl, ftpl_path)
            else:
                logging.error(f"  No ftpl output for state {label}")
                return None

        logging.info(f"  → {ftpl_path!s}")
        return str(ftpl_path)

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
        # Copy state PDB to output_dir, ensuring chain ID defaults to 'A'
        pdb_dest = str(Path(output_dir).joinpath(Path(pdb_path).name))
        _copy_pdb_with_chain_id(pdb_path, pdb_dest)

        ftpl = generate_ftpl_for_state(lig_id, pdb_path, label, output_dir)
        if ftpl:
            ftpls[label] = ftpl
        else:
            logging.error(f"  FAILED: No ftpl for state {label}")
    return ftpls


def merge_ftpl_files(per_state_ftpls, lig_id, output_path,
                     states=None, smiles="", warnings=None):
    """Merge per-state ftpl files into a single master .ftpl.

    Each state's ftpl has its own CONNECT, CHARGE, RADIUS, and CONFORMER
    blocks.  The merge:
      1. Writes header comment block with molecule metadata.
      2. Rebuilds CONFLIST to list ALL conformer states (not just neutral).
      3. Combines CONNECT, CHARGE, RADIUS, and CONFORMER blocks from every
         state in sorted label order.
      4. Appends any warnings as comments at the bottom.

    Args:
        per_state_ftpls: {label: ftpl_path} for each state.
        lig_id: 3-letter ligand code.
        output_path: Path for merged output .ftpl.
        states: List of ConformerState dicts (for header metadata).
        smiles: Canonical SMILES string (for header).
        warnings: List of warning strings to append as comments.

    Returns True on success.
    """
    if not per_state_ftpls:
        logging.error("No ftpl files to merge")
        return False

    from ..models import sort_conformer_labels
    from datetime import datetime

    # neutral_label = "01" if "01" in per_state_ftpls else list(per_state_ftpls.keys())[0]
    all_labels = sort_conformer_labels(per_state_ftpls.keys())

    parsed = {}
    for label, ftpl_path in per_state_ftpls.items():
        parsed[label] = _parse_ftpl_sections(ftpl_path, lig_id)

    # base = parsed[neutral_label]
    merged = []

    # ── Header comment block ──
    merged.append(">>>START of original comments, this file was generated by MCCE4 Topology Agent\n")
    merged.append(f"# Ligand: {lig_id}\n")
    if smiles:
        merged.append(f"# SMILES: {smiles}\n")
    # Per-state info
    if states:
        for s in states:
            sd = s if isinstance(s, dict) else s.to_dict() if hasattr(s, 'to_dict') else {}
            label = sd.get('label', '?')
            charge = sd.get('charge', 0)
            smi = sd.get('smiles', '')
            src = sd.get('source', '')
            rat = sd.get('rationale', '')
            merged.append(f"# {lig_id}{label}: charge={charge:+d}, source={src}\n")
            if smi:
                merged.append(f"#   SMILES: {smi}\n")
            if rat:
                merged.append(f"#   Rationale: {rat}\n")
    merged.append(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    merged.append("#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I\n")
    merged.append("#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn\n")
    merged.append("#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|  \n")
    merged.append("<<<END of original comments\n\n")

    # ── CONFLIST ──
    merged.append("# Values of the same key are appended and separated by \",\"\n")
    conftypes = ", ".join(f"{lig_id}{lb}" for lb in all_labels)
    merged.append(f"CONFLIST, {lig_id}: {lig_id}BK, {conftypes}\n\n")

    # ── CONNECT blocks from all states ──
    merged.append("# Atom definition\n")
    for label in all_labels:
        if parsed[label]["connect"]:
            merged.extend(parsed[label]["connect"])
            merged.append("\n")

    # ── CHARGE blocks from all states ──
    merged.append("# Atom charges\n")
    for label in all_labels:
        if parsed[label]["charge"]:
            merged.extend(parsed[label]["charge"])
            merged.append("\n")

    # ── RADIUS blocks from all states ──
    merged.append("# Atom radius, dielectric boundary radius, VDW radius, and energy well depth\n")
    for label in all_labels:
        if parsed[label].get("radius"):
            merged.extend(parsed[label]["radius"])
            merged.append("\n")

    # ── CONFORMER lines from all states ──
    merged.append("# Conformer parameters that appear in head3.lst: ne, Em0, nH, pKa0, rxn\n")
    for label in all_labels:
        merged.extend(parsed[label]["conformer"])

    # ── Warnings as comments at bottom ──
    if warnings:
        merged.append(f"\n# {'='*60}\n")
        merged.append("# WARNINGS from MCCE4 Topology Agent\n")
        merged.append(f"# {'='*60}\n")
        for w in warnings:
            # Wrap long warnings
            for i in range(0, len(w), 70):
                merged.append(f"# {w[i:i+70]}\n")

    with open(output_path, "w") as f:
        f.writelines(merged)

    logging.info(f"  Merged ftpl: {output_path}")
    logging.info(f"    States: {', '.join(all_labels)}")
    return True


VALID_HYBRIDIZATIONS = {"s", "sp", "sp2", "sp3", "d2sp3"}


def ensure_connect_hybridizations(ftpl_path, pdb_path=None, lig_id=None):
    """Ensure every CONNECT line in the ftpl has a valid hybridization.

    Reads the ftpl, checks each CONNECT line for a valid hybridization
    (s, sp, sp2, sp3, d2sp3).  Lines with missing or invalid hybridization
    (e.g. "ion", "udf", or blank) are corrected:
      - Using RDKit-perceived hybridization when available
      - Otherwise inferred from the number of bonded atoms listed on the line

    Rewrites the file in-place.  Returns the number of lines corrected.
    """
    # Try to build an RDKit hybridization map from the PDB
    rdkit_hybs = {}
    if pdb_path:
        try:
            from .rdkit_tools import get_mol_from_pdb, RDKIT_TO_MCCE_HYB
            mol = get_mol_from_pdb(pdb_path, remove_hs=False)
            if mol:
                for atom in mol.GetAtoms():
                    pdb_info = atom.GetPDBResidueInfo()
                    if pdb_info:
                        name = pdb_info.GetName().strip()
                        hyb = str(atom.GetHybridization()).split(".")[-1]
                        rdkit_hybs[name] = RDKIT_TO_MCCE_HYB.get(hyb, hyb.lower())
        except Exception as e:
            logging.debug(f"  RDKit hybridization lookup skipped: {e}")

    # Pattern to parse CONNECT lines:
    #   CONNECT, " C1 ", EMH01:  sp2, " N1 "," C2 "
    # Group 1: prefix up to and including the colon
    # Group 2: hybridization field (may be invalid)
    # Group 3: rest of line (bonded atom list)
    connect_pat = re.compile(
        r'^(CONNECT,\s*"(.{4})",\s*\S+:\s*)'   # prefix + atom name (grp2)
        r'(\S+)'                                 # hybridization field
        r'(,\s*.*)$'                             # bonded atoms
    )
    # Also catch CONNECT lines where hybridization might be entirely missing
    # e.g. CONNECT, " X  ", EMH01: , " C1 "
    connect_nohyb_pat = re.compile(
        r'^(CONNECT,\s*"(.{4})",\s*\S+:\s*)'   # prefix + atom name
        r'(,\s*".*)$'                            # bonded atoms (starts with comma)
    )

    with open(ftpl_path) as f:
        lines = f.readlines()

    fixed = 0
    new_lines = []
    for line in lines:
        m = connect_pat.match(line)
        if m:
            prefix, atom_name_raw, hyb_field, bonded = m.groups()
            atom_name = atom_name_raw.strip()
            hyb_clean = hyb_field.strip().lower()

            if hyb_clean not in VALID_HYBRIDIZATIONS:
                new_hyb = _infer_hybridization(atom_name, bonded, rdkit_hybs)
                # Preserve original field width (4 chars right-aligned)
                new_lines.append(f"{prefix}{new_hyb:>4s}{bonded}\n")
                logging.info(f"  Fixed hybridization: {atom_name} '{hyb_field}' → '{new_hyb}'")
                fixed += 1
                continue

        m2 = connect_nohyb_pat.match(line)
        if m2 and not connect_pat.match(line):
            prefix, atom_name_raw, bonded = m2.groups()
            atom_name = atom_name_raw.strip()
            new_hyb = _infer_hybridization(atom_name, bonded, rdkit_hybs)
            new_lines.append(f"{prefix}{new_hyb:>4s}{bonded}\n")
            logging.info(f"  Fixed hybridization: {atom_name} (missing) → '{new_hyb}'")
            fixed += 1
            continue

        new_lines.append(line)

    if fixed:
        with open(ftpl_path, "w") as f:
            f.writelines(new_lines)
        logging.info(f"  Corrected {fixed} CONNECT hybridization(s) in {ftpl_path}")

    return fixed


def _infer_hybridization(atom_name, bonded_str, rdkit_hybs):
    """Infer hybridization for an atom from RDKit data or bond count."""
    # Prefer RDKit if available
    if atom_name in rdkit_hybs:
        hyb = rdkit_hybs[atom_name]
        if hyb in VALID_HYBRIDIZATIONS:
            return hyb

    # Count bonded atoms from the CONNECT line
    num_bonds = bonded_str.count('"') // 2

    # Hydrogen is always s
    if atom_name.lstrip().startswith("H"):
        return "s"

    if num_bonds >= 4:
        return "sp3"
    elif num_bonds == 3:
        return "sp2"
    elif num_bonds == 2:
        return "sp2"
    elif num_bonds == 1:
        return "sp3"
    else:
        return "sp3"  # isolated atom fallback


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
                filled += 1
                matched = True

            # Try stripped match
            if not matched:
                stripped = atom_name.strip()
                mcce = to_mcce_atom_name(stripped)
                for cn, cv in charges.items():
                    if cn.strip() == stripped or to_mcce_atom_name(cn) == mcce:
                        new_lines.append(f"{prefix}{cv:7.3f} # auto-filled{suffix}\n")
                        filled += 1
                        matched = True
                        break

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
    from ..models import sort_conformer_labels
    labels = sort_conformer_labels(state_pdbs.keys())
    spacing = 25.0  # Angstroms between states

    lines = [
        f"# PyMOL script: {lig_id} protonation state comparison",
        "# Generated by mcce_ftpl",
        f"# Usage: pymol {Path(output_path).name}",
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
        # h_removed = diff.get("removed", [])

        lines.append(f"# ── State {label} ──")
        lines.append(f"load {Path(pdb_path).resolve()!s}, {obj_name}")
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
        # charge_sign = "+" if "p" in obj_name or "+" in label else ("-" if "m" in obj_name or "-" in label else "")
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
        "# Green H = added vs neutral, Red site = H removed vs neutral",
    ])

    with open(output_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    logging.info(f"  📄 PyMOL script: {output_path}")
    return output_path
