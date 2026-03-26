"""
Protonation state enumeration using Dimorphite-DL.
"""

import logging
from typing import List
from ..models import ConformerState
from ..config import CHARGE_TO_CONF


def enumerate_protonation_states(smiles: str, ph: float = 7.4,
                                  ph_tolerance: float = 1.0) -> List[ConformerState]:
    """Enumerate protonation states at target pH using Dimorphite-DL.

    Args:
        smiles: Input SMILES string.
        ph: Target pH.
        ph_tolerance: pH range tolerance (±).

    Returns:
        List of ConformerState objects for each unique charge state.
    """
    logging.info(f"🧪 Enumerating protonation states at pH {ph} (± {ph_tolerance})...")

    state_smiles = [smiles]  # fallback

    try:
        from dimorphite_dl.protonate import protonate_smiles
        results = protonate_smiles(
            smiles,
            min_ph=ph - ph_tolerance,
            max_ph=ph + ph_tolerance,
        )
        if results:
            state_smiles = list(set(results))
        logging.info(f"  ✓ Dimorphite-DL found {len(state_smiles)} state(s)")
    except ImportError:
        logging.warning("  Dimorphite-DL not installed (pip install dimorphite-dl)")
    except Exception as e:
        logging.warning(f"  Dimorphite-DL error: {e}")

    # Compute formal charges and build ConformerState objects
    states = _smiles_to_states(state_smiles, ph)

    for s in states:
        logging.info(f"  📋 {s.label}: charge={s.charge:+d}, nH={s.nH:+d}")

    return states


def _smiles_to_states(smiles_list: list, ph: float) -> List[ConformerState]:
    """Convert SMILES list to ConformerState objects with unique charges."""
    states = []
    seen_charges = set()

    try:
        from rdkit import Chem

        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            charge = Chem.GetFormalCharge(mol) if mol else 0

            if charge in seen_charges:
                continue
            seen_charges.add(charge)

            label = CHARGE_TO_CONF.get(charge, f"{charge:+d}"[-2:])
            nH = charge  # relative to neutral

            states.append(ConformerState(
                label=label, charge=charge,
                nH=nH if charge != 0 else 0,
                smiles=smi, source="dimorphite",
                rationale=f"Dimorphite-DL at pH {ph}"
            ))

    except ImportError:
        # No RDKit fallback
        states.append(ConformerState(
            label="01", charge=0, nH=0, smiles=smiles_list[0],
            source="fallback", rationale="No RDKit — assumed neutral"
        ))

    return states
