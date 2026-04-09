"""
mcce_dipole - Post-processing tool for MCCE4 dipole and quadrupole moment analysis.

Computes pH-dependent dipole moments from the MCCE4 microstate ensemble:
  1. Backbone-only dipole moment
  2. Quadrupole moment tensor
  3. Ionizable residue dipole moment
  4. Full protein dipole moment

Generates PyMOL CGO scripts for 3D visualization of dipole vectors.
"""

__version__ = "1.0.0"

E_ANG_TO_DEBYE = 4.80321  # 1 e·Å = 4.80321 Debye

IONIZABLE_RES = {"ASP", "GLU", "LYS", "ARG", "HIS", "CYS", "TYR", "NTR", "CTR"}
BACKBONE_ATOMS = {"N", "CA", "C", "O", "H", "HA"}
