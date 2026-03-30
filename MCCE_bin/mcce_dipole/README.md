# mcce_dipole — MCCE4 Dipole & Quadrupole Moment Analysis

Post-processing tool for MCCE4 (steps 1-4). Computes four electrostatic
moments from the conformer ensemble, Boltzmann-weighted at each pH.

## What It Computes

| # | Quantity | Atoms | Description |
|---|----------|-------|-------------|
| 1 | **Backbone Dipole** | N, CA, C, O, H, HA | Intrinsic polarity from aligned peptide bonds. Largely pH-independent. |
| 2 | **Ionizable Dipole** | ASP, GLU, LYS, ARG, HIS, CYS, TYR, NTR, CTR | pH-dependent — changes as titratable groups gain/lose protons. |
| 3 | **Full Protein Dipole** | All atoms | Total charge asymmetry at each pH. |
| 4 | **Quadrupole Tensor** | All atoms | Shape of charge distribution beyond dipole. Traceless 3x3, eigenvalues in e*A^2. |

## Quick Start

```bash
# From a completed MCCE4 run directory — just run it:
cd /path/to/mcce_run
ms_dipole.py

# Single microstate PQR:
ms_dipole.py --pqr state_0001.pqr
```

Required files (auto-detected in run directory):
- `step2_out.pdb` — atom coordinates + embedded partial charges
- `head3.lst` — conformer metadata
- `fort.38` — Monte Carlo occupancies at each pH

No `param/` directory needed — charges are read directly from step2_out.pdb.

## Output (always generated)

| File | Contents |
|------|----------|
| `dipole_ph_scan.csv` | Dipole magnitudes + quadrupole eigenvalues at every pH |
| `dipole_pH7_pymol.pml` | PyMOL script with all 4 moments as color-coded arrows |

## PyMOL Visualization Controls

All 4 moments are visible on load. Toggle interactively:

```
dipole_help              # show full menu

# Individual toggle
hide_backbone            show_backbone
hide_ionizable           show_ionizable
hide_full                show_full
hide_quadrupole          show_quadrupole

# Solo mode — show ONLY one
solo_backbone            solo_ionizable
solo_full                solo_quadrupole

# Bulk
show_all_dipoles         hide_all_dipoles
```

**Color code:** Blue = backbone, Red = ionizable, Green = full, Orange = quadrupole axes.

## Optional Flags

| Flag | Default | Description |
|------|---------|-------------|
| `--pqr FILE` | *(ensemble)* | Use a microstate PQR instead |
| `--dir DIR` | `.` | MCCE4 run directory |
| `--ph` | `7.0` | pH for PyMOL snapshot |
| `--pdb FILE` | step2_out.pdb | PDB for PyMOL display |
| `--arrow_scale` | `0.05` | Arrow length scaling |
| `-o PREFIX` | `dipole` | Output file prefix |

## Installation

```bash
cp ms_dipole.py /path/to/MCCE4-Alpha/MCCE_bin/
cp -r mcce_dipole/ /path/to/MCCE4-Alpha/MCCE_bin/
```

## File Structure

```
MCCE_bin/
├── ms_dipole.py          # control script (entry point)
├── mcce_dipole/           # package
│   ├── __init__.py        # constants
│   ├── __main__.py        # python -m mcce_dipole support
│   ├── parsers.py         # step2_out.pdb, head3.lst, fort.38, PQR readers
│   ├── compute.py         # dipole & quadrupole math
│   ├── visualize.py       # PyMOL CGO + toggle commands
│   └── README.md
```

## Dependencies

Python 3.8+, NumPy. No PyMOL needed to generate scripts.
