"""
visualize.py - Generate PyMOL .pml scripts for dipole visualization.

CGO arrow = CYLINDER (shaft) + SPHERE (tip marker).
Uses only the most basic CGO primitives — guaranteed to work in
every PyMOL version including Educational and Open-Source builds.

Color scheme:
  - Backbone dipole:   blue    (0.3, 0.5, 1.0)
  - Ionizable dipole:  red     (1.0, 0.2, 0.2)
  - Full dipole:       green   (0.2, 0.9, 0.3)
  - Quadrupole axes:   orange  (1.0, 0.7, 0.1)
"""

import numpy as np
from . import E_ANG_TO_DEBYE

COLORS = {
    "backbone":   (0.3, 0.5, 1.0),
    "ionizable":  (1.0, 0.2, 0.2),
    "full":       (0.2, 0.9, 0.3),
    "quadrupole": (1.0, 0.7, 0.1),
}


def _normalize(v):
    mag = np.linalg.norm(v)
    if mag < 1e-10:
        return np.zeros(3)
    return v / mag


def _cgo_arrow(origin, direction, magnitude, color,
               shaft_radius=0.6, tip_radius=1.4,
               scale=0.1, min_length=8.0):
    """
    CGO arrow: CYLINDER shaft + SPHERE tip.
    Works in every PyMOL version.
    """
    length = max(magnitude * scale, min_length)
    tip = origin + direction * length
    r, g, b = color

    cgo = "[\n"
    cgo += f"    CYLINDER,\n"
    cgo += f"    {origin[0]:.3f}, {origin[1]:.3f}, {origin[2]:.3f},\n"
    cgo += f"    {tip[0]:.3f}, {tip[1]:.3f}, {tip[2]:.3f},\n"
    cgo += f"    {shaft_radius:.3f},\n"
    cgo += f"    {r:.3f}, {g:.3f}, {b:.3f},\n"
    cgo += f"    {r:.3f}, {g:.3f}, {b:.3f},\n"
    cgo += f"    COLOR, {r:.3f}, {g:.3f}, {b:.3f},\n"
    cgo += f"    SPHERE, {tip[0]:.3f}, {tip[1]:.3f}, {tip[2]:.3f}, {tip_radius:.3f},\n"
    cgo += "]"
    return cgo


# ---------------------------------------------------------------------------
# PyMOL toggle commands
# ---------------------------------------------------------------------------
_PYMOL_TOGGLE_COMMANDS = r'''
_dipole_objects = []
_quad_objects = []

def _register_dipole(name):
    _dipole_objects.append(name)

def _register_quad(name):
    _quad_objects.append(name)

def show_backbone(quiet=0):
    cmd.enable("backbone_dipole")
    if not int(quiet): print("  Showing backbone dipole (blue)")
cmd.extend("show_backbone", show_backbone)

def hide_backbone(quiet=0):
    cmd.disable("backbone_dipole")
    if not int(quiet): print("  Hidden backbone dipole")
cmd.extend("hide_backbone", hide_backbone)

def show_ionizable(quiet=0):
    cmd.enable("ionizable_dipole")
    if not int(quiet): print("  Showing ionizable dipole (red)")
cmd.extend("show_ionizable", show_ionizable)

def hide_ionizable(quiet=0):
    cmd.disable("ionizable_dipole")
    if not int(quiet): print("  Hidden ionizable dipole")
cmd.extend("hide_ionizable", hide_ionizable)

def show_full(quiet=0):
    cmd.enable("full_dipole")
    if not int(quiet): print("  Showing full dipole (green)")
cmd.extend("show_full", show_full)

def hide_full(quiet=0):
    cmd.disable("full_dipole")
    if not int(quiet): print("  Hidden full dipole")
cmd.extend("hide_full", hide_full)

def show_quadrupole(quiet=0):
    for obj in _quad_objects:
        cmd.enable(obj)
    if not int(quiet): print("  Showing quadrupole axes (orange)")
cmd.extend("show_quadrupole", show_quadrupole)

def hide_quadrupole(quiet=0):
    for obj in _quad_objects:
        cmd.disable(obj)
    if not int(quiet): print("  Hidden quadrupole axes")
cmd.extend("hide_quadrupole", hide_quadrupole)

def show_all_dipoles():
    for obj in _dipole_objects + _quad_objects:
        cmd.enable(obj)
    print("  Showing all dipoles and quadrupole axes")
cmd.extend("show_all_dipoles", show_all_dipoles)

def hide_all_dipoles():
    for obj in _dipole_objects + _quad_objects:
        cmd.disable(obj)
    print("  Hidden all dipoles and quadrupole axes")
cmd.extend("hide_all_dipoles", hide_all_dipoles)

def solo_backbone():
    hide_all_dipoles()
    show_backbone(quiet=1)
    print("  Solo: backbone dipole (blue)")
cmd.extend("solo_backbone", solo_backbone)

def solo_ionizable():
    hide_all_dipoles()
    show_ionizable(quiet=1)
    print("  Solo: ionizable dipole (red)")
cmd.extend("solo_ionizable", solo_ionizable)

def solo_full():
    hide_all_dipoles()
    show_full(quiet=1)
    print("  Solo: full protein dipole (green)")
cmd.extend("solo_full", solo_full)

def solo_quadrupole():
    hide_all_dipoles()
    show_quadrupole(quiet=1)
    print("  Solo: quadrupole axes (orange)")
cmd.extend("solo_quadrupole", solo_quadrupole)

def dipole_help():
    print("")
    print("  +------------------------------------------------------+")
    print("  |          MCCE4 Dipole Visualization Controls          |")
    print("  +------------------------------------------------------+")
    print("  |  show_backbone  / hide_backbone    (blue)             |")
    print("  |  show_ionizable / hide_ionizable   (red)              |")
    print("  |  show_full      / hide_full        (green)            |")
    print("  |  show_quadrupole/ hide_quadrupole  (orange)           |")
    print("  |  show_all_dipoles / hide_all_dipoles                  |")
    print("  +------------------------------------------------------+")
    print("  |  solo_backbone  / solo_ionizable / solo_full          |")
    print("  |  solo_quadrupole                                      |")
    print("  +------------------------------------------------------+")
    print("  |  dipole_help      -- show this menu                   |")
    print("  +------------------------------------------------------+")
    print("")
cmd.extend("dipole_help", dipole_help)
'''


def generate_pymol_script(protein_pdb, results, output_path, ph_index=None,
                          arrow_scale=0.1):
    """
    Generate a .pml script for PyMOL visualization.

    All 4 dipole/quadrupole moments are created and visible by default.
    Toggle in PyMOL with 'dipole_help'.
    """
    center = results["center"]

    if ph_index is not None:
        mu_bb = results["backbone_dipole"][ph_index]
        mu_ion = results["ionizable_dipole"][ph_index]
        mu_full = results["full_dipole"][ph_index]
        Q_eigvals = results["quadrupole_eigenvalues"][ph_index]
        if "quadrupole_eigenvectors" in results:
            Q_eigvecs = results["quadrupole_eigenvectors"]
        else:
            Q = results["quadrupole_tensor"][ph_index]
            Q_eigvals, Q_eigvecs = np.linalg.eigh(Q)
        ph_label = f"pH {results['ph_values'][ph_index]:.1f}"
        net_q = results["net_charge"][ph_index]
    else:
        mu_bb = results["backbone_dipole"]
        mu_ion = results["ionizable_dipole"]
        mu_full = results["full_dipole"]
        Q_eigvals = results["quadrupole_eigenvalues"]
        Q_eigvecs = results["quadrupole_eigenvectors"]
        ph_label = "single state"
        net_q = results["net_charge"]

    mag_bb = np.linalg.norm(mu_bb)
    mag_ion = np.linalg.norm(mu_ion)
    mag_full = np.linalg.norm(mu_full)

    L = []
    L.append("# ============================================================")
    L.append(f"# MCCE4 Dipole Moment Visualization -- {ph_label}")
    L.append(f"# Generated by mcce_dipole v1.0")
    L.append("# Type 'dipole_help' in PyMOL for toggle commands")
    L.append("# ============================================================")
    L.append("")
    L.append(f"# Backbone dipole:   {mag_bb:.1f} D")
    L.append(f"# Ionizable dipole:  {mag_ion:.1f} D")
    L.append(f"# Full dipole:       {mag_full:.1f} D")
    L.append(f"# Net charge:        {net_q:.2f} e")
    L.append("")
    L.append("from pymol.cgo import *")
    L.append("from pymol import cmd")
    L.append("")

    # Toggle commands
    L.append(_PYMOL_TOGGLE_COMMANDS)

    # Load protein — show cartoon only, hide everything else
    L.append(f'cmd.load("{protein_pdb}", "protein")')
    L.append('cmd.hide("everything", "protein")')
    L.append('cmd.show("cartoon", "protein")')
    L.append('cmd.set("cartoon_fancy_helices", 1)')
    L.append('cmd.color("palegreen", "protein")')
    L.append('cmd.color("paleyellow", "protein and ss s")')
    L.append('cmd.set("cartoon_transparency", 0.15)')
    L.append('cmd.bg_color("black")')
    L.append("")

    # Center marker
    L.append(f'cmd.pseudoatom("prot_center", pos=[{center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}])')
    L.append('cmd.show("spheres", "prot_center")')
    L.append('cmd.set("sphere_scale", 0.8, "prot_center")')
    L.append('cmd.color("white", "prot_center")')
    L.append("")

    # Dipole arrows
    dipoles = [
        ("backbone_dipole",  mu_bb,  mag_bb,  COLORS["backbone"]),
        ("ionizable_dipole", mu_ion, mag_ion, COLORS["ionizable"]),
        ("full_dipole",      mu_full, mag_full, COLORS["full"]),
    ]

    for obj_name, mu_vec, mag, color in dipoles:
        if mag < 0.1:
            L.append(f"# {obj_name}: magnitude too small ({mag:.2f} D), skipping")
            continue

        direction = _normalize(mu_vec)
        cgo_str = _cgo_arrow(center, direction, mag, color, scale=arrow_scale)

        L.append(f"# {obj_name}: {mag:.1f} D")
        L.append(f"{obj_name}_cgo = {cgo_str}")
        L.append(f'cmd.load_cgo({obj_name}_cgo, "{obj_name}")')
        L.append(f'_register_dipole("{obj_name}")')
        L.append("")

    # Quadrupole axes (cylinders only)
    if Q_eigvecs is not None:
        L.append("# Quadrupole principal axes")
        for i in range(3):
            eigval = Q_eigvals[i]
            eigvec = Q_eigvecs[:, i]
            axis_length = abs(eigval) * 0.02
            if axis_length < 3.0:
                axis_length = 3.0

            direction = _normalize(eigvec)
            r, g, b = COLORS["quadrupole"]

            start = center - direction * axis_length
            end = center + direction * axis_length

            L.append(f"quad_axis_{i}_cgo = [")
            L.append(f"    CYLINDER,")
            L.append(f"    {start[0]:.3f}, {start[1]:.3f}, {start[2]:.3f},")
            L.append(f"    {end[0]:.3f}, {end[1]:.3f}, {end[2]:.3f},")
            L.append(f"    0.35,")
            L.append(f"    {r:.3f}, {g:.3f}, {b:.3f},")
            L.append(f"    {r:.3f}, {g:.3f}, {b:.3f},")
            L.append(f"]")
            L.append(f'cmd.load_cgo(quad_axis_{i}_cgo, "quad_axis_{i}")')
            L.append(f'_register_quad("quad_axis_{i}")')

        L.append("")

    # Group
    all_obj = []
    for obj_name, _, mag, _ in dipoles:
        if mag >= 0.1:
            all_obj.append(obj_name)
    all_obj.extend([f"quad_axis_{i}" for i in range(3)])
    if all_obj:
        L.append(f'cmd.group("dipoles", "{" ".join(all_obj)}")')
        L.append("")

    # Rendering
    L.append('cmd.set("ray_shadow", "on")')
    L.append('cmd.set("antialias", 2)')
    L.append('cmd.set("specular", 0.3)')
    L.append('cmd.center("protein")')
    L.append('cmd.zoom("protein", buffer=10)')
    L.append("")

    # Print summary
    L.append(f'print("")')
    L.append(f'print("  MCCE4 Dipole Visualization -- {ph_label}")')
    L.append(f'print("  -----------------------------------------")')
    L.append(f'print("  Backbone:   {mag_bb:7.1f} D  (blue)")')
    L.append(f'print("  Ionizable:  {mag_ion:7.1f} D  (red)")')
    L.append(f'print("  Full:       {mag_full:7.1f} D  (green)")')
    L.append(f'print("  Net charge: {net_q:7.2f} e")')
    L.append(f'print("")')
    L.append(f'print("  Type dipole_help for toggle commands")')
    L.append(f'print("")')

    with open(output_path, "w") as fh:
        fh.write("\n".join(L) + "\n")
    return output_path


def generate_ph_scan_csv(results, output_path):
    """Write pH scan CSV."""
    ph = results["ph_values"]
    with open(output_path, "w") as fh:
        fh.write("pH,backbone_D,ionizable_D,full_D,net_charge,Q_eig1,Q_eig2,Q_eig3\n")
        for i in range(len(ph)):
            mag_bb = np.linalg.norm(results["backbone_dipole"][i])
            mag_ion = np.linalg.norm(results["ionizable_dipole"][i])
            mag_full = np.linalg.norm(results["full_dipole"][i])
            nq = results["net_charge"][i]
            qe = results["quadrupole_eigenvalues"][i]
            fh.write(f"{ph[i]:.1f},{mag_bb:.2f},{mag_ion:.2f},{mag_full:.2f},"
                     f"{nq:.3f},{qe[0]:.2f},{qe[1]:.2f},{qe[2]:.2f}\n")
    return output_path
