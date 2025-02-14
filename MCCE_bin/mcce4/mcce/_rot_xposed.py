import math
import logging
from ._rot_swing import swing
from ._strip_cofactors import DEFAULT_RAD
from ._strip_cofactors import probe_rad
from ._strip_cofactors import area_k
from ._strip_cofactors import radius
from ._strip_cofactors import fibonacci_sphere

Xposed_threshold = 0.20  # only optimize those more exposed than this cut off.
Non_polar_res = {"ALA", "VAL", "LEU", "ILE", "GLY", "MET", "MEL", "TRY", "PHE", "HIL", "CYD", "CYL"}  # non polar residue and ligands are in exclude list


def conf_sas(conf, reference, include_backbone=True):
    sas = 0
    sas_atoms = [a for a in conf.atom]  # avoid passing conf.atom to sas_atoms, otherwise conf.atom will be permenantly changed)
    outside_atoms = []
    if sas_atoms and conf.resID[0] not in Non_polar_res:  # only do the calculation for a non-empty object sas
        point_preset = fibonacci_sphere(122)
        # group backbone atoms of the residue to sas_atoms
        bk_confID = conf.confID[:3] + "BK" + conf.confID[5:-3]
        for atom in reference:
            if conf.resID == atom.resID:
                if include_backbone:
                    if atom.confID[:-3] == bk_confID:
                        sas_atoms.append(atom)
                    else:
                        outside_atoms.append(atom)
            else: # exclude atoms in the same residue of sas_atoms
                outside_atoms.append(atom)
                

        x_lo = x_hi = sas_atoms[0].xyz[0]
        y_lo = y_hi = sas_atoms[0].xyz[1]
        z_lo = z_hi = sas_atoms[0].xyz[2]
        if len(sas_atoms) > 1:
            for atom in sas_atoms[1:]:
                if atom.xyz[0] < x_lo:
                    x_lo = atom.xyz[0]
                elif atom.xyz[0] > x_hi:
                    x_hi = atom.xyz[0]
                if atom.xyz[1] < y_lo:
                    y_lo = atom.xyz[1]
                elif atom.xyz[1] > y_hi:
                    y_hi = atom.xyz[1]
                if atom.xyz[2] < z_lo:
                    z_lo = atom.xyz[2]
                elif atom.xyz[2] > z_hi:
                    z_hi = atom.xyz[2]

        # filter those reference atoms out of range, expanding by twice probe radius and atom radii
        atoms_inrange = []
        for atom in sas_atoms+outside_atoms:
            if atom.element in radius:
                atom.r = radius[atom.element]
            else:
                atom.r = DEFAULT_RAD

            if x_lo - 2*probe_rad - atom.r - 2 < atom.xyz[0] < x_hi + 2*probe_rad + atom.r + 2 and\
               y_lo - 2*probe_rad - atom.r - 2 < atom.xyz[1] < y_hi + 2*probe_rad + atom.r + 2 and\
               z_lo - 2*probe_rad - atom.r - 2 < atom.xyz[2] < z_hi + 2*probe_rad + atom.r + 2:
                atoms_inrange.append(atom)
        # Test the reduction
        # print("%s: All = %d, In-range = %d" % (conf.confID, len(reference), len(atoms_inrange)))
        # NTR01A0001_001: All = 1010, In-range = 69
        # LYS01A0001_001: All = 1010, In-range = 112
        # VAL01A0002_001: All = 1010, In-range = 93
        # PHE01A0003_001: All = 1010, In-range = 168
        # ARG01A0005_001: All = 1010, In-range = 159
        # CYD01A0006_001: All = 1010, In-range = 123
        # GLU01A0007_001: All = 1010, In-range = 130
        # LEU01A0008_001: All = 1010, In-range = 193
        # ALA01A0009_001: All = 1010, In-range = 163
        # ALA01A0010_001: All = 1010, In-range = 121
        # ALA01A0011_001: All = 1010, In-range = 128
        # MET01A0012_001: All = 1010, In-range = 212
        # LYS01A0013_001: All = 1010, In-range = 125
        # ARG01A0014_001: All = 1010, In-range = 127
        # HIS01A0015_001: All = 1010, In-range = 134
        # LEU01A0017_001: All = 1010, In-range = 191
        # ASP01A0018_001: All = 1010, In-range = 124

        # Be noted the atoms_inrange contains the test conformer atoms
        # sas of sas_atoms
        sas_atoms_inprotein = 0.0
        for atom in sas_atoms:
            n_points = len(point_preset)
            counter = n_points
            rad_ext = atom.r + probe_rad
            for p_raw in point_preset:
                point = (p_raw[0]*rad_ext + atom.xyz[0], p_raw[1]*rad_ext + atom.xyz[1], p_raw[2]*rad_ext + atom.xyz[2])
                for atom2 in atoms_inrange:
                    if atom2 != atom:
                        dx = point[0] - atom2.xyz[0]
                        dy = point[1] - atom2.xyz[1]
                        dz = point[2] - atom2.xyz[2]
                        dd = dx*dx + dy*dy + dz*dz
                        rad2_ext = atom2.r + probe_rad
                        if dd < rad2_ext * rad2_ext:
                            counter -= 1
                            break

            #atom.sas = area_k * rad_ext * rad_ext * counter / n_points
            sas_atoms_inprotein += area_k * rad_ext * rad_ext * counter / n_points

        # sas of naked sas_atoms
        sas_atoms_naked = 0.0
        for atom in sas_atoms:
            n_points = len(point_preset)
            counter = n_points
            rad_ext = atom.r + probe_rad
            for p_raw in point_preset:
                point = (p_raw[0]*rad_ext + atom.xyz[0], p_raw[1]*rad_ext + atom.xyz[1], p_raw[2]*rad_ext + atom.xyz[2])
                for atom2 in sas_atoms:
                    if atom2 != atom:
                        dx = point[0] - atom2.xyz[0]
                        dy = point[1] - atom2.xyz[1]
                        dz = point[2] - atom2.xyz[2]
                        dd = dx*dx + dy*dy + dz*dz
                        rad2_ext = atom2.r + probe_rad
                        if dd < rad2_ext * rad2_ext:
                            counter -= 1
                            break

            #atom.sas = area_k * rad_ext * rad_ext * counter / n_points
            sas_atoms_naked += area_k * rad_ext * rad_ext * counter / n_points

        sas = sas_atoms_inprotein / sas_atoms_naked
        # print(conf.confID, sas)

    return sas

def rot_xposed(self):
    """ Make a most exposed conformer for surface residues
    Basic algorithm:
    - select the backbone and the first conformer of each residue as the staring structure
    - loop over each residue
        - calculate sas% or the residue
        - if sas % > exposed_cutoff:
            - swing() and find a better sas
            - if a better one is found, add this conformer 
    """
    # Load atoms of selected structure as SAS reference
    sas_reference = []
    for res in self.protein.residue:
        sas_reference += res.conf[0].atom
        if len(res.conf) > 1:
            sas_reference += res.conf[1].atom

    # Assign radius
    for atom in sas_reference:
        if atom.element in radius:
            atom.r = radius[atom.element]
        else:
            atom.r = DEFAULT_RAD

    for res in self.protein.residue:
        if len(res.conf) > 1:
            start_conf = res.conf[1]
            sas_original = sas_max = conf_sas(start_conf, sas_reference)
            # print(start_conf.confID, sas_max)
            if sas_original > Xposed_threshold:
                improved = True
                phi = 60/180*math.pi
                while improved:
                    improved = False
                    new_confs = self.swing(start_conf, phi)
                    for conf in new_confs:
                        sas = conf_sas(conf, sas_reference)
                        if sas > sas_max:
                            improved = True
                            sas_max = sas
                            start_conf = conf

                improved = True
                phi = 15/180*math.pi
                while improved:
                    improved = False
                    new_confs = self.swing(start_conf, phi)
                    for conf in new_confs:
                        sas = conf_sas(conf, sas_reference)
                        if sas > sas_max:
                            improved = True
                            sas_max = sas
                            start_conf = conf

                improved = True
                phi = 3/180*math.pi
                while improved:
                    improved = False
                    new_confs = self.swing(start_conf, phi)
                    for conf in new_confs:
                        sas = conf_sas(conf, sas_reference)
                        if sas > sas_max:
                            improved = True
                            sas_max = sas
                            start_conf = conf

                improved = True
                phi = 1/180*math.pi
                while improved:
                    improved = False
                    new_confs = self.swing(start_conf, phi)
                    for conf in new_confs:
                        sas = conf_sas(conf, sas_reference)
                        if sas > sas_max:
                            improved = True
                            sas_max = sas
                            start_conf = conf


                if sas_max > sas_original + 0.001:  # significantly more exposed conformer found
                    start_conf.history = start_conf.history[:2] + "X" + start_conf.history[3:]
                    start_conf.sas = sas_max
                    res.conf.append(start_conf)
                    logging.info("   %s: sas = %.3f -> %.3f" % (res.resID, sas_original, sas_max))
