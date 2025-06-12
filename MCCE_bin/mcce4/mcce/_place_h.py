# Place missing H atoms based on ftpl file CONNECT, use 1-4 vdw to place the first H at torsion minimum
import logging
from mcce4.geom import *
from mcce4.pdbio import *

bond_length = 1.09  # H-C bond length for all added H atoms
bond_angle_sp3 = 109.5  # H-C-H bond angle for sp3 hybridized atom
bond_angle_sp2 = 120.0  # H-C-H bond angle for sp2 hybridized atom
bond_angle_sp = 180.0  # H-C-H bond angle for sp hybridized atom


atom_mass = {" H" : 1,
             " C" : 12,
                " N" : 14,
                " O" : 16,
                " S" : 32,
                " P" : 31,
                " F" : 19,
                "CL": 35.5,
                "BR": 80,
                "I" : 127,
                "ZN": 65,
                "CA": 40,
                "MG": 24,
                "MN": 55,
                "FE": 56,
                "CU": 63,
                "CO": 59,
                "NI": 59,
                " B": 10,
                "SI": 28,
                "LI": 7,
                "BE": 9,
                "NA": 23,
                "AL": 27,
                " K": 39}


def sp3_2known(r0, r1, r2):
    """
    Place 2 H (r3 and r4) atoms to r0 when r1 and r2 are known
    """
    #	r3 r4
    #	  \/
    #	   r0 - r2
    #	   |
    #	   r1
    #
    # r0 is sp3 type.
    # r0, r1, r2's coordinates are known.
    # r3, r4's coordinates are to be determined.
    #
    # Solution: r3, r4 are placed at the bisector of r1 and r2
    n01 = vector_normalize(vector_vminusv(r1, r0))
    n02 = vector_normalize(vector_vminusv(r2, r0))
    bisect = vector_normalize(vector_scale(vector_vplusv(n01, n02), -1))  # bisector of n01 and n02, opposite direction
    norm102 = vector_normalize(vector_vxv(n01, n02))
    half_angle = math.radians(bond_angle_sp3 / 2)

    cos_half_angle = math.cos(half_angle)
    sin_half_angle = math.sin(half_angle)

    # r3 = r0 + bond_len * (bisect * cos(half_angle) + norm102 * sin(half_angle));
    # r4 = r0 + bond_len * (bisect * cos(half_angle) - norm102 * sin(half_angle));
    r3 = vector_vplusv(r0, vector_scale(vector_vplusv(vector_scale(bisect, cos_half_angle), vector_scale(norm102, sin_half_angle)), bond_length))
    r4 = vector_vplusv(r0, vector_scale(vector_vminusv(vector_scale(bisect, cos_half_angle), vector_scale(norm102, sin_half_angle)), bond_length))

    return (r3, r4)




def sp3_1known(r0, r1, r2):
    """
    Place 3 H (r3, r4, and r5) atoms to r0
    """
    #	r3 r4 r5             y
    #	  \|/               /
    #	   r0			x--0
    #	   |			   |
    #	   r1			   z (v10)
    #	   \			    \  
    #       r2               r2
    # r0 is sp3 type.
    # r0, r1, r2's coordinates are known.
    # r3, r4, r5's coordinates are to be determined.
    #
    # Solution:
    # 1. Calculate v01, this is the direction of r1 from r0
    # 2. Calculate uz = n01/|n01|, this is the unit vector of n01
    # 3a. If r2 is not defined, pick an arbitrary orthogonal vector to uz as uy
    # 3b. Calculate v012, this is the direction of plane r0, r1, r2
    #     Calculate uy = v012/|v012|, this is the unit vector of v012
    # 4. Calculate ux = uy x uz, this is the unit vector of uy cross uz
    # 5. ux, uy, uz is the new coordinate system, with uz (r01's unit vector) as z axis and other two as x and y axes
    # 6. Calculate r3, r4, r5 in the new coordinate system
    # 7. return r3, r4, r5 in the original coordinate system

    # 1. Calculate v01, this is the direction of r1 from r0
    v01 = vector_vminusv(r1, r0)
    # 2. Calculate uz = n01/|n01|, this is the unit vector of n01
    uz = vector_normalize(v01)

    if r2 is None:    # 3a. If r2 is not defined, pick an arbitrary orthogonal vector to uz
        uy = vector_normalize(vector_orthogonal(uz))   # pick an arbitrary orthogonal vector to uz
    else:            # 3b. Calculate v012, this is the normal direction of plane r0, r1, r2
        v012 = vector_vxv(vector_vminusv(r2, r1), vector_vminusv(r1, r0))
        uy = vector_normalize(v012)
    # 4. Calculate ux = uy x uz, this is the unit vector of uy cross uz
    ux = vector_normalize(vector_vxv(uy, uz))
    # 5. ux, uy, uz is the new coordinate system, with uz (r01's unit vector) as z axis and other two as x and y axes
    # 6. Calculate r3, r4, r5 in the new coordinate system
    cos_theta = math.cos(math.radians(bond_angle_sp3))
    sin_theta = math.sin(math.radians(bond_angle_sp3))
    u3 = vector_normalize(vector_vplusv(vector_scale(uz, cos_theta), vector_scale(ux, sin_theta)))
    r3 = vector_vplusv(r0, vector_scale(u3, bond_length))
    # r4
    x_component = vector_scale(ux, sin_theta * cos_theta)
    y_component = vector_scale(uy, sin_theta * sin_theta)
    z_component = vector_scale(uz, cos_theta)
    r4 = vector_vplusv(r0, vector_scale(vector_sum3v(x_component, y_component, z_component), bond_length))
    # r5
    x_component = vector_scale(ux, sin_theta * cos_theta)
    y_component = vector_scale(uy, -sin_theta * sin_theta)
    # z_component = vector_scale(uz, cos_theta)
    r5 = vector_vplusv(r0, vector_scale(vector_sum3v(x_component, y_component, z_component), bond_length))

    return r3, r4, r5


def sp3_0known(r0):
    """
    Place 3 H (r1, r2, r3) atoms to r0
    """
    #	r1  r2
    #	  \ /
    #	   r0--r3
    #	   |
    #	  r4
    #
    # r0 is sp3 type.
    # r0's coordinates are known.
    # r1, r2, r3 and r4's coordinates are to be determined.
    #
    # Solution: r1, r2, r3, r4 are placed at the vertices of a regular tetrahedron with r0 as the center (0,0,0)
    v1 = (1, 1, 1)
    v2 = (1, -1, -1)
    v3 = (-1, 1, -1)
    v4 = (-1, -1, 1)
    r1 = vector_vplusv(r0, vector_scale(vector_normalize(v1), bond_length))
    r2 = vector_vplusv(r0, vector_scale(vector_normalize(v2), bond_length))
    r3 = vector_vplusv(r0, vector_scale(vector_normalize(v3), bond_length))
    r4 = vector_vplusv(r0, vector_scale(vector_normalize(v4), bond_length))
    return (r1, r2, r3, r4)


def sp2_2known(r0, r1, r2):
    """
    Place 1 H (r3) atom to r0 when r1 and r2 are known
    """
    #	r3
    #	  \
    #	   r0 - r2
    #	   |
    #	   r1
    #
    # r0 is sp2 type.
    # r0, r1, r2's coordinates are known.
    # r3's coordinates are to be determined.
    #
    # Solution: r3 is placed at the bisector of r1 and r2
    bisect = vector_normalize(vector_vplusv(vector_normalize(vector_vminusv(r1, r0)), vector_normalize(vector_vminusv(r2, r0))))
    r3 = vector_vplusv(r0, vector_scale(bisect, -bond_length))  # place H atom at the opposite side of 3 known atoms
    return r3


def sp2_1known_(r0, r1, r2):
    """
    Place 2 H (r3 and r4) atoms to r0 when r1 and r2 are known
    """
    #	r3 r4
    #	  \/
    #	   r0
    #	   |
    #	   r1
    #	   /
    #	  r2
    #
    # r0 is sp2 type.
    # r0, r1, r2's coordinates are known.
    # r3, r4's coordinates are to be determined.
    #
    # Solution: r3, r4 are placed on the same plane defined by r0, r1, and r2, with bond angle of 120 degrees
    # 1. Calculate v01, this is the direction of r1 from r0
    v01 = vector_vminusv(r1, r0)
    # 2. Calculate v12, this is the direction of r2 from r1
    if r2 is None:
        v12 = vector_orthogonal(v01)
    else:
        v12 = vector_vminusv(r2, r1)
    # 3. Calculate the normal of plane n012, this is the cross product of v01 and v12
    n012 = vector_normalize(vector_vxv(v01, v12))
    ux = vector_normalize(vector_vxv(n012, v01))    # ux is the unit vector of v01 cross n012
    uy = vector_normalize(vector_scale(v01, -1))    # uy is the opposite direction of v01
    # 4. Calculate r3, r4 in the new coordinate system
    cos_theta = math.cos(math.radians(bond_angle_sp2-90))   # bond angle is 120 degrees, but we need the angle between ux and u4
    sin_theta = math.sin(math.radians(bond_angle_sp2-90))
    u4 = vector_normalize(vector_vplusv(vector_scale(ux, cos_theta), vector_scale(uy, sin_theta)))
    u3 = vector_normalize(vector_vplusv(vector_scale(ux, -cos_theta), vector_scale(uy, sin_theta)))
    r4 = vector_vplusv(r0, vector_scale(u4, bond_length))
    r3 = vector_vplusv(r0, vector_scale(u3, bond_length))

    return r3, r4

def sp2_0known(r0):
    """
    Place 3 H (r1, r2 and r3) atoms to r0
    """
    #	r1  r2
    #	  \ /
    #	   r0
    #	   |
    #	  r3
    #
    # r0 is sp2 type.
    # r0's coordinates are known.
    # r1, r2, r3's coordinates are to be determined.
    #
    # Solution: r1, r2, r3 are placed at the vertices of an equilateral triangle with r0 as the center (0,0,0)
    v1 = (1, math.sqrt(3), 0)
    v2 = (1, -math.sqrt(3), 0)
    v3 = (-2, 0, 0)
    r1 = vector_vplusv(r0, vector_scale(vector_normalize(v1), bond_length))
    r2 = vector_vplusv(r0, vector_scale(vector_normalize(v2), bond_length))
    r3 = vector_vplusv(r0, vector_scale(vector_normalize(v3), bond_length))
    return (r1, r2, r3)


def create_h(atom, hname, xyz):
    h = Atom()
    h.inherit(atom)
    h.name = hname
    h.element = " H"
    h.xyz = xyz
    h.connect12.append(atom)
    atom.connect12.append(h)
    return h


def place_h(self):
    self.reset_connect()
    self.make_connect12()

    for res in self.protein.residue:
        for conf in res.conf:
            for atom in conf.atom:
                if atom.element != " H":
                    # find connected H atoms on the heavy atom
                    should_have_h = [atom_name for atom_name in self.tpl.db[("CONNECT", atom.name, conf.confType)].connected if is_H(atom_name)]
                    actually_have_h = [a.name for a in atom.connect12 if is_H(a.name)]
                    missing_h = set(should_have_h) - set(actually_have_h)
                    if missing_h:
                        # print("%s: %s - %s = %s" % (atom.name, str(should_have_h), str(actually_have_h), str(missing_h)))
                        # center atom orbital type
                        orbital = self.tpl.db[("CONNECT", atom.name, conf.confType)].orbital
                        # connected atoms
                        connected_heavy_atoms = [a for a in atom.connect12 if not is_H(a.name)]
                        # print("Connected heavy atoms by name: %s" % [a.atomID for a in connected_heavy_atoms])
                        # print("Connected heavy atoms by obj:", [a for a in connected_heavy_atoms])
                        # Here is an edge case where this atom connects to multiple conformers in NTR and CTR
                        connected_heavy_atoms_name = []
                        unique_connected_heavy_atoms = []
                        for a in connected_heavy_atoms:
                            if a.name not in connected_heavy_atoms_name:
                                unique_connected_heavy_atoms.append(a)
                                connected_heavy_atoms_name.append(a.name)
                        connected_heavy_atoms = unique_connected_heavy_atoms    # replace with the unique list

                        # decide placing methods
                        if orbital.lower() == "ion":
                            continue    # do nothing for ions
                        elif orbital.lower() == "sp3":  # sp3 hybridized atom 
                            n_known = len(connected_heavy_atoms)
                            if n_known == 3:
                                #    r4
                                #     |
                                #    r0
                                #    /|\
                                # r1 r2 r3
                                #
                                # r0 is sp3 type.
                                # r0, r1, r2, r3's coordinates are known.
                                # r4's coordinate is to be determined.
                                if len(missing_h) == 1:    # place the 4th H atom at a defined place
                                    # place H atom at the opposite side of 3 known atoms
                                    n01 = vector_normalize(vector_vminusv(connected_heavy_atoms[0].xyz, atom.xyz))
                                    n02 = vector_normalize(vector_vminusv(connected_heavy_atoms[1].xyz, atom.xyz))
                                    n03 = vector_normalize(vector_vminusv(connected_heavy_atoms[2].xyz, atom.xyz))
                                    n04 = vector_normalize(vector_sum3v(n01, n02, n03))
                                    xyz = vector_vplusv(atom.xyz, vector_scale(n04, -bond_length))  # place H atom at the opposite side of 3 known atoms
                                    h = create_h(atom, missing_h.pop(), xyz)
                                    conf.atom.append(h)
                                else:
                                    print("Missing h atoms %s on %s with connected heavy atoms %s" % (str(missing_h), atom.name, [a.name for a in connected_heavy_atoms]))
                                    logging.error("Not enough slots for H atoms on atom %s in conformer %s" % (atom.name, conf.confID))
                                    sys.exit()
                            elif n_known == 2:
                                # r3 r4
                                #   \/
                                #   r0
                                #   /|
                                # r1 r2
                                    
                                # r0 is sp3 type.
                                # r0, r1, r2's coordinates are known.
                                # r3, r4's coordinates are to be determined.
                                if len(missing_h) <= 2:     # place the 3rd and 4th H atoms at defined places
                                    r0 = atom.xyz
                                    r1 = connected_heavy_atoms[0].xyz
                                    r2 = connected_heavy_atoms[1].xyz
                                    r3, r4 = sp3_2known(r0, r1, r2)
                                    if len(missing_h) == 1:
                                        h = create_h(atom, missing_h.pop(), r3)
                                        conf.atom.append(h)
                                    elif len(missing_h) == 2:
                                        h = create_h(atom, missing_h.pop(), r3)
                                        conf.atom.append(h)
                                        h = create_h(atom, missing_h.pop(), r4)
                                        conf.atom.append(h)
                                else:
                                    print("Missing h atoms %s on %s with connected heavy atoms %s" % (str(missing_h), atom.name, [a.name for a in connected_heavy_atoms]))
                                    logging.error("Not enough slots for H atoms on atom %s in conformer %s" % (atom.name, conf.confID))
                                    sys.exit()

                            elif n_known == 1:
                                #	r3 r4 r5
                                #	  \|/
                                #	   r0
                                #	   |
                                #	   r1
                                #	   \
                                #	    r2
                                # r0 is sp3 type.
                                # r0, r1, r2's coordinates are known.
                                # r3, r4, r5's coordinates are to be determined.
                                torsion_minimum = False
                                if len(missing_h) <= 3:
                                    r0 = atom.xyz
                                    r1 = connected_heavy_atoms[0].xyz
                                    # pick r2, the heaviest atom connected to r1
                                    extended_atoms = list(set(connected_heavy_atoms[0].connect12) - {atom})  # remove r0
                                    if len(extended_atoms) == 1:  # only one connected
                                        r2 = extended_atoms[0].xyz
                                    elif len(extended_atoms) == 0:  # no extended atom, place r3, r4, r5 at arbitrary places
                                        r3, r4, r5 = sp3_1known(r0, r1, None)
                                        # print("Atom %s has no connected atom" % connected_heavy_atoms[0].name)
                                        # print("Place r3, r4, r5 at arbitrary places")
                                    else:
                                        # init atom mass
                                        for a in extended_atoms:
                                            if a.element:
                                                if a.element in atom_mass:
                                                    a.mass = atom_mass[a.element]
                                                else:
                                                    a.mass = 100  # unknown atom, assign a large mass
                                            else:
                                                print("Error: atom %s has no element" % a.name)
                                                sys.exit()
                                        r2_atom = extended_atoms[0]
                                        for a in extended_atoms[1:]:
                                            if a.mass > r2_atom.mass:
                                                r2_atom = a
                                        r2 = r2_atom.xyz
                                    # calculate r3, r4, r5 with r3 at torsion minimum
                                    torsion_minimum = True
                                    r3, r4, r5 = sp3_1known(r0, r1, r2)                                     

                                else:
                                    print("Missing h atoms %s on %s with connected heavy atoms %s" % (str(missing_h), atom.name, [a.name for a in connected_heavy_atoms]))
                                    logging.error("Not enough slots for H atoms on atom %s in conformer %s" % (atom.name, conf.confID))
                                    sys.exit()

                                # now we have r3, r4, r5, place H atoms
                                if len(missing_h) == 1:
                                    h = create_h(atom, missing_h.pop(), r3)
                                    conf.atom.append(h)
                                elif len(missing_h) == 2:
                                    h = create_h(atom, missing_h.pop(), r3)
                                    conf.atom.append(h)
                                    h = create_h(atom, missing_h.pop(), r4)
                                    conf.atom.append(h)
                                elif len(missing_h) == 3:
                                    h = create_h(atom, missing_h.pop(), r3)
                                    conf.atom.append(h)
                                    h = create_h(atom, missing_h.pop(), r4)
                                    conf.atom.append(h)
                                    h = create_h(atom, missing_h.pop(), r5)
                                    conf.atom.append(h)
                                
                                if torsion_minimum:
                                    conf.history = conf.history[:6] + "M" + conf.history[7:]
                                    torsion_minimum = False

                            elif n_known == 0:
                                # this applies to some special cases where the atom has no connected heavy atom, such as NH4+ and H2O
                                r0 = atom.xyz
                                r1, r2, r3, r4 = sp3_0known(r0)
                                if len(missing_h) == 1:
                                    h = create_h(atom, missing_h.pop(), r1)
                                    conf.atom.append(h)
                                elif len(missing_h) == 2:
                                    h = create_h(atom, missing_h.pop(), r1)
                                    conf.atom.append(h)
                                    h = create_h(atom, missing_h.pop(), r2)
                                    conf.atom.append(h)
                                elif len(missing_h) == 3:
                                    h = create_h(atom, missing_h.pop(), r1)
                                    conf.atom.append(h)
                                    h = create_h(atom, missing_h.pop(), r2)
                                    conf.atom.append(h)
                                    h = create_h(atom, missing_h.pop(), r3)
                                    conf.atom.append(h)
                                elif len(missing_h) == 4:
                                    h = create_h(atom, missing_h.pop(), r1)
                                    conf.atom.append(h)
                                    h = create_h(atom, missing_h.pop(), r2)
                                    conf.atom.append(h)
                                    h = create_h(atom, missing_h.pop(), r3)
                                    conf.atom.append(h)
                                    h = create_h(atom, missing_h.pop(), r4)
                                    conf.atom.append(h)
                                else:
                                    print("Missing h atoms %s on %s with connected heavy atoms %s" % (str(missing_h), atom.name, [a.name for a in connected_heavy_atoms]))
                                    logging.error("Not enough slots for H atoms on atom %s in conformer %s" % (atom.name, conf.confID))
                                    sys.exit()
                                    

                        elif orbital.lower() == "sp2":
                            n_known = len(connected_heavy_atoms)
                            if n_known == 2:   # place the 3rd H atom at a defined place
                                # r3
                                #  \
                                #   r0 - r2
                                #   |
                                #   r1
                                #
                                # r0 is sp2 type.
                                # r0, r1, r2's coordinates are known.
                                # r3's coordinates are to be determined.
                                if len(missing_h) == 1:
                                    r0 = atom.xyz
                                    r1 = connected_heavy_atoms[0].xyz
                                    r2 = connected_heavy_atoms[1].xyz
                                    r3 = sp2_2known(r0, r1, r2)
                                    h = create_h(atom, missing_h.pop(), r3)
                                    conf.atom.append(h)
                                else:
                                    print("Missing h atoms %s on %s with connected heavy atoms %s" % (str(missing_h), atom.name, [a.name for a in connected_heavy_atoms]))
                                    logging.error("Not enough slots for H atoms on atom %s in conformer %s" % (atom.name, conf.confID))
                                    sys.exit()
                            elif n_known == 1:  # place the 2nd and 3rd H atoms at defined places
                                # r3 r4
                                #  \/
                                #   r0
                                #   |
                                #   r1
                                #   /
                                #  r2
                                #
                                # r0 is sp2 type.
                                # r0, r1, r2's coordinates are known.
                                # r3, r4's coordinates are to be determined.
                                r0 = atom.xyz
                                r1 = connected_heavy_atoms[0].xyz
                                torsion_minimum = False
                                # pick r2, the heaviest atom connected to r1
                                extended_atoms = list(set(connected_heavy_atoms[0].connect12) - {atom})  # remove r0
                                if len(extended_atoms) == 1:  # only one connected
                                    r2 = extended_atoms[0].xyz
                                elif len(extended_atoms) == 0:  # no extended atom, place r3, r4 at arbitrary places
                                    r2 = None
                                else:
                                    # init atom mass
                                    for a in extended_atoms:
                                        if a.element:
                                            if a.element in atom_mass:
                                                a.mass = atom_mass[a.element]
                                            else:
                                                a.mass = 100  # unknown atom, assign a large mass
                                        else:
                                            print("Error: atom %s has no element" % a.name)
                                            sys.exit()
                                    r2_atom = extended_atoms[0]
                                    for a in extended_atoms[1:]:
                                        if a.mass > r2_atom.mass:
                                            r2_atom = a
                                    r2 = r2_atom.xyz
                                    torsion_minimum = True
                                # now that we have r2, calculate r3, r4
                                r3, r4 = sp2_1known_(r0, r1, r2)
                                if len(missing_h) == 1:
                                    h = create_h(atom, missing_h.pop(), r3)
                                    conf.atom.append(h)
                                elif len(missing_h) == 2:
                                    h = create_h(atom, missing_h.pop(), r3)
                                    conf.atom.append(h)
                                    h = create_h(atom, missing_h.pop(), r4)
                                    conf.atom.append(h)
                                else:
                                    print("Missing h atoms %s on %s with connected heavy atoms %s" % (str(missing_h), atom.name, [a.name for a in connected_heavy_atoms]))
                                    logging.error("Not enough slots for H atoms on atom %s in conformer %s" % (atom.name, conf.confID))
                                    sys.exit()
                                
                                if torsion_minimum:
                                    conf.history = conf.history[:6] + "M" + conf.history[7:]
                                    torsion_minimum = False
                            elif n_known == 0:
                                # r1  r2
                                #   \ /
                                #    r0
                                #    |
                                #    r3
                                #
                                # r0 is sp2 type.
                                # r0's coordinates are known.
                                # r1, r2, r3's coordinates are to be determined.
                                r0 = atom.xyz
                                r1, r2, r3 = sp2_0known(r0)
                                if len(missing_h) == 1:
                                    h = create_h(atom, missing_h.pop(), r1)
                                    conf.atom.append(h)
                                elif len(missing_h) == 2:
                                    h = create_h(atom, missing_h.pop(), r1)
                                    conf.atom.append(h)
                                    h = create_h(atom, missing_h.pop(), r2)
                                    conf.atom.append(h)
                                elif len(missing_h) == 3:
                                    h = create_h(atom, missing_h.pop(), r1)
                                    conf.atom.append(h)
                                    h = create_h(atom, missing_h.pop(), r2)
                                    conf.atom.append(h)
                                    h = create_h(atom, missing_h.pop(), r3)
                                    conf.atom.append(h)
                                else:
                                    print("Missing h atoms %s on %s with connected heavy atoms %s" % (str(missing_h), atom.name, [a.name for a in connected_heavy_atoms]))
                                    logging.error("Not enough slots for H atoms on atom %s in conformer %s" % (atom.name, conf.confID))
                                    sys.exit()

                        elif orbital.lower() == "sp":
                            n_known = len(connected_heavy_atoms)
                            if n_known == 1:  # place H on the oppsite side of the known atom
                                # r2
                                #  \
                                #   r0
                                #    \
                                #     r1
                                r0 = atom.xyz
                                r1 = connected_heavy_atoms[0].xyz
                                r2 = vector_vplusv(r0, vector_scale(vector_vminusv(r1, r0), -bond_length))
                                if len(missing_h) == 1:
                                    h = create_h(atom, missing_h.pop(), r2)
                                    conf.atom.append(h)
                            else:
                                print("Atom %s on %s must of exactly one connected heavy atom to place missing H, connected %s" % (atom.name, conf.confID, [heavy.name for heavy in connected_heavy_atoms]))
                                logging.error("Not enough infornation to place H on atom %s in conformer %s" % (atom.name, conf.confID))
                                sys.exit()
                            
                            
                        else:
                            logging.error("Unknown orbital type %s for atom %s in conformer" % (orbital, atom.name, conf.confType))
                            sys.exit()


