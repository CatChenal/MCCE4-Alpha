# Optimize H poisition for hydrogen bond
"""
Solution:
1. make 12 connectivity of atoms in conformer
2. examine heavy atom pairs between conformers including backbone to identify potential hydrogen bonds.
3. for each potential hydrogen bond, check if the donor has alternative H placement (sp3_1known, sp3_0known).
4. if the heavy atom has alternative H placement, place H to make D-H...A angle to be 180 degree. Mark this conformer as a hbond conformer.
"""

import logging
from mcce4.geom import *
from mcce4.pdbio import *
from ._place_h import *

DNEAR = 2.5     # Min Distance cutoff for hydrogen bond between heavy atoms
DFAR  = 3.6     # Max Distance cutoff for hydrogen bond between heavy atoms
MIN_NCRG = -0.2 # Minimum heavy atom charge for hydrogen bond donor/acceptor
MIN_HCRG = 0.2  # Minimum H atom charge for hydrogen bond donor/acceptor
EXCLUDED_ELEMENTS = {" C", " P", "CL", "BR", " I"}    # Exclude these elements from hydrogen bond
BLOCKING_ANGLE = 90  # Angle in degree to consider a blocking atom


def blocking_test(r2, acceptor):
    """Check if r2 -- acceptor - ? is blocked
    """
    # r2 (H) -- acceptor
    #               |
    #               ? (blocking atom)

    blocking = False
    for q in acceptor.conn12:
        Vaq = (q.xyz[0] - acceptor.xyz[0], q.xyz[1] - acceptor.xyz[1], q.xyz[2] - acceptor.xyz[2])
        Vah = (r2[0] - acceptor.xyz[0], r2[1] - acceptor.xyz[1], r2[2] - acceptor.xyz[2])
        angle = avv(Vaq, Vah) * 180 / np.pi
        if angle < BLOCKING_ANGLE:
            blocking = True
            break
    
    return blocking

def hbond_h(self):
    """Optimize H poisition for hydrogen bond
    """
    # initialize atom charges
    self.initialize_atom_qr()
    self.initialize_atom_id()

    # serialize heavy atom conformers
    for res in self.protein.residue:
        counter = {}
        for conf in res.conf:
            id = conf.history[:3]  # conf type + heavy atom making type
            if id in counter:
                counter[id] += 1
            else:
                counter[id] = 0  # count from 0

            conf.history = conf.history[:3] + "%03d"%(counter[conf.history[:3]]%1000) + conf.history[6:]


    # make 12 connectivity of atoms in conformer
    for res in self.protein.residue:
        for conf in res.conf:
            atoms_inconf = conf.atom
            
            # set two new properties, conn12 and parent_res
            for a in atoms_inconf:
                a.conn12 = []
            conf.parent_res = res

            # decide 12 connectivity
            for i in range(len(atoms_inconf)-1):
                for j in range(i+1, len(atoms_inconf)):
                    atom1 = atoms_inconf[i]
                    atom2 = atoms_inconf[j]
                    if dvv(atom1.xyz, atom2.xyz) < 2.0:
                        # check if these two atoms are connected in tpl file
                        if atom2.name in self.tpl.db[("CONNECT", atom1.name, conf.confType)].connected:
                            atom1.conn12.append(atom2)
                            atom2.conn12.append(atom1)

    

    # collect a list of donor/acceptor pairs based on charge and distance. This list is atom based.
    all_conformers = []
    for res in self.protein.residue:
        for conf in res.conf:
            all_conformers.append(conf)

    all_donors_acceptors = []
    for conf in all_conformers:
        for atom in conf.atom:
            # These are conditions for an atom to be a donor or acceptor (heavy atom)
            # 1. atom is not H
            # 2. atom charge is less than -0.2
            if not is_H(atom.name) and atom.charge < -0.2 and atom.element not in EXCLUDED_ELEMENTS:
                atom.parent_conf = conf
                all_donors_acceptors.append(atom)
    
    # Pick two atoms as potential hydrogen bond donors and acceptors
    cos_theta = math.cos(math.radians(bond_angle_sp3))
    sin_theta = math.sin(math.radians(bond_angle_sp3))
    for donor in all_donors_acceptors:
        for acceptor in all_donors_acceptors:
            if donor.parent_conf.parent_res != acceptor.parent_conf.parent_res:
                distance = dvv(donor.xyz, acceptor.xyz)  # within 2.0-4.0A h bond distance
                if DNEAR < distance < DFAR:
                    # examine donor H flexibility, only sp3_1known and sp3_0known are flexible
                    connected_heavy_atoms = [a for a in donor.conn12 if not is_H(a.name)]
                    orbital_type = self.tpl.db["CONNECT", donor.name, donor.parent_conf.confType].orbital
                    if len(connected_heavy_atoms) <= 1 and orbital_type.strip() == "sp3":  # qualified donor
                        should_have_h = [atom_name for atom_name in self.tpl.db[("CONNECT", donor.name, donor.parent_conf.confType)].connected if is_H(atom_name)]
                        if len(should_have_h) > 0:  # donor has H
                            # print(donor.atomID, acceptor.atomID)
                            new_conf = donor.parent_conf.clone()
                            new_conf.parent_res = donor.parent_conf.parent_res
                            # match atom by name
                            cloned_atom_by_name = {}
                            for a in new_conf.atom:
                                for b in donor.parent_conf.atom:
                                    if a.name == b.name:
                                        cloned_atom_by_name[a.name] = a

                            if len(connected_heavy_atoms) == 1:
                                # sp3_1known
                                # r3   r4
                                #  \  /
                                #   r0 - r2 (H) -- ra (Acceptor)
                                #    |
                                #    r1 (known)
                                r0 = donor.xyz
                                r1 = connected_heavy_atoms[0].xyz
                                ra = acceptor.xyz
                                # since r1 is known, r2 is constrained by the bond angle.
                                r01 = vector_vminusv(r1, r0)
                                r0a = vector_vminusv(ra, r0)
                                r10a = vector_vxv(r01, r0a)
                                ry = vector_vxv(r10a, r01)
                                uy = vector_normalize(ry)
                                ux = vector_normalize(r01)
                                x_component = vector_scale(ux, cos_theta)
                                y_component = vector_scale(uy, sin_theta)
                                r2 = vector_vplusv(r0, vector_scale(vector_sum3v(x_component, y_component, (0, 0, 0)), bond_length))

                                # check bloking r2 -- acceptor - ?, if this angle is less than blocking angle, return True
                                blocked = blocking_test(r2, acceptor)
                                if not blocked:
                                    r3, r4 = sp3_2known(r0, r1, r2)
                                    available_xyz = [r2, r3, r4]
                                    # print(available_xyz, [a.atomID for a in donor.conn12])
                                    for a in donor.conn12:
                                        if is_H(a.name):
                                            cloned_atom_by_name[a.name].xyz = available_xyz.pop(0)
                                    new_conf.history = donor.parent_conf.history[:6] + "H" + donor.parent_conf.history[7:]
                                    new_conf.parent_res.conf.append(new_conf)

                            elif len(connected_heavy_atoms) == 0:
                                # sp3_0known
                                # r3   r4
                                #  \  /
                                #   r0 - r1 (H) -- ra (Acceptor)
                                #    |
                                #    r2
                                r0 = donor.xyz
                                ra = acceptor.xyz
                                r1 = vector_vplusv(r0, vector_scale(vector_normalize(vector_vminusv(ra, r0)), bond_length))
                                blocked = blocking_test(r1, acceptor)
                                if not blocked:
                                    r2, r3, r4 = sp3_1known(r0, r1, None)
                                    available_xyz = [r1, r2, r3, r4]
                                    for a in donor.conn12:
                                        if is_H(a.name):
                                            cloned_atom_by_name[a.name].xyz = available_xyz.pop(0)
                                    new_conf.history = donor.parent_conf.history[:6] + "H" + donor.parent_conf.history[7:]
                                    new_conf.parent_res.conf.append(new_conf)

    # conformers sorted by the nature order in ftpl file
    for res in self.protein.residue:
        conf_types = self.tpl.db["CONFLIST", res.resID[0]]
        res.conf = sorted(res.conf, key=lambda x: conf_types.index(x.confType))  # sort by confType
        # serialize H conformers, history[7:10]
        counter = {}
        for conf in res.conf:
            id = conf.history[:7]
            if id in counter:
                counter[id] += 1
            else:
                counter[id] = 0
            conf.history = conf.history[:7] + "%03d"%(counter[id]%1000)
