import sys
import numpy as np
from ..pdbio import is_H
from ..geom import *
from ._vdw import vdw_conf


def rot_rotate(self):
    """ Make rotamers based on ROTATE rules in ftpl files
    """
    n_rotation_steps = int(self.prm.ROTATIONS.value)
    phi = 2*np.pi/n_rotation_steps

    for res in self.protein.residue:
        resName = res.resID[0]
        key = ("ROTATE", resName)
        if key in self.tpl.db:  # handle one rotate rule at a time
            # find the axis
            rotate_axis_byatom_names = self.tpl.db[key]
            for each_axis in rotate_axis_byatom_names:
                new_confs = []  # New conformers will be added by each rotation and multiplied by other rotation axes
                for conf in res.conf[1:]:  # only create rotamers for side chain conformers
                    # find the second atom of the axis, the second atom is always in the side chain
                    found = False
                    for atom in conf.atom:
                        if atom.name == each_axis[1]:
                            found = True
                            atom2 = atom
                            break
                    if not found:
                        print("Atom %s in ROTATE record %s %s was not found in this residue %s" % (each_axis[1], key, each_axis, res.resID))
                        sys.exit()

                    # find the first atom of the axis, the first atom must be in connect12 of the second atom
                    found = False
                    for atom in atom2.connect12:
                        if atom.name == each_axis[0]:
                            found = True
                            atom1 = atom
                            break
                    if not found:
                        print("Atom %s in ROTATE record %s %s was not found in this residue %s" % (each_axis[0], key, each_axis, res.resID))
                        sys.exit()

                    # find affected heavy atoms
                    affected_atoms = []  # initialize affected atom set, 
                    for atom in atom2.connect12:  # we start from atom2 and expand the affected atom list away from atom1
                        if atom != atom1 and not is_H(atom.name):
                            affected_atoms.append(atom)
                    for a_atom in affected_atoms:  # expand
                        for atom in a_atom.connect12:
                            if atom != atom1 and \
                                atom != atom2 and \
                                not is_H(atom.name) and \
                                not atom in affected_atoms and \
                                atom.resID == atom2.resID:  # avoid cross residue linkage just in case
                                affected_atoms.append(atom)
                    # print([a.name for a in affected_atoms])

                    # obtain rotate operation
                    axis = LINE()
                    axis.from2p(atom1.xyz, atom2.xyz)
                    op = OPERATION()
                    for i in range(1, n_rotation_steps):
                        op.reset()
                        op.roll(phi*i, axis)
                        new_xyz = {}  # atom name : xyz dictionary for affected atoms
                        for atom in affected_atoms:
                            new_xyz[atom.name] = op.apply(atom.xyz)
                        new_conf = conf.clone()     # this clone will copy all atoms except their 13 and 14 connectivity
                        new_conf.history = conf.history[:2] + "R" + conf.history[3:]
                        for atom in new_conf.atom:  # Update
                            if atom.name in new_xyz:
                                atom.xyz = new_xyz[atom.name]

                        new_confs.append(new_conf)
                
                res.conf += new_confs
