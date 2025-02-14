import sys
import numpy as np
from ..pdbio import is_H
from ..geom import *



def rot_swing(self):
    """ Make rotamers based on SWING rules in ftpl files
    """

    phi = float(self.prm.PHI_SWING.value)/180*np.pi
    #self.make_connect12()  # make connect12 every time new conformers are created

    for res in self.protein.residue:
        if len(res.conf) > 1:
            new_confs = []
            for conf in res.conf[1:]:
                new_confs += self.swing(conf, phi)
            res.conf += new_confs


                

def swing(self, conf, phi):
    """ Make a swing of phi in radiant on one conformer.
    This function returns a list of new conformers with correct connectivity.
    """

    # get rotatable bonds from ftpl file
    resName = conf.resID[0]
    key = ("ROTATE", resName)
    new_confs = []
    if key in self.tpl.db:  # handle one rotate rule at a time
        rotate_axis_byatom_names = self.tpl.db[key]
        for each_axis in rotate_axis_byatom_names:
            confs_base = [conf] + new_confs  # These are the confs to be operated on for the rotation axis
            for each_conf in confs_base:
                # find the second atom of the axis, the second atom is always in the side chain
                found = False
                for atom in each_conf.atom:
                    if atom.name == each_axis[1]:
                        found = True
                        atom2 = atom
                        break
                if not found:
                    print("Atom %s in ROTATE record %s %s was not found in this residue %s" % (each_axis[1], key, each_axis, each_conf.resID))
                    sys.exit()

                # find the first atom of the axis, the first atom must be in connect12 of the second atom
                found = False
                for atom in atom2.connect12:
                    if atom.name == each_axis[0]:
                        found = True
                        atom1 = atom
                        break
                if not found:
                    print([n.name for n in atom.connect12])
                    print("Atom %s in ROTATE record %s %s was not found in this residue %s" % (each_axis[0], key, each_axis, each_conf.resID))
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

                # obtain rotate operation
                axis = LINE()
                axis.from2p(atom1.xyz, atom2.xyz)
                op = OPERATION()
                for i in [-1, 1]:
                    op.reset()
                    op.roll(phi*i, axis)
                    new_xyz = {}  # atom name : xyz dictionary for affected atoms
                    for atom in affected_atoms:
                        new_xyz[atom.name] = op.apply(atom.xyz)
                    new_conf = each_conf.clone()     # this clone will copy all atoms and connect12
                    new_conf.history = each_conf.history[:2] + "R" + conf.history[3:]
                    for atom in new_conf.atom:  # Update
                        if atom.name in new_xyz:
                            atom.xyz = new_xyz[atom.name]
                    new_confs.append(new_conf)

    return new_confs