import sys
import numpy as np
from ..pdbio import *
from ..geom import *
from ._vdw import vdw_conf

Hbond_atom_elements = {" O", " N", " F"}  # Hbond atoms need to have this name
Hbond_atom_charge = -0.2  # Hbond atoms need to have chare more negative than this
Hbond_distance_ini2 = 12.0 * 12.0  # Optimize if the initial distance is within this distance, squared
Hbond_distance_end2 = 3.0 * 3.0  # Make the final distance as close as possible to this distance, squared


class HBOND_ELEMENT:    # one Hbond screen candidate is a pair of atom and the conformer it belongs to. For now, the conformer is always conf[1]
    def __init__(self, atom, conf, res):
        self.atom = atom         
        self.conf = conf
        self.res = res

def look_for_hbond(self, c1, c2):
    """ Based on candidates c1 and c2, find optimized conformers.
    """
    base_pair = (c1, c2)  # optimized pair
    d2 = ddvv(c1.atom.xyz, c2.atom.xyz)
    base_doff = abs(d2 - Hbond_distance_end2)

    if c1.conf != c2.conf and d2 < Hbond_distance_ini2: # atoms in the same conf or too far will be ignored
        #print("optimizing ", c1.conf.confID, c2.conf.confID, end=" ")
        improved = True
        phi = 60/180*np.pi
        while improved:
            improved = False
            new1_confs = self.swing(base_pair[0].conf, phi)
            new2_confs = self.swing(base_pair[1].conf, phi)
            for conf1 in new1_confs+[base_pair[0].conf]:
                atom1 = None
                for atom in conf1.atom:
                    if c1.atom.name == atom.name:
                        atom1 = atom
                        break
                for conf2 in new2_confs+[base_pair[1].conf]:
                    atom2 = None
                    for atom in conf2.atom:
                        if c2.atom.name == atom.name:
                            atom2 = atom
                            break
                    if atom1 and atom2:
                        d2 = ddvv(atom1.xyz, atom2.xyz)
                        if abs(d2 - Hbond_distance_end2) < base_doff:
                            base_doff = abs(d2 - Hbond_distance_end2)
                            base_pair = (HBOND_ELEMENT(atom1, conf1, c1.res), HBOND_ELEMENT(atom2, conf2, c2.res))
                            improved = True

        improved = True
        phi = 15/180*np.pi
        while improved:
            improved = False
            new1_confs = self.swing(base_pair[0].conf, phi)
            new2_confs = self.swing(base_pair[1].conf, phi)
            for conf1 in new1_confs+[base_pair[0].conf]:
                atom1 = None
                for atom in conf1.atom:
                    if c1.atom.name == atom.name:
                        atom1 = atom
                        break
                for conf2 in new2_confs+[base_pair[1].conf]:
                    atom2 = None
                    for atom in conf2.atom:
                        if c2.atom.name == atom.name:
                            atom2 = atom
                            break
                    if atom1 and atom2:
                        d2 = ddvv(atom1.xyz, atom2.xyz)
                        if abs(d2 - Hbond_distance_end2) < base_doff:
                            base_doff = abs(d2 - Hbond_distance_end2)
                            base_pair = (HBOND_ELEMENT(atom1, conf1, c1.res), HBOND_ELEMENT(atom2, conf2, c2.res))
                            improved = True

        improved = True
        phi = 3/180*np.pi
        while improved:
            improved = False
            new1_confs = self.swing(base_pair[0].conf, phi)
            new2_confs = self.swing(base_pair[1].conf, phi)
            for conf1 in new1_confs+[base_pair[0].conf]:
                atom1 = None
                for atom in conf1.atom:
                    if c1.atom.name == atom.name:
                        atom1 = atom
                        break
                for conf2 in new2_confs+[base_pair[1].conf]:
                    atom2 = None
                    for atom in conf2.atom:
                        if c2.atom.name == atom.name:
                            atom2 = atom
                            break
                    if atom1 and atom2:
                        d2 = ddvv(atom1.xyz, atom2.xyz)
                        if abs(d2 - Hbond_distance_end2) < base_doff:
                            base_doff = abs(d2 - Hbond_distance_end2)
                            base_pair = (HBOND_ELEMENT(atom1, conf1, c1.res), HBOND_ELEMENT(atom2, conf2, c2.res))
                            improved = True

        improved = True
        phi = 1/180*np.pi
        while improved:
            improved = False
            new1_confs = self.swing(base_pair[0].conf, phi)
            new2_confs = self.swing(base_pair[1].conf, phi)
            for conf1 in new1_confs+[base_pair[0].conf]:
                atom1 = None
                for atom in conf1.atom:
                    if c1.atom.name == atom.name:
                        atom1 = atom
                        break
                for conf2 in new2_confs+[base_pair[1].conf]:
                    atom2 = None
                    for atom in conf2.atom:
                        if c2.atom.name == atom.name:
                            atom2 = atom
                            break
                    if atom1 and atom2:
                        d2 = ddvv(atom1.xyz, atom2.xyz)
                        if abs(d2 - Hbond_distance_end2) < base_doff:
                            base_doff = abs(d2 - Hbond_distance_end2)
                            base_pair = (HBOND_ELEMENT(atom1, conf1, c1.res), HBOND_ELEMENT(atom2, conf2, c2.res))
                            improved = True


    if base_doff > 3.25 or base_pair == (c1, c2):  # do not return anything if the atoms are beyond H bond distance (3.5*3.5 - 3*3)
        return None
    else:
        return base_pair


def rot_hdirected(self):
    """ Make hydrogen bond directed rotamers, using swing to find optimum O, N distances
    Basic Algorithm:
    - Define atom pairs: Hbond donor and acceptors are defined a negatively charged O, N and F atoms
    - Optimum distance: The Hbond atoms will be optimized to be as close as possible to Hbond_distance_end
    - procedure:
        - Search atom pair candidates and make a list of (atom <- conformer[1])
        - Use swing() to optimize the atom pairs
    """
    # Identify Hbond donors and acceptors
    Hbond_candidates = []
    for res in self.protein.residue:
        if len(res.conf) > 1:
            conf = res.conf[1]
            for atom in conf.atom:
                if atom.element in Hbond_atom_elements and atom.charge < Hbond_atom_charge:
                    candidate = HBOND_ELEMENT(atom, conf, res)
                    Hbond_candidates.append(candidate)

    logging.info("   Found %d potential hydrogen bond donors and acceptors" % len(Hbond_candidates))
    # go over all pairs to optimize each pair
    for ie1 in range(len(Hbond_candidates)-1):
        candidate1 = Hbond_candidates[ie1]
        for ie2 in range(ie1+1, len(Hbond_candidates)):
            candidate2 = Hbond_candidates[ie2]
            optimized_confs = look_for_hbond(self, candidate1, candidate2)
            if optimized_confs:
                new_candidate1 = optimized_confs[0]
                if new_candidate1.atom != candidate1.atom:
                    new_candidate1.conf.history = new_candidate1.conf.history[:2] + "H" + new_candidate1.conf.history[3:]
                    new_candidate1.res.conf.append(new_candidate1.conf)
                new_candidate2 = optimized_confs[1]
                if new_candidate2.atom != candidate2.atom:
                    new_candidate2.conf.history = new_candidate2.conf.history[:2] + "H" + new_candidate2.conf.history[3:]
                    new_candidate2.res.conf.append(new_candidate2.conf)

                
