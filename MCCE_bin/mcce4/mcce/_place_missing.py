#!/usr/bin/env python

import glob
import logging
from mcce4.geom import *
from mcce4.pdbio import *
from mcce4.mcce import _ideal_structure

class AtomTemplate:
    """
    This class defines non-H atom position template for a residue. Missing atoms
    will be completed by looking up the relative positions in this template
    """
    def __init__(self) -> None:
        self.atoms = {}

    def add_atom(self, atom):
        if atom.element != " H":
            self.atoms[atom.name] = atom.xyz


def load_atom_templates(folder):
    """
    Load ideal structures from a folder to a dictionary.
    Residue name is the key to a templates (dict)
    Atom name is the key to template indexed by residue name
    Value is (x, y, z) coordinates

    return a dictionary of dictionary
    """
    templates = {}
    files = glob.glob("%s/*.pdb" % folder)
    for fn in files:
        lines = open(fn).readlines()
        for line in lines:
            if line[:6] == "HETATM" or line[:6] == "ATOM  ":
                atom = Atom()
                atom.load_pdbline(line)
                resname = atom.resName
                if resname in templates:
                    if atom.name in templates[resname]:
                        logging.warning("Duplicate atom %s %s found when loading atom templates." % (resname. atom.name))
                    templates[resname][atom.name] = atom.xyz
                else:
                    templates[resname] = {atom.name: atom.xyz}

    return templates

   
def place_missing_heavy(self):
    """
    Place missing heavy atoms
    """
    self.reset_connect()
    self.make_connect12()

    n_completed = 0
    # load templates, templates hold a dictionary resname as key and a dictionary of atoms as value
    templates = load_atom_templates(_ideal_structure)

    #print(templates)
    # check connected atoms on each atom
    for res in self.protein.residue:
        if len(res.conf) > 1:
            for conf in res.conf[1:]:  # will not check missing heavy atoms on backbone atoms
                completed_missing = []
                for atom in conf.atom:
                    should_be_connected = self.tpl.db[("CONNECT", atom.name, conf.confType)].connected
                    actually_connected = [a.name for a in atom.connect12]
                    # filter out H atoms
                    should_be_connected_heavy = set([a for a in should_be_connected if not is_H(a)])
                    actually_connected_heavy = set([a for a in actually_connected if not is_H(a)])
                    completed_missing_heavy = set([a.name for a in completed_missing if not is_H(a.name)])
                    missing_heavy = should_be_connected_heavy - actually_connected_heavy - {" ?  "} - completed_missing_heavy  
                    # ignore outside connected atom as well
                    if missing_heavy:
                        # print("%s: %s - %s = %s" % (atom.name, str(should_be_connected_heavy), str(actually_connected_heavy), str(missing_heavy)))                    
                        
                        # add missing atom to conf.atom, will be checked in the same loop
                        connected_heavy = [a for a in atom.connect12 if not is_H(a.name)]
                        n_connected = len(connected_heavy)
                        if n_connected >= 2:  # 2 known connected atoms, match
                            resName = atom.resName
                            c_atom_p = atom.xyz                         # center atom location
                            c_atom_t = templates[resName][atom.name]    # center atom in template
                            
                            check_atom = connected_heavy[0]
                            if check_atom.name not in templates[resName]:
                                logging.error("Connected atom \"%s\" not found in template %s/%s_ideal.pdb" % (check_atom.name, _ideal_structure, resName))
                                logging.error("Exiting")
                                sys.exit()
                            m1_atom_p = check_atom.xyz
                            m1_atom_t = templates[resName][check_atom.name]

                            check_atom = connected_heavy[1]
                            if check_atom.name not in templates[resName]:
                                logging.error("Connected atom \"%s\" not found in template %s/%s_ideal.pdb" % (check_atom.name, _ideal_structure, resName))
                                logging.error("Exiting")
                                sys.exit()
                            m2_atom_p = check_atom.xyz
                            m2_atom_t = templates[resName][check_atom.name]
                            # print(c_atom_p, m1_atom_p, m2_atom_p)
                            # print(c_atom_t, m1_atom_t, m2_atom_t)
                            op = geom_3v_onto_3v(c_atom_t, m1_atom_t, m2_atom_t, c_atom_p, m1_atom_p, m2_atom_p)  # align template to known atoms
                            for atom_name in list(missing_heavy):
                                t_atom_xyz = templates[resName][atom_name]
                                # create a new atom                                
                                new_atom = Atom()                                
                                new_atom.inherit(atom)
                                new_atom.name = atom_name
                                new_atom.xyz = op.apply(t_atom_xyz)
                                # decide to add this atom
                                conf_atom_names = [x.name for x in conf.atom ]
                                if new_atom.name not in conf_atom_names:
                                    completed_missing.append(new_atom)

                        elif n_connected == 1:  # 1 known connected atom, look for the secondary connected atom
                            resName = atom.resName
                            c_atom_p = atom.xyz                         # center atom location
                            c_atom_t = templates[resName][atom.name]    # center atom in template
                            
                            check_atom = connected_heavy[0]
                            if check_atom.name not in templates[resName]:
                                logging.error("Connected atom \"%s\" not found in template %s/%s_ideal.pdb" % (check_atom.name, _ideal_structure, resName))
                                logging.error("Exiting")
                                sys.exit()
                            m1_atom_p = check_atom.xyz
                            m1_atom_t = templates[resName][check_atom.name]

                            # look for a secondary connected atom                            
                            found_secondary = False
                            for a in [x for x in check_atom.connect12 if not is_H(x.name)]:
                                if a != atom:
                                    found_secondary = True
                                    m2_atom_p = a.xyz
                                    m2_atom_t = templates[resName][a.name]
                                    op = geom_3v_onto_3v(c_atom_t, m1_atom_t, m2_atom_t, c_atom_p, m1_atom_p, m2_atom_p)  # align template to known atoms
                                    for atom_name in list(missing_heavy):
                                        t_atom_xyz = templates[resName][atom_name]
                                        # create a new atom                                
                                        new_atom = Atom()                                
                                        new_atom.inherit(atom)
                                        new_atom.name = atom_name
                                        new_atom.xyz = op.apply(t_atom_xyz)
                                        # decide to add this atom
                                        conf_atom_names = [x.name for x in conf.atom ]
                                        if new_atom.name not in conf_atom_names:
                                            completed_missing.append(new_atom)

                                    break # only one secondary atom is required
                            
                            if not found_secondary:
                                logging.error("Could not find enough known atoms to place Missing - %s - %s - ?" % (atom.name, check_atom.name))
                                sys.exit()

                        elif n_connected == 0: # 0 known connected atom, point match, arbitray placement
                            print("0 known, to be done")
                if completed_missing:
                    n_completed += len(completed_missing)
                    conf.atom += completed_missing
                    for added_atom in completed_missing:
                        logging.info("Added \"%s\" for %s" % (added_atom.name, added_atom.confID))

    return n_completed