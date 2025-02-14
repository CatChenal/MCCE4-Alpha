from copy import deepcopy

def rot_swap(self):
    """ Make swap rotamers based on ROT_SWAP rules in ftpl files
    """
    self.reset_connect()
    for res in self.protein.residue:
        resName = res.resID[0]
        key = ("ROT_SWAP", resName)
        if key in self.tpl.db:
            for swap_pair in self.tpl.db[key]:
                atom1_name, atom2_name = swap_pair
                nconf = len(res.conf)
                for i_conf in range(1, nconf):  # ignore the backbone conformer
                    conf = res.conf[i_conf]
                    # find the first atom
                    atom1 = None
                    atom2 = None
                    for atom in conf.atom:
                        if atom.name == atom1_name:
                            atom1 = atom
                        elif atom.name == atom2_name:
                            atom2 = atom
                    if atom1 and atom2:  # found both, create a new conformer
                        new_conf = deepcopy(conf)
                        new_conf.history = conf.history[:2] + "W" + conf.history[3:]
                        for new_atom in new_conf.atom:
                            if new_atom.name == atom1_name:  # swap to atom2
                                new_atom.xyz = atom2.xyz
                                new_atom.r_bound = atom2.r_bound
                                new_atom.charge = atom2.charge
                            elif new_atom.name == atom2_name:  # swap to atom1
                                new_atom.xyz = atom1.xyz
                                new_atom.r_bound = atom1.r_bound
                                new_atom.charge = atom1.charge
                            new_atom.history = new_conf.history

                        # insert this new conf
                        res.conf.append(new_conf)
                        
