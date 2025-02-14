import logging

indention_1 = "   "

def atom_number_check(self, ignore_H = True):
    """ Check missing and duplicate atoms at the conformer level.
    Optionally ignore H atoms
    Log warnning with level 1 indent if missings and duplicates are found
    """


    # make a conformer - {atoms} database by looking up CONNECT records in self.tpl.db
    atoms_in_conf = {}
    for key in self.tpl.db.keys():
        if key[0] == "CONNECT":
            conf = key[2]
            atom = key[1].strip('"')
            if conf in atoms_in_conf:
                atoms_in_conf[conf].add(atom)
            else:
                atoms_in_conf[conf] = {atom}

    # reduce the database if H is ignored
    if ignore_H:
        for key, value in atoms_in_conf.items():
            h_atoms = set()
            for atom in value:
                element = atom[:2]
                if len(atom.strip()) == 4 and atom[0] == "H":
                    element = " H"
                if element == " H":
                    h_atoms.add(atom)
            new_value = value - h_atoms
            atoms_in_conf[key] = new_value

    # debug: check the content
    # for conf, atoms in atoms_in_conf.items():
    #     print(conf, atoms)
    for i_res in range(len(self.protein.residue)):
        res = self.protein.residue[i_res]
        for conf in res.conf:
            detected_atoms = set()
            #print(conf.confID, conf.confType)
            for atom in conf.atom:
                if atom.name not in detected_atoms:
                    detected_atoms.add(atom.name)
                else: # already in this conformer, duplicate
                    logging.warning("%s Duplicated atom \"%s\" in conformer %s" % (indention_1, atom.name, conf.confID))
            if conf.confType in atoms_in_conf:
                missing_atoms = atoms_in_conf[conf.confType] - detected_atoms
            else:  # conf[0] is backbone and may have 0 atom and not defined by CONNECT in ftpl file
                missing_atoms = set()
            if missing_atoms:
                for atom in list(missing_atoms):
                    # check missing atom exceptions
                    exception = False
                    if i_res > 0:  # next to NTR
                        if self.protein.residue[i_res-1].resID[0] == "NTR" or \
                        self.protein.residue[i_res-1].resID[0] == "NTG":
                            if atom in {" N  ", " CA "}:
                                exception = True
                    if i_res < (len(self.protein.residue)-1):  # next to CTR
                        if self.protein.residue[i_res+1].resID[0] == "CTR":
                            if atom in {" O  ", " C  "}:
                                exception = True
                    if res.resID[0] == "CTR":  # OXT case
                        if atom == " OXT":
                            exception = True
                    if not exception:
                        logging.warning("%s Missing atom \"%s\" in %s" % (indention_1, atom, conf.confID))