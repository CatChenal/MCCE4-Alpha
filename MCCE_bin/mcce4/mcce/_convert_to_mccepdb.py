from copy import deepcopy
from mcce4.pdbio import *


def initialize_atom_qr(self):
    "init charge and radius"
    for res in self.protein.residue:
        for conf in res.conf:
            for atom in conf.atom:
                key = ("CHARGE", conf.confType, atom.name)
                if key in self.tpl.db:
                    atom.charge = self.tpl.db[key]
                else:
                    logging.debug("Parameter key %s doesn't exist, defaul to value 0." % str(key))
                    atom.charge = 0.0
                key = ("RADIUS", conf.confType, atom.name)
                if key in self.tpl.db:
                    radius_param = self.tpl.db[key]
                    atom.r_bound = radius_param.r_bound
                    atom.r_vdw = radius_param.r_vdw
                    atom.e_vdw = radius_param.e_vdw
                else:
                    logging.debug("Parameter key %s doesn't exist, default to value 0." % str(key))
                    atom.r_bound = 0.0
                    atom.r_vdw = 0.0
                    atom.e_vdw = 0.0


def initialize_atom_id(self):
    """init serial, atomID, confNum, and confID"""
    serial = 0  # serial starts from 0
    for res in self.protein.residue:
        confNum = 0
        for conf in res.conf:
            confID = ""
            for atom in conf.atom:
                atom.element = atom.name[:2]  # element name defaults to the first two char
                if len(atom.name.strip()) == 4 and atom.name[0] == "H": # special case H, 4 char, H is the first char
                    atom.element = " H"
                atom.serial = serial
                atom.confNum = confNum
                confID = atom.confID = "%5s%c%04d%c%03d" % (atom.confType, atom.chainID, atom.resSeq, atom.iCode, atom.confNum)
                atom.atomID = "%4s%3s%04d%c%03d" % (atom.name, atom.resName, atom.resSeq, \
                                                    atom.chainID, atom.confNum)
                serial += 1
            confNum += 1
            if confID:  # only update when atoms exist, 0 atom conformer keeps its ID
                conf.confID = confID

def convert_to_mccepdb(self):
    """Convert from generic pdb structure to internal mccepdb format
    self.prm holds run.prm
    self.tpl holds topology records
    self.structure.lines holds the first structure as atom lines
    """
    self.protein = Protein()  # reset a protein object to receive the converted structure

    if self.structure:
        lines = [x.strip("\n") for x in self.structure.lines if x[:6] == "ATOM  " or x[:6] == "HETATM"]
        for atom_line in lines:
            # atom info from PDB format
            atom = Atom()
            atom.name = atom_line[12:16]
            atom.altLoc = atom_line[16]
            if atom.altLoc == " ":
                atom.altLoc = "_"
            atom.resName = atom_line[17:20]
            atom.chainID = atom_line[21]
            atom.resSeq = int(atom_line[22:26])
            atom.iCode = atom_line[26]
            if atom.iCode == " ":
                atom.iCode = "_"
            atom.xyz = (float(atom_line[30:38]), float(atom_line[38:46]), float(atom_line[46:54]))
            if len(atom.name.strip()) == 4 and atom.name[0] == "H":
                atom.element = " H"
            else:
                atom.element = atom.name[:2]

            if atom.element == " H":  # skip all hydrogen atoms
                continue

            #print("HERE", atom_line, atom.altLoc)
            # Fit the atom to the first available conformer
            atom.resID = (atom.resName, atom.chainID, atom.resSeq, atom.iCode)
            ## mcce internals:

            # find which conformer this heavy atom belongs to
            key = ("CONFLIST", atom.resName)
            possible_conftypes = self.tpl.db[key]

            found_conftype = ""
            for conftype in possible_conftypes:
                key = ("CONNECT", atom.name, conftype)
                if key in self.tpl.db:
                    found_conftype = conftype
                    break

            if found_conftype:
                atom.confType = found_conftype
            elif atom.name == " OXT" and not found_conftype:  # special case, add to CTR
                logging.info("Atom %s in residue %s is assigned to CTR" % (atom.name, atom.resID))
                atom.confType = atom.resName + "BK"   # keep in backbone, make_termini will only check backbone
            else:
                logging.warning("Atom %s of residue %s not defined in ftpl file, drop." % (atom.name, atom.resID))
                continue

            if atom.altLoc == "_":
                label = "O"
            else:
                label = atom.altLoc
            atom.history = "%s%c000%s" % (atom.confType[3:5], label, atom.history[6:])

            # assign this atom to the structure: self -> protein -> residue -> conf -> atom
            found_matched_res = False
            found_matched_conf = False
            for res in self.protein.residue:
                if res.resID == atom.resID:
                    found_matched_res = True
                    for conf in res.conf:
                        #print("HERE %c - %c" % (conf.altLoc, atom.altLoc))
                        if conf.confType == atom.confType and conf.altLoc == atom.altLoc:
                            found_matched_conf = True
                            conf.atom.append(atom)
                    if found_matched_res and not found_matched_conf:
                        new_conf = Conformer()
                        new_conf.init_by_atom(atom)
                        res.conf.append(new_conf)
            if not found_matched_res:
                new_conf = Conformer()
                new_conf.init_by_atom(atom)
                new_res = Residue()
                new_res.resID = new_conf.resID
                new_res.conf.append(new_conf)
                self.protein.residue.append(new_res)


        # clean up conformers, make sure BK conformer is added, and combine conformers to one conformer type (usually "01")
        ## Move BK to conf slot 0 and make sure only one BK
        for res in self.protein.residue:
            if len(res.conf) > 0:
                if res.conf[0].confType[3:5] != "BK":  # the first conf slot is not BK
                    # if BK is in other slots, switch
                    iconf_found = 0
                    for iconf in range(1, len(res.conf)):
                        if res.conf[iconf].confType[3:5] == "BK":
                            iconf_found = iconf
                            break
                    if iconf_found > 0:
                        res.conf[0], res.conf[iconf_found] = res.conf[iconf_found], res.conf[0]
                    # if not, insert a blank BK conformer
                    else:
                        new_conf = Conformer()
                        new_conf.resID = res.resID
                        new_conf.altLoc = res.conf[0].altLoc  # inherit altLoc from the first conformer
                        new_conf.confType =  res.conf[0].confType[:3] + "BK"
                        new_conf.confID =  "%5s%c%04d%c%03d" % (new_conf.confType,
                                                                res.conf[0].atom[0].chainID, 
                                                                res.conf[0].atom[0].resSeq, 
                                                                res.conf[0].atom[0].iCode, 
                                                                0)
                        new_conf.history =  "BK" + "_"*8
                        res.conf.insert(0, new_conf)
                else:
                    # if more than 1 BK conformer exists, delete other BK and conformers share the same altLoc
                    iconf_found = []
                    for iconf in range(1, len(res.conf)):
                        if res.conf[iconf].confType[3:5] == "BK":
                            iconf_found.append(iconf)
                    if iconf_found:  # found BK duplicates
                        logging.warning("AltLoc detected on backbone atoms, ignored." )
                        for ic in reversed(iconf_found):
                            logging.warning("   Drop conformer %s - %s with AltLoc %c" % (res.conf[ic].confType, res.conf[ic].history, res.conf[ic].altLoc))
                            res.conf.pop(ic)

        ## scan side chain conformers, one altLoc should only have one conformer, share common atoms
        for res in self.protein.residue:
            altLoc = []
            if len(res.conf) > 1:
                for conf in res.conf[1:]:
                    altLoc.append(conf.altLoc)
                #print(res.resID, altLoc)
                for aL in altLoc:  # scan aL in side chain conformers, warn if multiple side chain conf found
                    counter = 0
                    for conf in res.conf[1:]:
                        if conf.altLoc == aL:
                            counter += 1
                    if counter > 1:
                        logging.warning("Atoms sharing the same alternative location %s - %c were not in one conformer" % (res.resID, aL))
                        raise Exception("Heavy atoms not in the same conformer.")
                if len(altLoc) > 1:  # multiple altLoc 
                    i_aL = -1
                    if "_" in altLoc:
                        i_aL = altLoc.index("_")
                    if i_aL >= 0:  # common atoms found with altLoc "_"
                        common_atoms = res.conf[i_aL+1].atom
                    else:
                        common_atoms = []
                    for conf in res.conf[1:]:
                        if conf.altLoc != "_":
                            conf.atom = deepcopy(common_atoms) + conf.atom
                    if i_aL >= 0:
                        res.conf.pop(i_aL+1)  # +1 because conf counts from BK conformer, while i_aL starts from 0 for sidechains
            
        # assign confNum, confID to atom and conf
        logging.info("Initializing conformer and atom identifiers in mccepdb format")
        initialize_atom_id(self)

        # initialize atom properties: charge and radius
        logging.info("Assigning charge and radiues to atoms")
        initialize_atom_qr(self)
