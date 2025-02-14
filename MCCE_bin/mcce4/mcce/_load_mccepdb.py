from ..pdbio import *

def load_mccepdb(self, fname):

    rawlines = open(fname).readlines()
    lines = [x.strip("\n") for x in rawlines if x[:6] == "ATOM  " or x[:6] == "HETATM"]
    for line in lines:
        atom = Atom()
        atom.loadline(line, self.tpl)


        # reverse search if this atom belongs to an existing conformer and residue
        found_conf = False
        found_res = False
        for i_res in range(len(self.protein.residue) - 1, -1, -1):
            insert_res = i_res
            for i_conf in range(len(self.protein.residue[i_res].conf) - 1, -1, -1):
                if self.protein.residue[i_res].conf[i_conf].confID == atom.confID:
                    found_conf = True
                    self.protein.residue[i_res].conf[i_conf].atom.append(atom)
                    # altLoc atom overwrites the history
                    if atom.history[2:6] != "O000":
                        self.protein.residue[i_res].conf[i_conf].history = atom.history
                        # print(atom.history[2:6])
                        # print(self.protein.residue[i_res].conf[i_conf].confID, self.protein.residue[i_res].conf[i_conf].history)
                    break
            if found_conf:
                break
            elif self.protein.residue[i_res].resID == atom.resID:  # not in conf but residue ID matches
                conf = Conformer()
                conf.confType = atom.confType
                conf.confID = atom.confID
                conf.resID = atom.resID
                conf.resSeq = atom.resSeq
                conf.iCode = atom.iCode
                conf.history = atom.history
                conf.atom.append(atom)
                self.protein.residue[i_res].conf.append(conf)
                found_res = True

        if not found_conf and not found_res:  # new residue
            conf = Conformer()
            conf.confType = atom.confType
            conf.confID = atom.confID
            conf.resID = atom.resID
            conf.resSeq = atom.resSeq
            conf.iCode = atom.iCode
            conf.history = atom.history
            conf.atom.append(atom)
            res = Residue()
            res.resID = conf.resID
            res.conf.append(conf)
            self.protein.residue.append(res)

    # Insert an empty conformer for cofactors that do not have backbone
    for res in self.protein.residue:
        first_confID = res.conf[0].confID
        if first_confID[3:5] != "BK":
            conf = Conformer()
            conf.confID = "%sBK%s000" % (first_confID[:3], first_confID[5:11])
            conf.confType = first_confID[:3] + "BK"
            conf.history = "BK" + "_"*8

            conf.resID = res.resID
            res.conf.insert(0, conf)

    # Assign an index number to each conformer
    n_conf = 0
    for res in self.protein.residue:
        if len(res.conf) <= 1:   # backbone only
            continue
        for conf in res.conf[1:]:
            conf.i = n_conf
            n_conf += 1

