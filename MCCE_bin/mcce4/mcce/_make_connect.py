""" Make connectivity for MCCE object, assuming radius parameters have been assigned to atoms already
"""
import logging
import math
from mcce4.geom import *
from mcce4.pdbio import is_H
from ._identify_ligands import match_strs

_BONDDISTANCE_scaling = 0.54  # calibrated by 1akk

# Detect and rename ligands by predefined rules copied from param_PARSE/ligand_detect_rules.ftpl
# LIGAND_ID, CYS, CYS: " SG " - " SG "; 2.03 +- 0.90; CYD, CYD
# LIGAND_ID, CYS, HEC: " SG " - " CA*"; 1.90 +- 1.00; CYL, HEC
# LIGAND_ID, CYS, HEM: " SG " - " CA*"; 1.90 +- 1.00; CYL, HEM
# LIGAND_ID, HIS, HEM: " NE2" - "FE  "; 2.10 +- 0.70; HIL, HEM
# LIGAND_ID, HIS, HEA: " NE2" - "FE  "; 2.10 +- 0.70; HIL, HEA
# LIGAND_ID, HIS, HEB: " NE2" - "FE  "; 2.10 +- 0.70; HIL, HEB
# LIGAND_ID, HIS, HEC: " NE2" - "FE  "; 2.10 +- 0.70; HIL, HEC
# LIGAND_ID, MET, HEM: " SD " - "FE  "; 2.30 +- 0.50; MEL, HEM
# LIGAND_ID, MET, HEA: " SD " - "FE  "; 2.30 +- 0.50; MEL, HEA
# LIGAND_ID, MET, HEB: " SD " - "FE  "; 2.30 +- 0.50; MEL, HEB
# LIGAND_ID, MET, HEC: " SD " - "FE  "; 2.30 +- 0.50; MEL, HEC

ligand_rules = {("CYD", "CYD"): (" SG ", " SG ", 2.03, 0.90),
                ("CYL", "HEC"): (" SG ", " CA*", 1.90, 1.00),
                ("CYL", "HEM"): (" SG ", " CA*", 1.90, 1.00),
                ("HIL", "HEM"): (" NE2", "FE  ", 2.10, 0.70),
                ("HIL", "HEA"): (" NE2", "FE  ", 2.10, 0.70),
                ("HIL", "HEB"): (" NE2", "FE  ", 2.10, 0.70),
                ("HIL", "HEC"): (" NE2", "FE  ", 2.10, 0.70),
                ("MEL", "HEM"): (" SD ", "FE  ", 2.30, 0.50),
                ("MEL", "HEA"): (" SD ", "FE  ", 2.30, 0.50),
                ("MEL", "HEB"): (" SD ", "FE  ", 2.30, 0.50),
                ("MEL", "HEC"): (" SD ", "FE  ", 2.30, 0.50),
                ("PAA", "HEM"): (" CAA", " C2A", 1.60, 0.50),
                ("PDD", "HEM"): (" CAD", " C3D", 1.60, 0.50),
                ("PAA", "HEA"): (" CAA", " C2A", 1.60, 0.50),
                ("PDD", "HEA"): (" CAD", " C3D", 1.60, 0.50),
                ("PAA", "HEB"): (" CAA", " C2A", 1.60, 0.50),
                ("PDD", "HEB"): (" CAD", " C3D", 1.60, 0.50),
                ("PAA", "HEC"): (" CAA", " C2A", 1.60, 0.50),
                ("PDD", "HEC"): (" CAD", " C3D", 1.60, 0.50),
}



def make_connect12(self):
    """Make 1-2 connectivity for each atom in protein
    """
    self.reset_connect()

    for i_res in range(len(self.protein.residue)):
        res = self.protein.residue[i_res]
        for conf in res.conf:  # search all conformers including backbone
            for atom in conf.atom:
                key = ("CONNECT", atom.name, atom.confType)
                connected_atoms = self.tpl.db[key].connected
                for c_atom in connected_atoms:
                    found = False
                    if "?" in c_atom:  # ligated to an atom outside the residue
                        for res2 in self.protein.residue:
                            if res != res2:  # skip the same residue
                                for conf2 in res2.conf:
                                    for atom2 in conf2.atom:
                                        key2 = ("CONNECT", atom2.name, atom2.confType)
                                        connected_atoms2 = self.tpl.db[key2].connected
                                        ligated2 = False
                                        for ligated_atom2 in connected_atoms2:
                                            if "?" in ligated_atom2:
                                                ligated2 = True
                                                break
                                        if ligated2:
                                            r = (atom.r_vdw + atom2.r_vdw) * _BONDDISTANCE_scaling
                                            CUTOFF2 = r * r
                                            d2 = ddvv(atom.xyz, atom2.xyz)
                                            if d2 < CUTOFF2 and atom2 not in atom.connect12:
                                                atom.connect12.append(atom2)
                                                found = True
                                                #print(atom.atomID, atom2.atomID)
                                            # else:
                                            #     print("%s - %s: d=%.3f" % (atom.atomID, atom2.atomID, math.sqrt(d2)))
                                            else: # in addition to detecting a ligated atom by distance, need to use LIGAND parameters to capture more 
                                                res1_name = res.resID[0]
                                                res2_name = res2.resID[0]
                                                atom1_name = atom.name
                                                atom2_name = atom2.name
                                                key = (res1_name, res2_name)
                                                swapped_key = (res2_name, res1_name)
                                                found_key = False
                                                if key in ligand_rules:
                                                    atom1_name_inrule, atom2_name_inrule, distance, tolerance = ligand_rules[key]
                                                    found_key = True
                                                elif swapped_key in self.tpl.db:
                                                    atom2_name_inrule, atom1_name_inrule, distance, tolerance = ligand_rules[swapped_key]
                                                    found_key = True
                                                if found_key:
                                                    #print("HERE")
                                                    #print(atom1_name, atom2_name, atom1_name_inrule, atom2_name_inrule)
                                                    if match_strs(atom1_name, atom1_name_inrule) and match_strs(atom2_name, atom2_name_inrule):
                                                        #print("THERE")
                                                        d = math.sqrt(d2)
                                                        if distance - tolerance <= d <= distance + tolerance:
                                                            if atom2 not in atom.connect12:
                                                                atom.connect12.append(atom2)
                                                                found = True
                                                                #print("Checking LIGAND pair: %s.%s - %s.%s" % (res1_name, atom1_name, res2_name, atom2_name))
                                                            if atom not in atom2.connect12:
                                                                atom2.connect12.append(atom)
                                                                found = True

                                    if found:   # one "?" for one ligand
                                        break
                        
                        if not found and (res.resID[0] == "NTR" or res.resID[0] == "NTG") and atom.name == " CA ":
                            # NTR CA connects to CB of next residue
                            res2 = self.protein.residue[i_res+1]
                            # find " CB "
                            found = False
                            for conf2 in res2.conf[1:]:
                                for atom2 in conf2.atom:
                                    if atom2.name == " CB ":
                                        r = (atom.r_vdw + atom2.r_vdw) * _BONDDISTANCE_scaling
                                        CUTOFF2 = r * r
                                        if ddvv(atom.xyz, atom2.xyz) < CUTOFF2:
                                            if atom2 not in atom.connect12:
                                                atom.connect12.append(atom2)
                                                found = True

                        if not found and atom.name == " N  ":
                            # " N  " in the first residue of a chain, a de facto NTR
                            if i_res == 0 or self.protein.residue[i_res-1].resID[1] != atom.chainID:
                                found = True
                            

                        if not found and res.resID[0] == "CTR" and atom.name == " C  ":   # CTR C connects to CA of previous atom
                            res2 = self.protein.residue[i_res-1]
                            found = False
                            for atom2 in res2.conf[0].atom:
                                if atom2.name == " CA ":
                                    r = (atom.r_vdw + atom2.r_vdw) * _BONDDISTANCE_scaling
                                    CUTOFF2 = r * r
                                    if ddvv(atom.xyz, atom2.xyz) < CUTOFF2:
                                        if atom2 not in atom.connect12:
                                            atom.connect12.append(atom2)
                                            found = True

                        if not found:  # no actual CTR case
                            if atom.name == " C  ":
                                #print("No CTR after atom \"%s\"?" % atom.atomID)
                                found = True

                        # if not found:
                        #     if not ("CTR" in atom.atomID):  # ignore CTR due to CA not specified as ligand
                        #         logging.warning("Ligand atom bond to \"%s\" was not found" % atom.atomID)
                        #         print(",".join([a.atomID for a in atom.connect12]))

                    else:  # a named atom
                        # 1) backbone
                        for atom2 in res.conf[0].atom:
                            if atom2.name == c_atom:
                                atom.connect12.append(atom2)
                                found = True
                                break
                        # 2) own conformer
                        if not found:
                            for atom2 in conf.atom:
                                if atom2.name == c_atom:
                                    atom.connect12.append(atom2)
                                    found = True
                                    break
                        # 3) a backbone atom could connect to side chain conformer atoms
                        if not found:
                            if conf.confType[3:5] == "BK":
                                for conf2 in res.conf[1:]:
                                    for atom2 in conf2.atom:
                                        if atom2.name == c_atom:
                                            atom.connect12.append(atom2)
                                            found = True
                        # 4.1) NTR case, residue before
                        if not found:
                            if atom.name == " CB " and c_atom == " CA ":  # NTR separated from this res
                                for conf2 in self.protein.residue[i_res-1].conf[1:]:
                                    for atom2 in conf2.atom:
                                        if atom2.name == c_atom:
                                            atom.connect12.append(atom2)
                                            found = True
                                            break

                        # 4.2) NTR case, " C  " connects " CA " of NTR
                        if not found:
                            if atom.name == " C  " and c_atom == " CA ":  # NTR separated from this res
                                for conf2 in self.protein.residue[i_res-1].conf[1:]:
                                    for atom2 in conf2.atom:
                                        if atom2.name == c_atom:
                                            atom.connect12.append(atom2)
                                            found = True
                                            break

                        # 5) CTR case, " CA " connects " C  " of CTR
                        if not found:
                            if atom.name == " CA " and c_atom == " C  ":  # CTR is separated from this residue
                                if i_res + 1 >= len(self.protein.residue): # last residue
                                    found = True
                                else:
                                    for conf2 in self.protein.residue[i_res+1].conf[1:]:
                                        for atom2 in conf2.atom:
                                            if atom2.name == c_atom:
                                                atom.connect12.append(atom2)
                                                found = True
                                                break

                        # 6) " CB " connects to " CA " of NTR case
                        if not found:
                            if atom.name == " CB " and c_atom == " CA ":
                                res2 = self.protein.residue[i_res - 1]
                                for conf2 in res2.conf:
                                    for atom2 in conf2.atom:
                                        if atom2.name == " CA ":
                                            r = (atom.r_vdw + atom2.r_vdw) * _BONDDISTANCE_scaling
                                            CUTOFF2 = r * r
                                            if ddvv(atom.xyz, atom2.xyz) < CUTOFF2:
                                                if atom2 not in atom.connect12:
                                                    atom.connect12.append(atom2)
                                                    found = True
                                                    break

                        if not found:
                            # Ignore missing H
                            if is_H(c_atom) or is_H(atom.name):
                                pass
                            else:
                                logging.warning("Atom \"%s\" bond to \"%s\" was not found" % (c_atom, atom.atomID))

    return

def reset_connect(self):
    """Reset 1-2 connectivity for each atom in protein to []
    This function is useful when making conformers. If atom.connect12 is not reset, the connectivity will be messed up
    """
    for res in self.protein.residue:
        for conf in res.conf:
            for atom in conf.atom:
                atom.connect12 = []
                atom.connect13 = []
                atom.connect14 = []

def print_connect12(self):
    for res in self.protein.residue:
        for conf in res.conf:
            for atom in conf.atom:
                key = ("CONNECT", atom.name, atom.confType)
                connected_atoms = self.tpl.db[key].connected
                print(atom.atomID, str(connected_atoms))
                for atom2 in atom.connect12:
                    print("   -> %s" % atom2.atomID)

def make_connect13(self):
    for res in self.protein.residue:
        res.serialize()
        for conf in res.conf:
            for atom in conf.atom:
                for atom2 in atom.connect12:
                    if atom2 != atom:
                        for atom3 in atom2.connect12:
                            if (atom3 != atom) and (atom3 not in atom.connect12) and (atom3 not in atom.connect13):
                                if atom3.confNum == 0: # connectivity to backbone allowed
                                    atom.connect13.append(atom3)
                                else:
                                    if atom3.resID != atom.resID: # connectivity outside residue allowed
                                        atom.connect13.append(atom3)
                                    else:  # now atom3 is non backbone, within the residue, after above conditions 
                                        if atom3.confNum == atom.confNum:  # connectivity within sidechain allowed
                                            atom.connect13.append(atom3)

                    else:
                        print("Warning: Atom \"%s\" has itself in connect12." % atom.atomID)

def make_connect14(self):
    for res in self.protein.residue:
        res.serialize()
        for conf in res.conf:
            for atom in conf.atom:
                for atom3 in atom.connect13:
                    if atom3 == atom: continue
                    for atom4 in atom3.connect12:
                        if (atom4 != atom) \
                                and (atom4 not in atom.connect12) \
                                and (atom4 not in atom.connect13) \
                                and (atom4 not in atom.connect14):
                                if atom4.confNum == 0: # connectivity to backbone allowed
                                    atom.connect14.append(atom4)
                                else:
                                    if atom4.resID != atom.resID: # connectivity outside residue allowed
                                        atom.connect14.append(atom4)
                                    else:  # now atom3 is non backbone, within the residue, after above conditions 
                                        if atom4.confNum == atom.confNum:  # connectivity within sidechain allowed
                                            atom.connect14.append(atom4)


