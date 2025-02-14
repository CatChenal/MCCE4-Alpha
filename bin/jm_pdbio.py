#!/usr/bin/env python
"""Sub routines to input and write out MCCE PDB files."""

import math
import os
import logging
import glob
import time
import copy
import numpy as np

from geom import *


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


KCAL2KT = 1.688
# bond distance scaling factor: cutoff = k*(r_vdw1 + r_vdw2)
BONDDISTANCE_scaling = 0.54  # calibrated by 1akk


# scaling factor for vdw
VDW_SCALE14 = 0.5

# set any big conf vdw to 999
VDW_CUTOFF_FAR = 10     # set atom vdw to 0 if atoms are further than this value
VDW_CUTOFF_FAR2 = VDW_CUTOFF_FAR * VDW_CUTOFF_FAR     # set atom vdw to 0 if atoms are further than this value
VDW_CUTOFF_NEAR2 = 1      # set atom vdw to 999 if atoms are closer than this value
VDW_UPLIMIT = 999.0     # set conf vdw to 999 if bigger than this number


def ddvv(xyz1, xyz2):
    """Distance squared between two vectors."""
    dx=xyz1[0]-xyz2[0]
    dy=xyz1[1]-xyz2[1]
    dz=xyz1[2]-xyz2[2]
    return dx*dx+dy*dy+dz*dz

class Atom:
    def __init__(self):
        self.serial = 0
        self.name = "  X "
        self.altLoc = " "
        self.resName = "UNK"
        self.chainID = "A"
        self.resSeq = 0
        self.iCode = " "
        self.confNum = 0
        self.atomID = ""
        self.confID = ""
        self.confType = ""
        self.resID = ""
        self.xyz = (0.0, 0.0, 0.0)
        self.connectivity_param = ""
        self.r_bound = 0.0
        self.charge = 0.0
        self.r_vdw = 0.0
        self.e_vdw = 0.0
        self.connect12 = []
        self.connect13 = []
        self.connect14 = []
        self.history = ""
        return

    def loadline(self, line):
        self.serial = int(line[6:11])
        self.name = line[12:16]
        self.altLoc = line[16]
        self.resName = line[17:20]
        self.chainID = line[21]
        self.resSeq = int(line[22:26])
        self.iCode = line[26]
        self.confNum = int(line[27:30])
        self.xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        self.r_bound = float(line[54:62])
        self.charge = float(line[62:74])
        self.confType = "%3s%2s" % (self.resName, line[80:82])
        self.history = line[80:].strip()

        self.compose_atomID()
        self.confID = "%5s%c%04d%c%03d" % (self.confType, self.chainID, self.resSeq, self.iCode, self.confNum)
        self.resID = "%3s%04d%c" % (self.resName, self.resSeq, self.chainID)

        # extended records
        connect_key = ("CONNECT", self.name, self.confType)
        self.connectivity_param = env.param[connect_key]
        radius_key = ("RADIUS", self.confType, self.name)
        if radius_key in env.param:
            radius_values = env.param[radius_key]
            self.r_vdw = radius_values.r_vdw
            self.e_vdw = radius_values.e_vdw
        else: # Use default value as C
            print("Warning, parameter %s not found, using default values as C" % str(radius_key))
            self.r_vdw = 1.908
            self.e_vdw = 0.086
        return

    def compose_atomID(self):
        self.atomID = "%4s%3s%04d%c%03d" % (self.name, self.resName, self.resSeq, self.chainID, self.confNum)

    def writeline(self):
        line = "ATOM  %5d %4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%8.3f%12.3f\n" %\
               (self.serial, self.name, self.altLoc, self.resName, self.chainID,\
                self.resSeq, self.iCode, self.xyz[0], self.xyz[1],\
                self.xyz[2], self.r_bound, self.charge)
        return line

class Blob:
    def __init__(self, conf):
        self.center = (0.0, 0.0, 0.0)
        self.radius = 0.0
        x = 0.0
        y = 0.0
        z = 0.0
        for atom in conf.atom:
            x += atom.xyz[0]
            y += atom.xyz[1]
            z += atom.xyz[2]
        n = len(conf.atom)

        if n > 0:
            self.center = (x/n, y/n, z/n)

            r_max = max([x.r_vdw for x in conf.atom])

            d_far2 = 0.0
            atom_far = ""
            for atom in conf.atom:
                d2 = ddvv(self.center, atom.xyz)
                if d2 > d_far2:
                    d_far2 = d2
                    atom_far = atom
            d_far = math.sqrt(d_far2)

            self.radius = d_far + r_max

        return


class Conformer:
    def __init__(self):
        self.confID = ""
        self.resID = ""
        self.i = 0    # index in vdw matrix
        self.atom = []
        self.vdw0 = 0.0
        self.vdw1 = 0.0
        self.crg = 0.0
        self.tors = 0.0
        self.history = ""
        self.mark = ""
        return
    
    def update_crg(self):
        self.crg = 0.0
        for atom in self.atom:
            self.crg += atom.charge
        return

class Residue:
    def __init__(self):
        self.resID = ""
        self.conf = []
        return

class Protein:
    def __init__(self):
        self.residue = []
        self.vdw_pw = {}  # place holder, dictionary with confID pair as the key
        return

    def loadpdb(self, fname):
        rawlines = open(fname).readlines()
        lines = [x.strip("\n") for x in rawlines if x[:6] == "ATOM  " or x[:6] == "HETATM"]
        for line in lines:
            atom = Atom()
            atom.loadline(line)
            # reverse search if this atom belongs to an existing conformer and residue
            found_conf = False
            found_res = False
            for i_res in range(len(self.residue) - 1, -1, -1):
                insert_res = i_res
                for i_conf in range(len(self.residue[i_res].conf) - 1, -1, -1):
                    if self.residue[i_res].conf[i_conf].confID == atom.confID:
                        found_conf = True
                        self.residue[i_res].conf[i_conf].atom.append(atom)
                        break
                if found_conf:
                    break
                elif self.residue[i_res].resID == atom.resID:  # not in conf but residue ID matches
                    conf = Conformer()
                    conf.confID = atom.confID
                    conf.resID = atom.resID
                    conf.history = atom.history
                    conf.atom.append(atom)
                    self.residue[i_res].conf.append(conf)
                    found_res = True

            if not found_conf and not found_res:  # new residue
                conf = Conformer()
                conf.confID = atom.confID
                conf.resID = atom.resID
                conf.history = atom.history
                conf.atom.append(atom)
                res = Residue()
                res.resID = conf.resID
                res.conf.append(conf)
                self.residue.append(res)

        # Insert an empty conformer for cofactors that do not have backbone
        for res in self.residue:
            first_confID = res.conf[0].confID
            if first_confID[3:5] != "BK":
                conf = Conformer()
                conf.confID = "%sBK%s000" % (first_confID[:3], first_confID[5:11])
                conf.resID = res.resID
                res.conf.insert(0, conf)

        # Assign an index number to each conformer
        n_conf = 0
        for res in self.residue:
            if len(res.conf) <= 1:   # backbone only
                continue
            for conf in res.conf[1:]:
                conf.i = n_conf
                n_conf += 1

        # Create a blob of conformer to screen vdw calculation
        for res in self.residue:
            for conf in res.conf:
                conf.blob = Blob(conf)
        return

    def make_connect12(self):
        # make connect table, need to start from connect records!
        # an optimized version should search connect12 for backbone and side chain separately.
        # backbone only connects the residue before and after, while side chain may connect with all other side chains
        # Special rules for backbone atoms:
        #    if ligated, search only residue before and after
        #    if connected atom is missing (CA), search residue before and after in case it is terminal residue
        #    all side chain conformers

        for i_res in range(len(self.residue)):
            res = self.residue[i_res]
            for conf in res.conf:  # search all conformers including backbone
                for atom in conf.atom:
                    connect_key = ("CONNECT", atom.name, atom.confType)
                    connected_atoms = env.param[connect_key].connected
                    for c_atom in connected_atoms:
                        found = False
                        if "?" in c_atom:  # ligated
                            for res2 in self.residue:
                                if res != res2:
                                    for conf2 in res2.conf:
                                        for atom2 in conf2.atom:
                                            connect_key2 = ("CONNECT", atom2.name, atom2.confType)
                                            connected_atoms2 = env.param[connect_key2].connected
                                            ligated2 = False
                                            for ligated_atom2 in connected_atoms2:
                                                if "?" in ligated_atom2:
                                                    ligated2 = True
                                                    break
                                            if ligated2:
                                                r = (atom.r_vdw + atom2.r_vdw) * BONDDISTANCE_scaling
                                                CUTOFF2 = r * r
                                                # if " C  " in atom.atomID:
                                                #     print("%s <-> %s: d=%.3f cut=%.3f" % (atom.atomID, atom2.atomID, math.sqrt(ddvv(atom.xyz, atom2.xyz)), math.sqrt(CUTOFF2)))
                                                #     print(atom2.atomID, [x.atomID for x in atom.connect12])
                                                #     print("Here", found)
                                                if ddvv(atom.xyz, atom2.xyz) < CUTOFF2:
                                                    if atom2 not in atom.connect12:
                                                        atom.connect12.append(atom2)
                                                        found = True  # after ligand found, do not break, continue to search other conformers within residue

                                    if found:   # one "?" for one ligand
                                        # print(conf.confID, atom.atomID, c_atom, res2.resID)
                                        break

                            if not found and res.resID[:3] == "NTR" and atom.name == " CA ":   # NTR CA connects to CB of next residue
                                res2 = self.residue[i_res+1]
                                # find " CB "
                                found = False
                                for conf2 in res2.conf[1:]:
                                    for atom2 in conf2.atom:
                                        if atom2.name == " CB ":
                                            r = (atom.r_vdw + atom2.r_vdw) * BONDDISTANCE_scaling
                                            CUTOFF2 = r * r
                                            if ddvv(atom.xyz, atom2.xyz) < CUTOFF2:
                                                if atom2 not in atom.connect12:
                                                    atom.connect12.append(atom2)
                                                    found = True

                            if not found and res.resID[:3] == "CTR" and atom.name == " C  ":   # CTR C connects to CA of previous atom
                                res2 = self.residue[i_res-1]
                                # find " CA "
                                found = False
                                for atom2 in res2.conf[0].atom:
                                    if atom2.name == " CA ":
                                        r = (atom.r_vdw + atom2.r_vdw) * BONDDISTANCE_scaling
                                        CUTOFF2 = r * r
                                        if ddvv(atom.xyz, atom2.xyz) < CUTOFF2:
                                            if atom2 not in atom.connect12:
                                                atom.connect12.append(atom2)
                                                found = True

                            if not found:  # no actual CTR case
                                if atom.name == " C  ":
                                    #print("No CTR after atom \"%s\"?" % atom.atomID)
                                    found = True

                            if not found:
                                if not "CTR" in atom.atomID:  # ignore CTR due to CA not specified as ligand
                                    print("Warning: Ligand atom bond to \"%s\" was not found" % atom.atomID)

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
                                if conf.confID[3:5] == "BK":
                                    for conf2 in res.conf[1:]:
                                        for atom2 in conf2.atom:
                                            if atom2.name == c_atom:
                                                atom.connect12.append(atom2)
                                                found = True
                            # 4) NTR case, residue before
                            if not found:
                                if atom.name == " C  " and c_atom == " CA ":  # NTR separated from this res
                                    for conf2 in self.residue[i_res-1].conf[1:]:
                                        for atom2 in conf2.atom:
                                            if atom2.name == c_atom:
                                                atom.connect12.append(atom2)
                                                found = True
                                                break

                            # 5) CTR case, " C  " bond to " CA " not found
                            if not found:
                                if atom.name == " CA " and c_atom == " C  ":  # CTR is separated from this residue
                                    if i_res + 1 >= len(self.residue): # last residue
                                        found = True
                                    else:
                                        for conf2 in self.residue[i_res+1].conf[1:]:
                                            for atom2 in conf2.atom:
                                                if atom2.name == c_atom:
                                                    atom.connect12.append(atom2)
                                                    found = True
                                                    break

                            # 6) " CB " connects to " CA " of NTR case
                            if not found:
                                if atom.name == " CB " and c_atom == " CA ":
                                    res2 = self.residue[i_res - 1]
                                    for conf2 in res2.conf:
                                        for atom2 in conf2.atom:
                                            if atom2.name == " CA ":
                                                r = (atom.r_vdw + atom2.r_vdw) * BONDDISTANCE_scaling
                                                CUTOFF2 = r * r
                                                if ddvv(atom.xyz, atom2.xyz) < CUTOFF2:
                                                    if atom2 not in atom.connect12:
                                                        atom.connect12.append(atom2)
                                                        found = True
                                                        break

                            if not found:
                                print("Warning: Atom \"%s\" bond to \"%s\" was not found" % (c_atom, atom.atomID))

        return

    def print_connect12(self):
        for res in self.residue:
            for conf in res.conf:
                for atom in conf.atom:
                    print(atom.atomID)
                    for atom2 in atom.connect12:
                        print("   -> %s" % atom2.atomID)
        return

    def make_connect13(self):
        for res in self.residue:
            for conf in res.conf:
                for atom in conf.atom:
                    for atom2 in atom.connect12:
                        if atom2 != atom:
                            for atom3 in atom2.connect12:
                                if (atom3 != atom) and (atom3 not in atom.connect12) and (atom3 not in atom.connect13):
                                    atom.connect13.append(atom3)
                        else:
                            print("Warning: Atom \"%s\" has itself in connect12." % atom.atomID)
        return

    def make_connect14(self):
        for res in self.residue:
            for conf in res.conf:
                for atom in conf.atom:
                    for atom3 in atom.connect13:
                        if atom3 == atom: continue
                        for atom4 in atom3.connect12:
                            if (atom4 != atom) \
                                    and (atom4 not in atom.connect12) \
                                    and (atom4 not in atom.connect13) \
                                    and (atom4 not in atom.connect14):
                                atom.connect14.append(atom4)
        return

    def print_connect13(self):
        for res in self.residue:
            for conf in res.conf:
                for atom in conf.atom:
                    print(atom.atomID)
                    for atom2 in atom.connect13:
                        print("   -X-> %s" % atom2.atomID)
        return

    def print_connect14(self):
        for res in self.residue:
            for conf in res.conf:
                for atom in conf.atom:
                    print(atom.atomID)
                    for atom2 in atom.connect14:
                        print("   -X-X-> %s" % atom2.atomID)
        return

    def exportpdb(self, fname):
        lines = []
        for res in self.residue:
            for conf in res.conf:
                for atom in conf.atom:
                    lines.append(atom.writeline())
        open(fname, "w").writelines(lines)
        return

    def print_atom_structure(self):
        for res in self.residue:
            print("Residue %s" % res.resID)
            for conf in res.conf:
                print("-->Conformer %s" % conf.confID)
                for atom in conf.atom:
                    print("---->Atom %s" % atom.atomID)
        return

    def calc_vdw(self, verbose=False):
        # do it on two sides so the two-way interaction numbers can be checked
        for res1 in self.residue:

            if len(res1.conf) <= 1:
                continue    # only backbone
            for conf1 in res1.conf[1:]:
                if verbose:
                    print("   vdw - %s ..." % conf1.confID)
                saved_time = time.time()

                # compute vdw0
                conf1.vdw0 = vdw_conf(conf1, conf1)

                for res2 in self.residue:
                    if len(res2.conf) == 1:
                        continue  # only backbone
                    if res1 == res2: # we need to do self to self vdw - vdw0
                        vdw = vdw_conf(conf1, conf1)
                        if abs(vdw) > 0.001:
                            self.vdw_pw[(conf1.confID, conf1.confID)] = vdw
                        continue
                    for conf2 in res2.conf[1:]:
                        vdw = vdw_conf(conf1, conf2)
                        if abs(vdw) > 0.001:
                            #print("%s - %s: %.3f" % (conf1.confID, conf2.confID, vdw))
                            self.vdw_pw[(conf1.confID, conf2.confID)] = vdw


                # compute vdw1, vdw to all backbone
                vdw1 = 0.0
                for res2 in self.residue:
                    conf2 = res2.conf[0]
                    vdw = vdw_conf(conf1, conf2)
                    vdw1 += vdw
                conf1.vdw1 = vdw1
                elapsed = time.time() - saved_time
                #print("   %s ... %.3f seconds" % (conf1.confID, elapsed))

    def calc_vdw_virtual(self, delta="", verbose=False):
        # Create virtual conformers
        if delta:
            for res in self.residue:
                conflist = env.param[("CONFLIST", res.resID[:3])]
                # detect dummy conformer in the list
                conf_types = set([c[-2:] for c in conflist])
                if "DM" in conf_types:
                    if len(res.conf) > 1:  # Only multiply none backbone conformer
                        new_confs = [res.conf[0]]
                        for conf in res.conf[1:]:
                            new_confs.append(conf)
                            new_conf = copy.deepcopy(conf)
                            new_conf.confID = new_conf.confID[:10]+"-"+new_conf.confID[11:]
                            new_conf.history = "VIRTUAL-"
                            for atom in new_conf.atom:
                                atom.r_vdw -= delta
                                if atom.r_vdw < 0.0:
                                    atom.r_vdw = 0.0
                            new_confs.append(new_conf)

                            new_conf = copy.deepcopy(conf)
                            new_conf.confID = new_conf.confID[:10]+"+"+new_conf.confID[11:]
                            new_conf.history = "VIRTUAL+"
                            for atom in new_conf.atom:
                                atom.r_vdw += delta
                                if atom.r_vdw < 0.0:
                                    atom.r_vdw = 0.0
                            new_confs.append(new_conf)
                        
                        res.conf = new_confs

        # Recalculate vdw
        self.calc_vdw(verbose=verbose)        


    def calc_tors(self) -> None:
        """Calculate torsion energy for all conformers, update conformer tors in-place.
        """
        for res in self.residue:
            if len(res.conf) <= 1: continue    # only backbone
            for conf in res.conf[1:]:
                conf.tors = torsion_conf(conf)

    def calc_tors_virtual(self):
        """
        Calculate torsion energy for vdw virtual conformers by inheriting from their parent conformers.
        """
        tors_db = {}  # store torsion energy in a dictionary
        for res in self.residue:
            if len(res.conf) > 1:  # no need to do bk conf
                for conf in res.conf[1:]:
                    if conf.confID[10] not in {"+", "-"}:   # not a virtual conformer
                        tors_db[conf.confID] = conf.tors
        # check vdw virtual conformers and assign their parent tors values
        for res in self.residue:
            if len(res.conf) > 1:  # no need to do bk conf
                for conf in res.conf[1:]:
                    if conf.confID[10] in {"+", "-"}:   # is virtual conformer
                        parent_confID = conf.confID[:10] + "_" + conf.confID[11:]
                        conf.tors = tors_db[parent_confID]

    
    def connect_reciprocity_check(self):
        # connectivity should be reciprocal except backbone atoms
        for res in self.residue:
            for conf in res.conf:
                for atom in conf.atom:
                    for atom2 in atom.connect12:
                        if atom not in atom2.connect12:
                            print("Atom %s in connect12 of atom %s but the other way is not true" % (atom2.atomID, atom.atomID))
                    for atom2 in atom.connect13:
                        if atom2.confType[-2:] == "BK":
                            continue
                        if atom not in atom2.connect13:
                            print("Atom %s in connect13 of atom %s but the other way is not true" % (atom2.atomID, atom.atomID))
                    for atom2 in atom.connect14:
                        if atom2.confType[-2:] == "BK":
                            continue
                        if atom not in atom2.connect14:
                            print("Atom %s in connect14 of atom %s but the other way is not true" % (atom2.atomID, atom.atomID))

    def vdw_reciprocity_check(self):
        for key, value in self.vdw_pw.items():
            r_key = (key[1], key[0])
            if r_key in self.vdw_pw:
                r_value = self.vdw_pw[r_key]
                if abs(value - r_value) > 0.001:
                    print("vdw(%s<->%s = %.3f <-> %.3f" % (key[0], key[1], value, r_value))
            else:
                print("vdw(%s->%s = %.3f but other way not reported" % (key[0], key[1], value))

        print("VDW two sides checked.")
        return

    def print_confindex(self):
        for res in self.residue:
            if len(res.conf) > 1:
                for conf in res.conf[1:]:
                    print("%s -> %4d" % (conf.confID, conf.i))

        return
    
    def update_confcrg(self):
        for res in self.residue:
            for conf in res.conf:
                conf.update_crg()
        return

class CONNECT_param:
    def __init__(self, value_str):
        fields = value_str.split(",")
        self.orbital = fields[0].strip()
        self.connected = [x.strip().strip("\"") for x in fields[1:]]


class RADIUS_param:
    def __init__(self, value_str):
        fields = value_str.split(",")
        self.r_bound = float(fields[0])
        self.r_vdw = float(fields[1])
        self.e_vdw = float(fields[2])

class CONFORMER_param:
    def __init__(self, value_str):
        self.param = {}
        fields = value_str.split(",")
        for f in fields:
            sf = f.split("=")
            key = sf[0].strip().lower()
            value = float(sf[1])
            self.param[key] = value


class TORS:
    """Internal torsion data structure to hold v2, nfold and gamma"""
    def __init__(self) -> None:
        """
        Construct elements v2, nfold, and gamma for one torsion term
        :return: None
        """
        self.v2 = 0.0
        self.n_fold = 0
        self.gamma = 0.0


class TORSION_param:
    """This class defines the data structure of TORSION parameters."""
    def __init__(self, value_str) -> None:
        """Construct the TORSION paramter and assign initial values from value_str.

        :param value_str: A string read from parameter file that defines the torsion parameter.
        :type value_str: string
        :return: None
        """
        d2r = math.pi/180.0  # cofficient to convert degrees to radians

        self.atom1 = value_str[:4]
        self.atom2 = value_str[5:9]
        self.atom3 = value_str[10:14]
        self.tors_terms = []
        fields = value_str[15:].split()
        relx = fields.pop(0).strip().upper()
        if relx == "F":
            self.relx = False
        else:
            self.relx = True

        n_term = int(len(fields)/3)
        for i in range(n_term):
            tors = TORS()
            tors.v2 = float(fields[i*3+0])
            tors.n_fold = int(fields[i*3+1])
            tors.gamma = float(fields[i*3+2]) * d2r
            self.tors_terms.append(tors)

class ENV:
    def __init__(self):
        self.runprm = {}
        self.param = {}
        #self.load_runprm()
        #self.print_runprm()
        #self.load_ftpl()

    def load_runprm(self):
        filename = "run.prm"
        if os.path.isfile(filename):
            lines = open(filename).readlines()
            for line in lines:
                entry_str = line.strip().split("#")[0]
                fields = entry_str.split()
                if len(fields) > 1:
                    key_str = fields[-1]
                    if key_str[0] == "(" and key_str[-1] == ")":
                        key = key_str.strip("()").strip()
                        value = fields[0]
                    self.runprm[key] = value
        else:  # assign key runprm with default values
            dist_folder = os.path.dirname(os.path.dirname(__file__))
            self.runprm["MCCE_HOME"] = dist_folder

    def print_runprm(self):
        for key, value in self.runprm.items():
            print("%s:%s" % (key, value))

    def read_ftpl_file(self, fname):
        lines = open(fname).readlines()
        for line in lines:
            end = line.find("#")
            line = line[:end]
            fields = line.split(":")
            if len(fields) != 2:
                continue

            key_string = fields[0].strip()
            keys = key_string.split(",")
            key1 = keys[0].strip().strip("\"")
            if len(keys) > 1:
                key2 = keys[1].strip().strip("\"")
            else:
                key2 = ""
            if len(keys) > 2:
                key3 = keys[2].strip().strip("\"")
            else:
                key3 = ""

            value_string = fields[1].strip()

            # Connectivity records
            if key1 == "CONFLIST":
                self.param[(key1, key2)] = [x.strip() for x in value_string.strip().split(",")]
            elif key1 == "CONNECT":
                self.param[(key1, key2, key3)] = CONNECT_param(value_string)
            # VDW parameters, for now use 00always_needed.tpl for vdw parameters
            elif key1 == "RADIUS":
                self.param[(key1, key2, key3)] = RADIUS_param(value_string)
            elif key1 == "CONFORMER":
                self.param[(key1, key2)] = CONFORMER_param(value_string)

    def load_ftpl(self):
        if "FTPLDIR" in self.runprm:
            ftpldir = self.runprm["FTPLDIR"]
        else:
            ftpldir = self.runprm["MCCE_HOME"]+"/param"
        cwd = os.getcwd()
        os.chdir(ftpldir)

        files = glob.glob("*.ftpl")
        files.sort()
        logger.info("Reading parameters from %s" % ftpldir)
        for fname in files:
            self.read_ftpl_file(fname)


        # Update vdw with 00always_needed.tpl
        fname = "00always_needed.tpl"
        lines = open(fname).readlines()
        logger.info("Updating vdw parameters from %s" % os.path.abspath(fname))
        for line in lines:
            end = line.find("#")
            line = line[:end].strip()
            if len(line) < 20:
                continue
            key1 = line[:9].strip()
            key2 = line[9:15].strip()
            key3 = line[15:19]
            if key1=="VDW_RAD" or key1=="VDW_EPS":
                value = float(line[20:].strip())
                # print("%s %s \"%s\" : %.4f" % (key1, key2, key3, value))
                new_key = ("RADIUS", key2, key3)
                if new_key in self.param:
                    param_value = self.param[new_key]
                else:
                    print("Warning: %s %s \"%s\" not defined in RADIUS parameter" % (key1, key2, key3))
                    param_value = RADIUS_param("2,0,0")
                if key1 == "VDW_RAD":
                    param_value.r_vdw = value
                elif key1 == "VDW_EPS":
                    param_value.e_vdw = value
                self.param[new_key] = param_value

        logger.info("Loading TORSION parameters from %s" % os.path.abspath(fname))
        for line in lines:
            end = line.find("#")
            line = line[:end].strip()
            if len(line) < 20:
                continue
            key1 = line[:9].strip()
            key2 = line[9:15].strip()
            key3 = line[15:19]
            if key1 == "TORSION":
                new_key = (key1, key2, key3)  # key1: TORSION, key2: residue name, key3: 4-char atom name
                param_value = TORSION_param(line[20:])
                self.param[new_key] = param_value


        os.chdir(cwd)

        # read from user_param
        ftpldir = "user_param"
        if os.path.isdir(ftpldir):
            print("Reading parameters from %s" % ftpldir)
            os.chdir(ftpldir)
            files = glob.glob("*.ftpl")
            files.sort()
            for fname in files:
                self.read_ftpl_file(fname)
            os.chdir(cwd)

        # read and convert from new.tpl
        # Only CONNECT records are read and converted from new.tpl
        newtpl = "new.tpl"
        if os.path.exists(newtpl):
            lines = open(newtpl).readlines()

            print("Reading CONFLIST parameters from %s" % newtpl)
            for line in lines:
                if len(line) > 10:
                    key1 = line[:9].strip()
                    if key1 == "CONFLIST":
                        key2 = line[9:12]
                        value = line[20:].strip().split()
                        self.param[(key1, key2)] = value

            print("Reading CONNECT parameters from %s" % newtpl)
            for line in lines:
                if len(line) > 10:
                    key1 = line[:9].strip()
                    if key1 == "CONNECT":
                        key3 = line[9:14]
                        key2 = line[15:19]
                        orbital = line[20:29].strip()
                        atoms_str = line[30:].rstrip()+"   "
                        n_atoms = len(atoms_str)//10
                        atoms = []
                        for i in range(n_atoms):
                            a = atoms_str[:10]
                            ires = a[:5].strip()
                            atom_name = a[5:9]
                            if ires != "0":
                                atom_name = " ?  "
                            atoms.append(atom_name)
                            atoms_str = atoms_str[10:]

                        new_atoms_str = ",".join(["\"%s\"" % x for x in atoms])
                        value_string = "%s, %s" % (orbital, new_atoms_str)
                        self.param[(key1, key2, key3)] = CONNECT_param(value_string)
                        #print("(%s, %s, %s): %s" % (key1, key2, key3, value_string))

        # read from extra.tpl
        # self.print_runprm()
        # print(self.runprm)
        if "EXTRA" in self.runprm:
            extratpl = self.runprm["EXTRA"]
        else:
            extratpl = self.runprm["MCCE_HOME"]+"/extra.tpl"
        lines = open(extratpl).readlines()
        for line in lines:
            if len(line) > 10:
                key1 = line[:9].strip()
                key2 = line[9:14].strip()
                value = float(line[20:].strip())
                self.param[(key1, key2)] = value

        return

    def print_param(self):
        for key, value in self.param.items():
            if len(key) == 3:
                key1, key2, key3 = key
            elif len(key) == 2:
                key1, key2 = key
            if key1 == "CONNECT":
                print("%s:%s, %s" % (key, value.orbital, value.connected))
            elif key1 == "RADIUS":
                print("%s: %6.3f, %6.3f %6.3f" % (key, value.r_bound, value.r_vdw, value.e_vdw))
            elif key1 == "TORSION":
                print("%s: %s, %s, %s" % (key, value.atom1, value.atom2, value.atom3))
                for term in value.tors_terms:
                    print("            V2=%.3f N=%d, Gamma=%.3f" % (term.v2, term.n_fold, term.gamma))
            else:
                print(key, value)

def vdw_conf(conf1, conf2, cutoff=0.001, verbose=False, display=False):
    vdw = 0.0

    d = math.sqrt(ddvv(conf1.blob.center, conf2.blob.center))
    #print(d, conf1.blob.radius + 6 + conf2.blob.radius)
    if d > conf1.blob.radius + 6 + conf2.blob.radius:
        if display:
            print("%s at (%.3f, %.3f, %.3f) <-> %s at (%.3f, %.3f, %.3f): d = %.3f" %
                  (conf1.confID, conf1.blob.center[0], conf1.blob.center[1], conf1.blob.center[2],
                   conf2.confID, conf2.blob.center[0], conf2.blob.center[1], conf2.blob.center[2],
                   d))

    else:
        if display:
            print("%12s %16s     %10s %8s %6s   %6s %6s %6s %6s %6s %6s" % ("ATOM1",
                                                                 "ATOM2",
                                                                 "vdw",
                                                                 "dist",
                                                                 "cnct",
                                                                 "r1",
                                                                 "e1",
                                                                 "r2",
                                                                 "e2",
                                                                            "R_sum",
                                                                            "E_par"))
        for atom1 in conf1.atom:
            for atom2 in conf2.atom:
                vdw_a2a = vdw_atom(atom1, atom2)
                vdw += vdw_a2a
                if verbose and abs(vdw_a2a) >= cutoff:
                    d2 = ddvv(atom1.xyz, atom2.xyz)
                    if atom1 == atom2:
                        connect = "self"
                    elif atom1 in atom2.connect12:
                        connect = "1--2"
                    elif atom1 in atom2.connect13:
                        connect = "1--3"
                    elif atom1 in atom2.connect14:
                        connect = "1--4"
                    else:
                        connect = "none"
                    if display:  # display details
                        R_sum = atom1.r_vdw+atom2.r_vdw
                        E_par = math.sqrt(atom1.e_vdw*atom2.e_vdw)
                        print("%s -> %s: %8.3f %8.3f %6s   %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f" % (atom1.atomID,
                                                                                                   atom2.atomID,
                                                                                                   vdw_a2a,
                                                                                                   math.sqrt(d2),
                                                                                                   connect,
                                                                                                   atom1.r_vdw,
                                                                                                   atom1.e_vdw,
                                                                                                   atom2.r_vdw,
                                                                                                   atom2.e_vdw,
                                                                                                   R_sum,
                                                                                                   E_par))
                    else:  # display essential information
                        print(
                            "%s -> %s: %8.3f %8.3f %6s" % (atom1.atomID, atom2.atomID, vdw_a2a, math.sqrt(d2), connect))
        if vdw >= VDW_UPLIMIT:
            vdw = 999.0

        if conf1 == conf2:
            vdw = 0.5*vdw

    return vdw

def vdw_lj(d2, eps, r0, scale):
    sig_d2 = r0 * r0 / d2
    sig_d6 = sig_d2 * sig_d2 * sig_d2
    sig_d12 = sig_d6 * sig_d6
    return scale * eps * (sig_d12 - 2.0 * sig_d6)

def vdw_atom(atom1, atom2):
    # A good post: https://mattermodeling.stackexchange.com/questions/4845/how-to-create-a-lookup-table-of-%CF%B5-and-%CF%83-values-for-lennard-jones-potentials
    # Parameter source: http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs
    #   CHARM36-jul2022.ff/ffnonbonded.itp
    #   σ and ε values are in nm and kJ/mol in this file
    # Lorentz-Berthelot combining rules:
    #   σij = 0.5*(σi + σj)
    #   ϵij = sqrt(ϵi * ϵj)
    # p_lj = ϵij[(σij/r)^12 - 2(σij/r)^6]
    # σij is the distance where LJ potential reaches minimum: -ϵij
    # r is the atom distance

    p_lj = 0.0

    # Using the equation in vdw.c. Need to work on new parameter set.
    if atom1.xyz[0] - atom2.xyz[0] < VDW_CUTOFF_FAR and \
       atom1.xyz[1] - atom2.xyz[1] < VDW_CUTOFF_FAR and \
       atom1.xyz[2] - atom2.xyz[2] < VDW_CUTOFF_FAR:
        d2 = ddvv(atom1.xyz, atom2.xyz)
        r = math.sqrt(d2)

        if atom1 != atom2 and atom2 not in atom1.connect12 and atom2 not in atom1.connect13:
            if d2 > VDW_CUTOFF_FAR2:
                p_lj = 0.0
            elif d2 < VDW_CUTOFF_NEAR2:
                p_lj = 999.0
            else:
                r1 = atom1.r_vdw
                e1 = atom1.e_vdw
                r2 = atom2.r_vdw
                e2 = atom2.e_vdw
                if atom2 in atom1.connect14:
                    scale = VDW_SCALE14
                else:
                    scale = 1.0

                r0 = r1 + r2
                eps = math.sqrt(e1 * e2)

                width = 1.0
                r0_min = r0 - width/2.0
                r0_max = r0 + width/2.0
                #r0_max = r0

                if r < r0_min:
                   p_lj = vdw_lj(d2, eps=eps, r0=r0_min, scale=scale)
                elif r0_min <= r < r0_max:
                   p_lj = -eps
                else:
                   p_lj = vdw_lj(d2, eps=eps, r0=r0_max, scale=scale)

                #p_lj = (
                #       vdw_lj(d2, eps=eps, r0=r0_min, scale=scale) * (1 - np.heaviside(r - r0_min, 1)) +
                #       (-eps) * (np.heaviside(r - r0_min, 1) - np.heaviside(r - r0_max, 1)) +
                #       vdw_lj(d2, eps=eps, r0=r0_max, scale=scale) * np.heaviside(r - r0_max, 1)
                #)

                # #print("===%s===%s===" % (atom1.atomID, atom2.atomID))
                # if (atom1.atomID == " HB2ASP0018A001" and atom2.atomID == " OD2ASP0018A001") or \
                #    (atom1.atomID == " HB2ASP0018A003" and atom2.atomID == " OD2ASP0018A003"):
                #     print("===%s -> %s: %8.3f===" % (atom1.atomID, atom2.atomID, p_lj))
                #     print("%s: r_vdw=%8.3f, e_vdw=%8.3f" % (atom1.atomID, atom1.r_vdw, atom1.e_vdw))
                #     print("%s: r_vdw=%8.3f, e_vdw=%8.3f" % (atom2.atomID, atom2.r_vdw, atom2.e_vdw))

    #     else:
    #         p_lj = 0.0
    # else:
    #     p_lj = 0.0

    return p_lj

def torsion_emprical(conf):  # estimate torsion energy by 1-4 vdw
    vdw = 0.0
    calculated_pairs = []
    for atom1 in conf.atom:
        for atom2 in atom1.connect14:
            if (atom2.confNum == atom1.confNum or atom2.confNum == 0) and ((atom2, atom1) not in calculated_pairs):
                vdw += vdw_atom(atom1, atom2)
                calculated_pairs.append((atom2, atom1))   # ensure only calculate a pair once

    return vdw


def torsion_angle(v0, v1, v2, v3):
    """ Calculate torsion angle of 4 points
    :param v0: coordinates of the 1st point
    :type v0: tuple (x, y, z)
    :param v1: coordinates of the 2nd point
    :type v1: tuple (x, y, z)
    :param v2: coordinates of the 3rd point
    :type v2: tuple (x, y, z)
    :param v3: coordinates of the 4th point
    :type v3: tuple (x, y, z)
    :return: torsion angle in radians
    
    C Code:
    float angle,cos_theta;
    VECTOR i,j,k;
    VECTOR r21,r10,r10_p,r23;
    
    r21 = vector_vminusv(v1,v2);
    r10 = vector_vminusv(v0,v1);
    r23 = vector_vminusv(v3,v2);
    k   = vector_normalize(r21);
    i   = vector_normalize(vector_vminusv(r23, vector_rescale(k, vdotv(k,r23))));
    j   = vector_vxv(k,i);
    r10_p = vector_normalize(vector_vminusv(r10, vector_rescale(k, vdotv(k,r10))));
    
    cos_theta = vdotv(r10_p,i);
    //printf("%8.3f,%8.3f,%8.3f,%8.3f,%8.3f,%8.3f,%8.3f\n",cos_theta,i.x,i.y,i.z,k.x,k.y,k.z);
    if (cos_theta > 1.) cos_theta = 1.;
    else if (cos_theta < -1.) cos_theta = -1.;
    angle = acos(cos_theta);
    
    if (vdotv(r10_p,j) < 0.)
        angle = 2.*env.PI - angle;
    
    return angle;
    """
    
    r21 = np.array(vector_vminusv(v1, v2))
    r10 = np.array(vector_vminusv(v0, v1))
    r23 = np.array(vector_vminusv(v3, v2))
    k   = vector_normalize(r21);
    i   = vector_normalize(vector_vminusv(r23, k * np.dot(k, r23)));
    j   = np.cross(k, i)
    r10_p = vector_normalize(vector_vminusv(r10, k * np.dot(k, r10)))

    cos_theta = np.dot(r10_p, i)
    if cos_theta > 1.0:
        cos_theta = 1.0
    elif cos_theta < -1.0:
        cos_theta = -1.0
    angle = np.arccos(cos_theta)

    if np.dot(r10_p, j) < 0.0:
        angle = 2*np.pi -angle

    return angle

    



def torsion_conf(conf):
    """
    Calculate total torsion energy of a given conformer
    :param conf: Input conformer
    :return: Torsion energy value

    C code:
    float e = 0, phi;
    int i_atom, i_term;
    ATOM *atom0_p, *atom1_p, *atom2_p, *atom3_p;
    TORS tors;
    
    for (i_atom=0; i_atom<conf_p->n_atom; i_atom++) {
        float k_tor_e = 0.;
        if (!conf_p->atom[i_atom].on) continue;
        if (torsion_atoms(conf_p, i_atom, &atom0_p, &atom1_p, &atom2_p, &atom3_p, &tors, 1) == -1) 
            continue;
        phi = torsion_angle(atom0_p->xyz,atom1_p->xyz,atom2_p->xyz,atom3_p->xyz);
        for (i_term=0;i_term<tors.n_term; i_term++) {
            k_tor_e += torsion(phi, tors.V2[i_term], tors.n_fold[i_term], tors.gamma[i_term]);
        }
        e += k_tor_e;

        /* print large single bond torsion 
        if (k_tor_e > 0.5)
            printf("Large torsion: %s-%s-%s-%s, %8.3f\n",atom0_p->name,atom1_p->name,atom2_p->name,atom3_p->name,k_tor_e);
        */
    }
    return e;
        
    """
    e = 0.0   # default value
    for atom in conf.atom:
        key1 = "TORSION"
        key2 = atom.confType
        key3 = atom.name
        key = (key1, key2, key3)
        k_tor_e = 0.0
        
        found = False
        if key in env.param:
            found = True
        else:
            key2 = atom.resName
            key = (key1, key2, key3)
            if key in env.param:
                found = True

        if found:
            t_param = env.param[key]
            atom0 = atom
            atom1_name = t_param.atom1
            atom2_name = t_param.atom2
            atom3_name = t_param.atom3
            # find atom 1, 2, and 3
            found = False
            for a1 in atom0.connect12:
                if found: break
                if atom1_name == a1.name:
                    for a2 in a1.connect12:
                        if found: break
                        if atom2_name == a2.name:
                            for a3 in a2.connect12:
                                if found: break
                                if atom3_name == a3.name:
                                    found = True
                                    atom1 = a1
                                    atom2 = a2
                                    atom3 = a3
                                    break
            if found:
                #print(atom0.name, atom1.name, atom2.name, atom3.name)
                phi = torsion_angle(atom0.xyz, atom1.xyz, atom2.xyz, atom3.xyz)
                for term in t_param.tors_terms:
                    k_tor_e += torsion(phi, term)

        e += k_tor_e

    return e


def torsion(phi, term):
    """Calculate torsion energy given torsion angle and parameters
    :param phi: torsion angle
    :type phi: float, angle in radians
    :return: torsion energy value in kcal/mol
    :rtype: float

    C Code:
    float e;
    e = V2 * (1. + cos(n_fold*phi - gamma));
    return e;
    """
    e = term.v2 * (1.0 + np.cos(term.n_fold*phi - term.gamma))

    return e


env = ENV()


if __name__ == "__main__":

    env.load_runprm()
    env.load_ftpl()

    env.print_param()

    # pdbfile = "step2_out.pdb"
    # protein = Protein()
    # protein.loadpdb(pdbfile)
    #protein.print_confindex()
    #protein.make_connect12()
    #protein.make_connect13()
    #protein.make_connect14()

    # protein.print_connect12()
    # protein.print_connect13()
    # protein.print_connect14()
    # protein.exportpdb("a.pdb")
    # protein.print_atom_structure()

    #protein.calc_vdw()
    # protein.connect_reciprocity_check()
    # protein.vdw_reciprocity_check()

    #vdw_by_conf_pair(protein, "ASPBKA0002_000", "NTG01A0001_001", 0.001)

