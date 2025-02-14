#!/usr/bin/env python

"""
Module: pbe_base.py

Holds classes needed by PBE Solvers and step3.py.

NOTE:
  WIP: Not yet tested with either step3.py or pbs_interfaces.py.
"""

import json
import logging
from mccesteps import record_runprm


logging.basicConfig(level=logging.DEBUG)


global protein, run_options


class ElePW:
    def __init__(self):
        self.mark = ""
        self.multi = 0.0
        self.single = 0.0
        self.scaled = 0.0
        self.averaged = 0.0
        return


class ExchangeAtom:
    def __init__(self, atom):
        self.x = atom.xyz[0]
        self.y = atom.xyz[1]
        self.z = atom.xyz[2]
        self.r = atom.r_bound  # radius
        self.c = 0.0  # default to boundary defining atom, charge should be set as 0
        self.p = 0.0
        return


class Exchange:
    """
    Holds the data passed to the PB wrapper, together with runoptions.
    """
    # We have to abandon the mop on and off mechanism when modifying the dielectric boundary.
    # In the parallel for loop, we don't want the boundary revising to be dependent on previous step.
    # Each step in the loop should be an addition or deletion from a static starting point.
    # Therefore, we will compose
    # * backbone atoms
    # * index to match multiple atoms to boundary line number
    # * method to compose single side chain condition
    # * method to compose multi side chain condition

    def __init__(self, protein):
        # initilaize backbone:
        self.backbone_xyzrcp = []
        self.backbone_atom = []
        for res in protein.residue:
            if res.conf:
                for atom in res.conf[0].atom:
                    xyzrcp = ExchangeAtom(atom)
                    self.backbone_xyzrcp.append(xyzrcp)
                    self.backbone_atom.append([atom])
                    # the atom is in a list because it is allowed to have multiple
                    # atoms to match the same line in xyzrcp line

        self.float_bnd_xyzrcp = []
        self.float_bnd_atom = []
        self.single_bnd_xyzrcp = []
        self.single_bnd_atom = []
        self.multi_bnd_xyzrcp = []
        self.multi_bnd_atom = []
        return

    def compose_float(self, protein, ir, ic):
        """Compose a floating side chain boundary condition for rxn0 calculation.
        Only atoms in residue[ir], conformer[ic] are in this list
        """
        self.float_bnd_xyzrcp = []
        self.float_bnd_atom = []

        for atom in protein.residue[ir].conf[ic].atom:
            xyzrcp = ExchangeAtom(atom)
            xyzrcp.c = atom.charge
            self.float_bnd_xyzrcp.append(xyzrcp)
            self.float_bnd_atom.append([atom])

        return

    def compose_single(self, protein, ir, ic):
        """Compose a single side chain boundary condition.
        The atoms are added in addition to backbone.
        Atoms other than in residue[ir], conformer[ic] are then appended.
        Atoms in residue[ir], conformer[ic] are appended last
        """
        self.single_bnd_xyzrcp = self.backbone_xyzrcp.copy()
        self.single_bnd_atom = self.backbone_atom.copy()

        for ires, protres in enumerate(protein.residue):
            # print(protres.resID)
            if ires == ir:  # this is the residue we want to put desired side chain conf
                for atom in protres.conf[ic].atom:
                    xyzrcp = ExchangeAtom(atom)
                    xyzrcp.c = atom.charge
                    self.single_bnd_xyzrcp.append(xyzrcp)
                    self.single_bnd_atom.append([atom])
            else:
                # find the first charged conformer if any, otherwise use the first
                if len(protres.conf) > 1:  # skip dummy or backbone only residue
                    i_useconf = 1  # default the first conformer
                    for iconf, protres_conf in enumerate(protres.conf, start=1):
                        # print(protres_conf.confID, protres_conf.crg)
                        if abs(protres_conf.crg) > 0.0001:
                            i_useconf = iconf
                            break

                    # just for boundary
                    for atom in protres.conf[i_useconf].atom:
                        xyzrcp = ExchangeAtom(atom)
                        self.single_bnd_xyzrcp.append(xyzrcp)
                        self.single_bnd_atom.append([atom])

        # Error checking, basic
        # for atom_list in self.ibound2atoms:
        #     print(len(atom_list), atom_list[0].atomID)
        #     if len(atom_list) != 1:
        #         print("ERROR")
        #         break

        if len(self.single_bnd_xyzrcp) != len(self.single_bnd_atom):
            logging.critical(
                "%s Single-sidechain boundary record length should be equal: %d, %d"
                % (protein.residue[ir].conf[ic].confID,
                   len(self.single_bnd_xyzrcp),
                   len(self.single_bnd_atom),
                   )
                )
        
        logging.debug(
            "%s Single-sidechain boundary record length should be equal: %d, %d"
            % (
                protein.residue[ir].conf[ic].confID,
                len(self.single_bnd_xyzrcp),
                len(self.single_bnd_atom),
            )
        )

        return

    def compose_multi(self, protein, ir, ic):
        """Compose a multi side chain boundary condition.
        The atoms are added in addition to backbone.
        Atoms other than in residue[ir], conformer[ic] are appended next.
        Atoms in residue[ir], conformer[ic] are appended last.
        When appending an atom, this subroutine will check if the same atom
        (xyzrc identical) already exists. If yes, just update icound
        """
        self.multi_bnd_xyzrcp = self.backbone_xyzrcp.copy()
        self.multi_bnd_atom = self.backbone_atom.copy()

        for ires, protres in enumerate(protein.residue):
            # print(protres.resID)
            if ires == ir:  # this is the residue we want to put desired side chain conf
                for atom in protres.conf[ic].atom:
                    xyzrcp = ExchangeAtom(atom)
                    xyzrcp.c = atom.charge
                    self.multi_bnd_xyzrcp.append(xyzrcp)
                    self.multi_bnd_atom.append([atom])
            else:
                # other residues will have all conformers with 0 charge.
                # skip dummy or backbone only residue
                if len(protres.conf) > 1:
                    # this is the xyzrcp record for all side chain atoms of this residue
                    residue_bnd_xyzrcp = (
                        []
                    )  
                    # this points to the atom records of each line in residue_bnd_xyzrpc
                    residue_bnd_atom = (
                        []
                    )
                    
                    for iconf, protres_conf in enumerate(protres.conf, start=1):
                        for atom in protres_conf.atom:
                            xyzrcp = ExchangeAtom(atom)
                            # test if this atom existed within this residue already
                            found = False
                            for ib, resbnd_xyzrcp in enumerate(residue_bnd_xyzrcp):
                                if (
                                    abs(resbnd_xyzrcp.x - xyzrcp.x) < 0.001
                                    and abs(resbnd_xyzrcp.y - xyzrcp.y) < 0.001
                                    and abs(resbnd_xyzrcp.z - xyzrcp.z) < 0.001
                                    and abs(resbnd_xyzrcp.r - xyzrcp.r) < 0.001
                                ):  # identical atom
                                    residue_bnd_atom[ib].append(atom)
                                    found = True
                                    break

                            if not found:
                                residue_bnd_xyzrcp.append(xyzrcp)
                                residue_bnd_atom.append([atom])

                    # merge this residue to the multi-bnd
                    self.multi_bnd_xyzrcp += residue_bnd_xyzrcp
                    self.multi_bnd_atom += residue_bnd_atom

        # Basic error checking
        if len(self.multi_bnd_xyzrcp) != len(self.multi_bnd_atom):
            logging.critical(
                "%s Multi-sidechain boundary record length should be equal: %d, %d"
                % (protein.residue[ir].conf[ic].confID,
                    len(self.multi_bnd_xyzrcp),
                    len(self.multi_bnd_atom),
                    )
                )
            
        logging.debug(
            "%s Multi-sidechain boundary record length should be equal: %d, %d"
            % (
                protein.residue[ir].conf[ic].confID,
                len(self.multi_bnd_xyzrcp),
                len(self.multi_bnd_atom),
            )
        )

        return

    def write_boundary(self, fname: str, kind: str):
        """This writes out both the xyzrcp and atom index files, for error
        checking and potentially as data exchange with PB wrapper.
        Args:
         - fname: Name for stem of file (preceding the '.xyzrcp' or '.atoms' extension)
         - kind (str): One of ["float", "single", "multi"]
        """
        valid_kinds = ["float", "single", "multi"]
        if kind not in valid_kinds:
            raise ValueError(f"'kind' must be one of {valid_kinds}; given: {kind}.")

        if kind == "float":
            bnd_xyzrcp = self.float_bnd_xyzrcp
            bnd_atom = self.float_bnd_atom
        elif kind == "single":
            bnd_xyzrcp = self.single_bnd_xyzrcp
            bnd_atom = self.single_bnd_atom
        else:
            bnd_xyzrcp = self.multi_bnd_xyzrcp
            bnd_atom = self.multi_bnd_atom

        # write xyzrcp
        xyzrcp_frmt = "{:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f}\n"
        with open(fname + ".xyzrcp", "w") as ofh:
            for xp in bnd_xyzrcp:
                ofh.write(xyzrcp_frmt.format(xp.x, xp.y, xp.z, xp.r, xp.c, xp.p))

        # write index file
        xyzrc_frmt = "{:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f}"
        with open(fname + ".atoms", "w") as ofh:
            for matched in bnd_atom:
                xyzrc = xyzrc_frmt.format(
                    matched[0].xyz[0],
                    matched[0].xyz[1],
                    matched[0].xyz[2],
                    matched[0].r_bound,
                    matched[0].charge,
                )
                ofh.write(f"{xyzrc} {' '.join([atom.atomID for atom in matched])}\n")

        return


class RunOptions:
    def __init__(self, args):
        self.inputpdb = "step2_out.pdb"  # implicit input pdb file
        self.start = args.c[0]
        self.end = args.c[1]
        self.d = float(args.d)
        self.s = args.s
        self.p = args.p
        self.t = args.t
        self.ftpl = args.ftpl
        self.salt = args.salt
        self.vdw = args.vdw
        self.fly = args.fly
        self.debug = args.debug
        self.refresh = args.refresh
        if args.l:  # load options from the specified file
            lines = open(args.l).readlines()
            for line in lines:
                line = line.split("#")[0].strip()
                # Extract an option line and strip off the spaces:
                fields = [x.strip() for x in line.split()]
                if len(fields) > 0:
                    key = fields[0]
                    if key == "-c":
                        self.start = int(fields[1])
                        self.end = int(fields[2])
                    elif key == "-d":
                        self.d = float(fields[1])
                    elif key == "-s":
                        self.s = fields[1]
                    elif key == "-p":
                        self.p = int(fields[1])
                    elif key == "-t":
                        self.t = fields[1]
                    elif key == "--vdw":
                        self.vdw = True
                    elif key == "--fly":
                        self.fly = True
                    elif key == "--refresh":
                        self.refresh = True
                    elif key == "--debug":
                        self.debug = True
                    elif key == "-ftpl":
                        self.ftpl = fields[1]

        if args.load_runprm:  # load additional run.prm file
            with open(args.load_runprm) as lines:
                for line in lines:
                    entry_str = line.strip().split("#")[0]
                    fields = entry_str.split()
                    if len(fields) > 1:
                        key_str = fields[-1]
                        if key_str[0] == "(" and key_str[-1] == ")":
                            key = key_str.strip("()").strip()
                            value = fields[0]
                            setattr(self, key, value)

        # unconditionally force to run step 3 only in step3.py
        setattr(self, "DO_ENERGY", "t")

        # convert attributes to runprm dict
        runprm = {}
        attributes = vars(self)
        for key, entry in attributes.items():
            if not key.startswith("_"):
                runprm[key] = entry
        record_runprm(runprm, "#STEP3")

    def toJSON(self):
        return json.dumps(self.__dict__, indent=4)
