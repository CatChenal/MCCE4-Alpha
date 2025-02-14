#!/usr/bin/env python

"""
This module defines the mcce internal data structures.

Classes:
  - Atom:
  - Conformer:
  - CONFORMER_param:
  - CONNECT_param:
  - LIGAND_ID_param:
  - Protein:
    Hierarchical structure of protein -> Residue -> Conformer -> Atom

  - RADIUS_param:
  - Residue:
  - Structure:
  - TORS:
  - TORSION_param:
  - TPL:
    Molecular topology parameters

  - _Biounit:
  - _Model:

"""
# ...................................................................................................| 100
from collections import defaultdict, Counter
import glob
import logging
import os
from pathlib import Path
from pprint import pformat
import re
import sys
from textwrap import wrap
from typing import Tuple, Union
from mcce4.constants import IONIZABLE_RES
from .geom import *

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class Atom:
    """This class defines atom properties and operations.

    Attributes:
        record (str) : Either ATOM of HETATM
        serial (int): Atom serial number. This number is not a part of atom identification.
        name (str): 4-char atom name.

        # FIX:
        # Conformer class says: "altLoc is now used to create initial conformers";
        # _Model class says: "the altLoc code used for spliting the model";
        # This Atom class says:
        altLoc (str): Alternative location marker. Not used in mccepdb format.

        resName (str): 3-char residue name.
        chainID (str): 1-char chain ID.
        resSeq (int): Residue sequence number.
        iCode (str): Insertion code of the residue. Use "_" in place of " " for mccepdb format.
        xyz (tuple): x, y, z corrdinates of the atom.
        element (str): 2-char element name, right adjusted.
        confNum (int): Conformer number.
        atomID (str): Unique atom ID = name + resName + resSeq + chainId + confNum
        confID (str): Unique conformer ID = confType + chainID + resSeq + iCode + confNum
        confType (str): 5-char conformer type as defined in CONFLIST from ftpl file.
        connectivity_param (CONNECT_param): Connectivity parameter.
        r_bound (float): radius for dielectric boundary in angstroms.
        charge (float): Atom charge.
        r_vdw (float): van der Waals radius in angstroms.
        e_vdw (float): van der Waals energy well depth in kcal/mol.
        connect12 [list]: connected atoms at distance level 1-2.
        connect13 [list]: connected atoms at distance level 1-3.
        connect14 [list]: connected atoms at distance level 1-4.
        history (str): conformer history string.
    """

    def __init__(self) -> None:
        """
        Declare variable types and default values
        """
        self.record = "ATOM"  # _atom_site.group_PDB in cif
        self.serial = 0  # from cif/pdb/pqr
        self.name = "  X "  # from cif/pdb/pqr
        self.altLoc = " "  # from cif/pdb/pqr
        self.resName = "UNK"  # from cif/pdb/pqr
        self.chainID = "A"  # from cif/pdb/pqr
        self.resSeq = 0  # from cif/pdb/pqr
        self.iCode = "_"  # from cif/pdb/pqr
        self.xyz = (0.0, 0.0, 0.0)  # from cif/pdb/pqr
        self.element = "  "  # element name converted from atom name
        # mcce internals:
        self.confNum = 0  # conformer number
        self.atomID = ""  # atom ID given residue, sequence, conformer and atom name.
        self.confID = ""  # conformer ID in context of residue, sequence, and conf name.
        self.confType = ""  # 5-char conftype as in CONFLIST from ftpl file
        # defined in ftpl file:
        self.resID = ""  # residue ID given residue name, chainId, sequence number and insertion code.
        self.connectivity_param = ""  # connectivity parameter defined in ftpl file
        self.r_bound = 0.0  # radius for dielectric boundary
        self.charge = 0.0  # atom charge
        self.r_vdw = 0.0  # vdw radius
        self.e_vdw = 0.0  # vdw energy well depth
        self.connect12 = []  # connected atoms at 1-2 level
        self.connect13 = []  # connected atoms at 1-3 level
        self.connect14 = []  # connected atoms at 1-4 level
        self.history = "_" * 10  # conformer history string

    def load_pdbline(self, line):
        """Load an atom from a pdb file ATOM/HETATM line.

        Args:
            line (str): An atom line in pdb file. It may start with "ATOM  " or "HETATM".
        """
        self.record = line[:6]
        self.serial = int(line[6:11])
        self.name = line[12:16]
        self.altLoc = line[16]
        self.resName = line[17:20]
        self.chainID = line[21]
        self.resSeq = int(line[22:26])
        self.iCode = line[26]
        self.xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        self.resID = (self.resName, self.chainID, self.resSeq, self.iCode)

    def loadline(self, line, tpl):
        self.record = line[:6]
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
        self.element = self.name[:2]  # element name defaults to the first two char
        if len(self.name.strip()) == 4 and self.name[0] == "H": # special case H, 4 char, H is the first char
            self.element = " H"

        self.history = line[80:].strip()

        self.compose_atomID()

        self.confID = "%5s%c%04d%c%03d" % (
            self.confType,
            self.chainID,
            self.resSeq,
            self.iCode,
            self.confNum,
        )
        self.resID = (self.resName, self.chainID, self.resSeq, self.iCode)

        # extended records
        connect_key = ("CONNECT", self.name, self.confType)
        self.connectivity_param = tpl.db[connect_key]

        radius_key = ("RADIUS", self.confType, self.name)
        if radius_key in tpl.db:
            radius_values = tpl.db[radius_key]
            self.r_vdw = radius_values.r_vdw
            self.e_vdw = radius_values.e_vdw
        else:  # Use default value as C
            print(
                "Warning, parameter %s not found, using default values as C"
                % str(radius_key)
            )
            self.r_vdw = 1.908
            self.e_vdw = 0.086
        return

    def compose_atomID(self):
        self.atomID = "%4s%3s%04d%c%03d" % (
            self.name,
            self.resName,
            self.resSeq,
            self.chainID,
            self.confNum,
        )

    def inherit(self, atom):
        """Inherit atom property except connectivity from a template

        Args:
            atom (Atom class): template atom 
        """
        self.serial = atom.serial
        self.name = atom.name
        self.altLoc = atom.altLoc
        self.resName = atom.resName 
        self.chainID = atom.chainID
        self.resSeq = atom.resSeq
        self.iCode = atom.iCode
        self.confNum = atom.confNum
        self.xyz = atom.xyz
        self.r_bound = atom.r_bound
        self.charge = atom.charge
        self.r_vdw = atom.r_vdw
        self.e_vdw = atom.e_vdw
        self.confType = atom.confType
        self.element = atom.element
        self.history = atom.history
        self.compose_atomID()
        self.confID = atom.confID
        self.resID = atom.resID
        self.connect12 = []
        self.connect13 = []
        self.connect13 = []


    def as_ATOM_line(self) -> str:
        """Use an ATOM instance data to build an undifferentiated pdb ATOM line.
        Note:
          The constant 'ATOM' string is used instead of ATOM.record.
        """
        return "ATOM %6d %4s %3s %c%4d%c   %8.3f%8.3f%8.3f\n" % (
            self.serial,
            self.name,
            self.resName,
            self.chainID,
            self.resSeq,
            self.iCode,
            self.xyz[0],
            self.xyz[1],
            self.xyz[2],
        )

    def as_mcce_ATOM_line(self) -> str:
        """Convert an ATOM instance into an undifferentiated mccepdb ATOM record.
        Note:
          The constant 'ATOM' string is used instead of ATOM.record.
        """
        serial = self.serial % 100000
        return "ATOM  %5d %4s %3s %c%04d%c%03d%8.3f%8.3f%8.3f%8.3f%12.3f      %s\n" % (
            serial,
            self.name,
            self.resName,
            self.chainID,
            self.resSeq,
            self.iCode,
            self.confNum,
            self.xyz[0],
            self.xyz[1],
            self.xyz[2],
            self.r_bound,
            self.charge,
            self.history,
        )


class Conformer:
    """
    Define mcce conformer properties
    """

    def __init__(self) -> None:
        self.confID = ""  # conformer ID given residue, sequence, and conformer name.
        self.confType = ""  # 5-char conformer type as in ftpl
        self.resSeq = 0  # sequence number
        self.iCode = "_"  # Insertion code
        self.confNum = 0 # conf number
        self.chainID = "A"  # default chainID
        self.altLoc = " "  # altLoc is now used to create initial conformers
        self.resID = ""  # residue ID given residue name, and sequence number.
        self.i = 0  # index number of conformer in the whole molecule.
        self.atom = []  # list of atoms it contains
        self.vdw0 = 0.0  # internal vdw
        self.vdw1 = 0.0  # vdw to backbone
        self.crg = 0.0  # net charge
        self.history = ""  # conformer history string
        self.calculated = (
            False  # flag if all energy terms were calculated; reported in head3.lst.
        )

    def clone(self):
        """Clone self and return a confomer with updated atom 12 connectivity. This is useful when creating a new conformer."""
        new_conf = Conformer()
        new_conf.confType = self.confType
        new_conf.chainID = self.chainID
        new_conf.resSeq = self.resSeq
        new_conf.iCode = self.iCode
        new_conf.altLoc = self.altLoc
        new_conf.resID = self.resID
        new_conf.history = self.history
        for atom in self.atom:
            new_atom = Atom()
            new_atom.inherit(atom)
            new_atom.connect12_inherited = atom.connect12
            new_conf.atom.append(new_atom)

        in_conf_atomnames = [a.name for a in self.atom]

        # create name to atom mapping for new_conf
        name2atom = {}
        for atom in new_conf.atom:
            name2atom[atom.name] = atom

        # recreate connect12 for in-conf connected atoms
        for atom in new_conf.atom:  
            for atom2 in atom.connect12_inherited:
                if atom2.name in in_conf_atomnames:  # in the same side chain conformer, update with the atom in the new_conf
                    atom.connect12.append(name2atom[atom2.name])
                else:
                    atom.connect12.append(atom2)
                    atom2.connect12.append(atom)
            atom.connect12_inherited = []
        return new_conf


    def init_by_atom(self, atom):
        """Initialize this conformer by a valid atom record"""
        self.confType = atom.confType
        self.chainID = atom.chainID
        self.altLoc = atom.altLoc
        self.resID = atom.resID
        self.resSeq = atom.resSeq
        self.iCode = atom.iCode
        self.i = 0
        self.atom = [atom]
        self.history = atom.history

    def make_blob(self):
        self.blob = Blob()
        

class Blob:
    """
    Define a sphere with center and radius to cover all atoms in a cobnformer including their vdw_r
    """
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
            for atom in conf.atom:
                d2 = ddvv(self.center, atom.xyz)
                if d2 > d_far2:
                    d_far2 = d2
            d_far = np.math.sqrt(d_far2)

            self.radius = d_far + r_max



class Residue:
    """
    Define residue properties
    """

    def __init__(self) -> None:
        self.resID = ("", "", 0, "")  # unique residue ID in context of residue name, chain ID, sequence number, and insertion code.
        self.conf = []  # list of conformers

    def serialize(self):
        i_conf = 0
        for conf in self.conf:
            conf.confNum = i_conf
            i_conf += 1
            for atom in conf.atom:
                atom.confNum = conf.confNum
            conf.confID = "%5s%c%04d%c%03d" % (conf.confType, conf.chainID, conf.resSeq, conf.iCode, conf.confNum)
        # serialize history string: heavy atom
        counter = {}
        for conf in self.conf:
            id = conf.history[:3]  # conf type + heavy atom making type
            if id in counter:
                counter[id] += 1
            else:
                counter[id] = 0  # count from 0

            conf.history = conf.history[:3] + "%03d"%(counter[conf.history[:3]]%1000) + conf.history[6:]

class Protein:
    """
    Define protein property
    """

    def __init__(self) -> None:
        self.residue = []

    def serialize(self):
        """Set serial number to atoms and conformers"""
        i_atom = 0
        # FIX: unused
        # found_atoms = set()  # some atoms are shared by conformer,
        for res in self.residue:
            res.serialize()
            for conf in res.conf:
                for atom in conf.atom:
                    atom.confNum = conf.confNum
                    atom.serial = i_atom
                    i_atom += 1
                    atom.history = conf.history

    def dump(self, fname, prepend=[], append=[]) -> None:
        """Dump internal structure to a named file."""
        # tidy the serial numbers for atoms and conformers
        self.serialize()

        lines = prepend
        for res in self.residue:
            lines.append(
                "# Residue: %s %c%04d%c\n"
                % (res.resID[0], res.resID[1], int(res.resID[2]), res.resID[3])
            )
            for conf in res.conf:
                lines.append("## Conformer %s: %s %s\n" % (conf.confID, conf.confType, conf.history))
                for atom in conf.atom:
                    lines.append(atom.as_mcce_ATOM_line())
                lines.append("#%s\n" % ("-" * 89))
            lines.append("#%s\n" % ("=" * 89))

        lines = lines + append
        if fname:
            open(fname, "w").writelines(lines)
        else:
            sys.stdout.writelines(lines)


class CONNECT_param:
    """CONNECT record
    Defines ftpl value format."""

    def __init__(self, value_str):
        fields = value_str.split(",")
        self.orbital = fields[0].strip()
        self.connected = [x.strip().strip('"') for x in fields[1:]]


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
        """
        self.v2 = 0.0
        self.n_fold = 0
        self.gamma = 0.0


class TORSION_param:
    """This class defines the data structure of TORSION parameters."""

    def __init__(self, value_str) -> None:
        """Construct the TORSION parameter and assign initial values from value_str."""
        d2r = 0.017453  # degrees to radians conversion factor

        fields = value_str.strip().split(",")
        self.atom1 = fields[0].strip().strip('"')
        self.atom2 = fields[1].strip().strip('"')
        self.atom3 = fields[2].strip().strip('"')
        self.tors_terms = []
        term_str = ", ".join(fields[3:])
        terms = re.findall(r"\(.*?\)", term_str)
        for term in terms:
            tors = TORS()
            tors_value_strings = term.strip("()").split(",")
            tors.v2 = float(tors_value_strings[0])
            tors.n_fold = int(tors_value_strings[1])
            tors.gamma = float(tors_value_strings[2]) * d2r
            self.tors_terms.append(tors)


class LIGAND_ID_param:
    """This class defines the data structure of LIGAND_ID parameters.
    Examples:
    LIGAND_ID, CYS, CYS: " SG " - " SG "; 2.00 +- 0.20; CYL, CYL
    LIGAND_ID, HIS, HEM: " NE2" - "FE  "; 2.50 +- 0.25; HIL, HEM
    LIGAND_ID, HIS, HEA: " NE2" - "FE  "; 2.50 +- 0.25; HIL, HEA
    LIGAND_ID, HIS, HEB: " NE2" - "FE  "; 2.50 +- 0.25; HIL, HEB
    LIGAND_ID, HIS, HEC: " NE2" - "FE  "; 2.50 +- 0.25; HIL, HEC
    LIGAND_ID, MET, HEA: " SD " - "FE  "; 2.50 +- 0.25; HIL, HEA
    LIGAND_ID, MET, HEB: " SD " - "FE  "; 2.50 +- 0.25; HIL, HEB
    LIGAND_ID, MET, HEC: " SD " - "FE  "; 2.50 +- 0.25; HIL, HEC
    """

    def __init__(self, value_str) -> None:
        """Construct the LIGAND_ID parameter and assign initial values from value_str."""
        fields = value_str.strip().split(";")
        # atom pair
        atom1, atom2 = fields[0].strip().split("-")
        self.atom1 = atom1.strip().strip('"')
        self.atom2 = atom2.strip().strip('"')

        # ligand bond distance
        distance, tolerance = fields[1].strip().split("+-")
        self.distance = float(distance)
        self.tolerance = float(tolerance)

        # residue name to
        name1, name2 = fields[2].strip().split(",")
        self.res1_name = name1.strip()
        self.res2_name = name2.strip()


class TPL:
    def __init__(self):
        """Class constructor. Initialize internal database named "db"."""
        self.db = {}

    def read_ftpl_folder(self, tplfolder: str):
        """Read ftpl records from ftpl files under a folder and store in the internal dictionary.

        Parameters
        ----------
        tplfolder: str
        Absolute path of the ftpl file folder.

        See Also
        --------
        mcce4.pdbio.TPL.read_ftpl_file : function to read a single ftpl file

        Notes
        -----
        The ftpl entries are loaded into the object internal database db, named "db",
        which is a dictionnary.
        The key is a tuple composed of up to three key elements, and the value
        can be in various types, int, str, float, or predefined data structure.

        Examples
        --------
        >>> from mcce4.pdbio import TPL
        >>> tpl = TPL()
        >>> ftpl_folder = "/home/jmao/projects/MCCE4/param_PARSE"
        >>> tpl.read_ftpl_folder(ftpl_folder)
        >>> tpl.printme()
        """
        updated_tplfolder = update_dist_folder(tplfolder)
        files = glob.glob("%s/*.ftpl" % updated_tplfolder)
        files.sort()
        logging.info(
            'Reading parameters from *.ftpl files under folder "%s"' % updated_tplfolder
        )
        if not files:
            logging.error('No ftpl file in folder "%s". Quiting' % updated_tplfolder)
            sys.exit(1)
        else:
            for fname in files:
                self.read_ftpl_file(fname)

    def read_ftpl_file(self, fname):
        """Read ftpl records from a ftpl file and store in the internal dictionary.
        Parameters
        ----------
        fname: str
        Absolute path of the ftpl file intended to be loaded.

        See Also
        --------
        mcce4.pdbio.TPL.read_ftpl_folder : function to read all ftpl files in a folder

        Notes
        -----
        The ftpl entries are loaded into the object internal database db, named "db",
        which is a dictionnary.
        The key is a tuple composed of up to three key elements, and the value
        can be in various types, int, str, float, or predefined data structure.

        Examples
        --------
        >>> from mcce4.pdbio import TPL
        >>> tpl = TPL()
        >>> ftpl_file = "/home/jmao/projects/MCCE4/param_PARSE/asp.ftpl"
        >>> tpl.read_ftpl_file()
        >>> tpl.printme()
        """
        lines = open(fname).readlines()
        for line in lines:
            end = line.find("#")
            line = line[:end]
            fields = line.split(":")
            if len(fields) != 2:
                continue

            key_string = fields[0].strip()
            keys = key_string.split(",")
            key1 = keys[0].strip().strip('"')
            if len(keys) > 1:
                key2 = keys[1].strip().strip('"')
            else:
                key2 = ""
            if len(keys) > 2:
                # No add'l strip() after stripping quotes to perserve atom names with spaces
                key3 = keys[2].strip().strip('"')
            else:
                key3 = ""

            value_string = fields[1].strip()
            warn_msg_overwrite_key = "This key {} was already loaded, now overwriting"

            if key1 == "CONFLIST":
                key = (key1, key2)
                if key in self.db:
                    logging.warning(warn_msg_overwrite_key.format(key))
                self.db[key] = [x.strip() for x in value_string.strip().split(",")]
            elif key1 == "CONNECT":
                key = (key1, key2, key3)
                if key in self.db:
                    logging.warning(warn_msg_overwrite_key.format(key))
                self.db[key] = CONNECT_param(value_string)
            elif key1 == "RADIUS":
                key = (key1, key2, key3)
                if key in self.db:
                    logging.warning(warn_msg_overwrite_key.format(key))
                self.db[key] = RADIUS_param(value_string)
            elif key1 == "CONFORMER":
                key = (key1, key2)
                if key in self.db:
                    logging.warning(warn_msg_overwrite_key.format(key))
                self.db[key] = CONFORMER_param(value_string)
            elif key1 == "CHARGE":
                key = (key1, key2, key3)
                if key in self.db:
                    logging.warning(warn_msg_overwrite_key.format(key))
                self.db[key] = float(value_string)
            elif key1 == "ROTATE":
                fields = value_string.split(",")
                key = (key1, key2)
                if key in self.db:
                    logging.warning(warn_msg_overwrite_key.format(key))
                self.db[(key1, key2)] = []
                for atom_pair in fields:
                    splited_atoms = atom_pair.split("-")
                    if len(splited_atoms) == 2:
                        atom1 = splited_atoms[0].strip().strip('"')
                        atom2 = splited_atoms[1].strip().strip('"')
                        self.db[(key1, key2)].append((atom1, atom2))
            elif key1 == "ROT_SWAP":
                fields = value_string.split(",")
                key = (key1, key2)
                if key in self.db:
                    logging.warning(warn_msg_overwrite_key.format(key))
                self.db[(key1, key2)] = []
                for atom_pair in fields:
                    splited_atoms = atom_pair.split("-")
                    if len(splited_atoms) == 2:
                        atom1 = splited_atoms[0].strip().strip('"')
                        atom2 = splited_atoms[1].strip().strip('"')
                        self.db[(key1, key2)].append((atom1, atom2))
            elif key1 == "TORSION":
                key = (key1, key2, key3)
                if key in self.db:
                    logging.warning(warn_msg_overwrite_key.format(key))
                self.db[key] = TORSION_param(value_string)
            elif key1 == "LIGAND_ID":
                key = (key1, key2, key3)
                if key in self.db:
                    logging.warning(warn_msg_overwrite_key.format(key))
                self.db[key] = LIGAND_ID_param(value_string)
            elif key1 == "EXTRA" or key1 == "SCALING":
                key = (key1, key2)
                if key in self.db:
                    logging.warning(warn_msg_overwrite_key.format(key))
                self.db[key] = float(value_string)
            else:
                logging.warning(
                    "This entry (%s, %s, %s) in file %s was unrecognized therefore not processed."
                    % (key1, key2, key3, fname)
                )

    def tpl_lines(self):
        """Function to print the entries in db"""
        lines = []
        # FIX: key_str, value_str : not used
        for keys, value in self.db.items():
            key1 = keys[0]
            if len(keys) > 1:
                key2 = keys[1]
            if len(keys) > 2:
                key3 = keys[2]

            if key1 == "CONFLIST":
                key_str = ", ".join([key1, key2])
                value_str = ", ".join(value)
            elif key1 == "CONNECT":
                key_str = ", ".join([key1, '"%s"' % key2, key3])
                values = [value.orbital]
                for connected in value.connected:
                    values.append('"%s"' % connected)
                value_str = ", ".join(values)
            elif key1 == "RADIUS":
                key_str = ", ".join([key1, '"%s"' % key2, key3])
                value_str = "%8.3f, %8.3f, %8.3f" % (
                    value.r_bound,
                    value.r_vdw,
                    value.e_vdw,
                )
            elif key1 == "CONFORMER":
                key_str = ", ".join([key1, key2])
                values = []
                for p_key, p_value in value.param.items():
                    values.append("%s=%6.3f" % (p_key, p_value))
                value_str = ", ".join(values)
            elif key1 == "CHARGE":
                key_str = ", ".join([key1, '"%s"' % key2, key3])
                value_str = "%8.3f" % (value)
            elif key1 == "ROTATE":
                key_str = ", ".join([key1, key2])
                values = []
                for pair in value:
                    values.append('"%s" - "%s"' % pair)
                value_str = ", ".join(values)
            elif key1 == "ROT_SWAP":
                key_str = ", ".join([key1, key2])
                values = []
                for pair in value:
                    values.append('"%s" - "%s"' % pair)
                value_str = ", ".join(values)
            elif key1 == "TORSION":
                key_str = ", ".join([key1, key2, '"%s"' % key3])
                values = [
                    '"%s"' % value.atom1,
                    '"%s"' % value.atom2,
                    '"%s"' % value.atom3,
                ]
                for term in value.tors_terms:
                    values.append(
                        "(%.3f, %d, %.2f)" % (term.v2, term.n_fold, term.gamma)
                    )
                value_str = ", ".join(values)
            elif key1 == "LIGAND_ID":
                key_str = ", ".join([key1, key2, key3])
                value_str = '"%s" - "%s"; %.2f +- %.2f; %s, %s' % (
                    value.atom1,
                    value.atom2,
                    value.distance,
                    value.tolerance,
                    value.res1_name,
                    value.res2_name,
                )
            elif key1 == "EXTRA" or key1 == "SCALING":
                key_str = ", ".join([key1, key2])
                value_str = "%6.2f" % value
            else:
                print(
                    "Warning: Key %s does not have a pretty print instruction" % (key1)
                )
                continue  # This is necessary for skipping the print statement

            lines.append("%s: %s\n" % (key_str, value_str))
        return lines

    def print_me(self):
        lines = self.tpl_lines()
        sys.stdout.writelines(lines)

    def dump(self, fname="tpl.record"):
        lines = self.tpl_lines()
        with open(fname, "w") as fo:
            fo.writelines(lines)


class _Biounit:
    """
    Define biounit properties
    """

    def __init__(self) -> None:
        self.serial = ""  # biounit serial number
        self.author_determined = None  # biounit determination by author
        self.software_determined = None  # biounit determination by software
        self.chains = []  # biounit chains


class _Model:
    """
    Define Model
    """

    def __init__(self) -> None:
        self.serial = ""  # model serial number, 1 if single model without MODEL line
        self.altLoc = " "  # the altLoc code used for spliting the model.
        self.lines = []  # atom lines

    def _line2resid(self, line:str, keep_type:bool = False) -> Union[tuple, None]:
        if line[:6] not in ["ATOM  ", "HETATM"]:
            return None
        resName = line[17:20]
        chainID = line[21]
        seqNum = int(line[22:26])
        iCode = line[26]
        if not keep_type:
            return resName, chainID, seqNum, iCode
        else:
            return line[:6], resName, chainID, seqNum, iCode

    def _mk_newftpl(self, lines):
        elebd_radius = {
            " N": 1.5,
            " H": 1.0,
            " C": 1.7,
            " O": 1.4,
            " P": 1.85,
            " S": 1.85,
            " X": 1.85,
        }

        vdw_parm = {
            " C": (2.000, 0.150),
            " H": (1.000, 0.020),
            " O": (1.600, 0.200),
            " N": (1.750, 0.160),
            " S": (2.000, 0.200),
            " P": (2.000, 0.200),
            " X": (2.000, 0.173),
        }

        newftpl = []
        # CONFLIST
        newftpl.append("# Conformer definition\n")
        resName = lines[0][17:20]
        new_line = "CONFLIST, %s: %sBK\n" % (resName, resName)
        if new_line not in newftpl:
            newftpl.append(new_line)

        atom_names = [line[12:16] for line in lines]
        # ATOM name and bonds
        newftpl.append("\n# ATOM name and bonds\n")
        for atom_name in atom_names:
            new_line = 'CONNECT, "%s", %sBK: ion\n' % (atom_name, resName)
            if new_line not in newftpl:
                newftpl.append(new_line)

        # ATOM charges
        newftpl.append("\n# ATOM charges\n")
        for atom_name in atom_names:
            new_line = 'CHARGE, %sBK, "%s": 0.0\n' % (resName, atom_name)
            if new_line not in newftpl:
                newftpl.append(new_line)

        # ATOM radius
        newftpl.append(
            "\n# ATOM radius: dielelctric boundary radius, VDW radius, and energy well depth\n"
        )
        for atom_name in atom_names:
            element_name = atom_name[:2]
            if len(atom_name.strip()) == 4 and atom_name[0] == "H":
                element_name = " H"
            if element_name in elebd_radius:
                r_b = elebd_radius[element_name]
            else:
                r_b = elebd_radius[" X"]
            if element_name in vdw_parm:
                r_v, e_v = vdw_parm[element_name]
            else:
                r_v, e_v = vdw_parm[" X"]

            new_line = 'RADIUS, %sBK, "%s": %.3f, %.3f, %.3f\n' % (
                resName,
                atom_name,
                r_b,
                r_v,
                e_v,
            )
            if new_line not in newftpl:
                newftpl.append(new_line)

        return newftpl

    def detect_unknowns(self, tpl):
        # idetify unknown residues
        unknown_res_names = set()
        # temp fix, atomlines contaminated
        # cleaned_lines = [line for line in self.lines if line[:6]=="ATOM  " or line[:6]=="HETATM"]

        for line in self.lines:
            resName = line[17:20]
            key = ("CONFLIST", resName)
            if key not in tpl.db:
                unknown_res_names.add(resName)

        # collect the unknown cofactor with most number of atoms and convert to new.ftpl
        unknown_res_names = list(unknown_res_names)
        ftpl_lines = []
        for unknown_res_name in unknown_res_names:
            logging.info(
                '   Detected cofactor not defined in mcce: "%s"' % unknown_res_name
            )
            this_res = {}
            for atom_line in self.lines[1:]:
                resID = self._line2resid(atom_line)
                if resID is None:
                    continue
                if resID[0] == unknown_res_name:
                    if resID in this_res:
                        this_res[resID].append(atom_line)
                    else:
                        this_res[resID] = [atom_line]

            # get the resID with the max number of atoms
            n_max = 0
            # resID_max = ""
            # FIX: local variable 'resID_max' is assigned to but never used
            #
            for resID, atom_lines in this_res.items():
                if len(atom_lines) > n_max:
                    # resID_max = resID
                    n_max = len(atom_lines)

            # get atom lines of one residue
            atom_lines = this_res[resID]
            ftpl_lines += self._mk_newftpl(atom_lines)

        return ftpl_lines


class Structure:
    """
    Organize structure file records into structural data:
    Molecule
     |---Global properties such as resolution, biounits, and method
     |
     |---Model m
     |    |---ATOM i
     |    |---...
     |---...
    Notes:
      Model is an MCCE Model (in class _Model), meaning that ATOM & HETATM lines
      are converted to residues IF a matching topology is found.
      For this reason, when the Model's residues are split between 'ionizable' and 'other',
      e.g. via `Structure.get_model1_res_dict`, the 'other' group will include species from
      HETATM lines MINUS those for which a new.tpl file was created, while the 'ionizable'
      group will always be correct (in the context of a H+ titration).
    """

    def __init__(self, load_all_models: bool = False) -> None:
        self.pdbid = ""
        self.date = ""  # release date
        self.function = ""  # molecule function
        self.title = ""
        self.molecule = ""  # molecule name
        self.method = ""  # structure determination method, NMR, XRD, etc
        self.resolution = ""  # self reported structure resolution
        self.splits = []  # lists the other pdbids that make up a macroassembly
        self.load_all_models = load_all_models
        self.n_models = 1  # updated if multiple MODEL records found
        self.hetero_names = defaultdict(tuple)  # names of hetero speciess
        self.tot_waters = 0
        self.tot_cofactors = 0
        self.ter_res = []  # terminal res -> CTR in step1
        self.models = []  # list of lists of ATOM, HETATM lines for each model
        self.chains = []
        self.units_per_chain = []  # units per chain grouped by type (AAs or NTs)
        self.missing = None  # dict to store missing re, res atoms, ligand atoms
        self.biounits = []  # biounits, if no record of biounits, all placed
        self.ssbonds = None  # disulfide bonds; list of 3-tuples: CYS1, dist, CYS2
        self.links = (
            None  # hetatm/res connect.; list of 3-tuples: entity1, dist, entity2
        )
        # if sites found, self.sites :: dict with key=siteid, values=list with
        # site description as 1st item and other items :: each of the residues involved.
        self.sites = None  # sites identified by authors or sw + info from REMARK 800
        self.pdb_fp = None  # path of pdb passed to load_pdb()
        # DO_NOT_USE flag: set to True when pdb headers are processed (if present) &
        #                  if CAVEAT or Obsolete statement found.
        self.DO_NOT_USE = False
        self.has_header = True  # no info to parse if headers removed.

    def rename_lines(self, raw_lines, rename=""):
        # Get the renaming rules:
        renaming_rules = []  # a list of tuples in the form of (str_from, str_to)
        if rename:
            lines = open(rename).readlines()
            for line in lines:
                entry = line.split("#")[0]
                if len(entry) >= 30:
                    str_from = line[:14]
                    str_to = line[16:30]
                    renaming_rules.append((str_from, str_to))

        new_lines = []
        for line in raw_lines:
            new_line = line
            if line[:6] == "ATOM  " or line[:6] == "HETATM":
                to_be_matched = line[12:26]
                for rule in renaming_rules:
                    matched = True
                    str_from = rule[0]
                    str_to = rule[1]
                    for i in range(14):
                        if str_from[i] != "*" and str_from[i] != to_be_matched[i]:
                            matched = False
                            break
                    if matched:
                        # print(rule, line)
                        replaced = []
                        for i in range(14):
                            if str_to[i] == "*":
                                replaced.append(to_be_matched[i])
                            else:
                                replaced.append(str_to[i])
                        to_be_matched = "".join(replaced)

                new_line = line[:12] + to_be_matched + line[26:]
            # compose new line after all rules are processed
            new_lines.append(new_line)

        return new_lines

    def process_atomlines(self, atomlines: list):
        # collect models
        detected_model = False
        for line in atomlines:
            #print(line)
            if line.startswith("MODEL "):
                detected_model = True
                model = _Model()
                model.serial = line.split()[-1].strip()
                continue
            if detected_model:
                if line[:6] == "TER   ":
                    continue  # ignore this separator for ATOM and HETATM - jmao
                elif line[:6] in ["ATOM  ", "HETATM"]:
                    #print("---",line)
                    model.lines.append(line)
                else:  # encountered ENDMDL or END
                    self.models.append(model)
                    detected_model = False
        if detected_model:  # model ended with ATOM or HETATM, push this last model
            self.models.append(model)
            detected_model = False
        # some files, e.g. 1b16.pdb are missing a "MODEL" line:
        if not self.models:
            self.n_models = 0
            logger.warning(f"{self.pdb_fp!s} is missing a MODEL line.")
            return

        # set total model count
        self.n_models = len(self.models)

        # Use the first model to list unique polymeric chains
        self.chains = list(
            set([line[21] for line in self.models[0].lines if line.startswith("ATOM")])
        )
        self.chains.sort()
        self.ter_res = [
            line.split(maxsplit=2)[-1].strip()
            for line in self.models[0].lines
            if line.startswith("TER")
        ]
        return
    
    def process_headers(self, nonatomlines: list):
        """Process a list of pdb lines other than the coordinate lines.
        """
        nogo_kws = ["CAVEAT", "REMARK 5,"]
        for line in nonatomlines:
            if line.startswith(nogo_kws[0]) or line.startswith(nogo_kws[1]):
                self.DO_NOT_USE = True
                break

        if self.DO_NOT_USE:
            # get the title, pdb name or id & exit:
            keyword = "TITLE"
            for line in nonatomlines:
                if line.startswith(keyword):
                    contid = line[8:10].strip()
                    if not contid.isnumeric():
                        self.title += line[len(keyword):].strip() + " "
                    else:
                        self.title += line[10:].strip() + "; "
            if self.title:
                self.title = self.title.replace("  ", " ")
                self.title = self.title.strip()

            self.pdbid = self.pdb_fp.stem.upper()

            return

        # removed: processing of NUMMDL record:
        # not always present even if pdb has several models, e.g. 3kch

        keyword = "EXPDTA   "  # structure determination method
        for line in nonatomlines:
            if line.startswith(keyword):
                self.method = line[len(keyword):].strip()
                break

        keyword = "REMARK   2 RESOLUTION."
        for line in nonatomlines:
            if line.startswith(keyword):
                self.resolution = line[len(keyword):].strip()
                break

        keyword = "COMPND   "  # molecule name(s)
        molecules = []
        for line in nonatomlines:
            if line.startswith(keyword):
                _, data = line[len(keyword):].split(maxsplit=1)
                if data.startswith("MOLECULE") or data.startswith("OTHER_DETAILS"):
                    molecules.append(data.split(":")[-1].strip(" ;\n"))
        self.molecule = ", ".join(molecules)

        keyword = "HEADER    "
        for line in nonatomlines:
            if line.startswith(keyword):
                hdr = line[len(keyword):].rsplit(maxsplit=2)
                self.function, self.date, self.pdbid = hdr
                break

        keyword = "TITLE"
        for line in nonatomlines:
            if line.startswith(keyword):
                contid = line[8:10].strip()
                if not contid.isnumeric():
                    self.title += line[len(keyword):].strip() + " "
                else:
                    self.title += line[10:].strip() + "; "
        if self.title:
            self.title = self.title.replace("  ", " ")
            self.title = self.title.strip()

        keyword = "SPLIT "
        self.splits = [
            line[len(keyword):].strip()
            for line in nonatomlines
            if line.startswith(keyword)
        ]

        keyword = "SSBOND"
        # output 3-tuple: CYS1, dist, CYS2. No SymOP.
        # some old files may not have the distance
        self.ssbonds = [
            (
                line[11:25].strip().replace("  ", " "),
                line[73:].strip() or "?",
                line[25:59].strip().replace("  ", " "),
            )
            for line in nonatomlines
            if line.startswith(keyword)
        ]

        keyword = "LINK        "
        # output 3-tuple: entity1, dist, entity2. No SymOP.
        self.links = [
            (
                line[12:28].strip().replace("  ", " "),
                line[74:].strip() or "?",
                line[42:57].strip().replace("  ", " "),
            )
            for line in nonatomlines
            if line.startswith(keyword)
        ]

        # keyword = "REMARK 800 "
        kws = ["REMARK 800 SITE_IDENTIFIER", "REMARK 800 SITE_DESCRIPTION"]
        REM800 = []
        for i, line in enumerate(nonatomlines):
            if line.startswith(kws[0]) or line.startswith(kws[1]):
                REM800.append(line.split(":")[1].strip())
        # check:
        misformed_r800 = len(REM800) % 2 != 0
        if misformed_r800:
            print(
                "WARNING: length of REM800 list not even: should have paired site id & description."
            )

        if REM800 and not misformed_r800:
            # prep dict with rem 800 desc:
            sites_d = defaultdict(list)
            seen = set()
            for j, rem in enumerate(REM800):
                if j % 2 == 0:
                    if rem not in seen:
                        sid = rem
                        seen.add(sid)
                else:
                    sites_d[sid].append(rem)

            keyword = "SITE  "
            for line in nonatomlines:
                if line.startswith(keyword):
                    sid = line[11:14]
                    sites_d[sid].extend(wrap(line[18:].strip(), width=9))
            # assign regular dict to sites attribute:
            self.sites = dict(sites_d)

        keyword = "SEQRES"
        chain_res = defaultdict(int)
        chain_oxy = defaultdict(int)
        chain_doxy = defaultdict(int)
        for line in nonatomlines:
            if line.startswith(keyword):
                _, _, chainID, numRes, res1, *_ = line.split()
                len1 = len(res1)
                if len1 == 3:
                    chain_res[chainID] = numRes
                elif len1 == 2:
                    chain_doxy[chainID] = numRes
                else:
                    chain_oxy[chainID] = numRes
        if chain_res:
            tot = f"; Total: {sum(int(n) for n in chain_res.values())}"
            self.units_per_chain.append(
                "Residues: "
                + ", ".join(v for v in [f"{r}:{n}" for r, n in chain_res.items()])
                + tot
            )
        if chain_doxy:
            tot = f"; Total: {sum(int(n) for n in chain_doxy.values())}"
            self.units_per_chain.append(
                "dRNTs: "
                + ", ".join(v for v in [f"{d}:{n}" for d, n in chain_doxy.items()])
                + tot
            )
        if chain_oxy:
            tot = f"; Total: {sum(int(n) for n in chain_oxy.values())}"
            self.units_per_chain.append(
                "RNTs: "
                + ", ".join(v for v in [f"{o}:{n}" for o, n in chain_oxy.items()])
                + tot
            )

        # combine het names and their synonyms
        keywords = ["HETNAM", "HETSYN"]
        for line in nonatomlines:
            if line.startswith(keywords[0]) or line.startswith(keywords[1]):
                hetid, hetname = line[11:14], line[15:].strip()
                self.hetero_names[hetid] = (hetname, )

        keyword = "FORMUL"
        cofactors_counts = defaultdict(int)
        for line in nonatomlines:
            if line.startswith(keyword):
                try:
                    hetid, hetcount = line[12:15], int(line[19:].strip().split("(")[0])
                except ValueError:
                    # likely: single molecule formula: no N() grouping
                    hetid, hetcount = line[12:15], 1

                if line[18].strip() == "*":
                    self.tot_waters = hetcount
                else:
                    cofactors_counts[hetid] = hetcount
                    hval = self.hetero_names[hetid]
                    self.hetero_names[hetid] = (hval[0], hetcount)

        self.tot_cofactors = sum(cofactors_counts.values())

        # missing species ---------------------------------------------
        # res: rem 465; res atoms: rem 470; non-polymer atoms: rem 610:
        miss_d = dict()
        keyword = "REMARK 465"
        per_d = defaultdict(list)
        for line in nonatomlines:
            if line.startswith(keyword):
                if line.strip() == keyword:
                    continue
                if line[11:15] in ["MISS", "THE ", "EXPE", "IDEN", "SSSE"]:
                    continue
                if line[13] == "M" or "SSSEQI" in line:
                    continue
                _, _, res, chainID, seq = line.split()
                per_d[chainID].append((res, int(seq)))
        if per_d:
            for c in per_d:
                tot = len(per_d[c])
                per_d[c].append(("Count", tot))
            miss_d["Residues"] = dict(per_d)

        keyword = "REMARK 470"
        per_d = defaultdict(list)
        for line in nonatomlines:
            if line.startswith(keyword):
                if line.strip() == keyword:
                    continue
                if line[11:15] in ["MISS", "THE ", "RES=", "I=IN", "C=CH"]:
                    continue
                if "SSEQI" in line:
                    continue
                if line[13:19] == "MODELS":
                    continue
                if line.endswith("ATOMS"):
                    continue
                _, _, res, chainID, seq, *atms = line.split()
                per_d[chainID].append(f"{res}_{seq}=(" + ",".join(a for a in atms) + ")")
            if per_d:
                miss_d["ResAtoms"] = dict(per_d)

        keyword = "REMARK 610"
        per_d = defaultdict(list)
        for line in nonatomlines:
            if line.startswith(keyword):
                if line.strip() == keyword:
                    continue
                if line[11:15] in ["MISS", "THE ", "RES=", "I=IN"]:
                    continue
                if "SSEQI" in line:
                    continue
                _, _, res, *info = line.split()
                if len(info) == 1:
                    per_d["free"].append(f"{res} {info[0]}")
                else:
                    per_d[info[0]].append(f"{res} {info[1]}")
        if per_d:
            miss_d["LigandAtoms"] = [
                f"{k}: " + ", ".join([r for r in per_d[k]]) for k in per_d
            ]
        if miss_d:
            self.missing = miss_d
        # end missing species -----------------------------------------

        # load biounits
        keyword = "REMARK 350"
        biounit_lines = [
            line[len(keyword):].strip()
            for line in nonatomlines
            if line.startswith(keyword)
        ]
        detected_biounit = False
        for line in biounit_lines:
            if line.startswith("BIOMOLECULE:"):
                detected_biounit = True
                biounit = _Biounit()
                _, biounit.serial = line.split(":")
            if detected_biounit and "DETERMINED" in line:
                if "AUTHOR" in line:
                    _, biounit.author_determined = line.split(":")
                else:
                    _, biounit.software_determined = line.split(":")
            if detected_biounit and "TO CHAINS:" in line:
                _, chns = line.split(":")
                biounit.chains = chns
            if detected_biounit and "AND CHAINS:" in line:
                _, chns = line.split(":")
                biounit.chains += chns
                # done:
                self.biounits.append(biounit)
                detected_biounit = False

        if not self.pdbid or self.pdbid == "XXXX":
            self.pdbid = self.pdb_fp.stem.upper()

        return

    def load_pdb(self, fname: str, rename: str = "") -> None:
        """
        Load structure file in pdb format, all atoms stored in conformer 0.
        Args:
          - fname (str): pdb filepath
          - rename (str): rename rules filepath
        """
        self.pdb_fp = Path(fname)
        rawlines = self.pdb_fp.read_text().splitlines()

        if rename:
            rename_rule_file = update_dist_folder(rename)
            logging.info(f"Renamed atom lines following rules in {rename_rule_file!r}")
            renamed_lines = self.rename_lines(rawlines, rename=rename_rule_file)
        else:
            renamed_lines = rawlines

        # Divide raw lines into atom lines and non-atom lines for more efficient processing.
        nonatomlines = []
        atomlines = []
        for line in renamed_lines:
            if line[:6] in ["ATOM  ", "HETATM", "MODEL ", "TER   ", "ENDMDL", "END   "]:
                atomlines.append(line)
            else:
                nonatomlines.append(line)

        # process headers first as self.process_atomlines needs self.n_models
        if nonatomlines:
            self.process_headers(nonatomlines)
        else:
            self.has_header = False
            self.pdbid = self.pdb_fp.stem.upper()

        if not self.DO_NOT_USE:
            # ok to get the model(s):
            self.process_atomlines(atomlines)

        return

    def get_model1_heteros(self) -> Union[dict, None]:
        """Given a list of pdb coordinates lines from the first model, parse the hetero atoms
        to return a dictionary holding the count of free cofactors and waters, per chain (key).
        """
        if not self.n_models:
            return None
        # only use model 1:
        mdl1 = self.models[0]
        if not mdl1.lines:
            logger.error("Something went wrong: no model lines! Didn't use Structure.load_pdb(), perhaps?")
            return None

        # get unique hetero species per chain:
        hetero_set = defaultdict(set)
        for line in mdl1.lines:
            x = mdl1._line2resid(line, keep_type=True)
            # populate dict with chain as key:
            # x0: [ATOM | HETATM], x1: resName, x2: chainID, x3: seqNum, x4: iCode
            if x[0] == "HETATM":
                hetero_set[x[2]].add((x[1], x[3]))
        
        # get species count per chain:
        tots_per_chain = defaultdict(dict)
        for c in hetero_set:
            cntr = defaultdict(int)
            for val in hetero_set[c]:
                cntr[val[0]] += 1
            tots_per_chain[c] = dict(cntr)
            
        return dict(tots_per_chain)

    def get_model1_res(self) -> Union[Tuple[dict, int], None]:
        """
        Return a dict for the name & count of residues in the first model , per chain (key).
        """
        if not self.n_models:
            return None
        # only use model 1:
        mdl1 = self.models[0]
        if not mdl1.lines:
            logger.error("Something went wrong: no model lines! Didn't use Structure.load_pdb(), perhaps?")
            return None

        mdl_data = defaultdict(set)
        for line in mdl1.lines:
            x = mdl1._line2resid(line, keep_type=True)
            if x is None:
                continue
            # populate dict with chain as key:
            # x: [ATOM | HETATM], resName, chainID, seqNum, iCode
            mdl_data[x[2]].add((x[0], x[1], x[3]))

        for chn in mdl_data:
            toti, toto = 0, 0
            res_data = defaultdict(list)
            
            for val in mdl_data[chn]:
                if val[0] == "ATOM  ":
                    res_data[chn].append(val[1])

            res_count_per_chain = Counter(res_data[chn])
            tot_chn_res = sum(res_count_per_chain.values())

            residues = []
            for k in res_count_per_chain:
                cnt = res_count_per_chain[k]
                residues.append(f"{k}: {cnt}")
                if k in IONIZABLE_RES:
                    toti += cnt
                else:
                    toto += cnt

            chn_dict = {"RESIDUES": (", ".join(r for r in sorted(residues)),
                                     f"Total: {tot_chn_res}",
                                     f"Ionizable: {toti}",
                                     f"Ratio: {toti/tot_chn_res:.1%}",
                                     ),
                        }

        model_dict = defaultdict(dict)
        model_dict[chn] = chn_dict

        return dict(model_dict)

    def get_prerun_dict(self) -> dict:
        """Output info to dict for use in protinfo."""
        d = {}
        if self.title:
            name = f"{self.pdbid} :: {self.title.title()}"
        else:
            name = self.pdbid

        if self.DO_NOT_USE:
            d["UNUSABLE"] = "File headers contain CAVEAT or Obsolete keywords."
        elif not self.has_header:
            d["HEADERLESS"] = "File has no header records to parse."
        else:
            # fill info dict:
            d["Function"] = self.function
            d["First Release"] = self.date
            d["Method"] = self.method
            d["Resolution"] = self.resolution
            d["Molecule"] = self.molecule
            d["Seqres Species"] = "; ".join(self.units_per_chain)
            if self.hetero_names:
                d["Cofactors"] = self.hetero_names
                d["Total cofactors"] = self.tot_cofactors
            if self.tot_waters:
                d["Total waters"] = self.tot_waters

            if self.missing is not None:
                d["Missing"] = self.missing
            if self.biounits:
                d["Biounits"] = (
                    len(self.biounits),
                    "; ".join(f"{b.serial}: {b.chains}" for b in self.biounits),
                )
            d["Models"] = str(self.n_models)
            if self.splits:
                d["Splits"] = " ".join(s for s in self.splits)
            d["Chains"] = ", ".join(self.chains)
            if self.ter_res:
                d["TER Residues"] = "; ".join(self.ter_res)

            d["Model 1 Residues"] = self.get_model1_res()
            d["Model 1 Free Cofactors & Waters"] = self.get_model1_heteros()  # includes waters
            if self.ssbonds:
                d["Disulfides"] = [
                    f"{lnk[0]} -- {lnk[1]} \u212B --> {lnk[2]}" for lnk in self.ssbonds
                ]
            if self.links:
                d["Links"] = [
                    f"{lnk[0]} -- {lnk[1]} \u212B --> {lnk[2]}" for lnk in self.links
                ]
            if self.sites:
                d["Sites"] = self.sites

        return {"Name": name, "PDB.Structure": d}

    def print_prerun(self):
        pretty_prerun = pformat(self.get_prerun_dict(), sort_dicts=False)
        sys.stdout.write(pretty_prerun)

    def print_summary(self) -> None:
        lines = ["Structure Summary\n"]
        lines.append("%15s: %s\n" % ("PDBID", self.pdbid))
        lines.append("%15s: %s\n" % ("First Release", self.date))
        lines.append("%15s: %s\n" % ("Function", self.function))
        lines.append("%15s: %s\n" % ("Molecule", self.molecule))
        lines.append("%15s: %s\n" % ("Method", self.method))
        lines.append("%15s: %s\n" % ("Resolution", self.resolution))
        lines.append(
            "%15s: %s\n" % ("Models", ",".join([m.serial for m in self.models]))
        )
        lines.append("%15s: %s\n" % ("Chains", ",".join(self.chains)))
        lines.append("%15s: %d\n" % ("Biounits", len(self.biounits)))
        for unit in self.biounits:
            lines.append("%20s %s: %s\n" % ("Unit", unit.serial, str(unit.chains)))

        sys.stdout.writelines(lines)

    def printme(self) -> None:
        """pretty print myself"""
        fmt = "{:>15}: {:}\n"
        lines = ["Structure Summary\n"]
        lines.append(fmt.format("PDBID", self.pdbid))
        lines.append(fmt.format("First Release", self.date))
        lines.append(fmt.format("Function", self.function))
        lines.append(fmt.format("Molecule", self.molecule))
        lines.append(fmt.format("Method", self.method))
        lines.append(fmt.format("Resolution", self.resolution))
        lines.append(fmt.format("Models", ",".join([m.serial for m in self.models])))
        lines.append(fmt.format("Chains", ",".join(self.chains)))
        lines.append(fmt.format("Biounits", len(self.biounits)))
        for unit in self.biounits:
            lines.append("{:>20} {:}: {}\n".format("Unit", unit.serial, unit.chains))

        sys.stdout.writelines(lines)


def update_dist_folder(pathname):
    """Update a path name by replacing $DIST_FOLDER with actual distribution folder name"""
    # WHY?
    # While import statements remain the same in different distribution methods, the relative path
    # of other folders such as parameter folders depend on the distribution methods. We allow macro
    # $DIST_FOLDER in parameter files so that users won't have to figure this out. All distribution
    # dependent updates are kept in mcce/.__init__.py
    from .mcce import _dist_folder

    fields = pathname.split(os.path.sep)
    if fields[0] == "$DIST_FOLDER":
        fields[0] = _dist_folder
    updated_pathname = os.path.sep.join(fields)
    return updated_pathname


def is_H(atom_name):
    H = False
    atom = atom_name.strip()
    if len(atom) == 4:  # full 4 chars
        if atom_name[0] == "H":
            H = True
        else:
            H = False
    else:
        if atom_name[1] == "H":
            H = True
        else:
            H = False
    return H
