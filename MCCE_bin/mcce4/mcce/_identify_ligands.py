from ._strip_cofactors import Atom, Residue
from mcce4.geom import *
import math
import logging



# Detect and rename ligands by predefined rules

def _group_residues(lines):
    residues = []
    atoms = []
    
    # temp fix, atomlines contaminated
    # cleaned_lines = [line for line in lines if line[:6]=="ATOM  " or  line[:6]=="HETATM"]

    for line in lines:
        atom = Atom()
        atom.loadline(line)
        atoms.append(atom)

    resids = []
    for atom in atoms:
        if atom.resid in resids:
            ires = resids.index(atom.resid)
            residues[ires].atoms.append(atom)
        else:
            residue = Residue(atom.resid)
            residue.atoms = []
            residue.resid = atom.resid
            residue.atoms.append(atom)
            residues.append(residue)
            resids.append(atom.resid)

    return residues

def _rename_ligand(res, new_resName="XXX"):
    for atom in res.atoms:
        atom.resname = new_resName
        atom.resid = (new_resName, atom.resid[1:])

def match_strs(str1, str2):
    """Match two strings while allowing wildcard * to match any char
    The two strings are not order sensitive.
    Return False if not matched.
    Return True if matched.
    """
    matched = True
    nchar = len(str1)
    
    if nchar != len(str2):
        matched = False
    else:
        for i in range(nchar):
            if str1[i] != "*" and str2[i] != "*":  # only compare when no wildcard exists
                if str1[i] != str2[i]:
                    matched = False
                    break
    
    return matched


def identify_ligands(self):
    """Identify ligands in a structure and rename accordingly
    Return LINK lines
    """
    LINK_lines = []
    ssbond_serial = 1  # counts from 1

    residues = _group_residues(self.structure.lines)
    #print([res.resid for res in residues])
    n_residues = len(residues)
    for i_res in range(n_residues-1):
        res1_ID = residues[i_res].resid
        for j_res in range(i_res+1, n_residues):
            res2_ID = residues[j_res].resid
            key = ("LIGAND_ID", res1_ID[0], res2_ID[0])
            swapped_key = ("LIGAND_ID", res2_ID[0], res1_ID[0])
            found = False
            if key in self.tpl.db:
                found = True
                res1 = residues[i_res]
                res2 = residues[j_res]
                swapped = False
            elif swapped_key in self.tpl.db:
                found = True
                res2 = residues[i_res]
                res1 = residues[j_res]
                swapped = True
            #print(key, swapped_key, found)

            if found:
                ligand_param = self.tpl.db[("LIGAND_ID", res1.resid[0], res2.resid[0])]
                #print(res1.resid, res2.resid, vars(ligand_param))
                if not swapped:
                    atom1_name = ligand_param.atom1
                    atom2_name = ligand_param.atom2
                else:
                    atom1_name = ligand_param.atom2
                    atom2_name = ligand_param.atom1
                distance = ligand_param.distance
                tolerance = ligand_param.tolerance
                found_atom1 = False
                found_atom2 = False
                atom1 = []
                for atom in res1.atoms:
                    if match_strs(atom.name, atom1_name):
                        atom1.append(atom)
                        found_atom1 = True
                atom2 = []
                for atom in res2.atoms:
                    if match_strs(atom.name, atom2_name):
                        atom2.append(atom)
                        found_atom2 = True

                if found_atom1 and found_atom2:
                    d_min = 1.0E10  # a large number to start
                    # find the shortest distance in matched atoms
                    for a1 in atom1:
                        for a2 in atom2:
                            d = math.sqrt(ddvv(a1.xyz, a2.xyz))
                            if d < d_min:
                                d_min = d
                                a1_min = a1
                                a2_min = a2
                    if distance - tolerance <= d_min <= distance + tolerance:
                        logging.debug("Identified a ligand at d=%.2f between %s %s and %s %s" % (d_min,
                                                                                  atom1_name, 
                                                                                  res1.resid,
                                                                                  atom2_name,
                                                                                  res2.resid))
                        _rename_ligand(res1, new_resName=ligand_param.res1_name)
                        _rename_ligand(res2, new_resName=ligand_param.res2_name)
                        chainID_1 = res1_ID[1]
                        seqNum_1 = int(res1_ID[2])
                        chainID_2 = res2_ID[1]
                        seqNum_2 = int(res2_ID[2])
                        if res1_ID[0] == "CYS" and res2_ID[0] == "CYS":  # disulfur bond
                            fmt = "SSBOND %3d CYS %c %4d%c   CYS %c %4d%c %s  1555   1555 %5.2f\n"
                            line = fmt % (ssbond_serial, chainID_1, seqNum_1, " ",\
                                          chainID_2, seqNum_2, " ", " "*22, d_min)
                            LINK_lines.append(line)
                            ssbond_serial += 1
                        else:  # LINK
                            fmt = "LINK        %4s %3s %c%4d%c %s  %4s %3s %c%4d%c    1555   1555 %5.2f\n"
                            line = fmt % (a1_min.name, res1_ID[0], chainID_1, seqNum_1, " ", " "*12,\
                                          a2_min.name, res2_ID[0], chainID_2, seqNum_2, " ", d_min)
                            LINK_lines.append(line)
                            ssbond_serial += 1

                # Ignore when atom pair was not found.
                # else:
                #     logging.error("Atom pair %s of %s and %s of %s not found." % (atom1_name, 
                #                                                                   res1.resid,
                #                                                                   atom2_name,
                #                                                                   res2.resid))
                
    new_lines = []
    for res in residues:
        for atom in res.atoms:
            new_lines.append(atom.print_me())
    
    self.structure.lines = new_lines

    # clean up
    del(residues)

    return LINK_lines