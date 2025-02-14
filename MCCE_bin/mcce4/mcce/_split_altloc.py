from mcce4.pdbio import *
from mcce4.pdbio import _Model

import copy

# A backbone atom should have residue name in residue_names AND name in backbone_atoms
backbone_atoms = {" N  ", " CA ", " C  ", " O  "}
residue_names = {"GLY", "ALA", "VAL", "LEU", "ILE", "THR", "SER", "MET", "MEL",
                 "CYS", "CYD", "PRO", "PHE", "TYR", "TRP", "HIS", "HIL", "LYS",
                 "ARG", "ASP", "ASN", "GLN"}



def split_altloc(model):
    """Split the input structure if multile altLoc on backbone atoms are deteceted."""
    splited_models = []
    collected_altLoc = set()
    for line in model.lines:
        atom_name = line[12:16]
        res_name = line[17:20]
        altLoc = line[16]
        # print(line)
        # print("===%s===%s===%s===" % (atom_name, res_name, altLoc))
        if atom_name in backbone_atoms and \
            res_name in residue_names and \
            altLoc != " ":
            collected_altLoc.add(altLoc)

    altLocs = list(collected_altLoc)
    altLocs.sort()

    #print("HERE: %s" % str(altLocs))
    if len(altLocs) > 1:
        for altLoc in altLocs:
            # Due to the multi-sites backbone altLoc, the correct way to handle them is to identify 
            # segments first and group altLoc for each segments. However, this may results in 
            # combinatorial explosion when segment number goes up. Here I am going to use the matching 
            # altLoc from all residue when available, and use the first altLoc when not matched.
            new_model = _Model()
            new_model.altLoc = altLoc
            line = model.lines[0]
            previous_line = line  # intial value
            previous_atom = line[12:16] + line[17:20] + line[21:26]
            #print(altLoc)
            for this_line in model.lines[1:]:
                altLoc_this_atom = this_line[16]
                #print(altLoc_this_atom)
                this_atom = this_line[12:16] + this_line[17:20] + this_line[21:26]
                #print(this_atom, previous_atom)
                if this_atom != previous_atom:  # new atom found, save previous line
                    new_model.lines.append(previous_line)
                    previous_line = this_line
                    previous_atom = this_atom
                elif altLoc == altLoc_this_atom:  # exact match, use this line, only save when a new atom is encountered
                    # print(previous_line)
                    # print(this_line)
                    # print("-"*80)
                    previous_line = this_line

                elif altLoc_this_atom == " ":  # empty altLoc means shared stoms
                    previous_line = this_line
                elif altLoc_this_atom == "A":  # when above condition not satisfied, use this as default value
                    previous_line = this_line
            # when exit, save the last line
            
            new_model.lines.append(previous_line)
            new_model.lines = [l for l in new_model.lines if l[:6]=="ATOM  "or l[:6]=="HETATM"]
            splited_models.append(new_model)

    else:  # only one altLoc
        new_model = model
        new_model.lines = [l for l in new_model.lines if l[:6]=="ATOM  "or l[:6]=="HETATM"]
        splited_models.append(new_model)
    
    # for x in splited_models:
    #     for l in x.lines:
    #         print(l)

    return splited_models

