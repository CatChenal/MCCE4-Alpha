from mcce4.pdbio import *
from ._make_connect import make_connect12, print_connect12


ntr_atoms = {" N  ", " CA "}  # heavy atoms only
ctr_atoms = {" C  ", " O  ", " OXT"}  # heavy atoms only
residue_names = {"GLY", "ALA", "VAL", "LEU", "ILE", "THR", "SER", "MET", "MEL",
                 "CYS", "CYD", "PRO", "PHE", "TYR", "TRP", "HIS", "HIL", "LYS",
                 "ARG", "ASP", "ASN", "GLN", "GLU", "NTR", "NTG", "CTR"}

def make_termini(self):
    """ Make the terminus for the protein in MCCE object
    """
    # make 1-2 connectivity for atoms
    # make_connect12(self)
    # print_connect12(self)

    # find the first and last amino acids on each chain
    ntr_res_ids = []
    ctr_res_ids = []
    amino_acids = [res for res in self.protein.residue if res.resID[0] in residue_names]
    chains = set()
    for res in amino_acids:
        chain = res.resID[1] 
        if chain not in chains:
            chains.add(chain)
    chains = list(chains)
    for chain in chains:
        peptide = [res for res in amino_acids if res.resID[1] == chain]
        #print(len(peptide))
        if peptide[0].resID[0] not in {"NTR", "NTG"}:  # identified a NTR residue
            ntr_res_ids.append(peptide[0].resID)
        if peptide[-1].resID[0] != "CTR":  # identified a CTR residue
            ctr_res_ids.append(peptide[-1].resID)
        

    new_residues = []
    for res in self.protein.residue:
        #print(res.resID)
        if res.resID in ntr_res_ids:
            new_res = Residue()
            new_bkconf = Conformer()  # BK conformer for NTR, 0 atoms
            new_conf = Conformer()  # the first conformer that stores terminus atoms
            for atom in res.conf[0].atom:
                if atom.name in ntr_atoms:
                    new_conf.atom.append(atom)
            # clean old conf
            res.conf[0].atom = [atom for atom in res.conf[0].atom if atom.name not in ntr_atoms]
            # insert BK and the 1st conformers
            new_res.conf.append(new_bkconf)
            new_res.conf.append(new_conf)

            # refresh the conformer, residue related info in the newly created conformer
            if res.resID[0] == "GLY":  # -> NTG
                new_resName = "NTG"
            else:
                new_resName = "NTR"
            for atom in new_res.conf[1].atom:
                atom.resName = new_resName
                atom.confType = new_resName+"01"
                atom.resID = (atom.resName, atom.chainID, atom.resSeq, atom.iCode)
                atom.history = "01"+atom.history[2:]
                atom.confNum = 1
            
            new_res.conf[0].confType = new_resName+"BK"
            new_res.conf[0].confID = "%5s%c%04d%c%03d" % (new_resName+"BK", 
                                                          new_res.conf[1].atom[0].chainID, 
                                                          new_res.conf[1].atom[0].resSeq, 
                                                          new_res.conf[1].atom[0].iCode, 
                                                          0)
            new_res.conf[1].confType = new_resName+"01"
            new_res.conf[1].confID = "%5s%c%04d%c%03d" % (new_resName+"01", 
                                                          new_res.conf[1].atom[0].chainID, 
                                                          new_res.conf[1].atom[0].resSeq, 
                                                          new_res.conf[1].atom[0].iCode, 
                                                          new_res.conf[1].atom[0].confNum)

            new_res.conf[0].altLoc = new_res.conf[1].altLoc = new_res.conf[1].atom[0].altLoc
            new_res.conf[0].resID = new_res.conf[1].resID = new_res.conf[1].atom[0].resID
            new_res.conf[0].history = new_res.conf[1].history = new_res.conf[1].atom[0].history

            new_res.resID = new_res.conf[0].resID

            # add this ntr BEFORE the parent residue
            new_residues.append(new_res)
            new_residues.append(res)

        elif res.resID in ctr_res_ids:
            new_res = Residue()
            new_bkconf = Conformer()  # BK conformer for CTR, 0 atoms
            new_conf = Conformer()  # the first conformer that stores terminus atoms
            for atom in res.conf[0].atom:
                if atom.name in ctr_atoms:
                    new_conf.atom.append(atom)
            # clean old conf
            res.conf[0].atom = [atom for atom in res.conf[0].atom if atom.name not in ctr_atoms]

            # insert BK and the 1st conformers
            new_res.conf.append(new_bkconf)
            new_res.conf.append(new_conf)

            # refresh the conformer, residue related info in the newly created conformer
            new_resName = "CTR"
            for atom in new_res.conf[1].atom:
                atom.resName = new_resName
                atom.confType = new_resName+"01"
                atom.resID = (atom.resName, atom.chainID, atom.resSeq, atom.iCode)
                atom.history = "01"+atom.history[2:]
                atom.confNum = 1
            
            new_res.conf[0].confType = new_resName+"BK"
            new_res.conf[0].resSeq = new_res.conf[1].atom[0].resSeq
            new_res.conf[0].confID = "%5s%c%04d%c%03d" % (new_resName+"BK", 
                                                          new_res.conf[1].atom[0].chainID, 
                                                          new_res.conf[1].atom[0].resSeq, 
                                                          new_res.conf[1].atom[0].iCode, 
                                                          0)
            new_res.conf[1].confType = new_resName+"01"
            new_res.conf[1].confID = "%5s%c%04d%c%03d" % (new_resName+"01", 
                                                          new_res.conf[1].atom[0].chainID, 
                                                          new_res.conf[1].atom[0].resSeq, 
                                                          new_res.conf[1].atom[0].iCode, 
                                                          new_res.conf[1].atom[0].confNum)
            new_res.conf[1].resSeq = new_res.conf[1].atom[0].resSeq
            
            new_res.conf[0].altLoc = new_res.conf[1].altLoc = new_res.conf[1].atom[0].altLoc
            new_res.conf[0].resID = new_res.conf[1].resID = new_res.conf[1].atom[0].resID
            new_res.conf[0].history = new_res.conf[1].history = new_res.conf[1].atom[0].history

            new_res.resID = new_res.conf[0].resID


            # add this ntr AFTER the parent residue
            new_residues.append(res)
            new_residues.append(new_res)

        else:
            new_residues.append(res)
    
    self.protein.residue = new_residues


        


        # double check the sequence number and chain ID to find possible missing connectivity

    #print(self.protein.residue)
    pass