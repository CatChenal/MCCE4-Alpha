#!/usr/bin/env python

from mcce4.pdbio import *

class ROT_STAT_item:
    def __init__(self) -> None:
        """Rotamer making recording structure for each residue"""
        self.start = 0       # conformer count at the start
        self.swap = 0        # conformer count after swap
        self.rotate = 0      # conformer count after heavy atom rotamer
        self.swing = 0       # conformer count after self energy cleaning
        self.hbond = 0       # conformer count after hydrogen bond directed rotamer making
        self.repack = 0      # conformer count after repacking, GA repacking included
        self.xposed = 0      # conformer count after adding the most exposed conformer
        self.ioni = 0        # conformer count after ionization confomers are made, proton added in this step
        self.torh = 0        # conformer count after torsion minimum H is created
        self.oh = 0          # conformer count after H bond hydrogen is created
        self.prune = 0       # conformer count after clustering and election are made (limit conformer to <= 999)


class ROT_STAT:
    def __init__(self, prot) -> None:
        self.res_rot_stat = [ROT_STAT_item() for x in prot.residue]  # a list of ROT_STAT_item of residues

    def count_stat(self, prot, step=None):
        """Count the residue conformers and record to the step
        """
        attributes = [x for x in self.res_rot_stat[0].__dir__() if x[0] != "_"]
        if step and step in attributes:
            for i_res in range(len(prot.residue)):
                n_conf = len(prot.residue[i_res].conf)-1  # side chain conformers
                if n_conf < 0:
                    n_conf = 0
                    logging.warning("Detected 0 conformer residue %3s%c%04d%c" % prot.residue[i_res].resID)
                setattr(self.res_rot_stat[i_res], step, n_conf)
        else:
            logging.error("\"%s\" is not a valid rotamer making step name." % step)
            logging.error("Valid names are %s." % str(attributes))

    def write_stat(self, prot):
        header = "  Residue  Start   Swap Rotate  Swing Repack  Hbond Xposed   Ioni   TorH     OH  Prune\n"
        # items_in_header = len(header.strip().split()) - 1
        total_conf = ROT_STAT_item()
        attributes = [x for x in total_conf.__dir__() if x[0] != "_"]
        # print(attributes)

        lines = [header]

        for i_res in range(len(prot.residue)):
            res = prot.residue[i_res]
            res_name = "%3s%c%04d%c" % res.resID
            stat = self.res_rot_stat[i_res]
            line = "%9s %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d\n" % (res_name,
                                                                              stat.start,
                                                                              stat.swap,
                                                                              stat.rotate,
                                                                              stat.swing,
                                                                              stat.repack,
                                                                              stat.hbond,
                                                                              stat.xposed,
                                                                              stat.ioni,
                                                                              stat.torh,
                                                                              stat.oh,
                                                                              stat.prune)
            for attr in attributes:
                setattr(total_conf, attr, getattr(total_conf, attr) + getattr(stat, attr))
            lines.append(line)
        lines.append("-"*87 + "\n")

        line = "%9s %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d\n" % ("Total",
                                                                          total_conf.start,
                                                                          total_conf.swap,
                                                                          total_conf.rotate,
                                                                          total_conf.swing,
                                                                          total_conf.repack,
                                                                          total_conf.hbond,
                                                                          total_conf.xposed,
                                                                          total_conf.ioni,
                                                                          total_conf.torh,
                                                                          total_conf.oh,
                                                                          total_conf.prune)
        lines.append(line)

        return lines


class MCCE:
    def __init__(self, prm=None, tpl=None, structure=None):
        self.prm = prm
        self.tpl = tpl
        self.structure = structure
        self.protein = Protein()

    from .mcce._make_termini import make_termini
    from .mcce._convert_to_mccepdb import convert_to_mccepdb
    from .mcce._convert_to_mccepdb import initialize_atom_qr
    from .mcce._convert_to_mccepdb import initialize_atom_id
    from .mcce._make_connect import make_connect12
    from .mcce._make_connect import make_connect13
    from .mcce._make_connect import make_connect14
    from .mcce._make_connect import reset_connect
    from .mcce._make_connect import print_connect12
    from .mcce._center_protein import center_protein
    from .mcce._identify_ligands import identify_ligands
    from .mcce._atom_number_check import atom_number_check
    from .mcce._write_head1 import write_head1
    from .mcce._load_mccepdb import load_mccepdb
    from .mcce._rot_swap import rot_swap
    from .mcce._rot_rotate import rot_rotate
    from .mcce._rot_swing import rot_swing
    from .mcce._rot_swing import swing
    # from .mcce._rot_repack import rot_repack
    from .mcce._rot_repack_optimized import rot_repack
    from .mcce._rot_hdirected import rot_hdirected
    from .mcce._rot_xposed import rot_xposed
    from .mcce._clean_hvrot import clean_hvrot
    from .mcce._rot_ionization import rot_ionization
    from .mcce._place_missing import place_missing_heavy
    from .mcce._place_h import place_h
    from .mcce._hbond_h import hbond_h
    from .mcce._vdw import assign_vdw_param
    from .mcce._rot_prune import prune_conf
    from mcce4.mcce._vdw import make_blob

