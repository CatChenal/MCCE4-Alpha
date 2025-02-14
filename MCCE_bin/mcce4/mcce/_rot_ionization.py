import logging
from ..pdbio import *

def rot_ionization(self):
    """ Propogate conformers based on CONFLIST
    Basic algorithm:
    - For each residue, get CONFLIST record from ftpl database, and examine the number of CONNECT of that confType.
    - Loop over each non-backbone conformer, multiply the confomer with the confType it didn't have, update confType and history properly
    - Place missing atom at the end, in case some confType has different heavy atom connectivity
    """

    # make a dictionary to store confomer type and connect record count
    confType_counts = {}
    for key in self.tpl.db.keys():
        if key[0] == "CONNECT":
            confType = key[2]
            if confType in confType_counts:
                confType_counts[confType] += 1
            else:
                confType_counts[confType] = 1

    for res in self.protein.residue:
        key = ("CONFLIST", res.resID[0])
        if key in self.tpl.db:
            new_confs = []
            all_conf_types = self.tpl.db[key]
            nondummy_conf_types = set([x for x in all_conf_types if x in confType_counts and x[-2:] != "BK"])

            # duplicate to other conformer types
            if len(res.conf) > 1:
                for parent_conf in res.conf[1:]:
                    # print(nondummy_conf_types, parent_conf.confType)
                    other_conftypes = nondummy_conf_types - set([parent_conf.confType])
                    if other_conftypes:
                        for other_conftype in list(other_conftypes):
                            new_conf = parent_conf.clone()
                            new_conf.confType = other_conftype
                            new_conf.history = other_conftype[-2:] + parent_conf.history[2:]
                            new_confs.append(new_conf)

                # merge the new conformers and sort by confType
                all_confs = res.conf[1:] + new_confs
                ordered_confs = []
                for conf_type in all_conf_types:
                    ordered_confs += [c for c in all_confs if c.confType == conf_type]

                res.conf = [res.conf[0]] + ordered_confs



        else:
            logging.warning("   Residue %s doesn't have a CONFLIST entry" % res.resName)
