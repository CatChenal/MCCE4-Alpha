import sys
import logging
import numpy as np
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd
from ._vdw import vdw_conf
from ..constants import PH2KCAL

NCONF_LIMIT = 999               # Maximum number of conformers allowed for a residue
PRUNE_VDW = 3 * PH2KCAL * 8     # Van der Waals energy cutoff, in kcal/mol, those get cut will have < 0.001 occupancy even at 1/8 scaling
PRUNE_RMSD = 0.5                # Root mean square deviation cutoff for clustering, in Angstrom


def rmsd_conf(conf1, conf2):
    """Calculate the root mean square deviation between two conformers.
    """
    atoms1 = conf1.atom
    atoms2 = conf2.atom
    if len(atoms1) != len(atoms2):
        logging.error("Conformers have different number of atoms. %s vs %s" % conf1.confID, conf2.confID)
        sys.exit()

    matched_coord1 = []
    matched_coord2 = []
    for atom1 in atoms1:
        matched = False
        for atom2 in atoms2:
            if atom1.name == atom2.name:
                matched_coord1.append(atom1.xyz)
                matched_coord2.append(atom2.xyz)
                matched = True
                break
        if not matched:
            print("Atom %s in conf1 not found in conf2" % atom1.name)
            sys.exit()

    if len(matched_coord1) == len(atoms1):
        matched_coord1 = np.array(matched_coord1)
        matched_coord2 = np.array(matched_coord2)
        rmsd = np.sqrt(np.sum((matched_coord1 - matched_coord2)**2)/len(matched_coord1))
    else:
        logging.error("Not all atoms matched between conformers %s and %s" % (conf1.confID, conf2.confID))
        sys.exit()

    return rmsd


def prune_conf(self):
    """Prune conformers based on the number of conformers allowed for a residue.
    """
    # Three steps of pruning
    # 1. Prune by self vdw energy
    # 2. Cluster conformers and select representative structures
    # 3. Prune by the number of conformers allowed for a residue

    if not hasattr(self.prm, "PRUNE_VDW"):
        prune_vdw = PRUNE_VDW
        logging.warning("   PRUNE_VDW is not set in runprm, use the default value %.3f kcal/mol." % prune_vdw)
    else:
        prune_vdw = float(self.prm.PRUNE_VDW.value)
        logging.info("   PRUNE_VDW read from run.prm.default %.3f kcal/mol." % prune_vdw)

    if not hasattr(self.prm, "PRUNE_RMSD"):
        prune_rmsd = PRUNE_RMSD
        logging.warning("   PRUNE_RMSD is not set in runprm, use the default value %.3f Angstroms." % prune_rmsd)   
    else:
        prune_rmsd = float(self.prm.PRUNE_RMSD.value)
        logging.info("   PRUNE_RMSD read from run.prm.default %.3f Angstroms." % prune_rmsd)

    logging.info("   Preparing connectivity and vdw parameters for pruning")
    self.reset_connect()
    self.make_connect12()
    self.make_connect13()
    self.make_connect14()
    self.make_blob()
    self.assign_vdw_param()

    
    # Prune by self vdw energy
    backbone_confs = [res.conf[0] for res in self.protein.residue]
    for res in self.protein.residue:
        if len(res.conf) > 1:
            for conf in res.conf[1:]:
                vdw0 = vdw_conf(conf, conf)
                vdw1 = 0.0
                for conf2 in backbone_confs:
                    vdw1 += vdw_conf(conf, conf2)
                conf.vdw = vdw0 + vdw1
            vdw_min = min([conf.vdw for conf in res.conf[1:]])
            for conf in res.conf[1:]:
                if conf.history[2] != "O" and conf.vdw > vdw_min + prune_vdw:   # don't prune original conformer
                    logging.debug("   Conformer %s self+backbone vdw = %.2f, > %.2f cut off of min = %.2f, prune." % (conf.confID, conf.vdw, prune_vdw, vdw_min))
                    res.conf.remove(conf)


    # Prune by clustering
        # clustering is confined within each conformer type
    # Group by conformer types
    for res in self.protein.residue:
        confs_to_keep = []
        conf_Types = self.tpl.db["CONFLIST", res.resID[0]]
        confs_grouped = {}
        for t in conf_Types:
            confs_grouped[t] = []
        logging.debug("Backbone conf (keep): %s" % res.conf[0].confID)
        confs_to_keep.append(res.conf[0])    # keep the backbone conf
        for conf in res.conf[1:]:    # skip the backbone conf
            confs_grouped[conf.confType].append(conf)
        # Cluster within each group
        for t in conf_Types:
            confs = confs_grouped[t]
            logging.debug("%s, %s" % (t, [conf.confID for conf in confs]))
            if len(confs) > 1:
                # 1. Generate the distance matrix
                n_confs = len(confs)
                d = np.zeros((n_confs, n_confs))
                for i in range(n_confs):
                    for j in range(i+1, n_confs):
                        d[i][j] = rmsd_conf(confs[i], confs[j]) 
                        d[j][i] = d[i][j]
                
                # 2. Perform hierarchical clustering
                linkage_matrix = sch.linkage(ssd.squareform(d), method="average")
                
                # 3. Find the representative structure of each cluster
                clusters = sch.fcluster(linkage_matrix, t=prune_rmsd, criterion="distance")
                for cluster in set(clusters):
                    cluster_indices = [i for i, x in enumerate(clusters) if x == cluster]
                    logging.debug(str(cluster_indices))    # the indices of conformers in the group, not the cluster
                    # Find the center of the cluster
                    cluster_confs = [confs[i] for i in cluster_indices]
                    center_index = cluster_indices[0]
                    min_rmsd_sum = float('inf')
                    for i in cluster_indices:
                        rmsd_sum = sum(rmsd_conf(confs[i], confs[j]) for j in cluster_indices)
                        if rmsd_sum < min_rmsd_sum:
                            min_rmsd_sum = rmsd_sum
                            center_index = i
                    logging.debug("Center of cluster (keep): %s" % confs[center_index].confID)
                    confs_to_keep.append(confs[center_index])
            elif len(confs) == 1:
                center_index = 0
                logging.debug("Single conformer in this conformer type (keep): %s" % confs[center_index].confID)
                confs_to_keep.append(confs[center_index])
        
        res.conf = confs_to_keep


    # Prune by the conformer number limit, last resort, unlikely to be used
    for res in self.protein.residue:
        if len(res.conf) > NCONF_LIMIT:
            logging.warning("Residue %s has %d conformers, throw away confomers after %dth conformer." % (res.resID, len(res.conf), NCONF_LIMIT))
            res.conf = res.conf[:NCONF_LIMIT]


        