import numpy as np
from ._vdw import vdw_conf
import random
import time
import logging
from scipy.sparse import csr_matrix 

Max_vdw_value = 999.0  # max possible vdw value allowed
Max_repack_steps = 10  # max number of repacks in a cycle. The repack may terminate early if the ms is converged
Max_cycle_repeats = 20 # max number of no-discovery cycle repeats. If cycles didn't find new microstate consecutively for this number, exit early
Use_dirty_vdw = False


# repack at "step2 -level 3" is slow and the following measures will be tried:
# - Construct a dirty vdw function to eliminate sqrt and use universal r0 = 4.0 and eps = 0.15
# - Adjust the order of the condition in vdw function to improve efficiency
# - optimize the vdw calculation by making a precalculated table
# 
# Original (max repack reduced to 20): 132m50s
# Original optimized: 27m32s at 20 cycles, 42m32s at 500 cycles
# Dirty vdw: 72m31s at 20 cycles
# Dirty vdw optimized:

def vdw_total(conf, microstate, vdw_lookup_table, mcce):
    vdw_backbone, vdw_pairwise = vdw_lookup_table
    # backbone (vdw0+vdw1 included in the lookup table)
    vdw = vdw_backbone[conf.i]
    # sidechain of microstate
    for ires in range(len(microstate)):
        res = mcce.protein.residue[ires]
        iconf = microstate[ires]
        conf2 = res.conf[iconf]

        if conf.resID != res.resID and iconf != 0:
            vdw += vdw_pairwise[(conf.i, conf2.i)]
    return vdw


def precalculate_vdw(self):
    """
    Precalculate vdw table for conformer to conformer pairs
    The backbone is named as "BK", conformer key is confID
    """
    #self.make_connect12()
    self.make_connect13()
    self.make_connect14()
    self.make_blob()
    # init conf serial number
    counter = 0
    for res in self.protein.residue:
        if len(res.conf) > 1:  # exclude backbone to save memory 
            for conf in res.conf[1:]:
                conf.i = counter
                counter += 1

    vdw_backbone = []
    vdw_pairwise = csr_matrix((counter, counter), dtype=float).toarray()
    # backbone and internal
    backbone_confs = [res.conf[0] for res in self.protein.residue]
    for res in self.protein.residue:
        if len(res.conf) > 1:
            for conf in res.conf[1:]:
                vdw_bk = vdw_conf(conf, conf, dirty=Use_dirty_vdw)      # vdw0
                for conf2 in backbone_confs:  # vdw1
                    vdw_bk += vdw_conf(conf, conf2, dirty=Use_dirty_vdw)
                vdw_backbone.append(vdw_bk)

    # pairwise
    n_res = len(self.protein.residue)
    for ir1 in range(n_res - 1):
        res1 = self.protein.residue[ir1]
        if len(res1.conf) > 1:
            for conf1 in res1.conf[1:]:
                for ir2 in range(ir1+1, n_res):
                    res2 = self.protein.residue[ir2]
                    if len(res2.conf) > 1:
                        for conf2 in res2.conf[1:]:
                            vdw = vdw_conf(conf1, conf2, dirty=Use_dirty_vdw)
                            if abs(vdw) > 0.00001:
                                vdw_pairwise[conf1.i][conf2.i] = vdw
                                vdw_pairwise[conf2.i][conf1.i] = vdw 

    # check mutual equality, for debug only
    # print("Checking mutual equality of pairwise interaction")
    # keys = vdw_table.keys()
    # for key in keys:
    #     if key[1] == "BK":  # we only calculated side chain to backbone, no BK to side chain, they are supposed to be different
    #         continue
    #     reversed_key = (key[1], key[0])
    #     vdw1 = vdw_table[key]
    #     if reversed_key in vdw_table:
    #         vdw2 = vdw_table[reversed_key]
    #     else:
    #         vdw2 = 0.0
    #     if abs(vdw1-vdw2) > 0.0001:
    #         print("Error: %s -> %s = %.3f; %s <- %s = %.3f" % (key[0], key[1], vdw1, key[0], key[1], vdw2))
    

    return (vdw_backbone, vdw_pairwise)


def rot_repack(self):
    """ Repack rotamers to select low enegy rotamers
    """
    # connect12 is inherited, clean_hvrot did connect 13 and 14 already, we do this just make sure
    # current_time = time.time()
    occ_cutoff = float(self.prm.REPACK_CUTOFF.value)
    max_repacks = int(self.prm.REPACKS.value)

    logging.info("   Prepare vdw energy lookup table. This may take a while...")
    vdw_lookup_table = precalculate_vdw(self)
    logging.info("   Done calculating vdw energy lookup table.")
    
    # print(time.time() - current_time)
    random.seed(time.time())

    converged_ms = []  # stats of converged microstates

    # loop over a pre-defined number of initial microstates, as in (REPACKS)
    for ipack in range(max_repacks):
        # generate a random initial state
        microstate = []
        for res in self.protein.residue:
            if len(res.conf) <= 1:
                iconf = 0  # no side chain
            else:
                iconf = random.choice(list(range(1, len(res.conf))))
            microstate.append(iconf)
        previous_microstate = microstate.copy()

        for istep in range(Max_repack_steps):
            # get the order of residues that repack will optimize
            optimize_res = list(range(len(self.protein.residue)))
            random.shuffle(optimize_res)


            # loop over residues to find the lowest energy conformer
            vdw_sum = 0.0
            for ires in optimize_res:
                res = self.protein.residue[ires]
                if len(res.conf) > 2:  # only need to search through reside with more than 1 side chain conformers
                    vdw_min = vdw_total(res.conf[1], microstate, vdw_lookup_table, self)  # use the first conf as inital vdw_min
                    iconf_min = 1
                    for iconf in range(2, len(res.conf)):
                        vdw = vdw_total(res.conf[iconf], microstate, vdw_lookup_table, self)
                        if vdw < vdw_min:
                            vdw_min = vdw
                            iconf_min = iconf
                    microstate[ires] = iconf_min
                    vdw_sum += vdw_min

            # test if the microstate has converged, quit repacking if converged
            #print(previous_microstate, "--->", microstate)
            break_step = istep + 1  # record current break step number for log
            if microstate != previous_microstate:
                previous_microstate = microstate.copy()
            else:
                break  # converged
            
        converged_ms.append(microstate)
        logging.debug("   Repacking cycle %4d of %4d: exit vdw = %.3f kcal/mol at step %d" % (ipack+1, max_repacks, vdw_sum, break_step))
        
    # convert converged_ms to conf_occ, and select only the conformers pass the cutoff threshold (REPACK_CUTOFF)
    n_ms = len(converged_ms)
    for ires in range(len(self.protein.residue)):
        res = self.protein.residue[ires]
        iconf_stat = [ms[ires] for ms in converged_ms]
        for iconf in range(len(res.conf)):
            conf = res.conf[iconf]
            conf.occ = iconf_stat.count(iconf)/n_ms
            if conf.occ > occ_cutoff or iconf < 2:
                conf.keep = True
            else:
                conf.keep = False
        res.conf = [conf for conf in res.conf if conf.keep]


