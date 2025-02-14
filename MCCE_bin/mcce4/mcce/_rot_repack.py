import numpy as np
from ._vdw import vdw_conf
import random
import time
import logging

Max_repack_steps = 10  # max number of repacks in a cycle. The repack may terminate early if the ms is converged
Max_cycle_repeats = 20 # max number of no-discovery cycle repeats. If cycles didn't find new microstate consecutively for this number, exit early
Use_dirty_vdw = True


# repack at "step2 -level 3" is slow and the following measures will be tried:
# - Construct a dirty vdw function to eliminate sqrt and use universal r0 = 4.0 and eps = 0.15
# - Adjust the order of the condition in vdw function to improve efficiency
# 
# Original (max repack reduced to 20): 132m50s
# Original optimized:
# Dirty vdw: 72m31s
# Dirty vdw optimized:

def vdw_total(conf, microstate, mcce):
    # internal vdw0 
    vdw = vdw_conf(conf, conf, dirty=Use_dirty_vdw)
    # backbone vdw1
    for res in mcce.protein.residue:
        vdw += vdw_conf(conf, res.conf[0], dirty=Use_dirty_vdw)
    # sidechain of microstate
    for ires in range(len(microstate)):
        res = mcce.protein.residue[ires]
        iconf = microstate[ires]
        conf2 = res.conf[iconf]

        if conf.resID != res.resID and iconf != 0:
            vdw += vdw_conf(conf, conf2, dirty=Use_dirty_vdw)
    return vdw

def rot_repack(self):
    """ Repack rotamers to select low enegy rotamers
    """
    # connect12 is inherited, clean_hvrot did connect 13 and 14 already, we do this just make sure
    # current_time = time.time()
    occ_cutoff = float(self.prm.REPACK_CUTOFF.value)
    max_repacks = int(self.prm.REPACKS.value)
    max_repacks = 20

    self.make_connect13()
    self.make_connect14()
    self.make_blob()
    # print(time.time() - current_time)
    random.seed(time.time())

    converged_ms = []  # stats of converged microstates

    # loop over a pre-defined number of initial microstates, as in (REPACKS)
    cycle_repeats = 0
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
                    vdw_min = vdw_total(res.conf[1], microstate, self)  # use the first conf as inital vdw_min
                    iconf_min = 1
                    for iconf in range(2, len(res.conf)):
                        vdw = vdw_total(res.conf[iconf], microstate, self)
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
            
        # test if this is a newly discovered microstate
        if microstate in converged_ms:
            cycle_repeats += 1
        else:
            cycle_repeats = 0
        converged_ms.append(microstate)
        logging.debug("   Repacking cycle %4d of %4d: exit vdw = %.3f kcal/mol at step %d" % (ipack+1, max_repacks, vdw_sum, break_step))
        

        if cycle_repeats >= Max_cycle_repeats:
            logging.info("   The last %d cycles didn't find new microstate, optimization can exit early." % cycle_repeats)
            break

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


