#!/usr/bin/env python
"""
Test swing function
"""
import logging

from mcce4.runprm import *
from mcce4.pdbio import *
from mcce4.mcce import *
from mcce4 import __version__
from mcce4.mcce import _step1_out
from mcce4.mcce import _step2_out
from mcce4.mcce import _rot_stat
from mcce4.main import *

newtpl_name = "new"

if __name__ == "__main__":
    print("MCCE4 version: %s\n" % __version__)
    logging.info("Running SWING TEST")
    # 3 key elements for a calculation
    # run.prm - configuration of a run
    # tpl files - parameter files (optional in step 4)
    # input structure - (optional in step 2, 3, and 4)

    # Get runprm
    
    runprm = RunPrm()                   # create a runprm object and load default values
  
    
    # Get tpl
    tpl = TPL()
    tpl.read_ftpl_folder(runprm.FTPL_FOLDER.value)
    tpl.dump(fname="ftpl.record")
        
    newftpl_file = "%s.ftpl" % (newtpl_name)        
    if os.path.isfile(newftpl_file):
        tpl.read_ftpl_file(newftpl_file)


    # Here we don't pass structure to MCCE class, make an empty mcce.protein object instead
    mcce = MCCE(prm=runprm, tpl=tpl)
    mcce.convert_to_mccepdb()  # create an empty mcce.protein object
 
    # Load mccepdb to mcce.protein object
    logging.info("Loading structure from %s." % _step1_out)
    mcce.load_mccepdb(_step1_out)

    # Place missing heavy atoms
    logging.info("Place missing heavy atoms ...")
    while True:
        if mcce.place_missing_heavy() == 0:
            break

    logging.info("Make atom conenctivity ...")
    mcce.make_connect12()  # make connect12 every time new conformers are created

    rot_stat = ROT_STAT(mcce.protein)
    # test
    for res in mcce.protein.residue:
        if len(res.conf) > 1:
            new_confs = []
            for conf in res.conf[1:]:
                new_confs += mcce.swing(conf, 60/180*3.1416)
            res.conf += new_confs
    rot_stat.count_stat(mcce.protein, step="swing")
    lines = rot_stat.write_stat(mcce.protein)
    open(_rot_stat, "w").writelines(lines)

    mcce.protein.dump(_step2_out)

