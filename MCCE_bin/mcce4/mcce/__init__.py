import logging
import os

# Constants for mcce main program



## Input folder locations
_distribution = "github"
#_distribution = "conda"


_mod_folder = str(os.path.dirname(os.path.abspath(__file__)))

if _distribution == "github":
    _dist_folder = "/".join(_mod_folder.split("/")[:-3])
elif _distribution == "conda":
    _dist_folder = _mod_folder
else:
    _dist_folder = ""

_ideal_structure = os.path.join(_dist_folder, "ideal_structure") 

## Output file names
_runprm_record_fname = "run.prm.record"
_step1_out = "step1_out.pdb"  # used by step 2
_head1 = "head1.lst"  # may be used by step 2

_step2_out = "step2_out.pdb" # used by step 3
_head2 = "head2.lst"  # information only
_rot_stat = "rot_stat"  # information only

_step3_out = "step3_out.pdb" # information only, identical to step2_out.pdb
_energies = "energies"  # used by step 4
_head3 = "head3.lst"  # used by step 4



## Logging config
_format = "%(asctime)s: %(levelname)s: %(message)s"
logging.basicConfig(format=_format, level=logging.DEBUG, datefmt="%D %H:%M:%S")

