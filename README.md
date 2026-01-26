# Multi-Conformation Continuum Electrostatics

<p align="center">
  <img src="docs/images/mcce_logo1.png" alt="MCCE Logo" style="max-width: 100%; height: auto;">
</p>

# Welcome to MCCE4-Alpha!  
Please see our CHANGELOG at the bottom for the latest updates!

# [__📖 MCCE4-Alpha Tutorial__](https://gunnerlab.github.io/mcce4_tutorial/) 

Comprehensive documentation covering:
- Installation
- Guide: Detailed explanations of all settings
- Example Projects 

## CHANGELOG:
<!--- NOTE TO EDITOR: Use tis line to indicate that the user must/should update their clone"
  - __Apply changes: cd to your clone, then run `git pull`__
-->
_This section will reflect important changes and will provide you with information on how to apply them; For example, if new python packages are added to the environment file (mc4.yml), then the entry pertaining to that change will list the command(s) to update your environment._ 
* 2026-01-26:
  - Added `numba` in env file:
  - __Apply changes: cd to your clone, then run `git pull`__
  - __Apply changes: run these commands:__
  ```
  CLONE=$(dirname $(dirname "$(python3 -c "import os, sys; print(os.path.realpath(sys.argv[1]))" "$(which mcce)")"));
  conda env update -n mc4 -f $CLONE/mc4.yml
  ```

* 2026-01-20:
  - Comprehensive update of the tutorial site
  - Minimized README file


* 2026-01-08:
  - Updated python dependencies in mc4.yml
  - __Apply changes: run these commands:__
  ``` 
  CLONE=$(dirname $(dirname "$(python3 -c "import os, sys; print(os.path.realpath(sys.argv[1]))" "$(which mcce)")"));
  conda env update -n mc4 -f "$CLONE/mc4.yml
  ```

* 2025-11-25:
  - step1.py: Added error trapping on atom.loadline call
  - mfe.py: Updated & moved to MCCE_bin
  - __Apply changes: cd to your clone, then run `git pull`__

* 2025-11-11:
  - Fixed deleterious typo in bin/pdbs_interfaces.py
  - __Apply changes: cd to your clone, then run `git pull`__

* 2025-10-30:
  - Updated README: Added CHANGELOG, link to sudo_install.txt
  - Added topologies for SO4 and PO4 in param/.
  - Updated bin/step3.py with longer timeout value
  - Updated MCCE_bin/quick_install.sh
  - __Apply changes: cd to your clone, then run `git pull`__

---

## Help us improve MCCE4
This is a testing version of MCCE4 development. 
Please let us know about questions, comments or report any issues you encounter [here](https://github.com/GunnerLab/MCCE4-Alpha/issues).
Thank You and we hope you enjoy using MCCE4!  

## MCCE Wiki
[Learn about MCCE, installation, available tools, and research done with MCCE.](https://mccewiki.levich.net) (under construction)

---

Copyright (C) 2024 GunnerLab
This software is distributed under the terms the terms of the MIT licence
