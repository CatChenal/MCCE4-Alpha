# MCCE4-Alpha

Testing version of MCCE4 development. Please let us know about questions or comments! 

MCCE4 is now using [NextGenPB (NGPB) from the Rocchia Lab (IIT)](https://github.com/concept-lab/NextGenPB) as its default PBE solver.  

There are two ways you can install MCCE4-Alpha, which differ on whether a script is used: 
 * Option A: Semi-automated setup using provided scripts that download a generic NGPB image.
 * Option B: Manual setup that includes creation of a NGPB image optimized for your platform.

---
## Installation Option A: Quick Installation with Scripts
### Scripts:
These scripts automate many steps & download a generic NGPB image:
  * `MCCE_bin/quick_install.sh` (Linux, bash shell)
  * `MCCE_bin/quick_install_zsh.sh` (MacOS)

### 1. Clone the repository to a desired place on your computer (referred to as "clone_dir"):
```
 git clone https://github.com/GunnerLab/MCCE4-Alpha.git
```

### 2. Run the "Quick Install" script:
#### i. Go to your MCCE4-Alpha clone:
```
 cd ~/clone_dir/MCCE4-Alpha
```
#### ii. Run the approriate script:
* On MacOS, run:
```
 sh ./MCCE_bin/quick_install_zsh.sh
```

* On Linux, run:
```
 bash ./MCCE_bin/quick_install.sh
```

  Note: Creating 'install.log' is not required but is recommended as you could copy its contents if your created an "installation issue", which could help us fix an unexpected problem.

### 3. Follow the instructions displayed by the script to test your installation.
---

## Installation Option B: Installation with NGPB Optimized Image Creation:

### 1. Clone the repository to a desired place on your computer (referred to as "clone_dir"):
  * Git clone MCCE4-Alpha to a desired place on your computer:
  ```
   git clone https://github.com/GunnerLab/MCCE4-Alpha.git
  ```
 
  * Add the clone's bin paths to your `.bashrc` (`.bash_profile`) file then save it.
  ```
   export PATH="clone_dir/MCCE4-Alpha/bin:$PATH"
   export PATH="clone_dir/MCCE4-Alpha/MCCE_bin:$PATH"
  ```

  * Then apply the changes to your PATH variable by sourcing your `.bashrc` (`.bash_profile`) file, depending on your system.

  * Check a tool's command correct path location (tools do not require compiling):
  ```
    which p_info
  ```
  The command should return \your [clone_dir\]/MCCE4-Alpha/MCCE_bin/p_info


### 2. Compile Executables 
MCCE4 contains C and C++ libraries that must be compiled prior to use. Two compilations are required: for the mcce executable and for NGPB, the default PBE solver.  
  * NOTE: We also provide the DelPhi PBE solver executable file (`delphi`), without guarantying its runnability on all systems.

  * To compile mcce and ngpb:
    - Ensure you have sudo access as it is necessary for the installation of the NGPB container (~15min+)
    - `cd` into your MCCE4-Alpha clone directory
    - Run these commands:
      1. First and fast compilation command:
      ```
        make clean
      ```
      2. Second and slow compilation (& container setup) command:
        - IMPORTANT: In case this is not the first time the command is run, make sure to delete its target file, i.e.:
        ```
          rm bin/NextGenPB_MCCE4.sif
        ```
        The screen output of this long compilation is extensive and not recoverable if not directed to a file, so there are two way to run the command:
        - Run `make` command, without logging:
        ```
          make
        ```
        - Run `make` command, with redirection to a log file:
        ```
          make > make.log 2>&1
        ```
  
  * NOTE: To use the Openeye Zap solver, see "Zap Installation" below.

### 3. Test Installation
  * Create and activate a conda environment using MCCE4-Alpha environment file `mca.yml`. Choose either Command 1 or 2 below to create the environment:
    1. Command 1: To use the default environment name of 'mca':
       ```
        conda env create -f mca.yml
       ```
    2. Command 2: If you want something else, e.g. 'new_env' to be the environment name instead of 'mca':
       ```
        conda env create -f mca.yml -n new_env
       ```

  * Activate the mca env:
    ```
     conda activate mca
    ```

  * Check that a tool is functional; Its usage message should display:
    ```
     p_info
    ```
  * Display a command's help, e.g:
    ```
     step1.py -h
    ```


## Help us improve MCCE4
Please, report any issues you encounter [here](https://github.com/GunnerLab/MCCE4-Alpha/issues).
Thank you.
Enjoy trying MCCE4!  


## Supported PBE Solvers used in the energies calculation step (step 3) of MCCE4:
  
🚀  __NextGenPB (NGPB) from the Rocchia Lab (IIT)__ is now our default PBE solver.  
The former default PBE solver __Delphi from the Honig Lab__ and __Zap TK from OpenEye Scientific__ are also available.

---
### 🔹 NextGenPB (NGPB) - Rocchia Lab (IIT)
Developed by [Vincenzo di Florio](https://github.com/vdiflorio) at the [Rocchia Lab](https://github.com/concept-lab), IIT Genova, Italy.  

### 🔹 Zap TK - OpenEye Scientific 
To use **Zap TK**, you must obtain an **OpenEye license**:  
🔗 [Request a License](https://www.eyesopen.com/contact)  

#### Zap Installation:
  * Follow the OpenEye Toolkit Setup Instructions: 🔗 [License Setup](https://docs.eyesopen.com/toolkits/python/quickstart-python/license.html)

  * Add the toolkit to your mca environment:
    ```
     conda activate mca
     conda install -c openeye openeye-toolkits
    ```

  * Test if your license and installation are working with:
    ```
      oecheminfo.py
    ```

  * When using Zap, we highly recommend setting the salt concentration to 0.05 for best results, i.e:
    ```
     step3.py -s zap -salt 0.05
    ```

---

# Tools Info
Many (but not all!), tools have a succint description; you can list them all with:
```
 cat MCCE4-Alpha/MCCE_bin/tools.info
```

# MCCE Wiki
[Learn about MCCE, installation, available tools, and research done with MCCE.](https://mccewiki.levich.net/shelves) (under construction)


## Running Batches of Proteins

Sometimes it is convenient to run multiple proteins at time. Try using p_batch, located in the MCCE_bin.

p_batch accepts a directory containing pdb files, and optionally a shell script with custom instructions (level of conformers, what steps, dielectric const, etc.). For example, if the current working directory has a directory named "protein_list" containing 4lzt.pdb and 1a2p.pdb, you could run the following command to begin a default run of those proteins.

```
 p_batch protein_list
```

If you want to change the shell instructions, edit the default_script.sh, or create your own similarly structured file. If your shell script was named "custom_script.sh", you would run it like

```
 p_batch protein_list custom_script.sh
```

Additional custom instructions can be included in a "run.prm.custom" file. For example, if you want to simulate a membrane slab, or change the distance limit for reporting clashes, you can create your own "run.prm.custom" file (place in the same directory where you run p_batch). The many options for "run.prm.custom" files can be found in "MCCE4-Alpha/runprms/run.prm.full". You only need to include the desired options- default instructions are found in "MCCE4-Alpha/runprms/run.prm.default".

![Image](https://github.com/user-attachments/assets/6226520b-c3bf-40b6-bb07-ae78ad0c6e73)
