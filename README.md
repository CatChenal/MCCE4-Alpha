# MCCE4-Alpha

Public version of the private MCCE4 development repository.

## Installation Guide
### 1. Clone the repository and add its bin folders to your system PATH variable
  * Git clone MCCE4-Alpha to a desired place on your computer:
    ```
    git clone https://github.com/GunnerLab/MCCE4-Alpha.git
    ```
  * Add the clone's bin paths to your `.bashrc` (`.zshrc`) file then save it.
  If you cloned MCCE4-Alpha in ~/gunnerlab/ directory, you would add the following lines:
  ```
  export PATH="~/gunnerlab/MCCE4-Alpha/bin:$PATH"
  export PATH="~/gunnerlab/MCCE4-Alpha/MCCE_bin:$PATH"
  ```
  * Then apply the changes by sourcing or 'dotting' your `.bashrc` (`.zshrc`) file, depending on your system.

  * Check a tool's command correct path location (tools do not require compiling):
    ```
     which p_info
    ```
    The command should return \[your/root/dir/name\]/gunnerlab/MCCE4-Alpha/MCCE_bin/p_info


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
