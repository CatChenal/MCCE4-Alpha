# Multi-Conformation Continuum Electrostatics (MCCE4-Alpha)

<p align="center">
  <img src="docs/images/mcce_logo1.png" alt="MCCE Logo" style="max-width: 100%; height: auto;">
</p>

## Welcome to the **MCCE4-Alpha**! 
- This is a testing version of MCCE4 development. Please let us know about questions or comments!

In this guide walks you through calculating electrostatic potentials from a biomolecular structure.

🚀 **Get Started Now:** [**MCCE4 GitHub**](https://github.com/GunnerLab/MCCE4-Alpha)  
📖 **Full documentation:** [**MCCE4 Tutorial & Guide Website**](https://gunnerlab.github.io/mcce4_tutorial/)

---

## **Quick Introduction**

**MCCE4** is a physics-based computational biophysics program for predicting:

- **pKₐ values**
- **Protonation states**
- **Electrostatic properties** of biomolecules

In this program, protein side chain motions are simulated explicitly while the dielectric effect of solvent and bulk protein material is modeled by continuum electrostatics.

---

## **Documentation Overview**

### **Tutorial**
Step-by-step practical examples designed to help new users quickly run simulations and understand key features of MCCE4.

### **Guide**
Comprehensive documentation covering:
- Installation
- Advanced features
- Developer notes
- Detailed explanations of all settings

All information about installation and usage can be found in the [**MCCE4 Tutorial & Guide Website**](https://gunnerlab.github.io/mcce4_tutorial/)

---
There are two ways you can install MCCE4-Alpha, which differ on whether a script is used: 
 * Option A: Semi-automated setup using provided scripts that download a generic NGPB image.
 * Option B: Manual setup that includes creation of a NGPB image optimized for your platform.
   
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

# Basic 4-Step Run

Once the path to MCCE4-Alpha/bin is established, MCCE can be run. First, create a directory for the desired protein. We do not recommend mixing different MCCE runs in one directory.

`
mkdir 4LZT
cd 4LZT
`

Copy/link your chosen protein into the protein directory. If the desired protein exists on RCSB.org, you can use the tool "getpdb" to download a copy directly from there. For this example, we'll use 4LZT:

`
getpdb 4LZT
`

Now we can run MCCE. The easiest way is with "run_mcce4", which runs through the first four steps and calculates pKas for each residue of the PDB file, saving them to "pK.out" upon successful completion of the fourth step.

`
run_mcce4 4lzt.pdb
`

Other files are created by MCCE in the process of creating "pK.out". Learn about [these output files here.](https://mccewiki.levich.net/books/results/page/mcce-output-files) Learn more about [the four individual steps that comprise run_mcce4 here.](https://mccewiki.levich.net/books/mcce-tutorial-4lzt/page/calculate-pkas-of-lysozyme-mcce-steps-1-4)

[How do I customize a run?](https://mccewiki.levich.net/books/p-batch-tutorial/page/custom-mcce-runs-and-submit-shell)

[How do I run MCCE on multiple PDB files at once?](https://mccewiki.levich.net/books/p-batch-tutorial/page/how-do-i-run-multiple-proteins-at-once-p-batch-and-pro-batch)

# Tools Info
Many tools have a succint description in tools.info. You can list them all with:
```
 cat MCCE4-Alpha/MCCE_bin/tools.info
```

For more information about a given tool, check the wiki, or use the "-h" flag, e.g.:

`step3.py -h`

# MCCE Wiki
[Learn about MCCE, installation, available tools, and research done with MCCE.](https://mccewiki.levich.net) (under construction)




