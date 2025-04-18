# MCCE4-Alpha

ALPHA MCCE version development. In this development cycle, code **will** be pure Python except PBE solver, for the benefit of
* flexible concept prototyping
* quick developing
* available functions already in NumPy and SciPy

## MCCE4 Project Essential Links

MCCE4 will stay as a private repository that only lab members and collaborators can see. A public version will be released separately. 

* [Repository](https://github.com/GunnerLab/MCCE4) - where you get and deposit code
* [Project](https://github.com/orgs/GunnerLab/projects/4) - automated project monitoring
* [Documentation](doc) - a place that makes your code actually usable
  * [Feature Requests](doc/Features.md)
  * [Developer Manual](doc/DevManual.md)
  * [User Manual](doc/UserManual.md)

* [Discussions](https://github.com/GunnerLab/MCCE4/discussions) - discussion board for anything

## Installation Guide

Git clone MCCE4-Alpha to a desired place on your computer. For ease of use, we recommend adding to the path in your ".bashrc" file. Use vi or nano to open your ".bashrc" file. If MCCE4-Alpha was located in /gunnerlab/, you would add the following lines:

```
export PATH=/home/gunnerlab/MCCE4-Alpha/bin:$PATH
export PATH=/home/gunnerlab/MCCE4-Alpha/MCCE_bin:$PATH
```

The important thing is that both the /bin/ and /MCCE_bin/ folders are referenced and accessible. Now, use command "source .bashrc". Commands like "step1.py", "p_info.py" should now be accessible- try using "which step1.py". If the command returns the path to MCCE4-Alpha/bin, you'll know it worked. 
Enjoy trying MCCE4-Alpha!

## MCCE4 Updates

### 🔹 Poisson-Boltzmann Solvers in MCCE4  
🚀  **NextGenPB (NGPB) from the Rocchia Lab (IIT)** is now our default PB solver 
The traditional PB solver **Delphi from the Honig Lab** and **Zap TK from OpenEye Scientific** are additional options.

---
### 🔹 NextGenPB (NGPB) - Rocchia Lab (IIT)
Developed by [Vincenzo di Florio](https://github.com/vdiflorio) at the [Rocchia Lab](https://github.com/concept-lab), IIT Genova, Italy.  

### 🔹 Delphi 
The standard PB solver used in MCCE4. Well-established for electrostatic calculations.  

### 🔹 Zap TK - OpenEye Scientific 
To use **Zap TK**, you must obtain an **OpenEye license**:  
🔗 [Request a License](https://www.eyesopen.com/contact)  

**Installation Instructions**   
Follow the OpenEye Toolkit installation guide:  
🔗 [Quickstart Guide](https://docs.eyesopen.com/toolkits/python/quickstart-python/linuxosx.html)  

We recommend using a dedicated Conda environment:  
```
conda create -n oepython -c openeye openeye-toolkits python=3.10
conda install numpy scipy matplotlib pygraphviz pandas xlrd openpyxl requests
```

Add your license file to your PATH in your .bashrc file
```
export OE_LICENSE=/home/gunnerlab/openeye/oe_license.txt
```

Test if your license and installation are working with:
```
oecheminfo.py
```

Zap is run at MCCE4 step3. We highly recommend setting the salt concentration to 0.05 for best results.
```
step3.py -s zap -salt 0.05
```


## Transition stage
Transition stage is when c code version mcce coexists with python mcce4. In this transition statge, you will need to compile the program:

```
cd your_MCCE4_path
make clean
make all
```
**Note to install the NGPB container you must have sudo access ~15min.**

Make sure to add both paths, bin and MCCE_bin, to you environment variable PATH. If a program are in both bin folders, bin/ has higher priority. Once we fully tested the newer version, the old one will be removed. 
```
export PATH=your_MCCE4_path/bin:your_MCCE4_path/MCCE_bin:$PATH
```
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

## Latest updates
- 03/17/2025 MCCE4-Alpha has been created for testing distribution
- 05/31/2024 Stable-MCCE is merged to this repository.
- 06/06/2004 Torsion energy is restored to C code algorithm.
