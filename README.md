# MCCE4-Alpha

ALPHA MCCE version development. In this development cycle, code **will** be pure Python except PBE solver, for the benefit of
* flexible concept prototyping
* quick developing
* available functions already in NumPy and SciPy

## Installation Guide

1. Git clone the MCCE4-Alpha repsoitory to a desired place on your computer.
2. Add the path of your cloned repository into your ".bashrc" file. Use vi or nano to open your ".bashrc" file. 

Example: If you cloned MCCE4-Alpha to the directory /home/user/gunnerlab, you would add the following lines to your .bashrc file:
```
export PATH=/home/user/gunnerlab/MCCE4-Alpha/bin:$PATH
export PATH=/home/user/gunnerlab/MCCE4-Alpha/MCCE_bin:$PATH
```

3. Apply these PATH updates by sourcing your .bashrc file:
```
 source ~/.bashrc
```
or
```
 . ~/.bashrc
```

### Compile Executables 
MCCE4 contains c-code programs that must be compiled by the user prior to use. The three compilable programs are mcce, delphi and ngpb.
In this respository, MCCE4-Alpha, we provide the compiled delphi excecutable.

4. To compile mcce and ngpb, please run:
```
cd your_MCCE4_path
make clean
make
```
**Note: To install the NGPB container you must have sudo access ~15min.** 


### Test Installation
Commands like "step1.py", "p_info" should now be accessible; check with:
```
which step1.py
```
The command should return /home/user/gunnerlab/MCCE4-Alpha/bin/step1.py.  

Please, report any issues you encounter [here](https://github.com/GunnerLab/MCCE4-Alpha/issues).
Enjoy trying MCCE4!  

# MCCE Wiki
[Learn about MCCE, installation, available tools, and research done with MCCE.](https://mccewiki.levich.net/shelves) (under construction)

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

We recommend using a dedicated conda environment:
```
conda create -n oepython python=3.10
conda activate oepython
conda install -c openeye openeye-toolkits
conda install numpy scipy matplotlib seaborn pandas pygraphviz networkx requests
```

Add your license file to your PATH in your .bashrc file
```
export OE_LICENSE=/home/user/gunnerlab/openeye/oe_license.txt
```

Test if your license and installation are working with:
```
oecheminfo.py
```

Zap is run at MCCE4 step3. We highly recommend setting the salt concentration to 0.05 for best results.
```
step3.py -s zap -salt 0.05
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
