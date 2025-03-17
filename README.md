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

<pre>
export PATH=/home/gunnerlab/MCCE4-Alpha/bin:$PATH
export PATH=/home/gunnerlab/MCCE4-Alpha/MCCE_bin:$PATH
</pre>

The important thing is that both the /bin/ and /MCCE_bin/ folders are referenced and accessible. Now, use command "source .bashrc". Commands like "step1.py", "p_info.py" should now be accessible- try using "which step1.py". If the command returns the path to MCCE4-Alpha/bin, you'll know it worked. 
Enjoy trying MCCE4-Alpha!

## MCCE4 Updates
MCCE4 now supports multiple Poisson-Boltzmann (PB) solvers, offering flexibility in electrostatic calculations. 
While **Delphi remains the default PB solver**, we have integrated **NextGenPB(NGPB)** and **Zap TK from OpenEye Scientific** as additional options for enhanced performance and accuracy.
In the near furture **NextGenPB(NGPB)** will be used as default MCCE4 PB solver.

To use the **Zap TK**, please obtain a OpenEye License (https://www.eyesopen.com/contact).
Follow instructions for installing OpenEye Toolkits at https://docs.eyesopen.com/toolkits/python/quickstart-python/linuxosx.html

We recommend using a dedicated conda enviorment:

<pre>
conda create -n oepython -c openeye openeye-toolkits python=3.10
conda install numpy scipy matplotlib pygraphviz pandas xlrd openpyxl requests
</pre>

Test if your license is working properly with 
<pre>
oecheminfo.py
</pre>pre> 


## Transition stage
Transition stage is when c code version mcce coexists with python mcce4. In this transition statge, you will need to compile the program:

<pre>
cd your_MCCE4_path
make clean
make
</pre>

Make sure to add both paths, bin and MCCE_bin, to you environment variable PATH. If a program are in both bin folders, bin/ has higher priority. Once we fully tested the newer version, the old one will be removed. 
<pre>
export PATH=your_MCCE4_path/bin:your_MCCE4_path/MCCE_bin:$PATH
</pre>

## Latest updates
- 03/17/2025 MCCE4-Alpha has been created for testing distribution
- 05/31/2024 Stable-MCCE is merged to this repository.
- 06/06/2004 Torsion energy is restored to C code algorithm.
