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
* 06/06/2004 Torsion energy is restored to C code algorithm.
- 05/31/2024 Stable-MCCE is merged to this repository.
