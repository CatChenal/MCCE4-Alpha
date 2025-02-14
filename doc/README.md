# MCCE4 Wiki

## MCCE4 goal
MCCE4 is a project of rewriting the core code to improve the usability and scalability of biophysics simulation program [Stable-MCCE](https://github.com/GunnerLab/Stable-MCCE). 

MCCE4 will also 
* provide a public web interface for users who would try out MCCE
* provide tool set to help mcce runs and analyze results 

## MCCE4 Documentation

### Feature requests
[This section](Features.md) lists the proposed features and fixes. MCCE4 will keep the 4 step workflow of Stable-MCCE so the features will be categorized into:
  * step 0: evaluate the input structure
  * step 1: pre-mcce run: strip off water, identify ligands, split terminal residues and big cofactors into independent ionizable groups.
  * step 2: make side chain conformers
  * step 3: calculate energy look up table
  * step 4: simulated pH or Eh titration to determine charge states of each residue/cofactor
  * tools: Microstates analysis, energy analysis, and other tools for MCCE
  * web front: web interface to a mcce server

### Developer Manual
[This section](DevManual.md) provides brief description of the data structure and function usage with examples so that other developers can use them to build their code.

### User Manual
[This section](UserManual.md) provides detailed and organized description to MCCE program and tools with working examples.


---
[Home](README.md) | [Feature Requests](Features.md) | [Developer Manual](DevManual.md) | [User Manual](UserManual.md)
