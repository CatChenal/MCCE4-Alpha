# Feature Requests

## Pre-mcce
Check input structure

### Read PDB CIF format
**Background:** Due to the limitation of PDB format, Protein Data Bank is moving to CIF format.\
**Wanted features:** MCCE supports CIF format reading.\
**Comments:** Doable. (Junjun) 

### Increase strip ion list: DONE: `striph2o.py`
**Background:** By default, a small number of cofactors are stripped off when on surface. We need to expand this list.\
**Wanted features:** Expand the strip-off cofactor list or let users customize this list.\
**Comments:** The script striph2o.py in Stable-MCCE does exactly this. Run "striph2o.py -h" for help. (Junjun)

### Remove all waters but leave buried ions: DONE: `striph2o.py`
**Background:** By default, a small number of cofactors are stripped off when on surface. We need to expand this list.\
**Wanted features:** Expand the strip-off cofactor list or let users customize this list.\
**Comments:** The script striph2o.py in Stable-MCCE does exactly this. Run "striph2o.py -h" for help. (Junjun)

### Report chain breaks - note which ones should be NTR/CTR: DONE: `ProtInfo`
**Background:** Right now if a chain breakage happens, either at termials or at site of missing residues, mcce puts on NTR and NTR.\
**Wanted features:** Use chain ID to report the breakage.\
**Comments:** This should go in to pre-mcce check as well as step 1 default treatment. (Junjun)\
**DONE**: See [ProtInfo](https://github.com/GunnerLab/MCCE_ProtInfo)

### Report unknown ligands: DONE: `ProtInfo`
**Background:** Unknown ligands/cofactors do not have parameter files. They will be treated as 0-charge atom clusters and could \
potentially be a problem in calculation.\
**Wanted features:** Identify these ligands and let users decide what to do with them.\
**Comments:** Also provide instruction on how to make parameter file. (Junjun)\
**SCHEDULED**: See [ProtInfo issue #1](https://github.com/GunnerLab/MCCE_ProtInfo/issues/1)

### Report number of residues; Maximum dimension: DONE: `ProtInfo`
**Background:** Protein size is a useful piece of information.\
**Wanted features:** Report chains, biounit, size etc.\
**Comments:** Cat already has this done. I suggest these functions be integrated in a premcce class. (Junjun)\
**DONE**: See [ProtInfo](https://github.com/GunnerLab/MCCE_ProtInfo)

### Write CONFLIST for each ligands :-1:
**Background:** CONFLIST gives a rough idea about a residue's conformation choices identified as conformer types.\
**Wanted features:** Report conformer types except amino acids.\
**Comments:** Where, in run.log? Can you imaging the size of the log when running a lipid binding protein (assuming we'd have their\
tpl files)?

### IPECE
**Background:** Membrane proteins require low dielectric constant material as membrane to define correct dielectric boundary.\
**Wanted features:** Detect and add membrane atoms.\
**Comments:** This could be a summer intern project. Combining solvent accessibility and atom charges, one can find patches of hydrophic region.\
Proteins can be categorized into soluble proteins, peripheral proteins, and integral proteins based on how hydrophilic and hydrophobic surface \
residues are clustered. (Junjun) \
**Comments:** A slab of carbon atoms around a protein does not make much difference in mcce results.\
A better builder would be more reliable and incorporate biolipids; see [Charmm Membrane Builder](https://charmm-gui.org/?doc=tutorial&project=membrane).


## Step 1
Convert input structure to MCCE format, label terminal residues and ligands.

### More tpl files
**Background:** Currently we only covered a small fraction of ligands.\
**Wanted features:** Identify and prepare tpl files for more ligands.\
**Comments:** This could be a summer intern project. First find the most frequently appeared ligands, then prepare tools to help developing\
parameter files for these ligands. (Junjun)

### fix missing side chain atoms
**Background:** Some side chain atoms are missing.\
**Wanted features:** Complete these missing atoms.\
**Comments:** This could be done with extending the side chain bonds from the connectivity records in ftpl file or copying ideal positions\
by super-imposing the beginning atoms. (Junjun)

### Delete surface ions - where is the file: DONE: `striph2o.py`
**Background:** Surface ions are part of solution and should not be explicitly calculated.\
**Wanted features:** Strip off selected ions on the surface.\
**Comments:** See striph2o.py (Junjun)

### Delete surface or all waters: DONE: `striph2o.py`
**Background:** Water molecules are part of solution and should not be explicitly calculated.\
**Wanted features:** Strip off surface water molecules.\
**Comments:** See striph2o.py (Junjun)

### What to do with termini (NTR/CTR or leave as broken, non-ionizable chain)
**Background:** \
**Wanted features:** \
**Comments:** I am not quite clear about this question. (Junjun)


## Step 2

---
[Home](README.md) | [Feature Requests](Features.md) | [Developer Manual](DevManual.md) | [User Manual](UserManual.md)
