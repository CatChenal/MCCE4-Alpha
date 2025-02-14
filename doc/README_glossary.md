# glossary

Search for parameter keys begining with the user-provided prefix.

Notes:
  - The search is case-sensitive: pass an all-caps (sub)string if searching for keys, else a capitalized string.
  - If searching for step information for a specific step, use StepN, e.g. Step1.
  - The list of keys may contain obsolete keys or missing ones as the list is based on the old run.prm.full file.

## Examples:

```
> glossary REL
Result of query: 'REL'

RELAX_WAT
  Description: Step2. Do relaxation on waters
  Default: f

RELAX_H
  Description: Step2. Do relaxation on hydrogens
  Default: t

RELAX_E_THR
  Description: Step2. Energy threshold for keeping a conformer
  Default: -5.0

RELAX_NSTATES
  Description: Step2. Loop over N local microstates
  Default: 200

RELAX_N_HYD
  Description: Step2. Default number of hydroxyl positions
  Default: 6

RELAX_NITER
  Description: Step2. Maximum number of relaxation steps
  Default: 300

RELAX_CLASH_THR
  Description: Step2. Do not relax hydrogen if vdw of two sidechain conformers bigger than this
  Default: 5.

RELAX_PHI
  Description: Step2. Phi for each step of relaxation
  Default: 1.0

RELAX_TORQ_THR
  Description: Step2. Torque threshold for keep relaxing
  Default: 0.5
```

```
> glossary Step1
Result of query: 'Step1'

Step1: basic structural checks & pdb format conversion
  Key: DO_PREMCCE
  Default: f

Step1. Label terminal residues?
  Key: TERMINALS
  Default: t

Step1. Remove waters with SAS% above this value
  Key: H2O_SASCUTOFF
  Default: 0.1

Step1. Distance limit of reporting clashes
  Key: CLASH_DISTANCE
  Default: 2.0

Step1. Ignore hydrogens in input structure
  Key: IGNORE_INPUT_H
  Default: t
```
