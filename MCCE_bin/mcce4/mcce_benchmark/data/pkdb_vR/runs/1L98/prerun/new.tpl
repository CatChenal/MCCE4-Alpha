### This is a temporary parameter file made for residue BME ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST BME        BMEBK 

NATOM    BMEBK      4

IATOM    BMEBK  C1     0
IATOM    BMEBK  C2     1
IATOM    BMEBK  O1     2
IATOM    BMEBK  S2     3

ATOMNAME BMEBK    0  C1 
ATOMNAME BMEBK    1  C2 
ATOMNAME BMEBK    2  O1 
ATOMNAME BMEBK    3  S2 

CONNECT  BMEBK  C1  ion        0    C2   0    O1 
CONNECT  BMEBK  C2  ion        0    C1   0    S2 
CONNECT  BMEBK  O1  ion        0    C1 
CONNECT  BMEBK  S2  ion        0    C2 

