### This is a temporary parameter file made for residue _CU ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST _CU        _CUBK 

NATOM    _CUBK      1

IATOM    _CUBK CU      0

ATOMNAME _CUBK    0 CU  

CONNECT  _CUBK CU   ion      

### This is a temporary parameter file made for residue MPD ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST MPD        MPDBK 

NATOM    MPDBK      8

IATOM    MPDBK  C1     0
IATOM    MPDBK  C2     1
IATOM    MPDBK  O2     2
IATOM    MPDBK  CM     3
IATOM    MPDBK  C3     4
IATOM    MPDBK  C4     5
IATOM    MPDBK  O4     6
IATOM    MPDBK  C5     7

ATOMNAME MPDBK    0  C1 
ATOMNAME MPDBK    1  C2 
ATOMNAME MPDBK    2  O2 
ATOMNAME MPDBK    3  CM 
ATOMNAME MPDBK    4  C3 
ATOMNAME MPDBK    5  C4 
ATOMNAME MPDBK    6  O4 
ATOMNAME MPDBK    7  C5 

CONNECT  MPDBK  C1  ion        0    C2 
CONNECT  MPDBK  C2  ion        0    C1   0    O2   0    CM   0    C3 
CONNECT  MPDBK  O2  ion        0    C2 
CONNECT  MPDBK  CM  ion        0    C2 
CONNECT  MPDBK  C3  ion        0    C2   0    C4 
CONNECT  MPDBK  C4  ion        0    C3   0    O4   0    C5 
CONNECT  MPDBK  O4  ion        0    C4 
CONNECT  MPDBK  C5  ion        0    C4 

