### This is a temporary parameter file made for residue EPE ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST EPE        EPEBK 

NATOM    EPEBK      15

IATOM    EPEBK  N1     0
IATOM    EPEBK  C2     1
IATOM    EPEBK  C3     2
IATOM    EPEBK  N4     3
IATOM    EPEBK  C5     4
IATOM    EPEBK  C6     5
IATOM    EPEBK  C7     6
IATOM    EPEBK  C8     7
IATOM    EPEBK  O8     8
IATOM    EPEBK  C9     9
IATOM    EPEBK  C10   10
IATOM    EPEBK  S     11
IATOM    EPEBK  O1S   12
IATOM    EPEBK  O2S   13
IATOM    EPEBK  O3S   14

ATOMNAME EPEBK    0  N1 
ATOMNAME EPEBK    1  C2 
ATOMNAME EPEBK    2  C3 
ATOMNAME EPEBK    3  N4 
ATOMNAME EPEBK    4  C5 
ATOMNAME EPEBK    5  C6 
ATOMNAME EPEBK    6  C7 
ATOMNAME EPEBK    7  C8 
ATOMNAME EPEBK    8  O8 
ATOMNAME EPEBK    9  C9 
ATOMNAME EPEBK   10  C10
ATOMNAME EPEBK   11  S  
ATOMNAME EPEBK   12  O1S
ATOMNAME EPEBK   13  O2S
ATOMNAME EPEBK   14  O3S

CONNECT  EPEBK  N1  ion        0    C2   0    C6   0    C9 
CONNECT  EPEBK  C2  ion        0    N1   0    C3 
CONNECT  EPEBK  C3  ion        0    C2   0    N4 
CONNECT  EPEBK  N4  ion        0    C3   0    C5   0    C7 
CONNECT  EPEBK  C5  ion        0    N4   0    C6 
CONNECT  EPEBK  C6  ion        0    N1   0    C5 
CONNECT  EPEBK  C7  ion        0    N4   0    C8 
CONNECT  EPEBK  C8  ion        0    C7   0    O8 
CONNECT  EPEBK  O8  ion        0    C8 
CONNECT  EPEBK  C9  ion        0    N1   0    C10
CONNECT  EPEBK  C10 ion        0    C9   0    S  
CONNECT  EPEBK  S   ion        0    C10  0    O1S  0    O2S  0    O3S
CONNECT  EPEBK  O1S ion        0    S  
CONNECT  EPEBK  O2S ion        0    S  
CONNECT  EPEBK  O3S ion        0    S  

