### This is a temporary parameter file made for residue MES ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST MES        MESBK 

NATOM    MESBK      12

IATOM    MESBK  O1     0
IATOM    MESBK  C2     1
IATOM    MESBK  C3     2
IATOM    MESBK  N4     3
IATOM    MESBK  C5     4
IATOM    MESBK  C6     5
IATOM    MESBK  C7     6
IATOM    MESBK  C8     7
IATOM    MESBK  S      8
IATOM    MESBK  O1S    9
IATOM    MESBK  O2S   10
IATOM    MESBK  O3S   11

ATOMNAME MESBK    0  O1 
ATOMNAME MESBK    1  C2 
ATOMNAME MESBK    2  C3 
ATOMNAME MESBK    3  N4 
ATOMNAME MESBK    4  C5 
ATOMNAME MESBK    5  C6 
ATOMNAME MESBK    6  C7 
ATOMNAME MESBK    7  C8 
ATOMNAME MESBK    8  S  
ATOMNAME MESBK    9  O1S
ATOMNAME MESBK   10  O2S
ATOMNAME MESBK   11  O3S

CONNECT  MESBK  O1  ion        0    C2   0    C6 
CONNECT  MESBK  C2  ion        0    O1   0    C3 
CONNECT  MESBK  C3  ion        0    C2   0    N4 
CONNECT  MESBK  N4  ion        0    C3   0    C5   0    C7 
CONNECT  MESBK  C5  ion        0    N4   0    C6 
CONNECT  MESBK  C6  ion        0    O1   0    C5 
CONNECT  MESBK  C7  ion        0    N4   0    C8 
CONNECT  MESBK  C8  ion        0    C7   0    S  
CONNECT  MESBK  S   ion        0    C8   0    O1S  0    O2S  0    O3S
CONNECT  MESBK  O1S ion        0    S  
CONNECT  MESBK  O2S ion        0    S  
CONNECT  MESBK  O3S ion        0    S  

