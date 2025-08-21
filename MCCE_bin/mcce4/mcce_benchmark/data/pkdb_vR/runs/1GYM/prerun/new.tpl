### This is a temporary parameter file made for residue MYG ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST MYG        MYGBK 

NATOM    MYGBK      23

IATOM    MYGBK  C1A    0
IATOM    MYGBK  C2A    1
IATOM    MYGBK  C3A    2
IATOM    MYGBK  C4A    3
IATOM    MYGBK  C5A    4
IATOM    MYGBK  C6A    5
IATOM    MYGBK  C1     6
IATOM    MYGBK  C2     7
IATOM    MYGBK  C3     8
IATOM    MYGBK  C4     9
IATOM    MYGBK  C5    10
IATOM    MYGBK  C6    11
IATOM    MYGBK  N2    12
IATOM    MYGBK  O1A   13
IATOM    MYGBK  O2    14
IATOM    MYGBK  O3A   15
IATOM    MYGBK  O4A   16
IATOM    MYGBK  O5A   17
IATOM    MYGBK  O1    18
IATOM    MYGBK  O3    19
IATOM    MYGBK  O4    20
IATOM    MYGBK  O5    21
IATOM    MYGBK  O6    22

ATOMNAME MYGBK    0  C1A
ATOMNAME MYGBK    1  C2A
ATOMNAME MYGBK    2  C3A
ATOMNAME MYGBK    3  C4A
ATOMNAME MYGBK    4  C5A
ATOMNAME MYGBK    5  C6A
ATOMNAME MYGBK    6  C1 
ATOMNAME MYGBK    7  C2 
ATOMNAME MYGBK    8  C3 
ATOMNAME MYGBK    9  C4 
ATOMNAME MYGBK   10  C5 
ATOMNAME MYGBK   11  C6 
ATOMNAME MYGBK   12  N2 
ATOMNAME MYGBK   13  O1A
ATOMNAME MYGBK   14  O2 
ATOMNAME MYGBK   15  O3A
ATOMNAME MYGBK   16  O4A
ATOMNAME MYGBK   17  O5A
ATOMNAME MYGBK   18  O1 
ATOMNAME MYGBK   19  O3 
ATOMNAME MYGBK   20  O4 
ATOMNAME MYGBK   21  O5 
ATOMNAME MYGBK   22  O6 

CONNECT  MYGBK  C1A ion        0    C2A  0    C6A  0    O1A
CONNECT  MYGBK  C2A ion        0    C1A  0    C3A  0    O2 
CONNECT  MYGBK  C3A ion        0    C2A  0    C4A  0    O3A
CONNECT  MYGBK  C4A ion        0    C3A  0    C5A  0    O4A
CONNECT  MYGBK  C5A ion        0    C4A  0    C6A  0    O5A
CONNECT  MYGBK  C6A ion        0    C1A  0    C5A  0    O1 
CONNECT  MYGBK  C1  ion        0    C2   0    O1   0    O5 
CONNECT  MYGBK  C2  ion        0    C1   0    C3   0    N2 
CONNECT  MYGBK  C3  ion        0    C2   0    C4   0    O3 
CONNECT  MYGBK  C4  ion        0    C3   0    C5   0    O4 
CONNECT  MYGBK  C5  ion        0    C4   0    C6   0    O5 
CONNECT  MYGBK  C6  ion        0    C5   0    O6 
CONNECT  MYGBK  N2  ion        0    C2 
CONNECT  MYGBK  O1A ion        0    C1A
CONNECT  MYGBK  O2  ion        0    C2A
CONNECT  MYGBK  O3A ion        0    C3A
CONNECT  MYGBK  O4A ion        0    C4A
CONNECT  MYGBK  O5A ion        0    C5A
CONNECT  MYGBK  O1  ion        0    C6A  0    C1 
CONNECT  MYGBK  O3  ion        0    C3 
CONNECT  MYGBK  O4  ion        0    C4 
CONNECT  MYGBK  O5  ion        0    C1   0    C5 
CONNECT  MYGBK  O6  ion        0    C6 

