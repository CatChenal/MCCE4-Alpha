### This is a temporary parameter file made for residue THP ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST THP        THPBK 

NATOM    THPBK      25

IATOM    THPBK  P2     0
IATOM    THPBK  O4P    1
IATOM    THPBK  O5P    2
IATOM    THPBK  O6P    3
IATOM    THPBK  P1     4
IATOM    THPBK  O1P    5
IATOM    THPBK  O2P    6
IATOM    THPBK  O3P    7
IATOM    THPBK  O5'    8
IATOM    THPBK  C5'    9
IATOM    THPBK  C4'   10
IATOM    THPBK  O4'   11
IATOM    THPBK  C3'   12
IATOM    THPBK  O3'   13
IATOM    THPBK  C2'   14
IATOM    THPBK  C1'   15
IATOM    THPBK  N1    16
IATOM    THPBK  C2    17
IATOM    THPBK  O2    18
IATOM    THPBK  N3    19
IATOM    THPBK  C4    20
IATOM    THPBK  O4    21
IATOM    THPBK  C5    22
IATOM    THPBK  C5M   23
IATOM    THPBK  C6    24

ATOMNAME THPBK    0  P2 
ATOMNAME THPBK    1  O4P
ATOMNAME THPBK    2  O5P
ATOMNAME THPBK    3  O6P
ATOMNAME THPBK    4  P1 
ATOMNAME THPBK    5  O1P
ATOMNAME THPBK    6  O2P
ATOMNAME THPBK    7  O3P
ATOMNAME THPBK    8  O5'
ATOMNAME THPBK    9  C5'
ATOMNAME THPBK   10  C4'
ATOMNAME THPBK   11  O4'
ATOMNAME THPBK   12  C3'
ATOMNAME THPBK   13  O3'
ATOMNAME THPBK   14  C2'
ATOMNAME THPBK   15  C1'
ATOMNAME THPBK   16  N1 
ATOMNAME THPBK   17  C2 
ATOMNAME THPBK   18  O2 
ATOMNAME THPBK   19  N3 
ATOMNAME THPBK   20  C4 
ATOMNAME THPBK   21  O4 
ATOMNAME THPBK   22  C5 
ATOMNAME THPBK   23  C5M
ATOMNAME THPBK   24  C6 

CONNECT  THPBK  P2  ion        0    O4P  0    O5P  0    O6P  0    O5'
CONNECT  THPBK  O4P ion        0    P2 
CONNECT  THPBK  O5P ion        0    P2 
CONNECT  THPBK  O6P ion        0    P2 
CONNECT  THPBK  P1  ion        0    O1P  0    O2P  0    O3P  0    O3'
CONNECT  THPBK  O1P ion        0    P1 
CONNECT  THPBK  O2P ion        0    P1 
CONNECT  THPBK  O3P ion        0    P1 
CONNECT  THPBK  O5' ion        0    P2   0    C5'
CONNECT  THPBK  C5' ion        0    O5'  0    C4'
CONNECT  THPBK  C4' ion        0    C5'  0    O4'  0    C3'
CONNECT  THPBK  O4' ion        0    C4'  0    C1'
CONNECT  THPBK  C3' ion        0    C4'  0    O3'  0    C2'
CONNECT  THPBK  O3' ion        0    P1   0    C3'
CONNECT  THPBK  C2' ion        0    C3'  0    C1'
CONNECT  THPBK  C1' ion        0    O4'  0    C2'  0    N1 
CONNECT  THPBK  N1  ion        0    C1'  0    C2   0    C6 
CONNECT  THPBK  C2  ion        0    N1   0    O2   0    N3 
CONNECT  THPBK  O2  ion        0    C2 
CONNECT  THPBK  N3  ion        0    C2   0    C4 
CONNECT  THPBK  C4  ion        0    N3   0    O4   0    C5 
CONNECT  THPBK  O4  ion        0    C4 
CONNECT  THPBK  C5  ion        0    C4   0    C5M  0    C6 
CONNECT  THPBK  C5M ion        0    C5 
CONNECT  THPBK  C6  ion        0    N1   0    C5 

