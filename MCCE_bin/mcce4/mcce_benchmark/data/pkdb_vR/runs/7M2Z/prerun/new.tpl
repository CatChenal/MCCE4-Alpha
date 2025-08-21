### This is a temporary parameter file made for residue GDP ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST GDP        GDPBK 

NATOM    GDPBK      28

IATOM    GDPBK  PB     0
IATOM    GDPBK  O1B    1
IATOM    GDPBK  O2B    2
IATOM    GDPBK  O3B    3
IATOM    GDPBK  O3A    4
IATOM    GDPBK  PA     5
IATOM    GDPBK  O1A    6
IATOM    GDPBK  O2A    7
IATOM    GDPBK  O5'    8
IATOM    GDPBK  C5'    9
IATOM    GDPBK  C4'   10
IATOM    GDPBK  O4'   11
IATOM    GDPBK  C3'   12
IATOM    GDPBK  O3'   13
IATOM    GDPBK  C2'   14
IATOM    GDPBK  O2'   15
IATOM    GDPBK  C1'   16
IATOM    GDPBK  N9    17
IATOM    GDPBK  C8    18
IATOM    GDPBK  N7    19
IATOM    GDPBK  C5    20
IATOM    GDPBK  C6    21
IATOM    GDPBK  O6    22
IATOM    GDPBK  N1    23
IATOM    GDPBK  C2    24
IATOM    GDPBK  N2    25
IATOM    GDPBK  N3    26
IATOM    GDPBK  C4    27

ATOMNAME GDPBK    0  PB 
ATOMNAME GDPBK    1  O1B
ATOMNAME GDPBK    2  O2B
ATOMNAME GDPBK    3  O3B
ATOMNAME GDPBK    4  O3A
ATOMNAME GDPBK    5  PA 
ATOMNAME GDPBK    6  O1A
ATOMNAME GDPBK    7  O2A
ATOMNAME GDPBK    8  O5'
ATOMNAME GDPBK    9  C5'
ATOMNAME GDPBK   10  C4'
ATOMNAME GDPBK   11  O4'
ATOMNAME GDPBK   12  C3'
ATOMNAME GDPBK   13  O3'
ATOMNAME GDPBK   14  C2'
ATOMNAME GDPBK   15  O2'
ATOMNAME GDPBK   16  C1'
ATOMNAME GDPBK   17  N9 
ATOMNAME GDPBK   18  C8 
ATOMNAME GDPBK   19  N7 
ATOMNAME GDPBK   20  C5 
ATOMNAME GDPBK   21  C6 
ATOMNAME GDPBK   22  O6 
ATOMNAME GDPBK   23  N1 
ATOMNAME GDPBK   24  C2 
ATOMNAME GDPBK   25  N2 
ATOMNAME GDPBK   26  N3 
ATOMNAME GDPBK   27  C4 

CONNECT  GDPBK  PB  ion        0    O1B  0    O2B  0    O3B  0    O3A
CONNECT  GDPBK  O1B ion        0    PB 
CONNECT  GDPBK  O2B ion        0    PB 
CONNECT  GDPBK  O3B ion        0    PB 
CONNECT  GDPBK  O3A ion        0    PB   0    PA 
CONNECT  GDPBK  PA  ion        0    O3A  0    O1A  0    O2A  0    O5'
CONNECT  GDPBK  O1A ion        0    PA 
CONNECT  GDPBK  O2A ion        0    PA 
CONNECT  GDPBK  O5' ion        0    PA   0    C5'
CONNECT  GDPBK  C5' ion        0    O5'  0    C4'
CONNECT  GDPBK  C4' ion        0    C5'  0    O4'  0    C3'
CONNECT  GDPBK  O4' ion        0    C4'  0    C1'
CONNECT  GDPBK  C3' ion        0    C4'  0    O3'  0    C2'
CONNECT  GDPBK  O3' ion        0    C3'
CONNECT  GDPBK  C2' ion        0    C3'  0    O2'  0    C1'
CONNECT  GDPBK  O2' ion        0    C2'
CONNECT  GDPBK  C1' ion        0    O4'  0    C2'  0    N9 
CONNECT  GDPBK  N9  ion        0    C1'  0    C8   0    C4 
CONNECT  GDPBK  C8  ion        0    N9   0    N7   0    C5   0    C4 
CONNECT  GDPBK  N7  ion        0    C8   0    C5 
CONNECT  GDPBK  C5  ion        0    C8   0    N7   0    C6   0    C4 
CONNECT  GDPBK  C6  ion        0    C5   0    O6   0    N1 
CONNECT  GDPBK  O6  ion        0    C6 
CONNECT  GDPBK  N1  ion        0    C6   0    C2 
CONNECT  GDPBK  C2  ion        0    N1   0    N2   0    N3 
CONNECT  GDPBK  N2  ion        0    C2 
CONNECT  GDPBK  N3  ion        0    C2   0    C4 
CONNECT  GDPBK  C4  ion        0    N9   0    C8   0    C5   0    N3 

