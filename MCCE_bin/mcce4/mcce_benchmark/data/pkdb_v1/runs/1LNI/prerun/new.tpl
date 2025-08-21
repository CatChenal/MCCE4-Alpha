### This is a temporary parameter file made for residue SO4 ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST SO4        SO4BK 

NATOM    SO4BK      5

IATOM    SO4BK  S      0
IATOM    SO4BK  O1     1
IATOM    SO4BK  O2     2
IATOM    SO4BK  O3     3
IATOM    SO4BK  O4     4

ATOMNAME SO4BK    0  S  
ATOMNAME SO4BK    1  O1 
ATOMNAME SO4BK    2  O2 
ATOMNAME SO4BK    3  O3 
ATOMNAME SO4BK    4  O4 

CONNECT  SO4BK  S   ion        0    O1   0    O2   0    O3   0    O4 
CONNECT  SO4BK  O1  ion        0    S  
CONNECT  SO4BK  O2  ion        0    S  
CONNECT  SO4BK  O3  ion        0    S  
CONNECT  SO4BK  O4  ion        0    S  

### This is a temporary parameter file made for residue GOL ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST GOL        GOLBK 

NATOM    GOLBK      6

IATOM    GOLBK  C1     0
IATOM    GOLBK  O1     1
IATOM    GOLBK  C2     2
IATOM    GOLBK  O2     3
IATOM    GOLBK  C3     4
IATOM    GOLBK  O3     5

ATOMNAME GOLBK    0  C1 
ATOMNAME GOLBK    1  O1 
ATOMNAME GOLBK    2  C2 
ATOMNAME GOLBK    3  O2 
ATOMNAME GOLBK    4  C3 
ATOMNAME GOLBK    5  O3 

CONNECT  GOLBK  C1  ion        0    O1   0    C2 
CONNECT  GOLBK  O1  ion        0    C1 
CONNECT  GOLBK  C2  ion        0    C1   0    O2   0    C3 
CONNECT  GOLBK  O2  ion        0    C2 
CONNECT  GOLBK  C3  ion        0    C2   0    O3 
CONNECT  GOLBK  O3  ion        0    C3 

