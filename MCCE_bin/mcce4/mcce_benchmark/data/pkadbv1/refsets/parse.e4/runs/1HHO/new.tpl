### This is a temporary parameter file made for residue PO4 ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST PO4        PO4BK 

NATOM    PO4BK      5

IATOM    PO4BK  P      0
IATOM    PO4BK  O1     1
IATOM    PO4BK  O2     2
IATOM    PO4BK  O3     3
IATOM    PO4BK  O4     4

ATOMNAME PO4BK    0  P  
ATOMNAME PO4BK    1  O1 
ATOMNAME PO4BK    2  O2 
ATOMNAME PO4BK    3  O3 
ATOMNAME PO4BK    4  O4 

CONNECT  PO4BK  P   ion        0    O1   0    O2   0    O3   0    O4 
CONNECT  PO4BK  O1  ion        0    P  
CONNECT  PO4BK  O2  ion        0    P  
CONNECT  PO4BK  O3  ion        0    P  
CONNECT  PO4BK  O4  ion        0    P  

### This is a temporary parameter file made for residue OXY ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST OXY        OXYBK 

NATOM    OXYBK      2

IATOM    OXYBK  O1     0
IATOM    OXYBK  O2     1

ATOMNAME OXYBK    0  O1 
ATOMNAME OXYBK    1  O2 

CONNECT  OXYBK  O1  ion        0    O2 
CONNECT  OXYBK  O2  ion        0    O1 

