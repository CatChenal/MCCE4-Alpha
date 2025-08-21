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

