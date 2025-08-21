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

### This is a temporary parameter file made for residue MRD ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST MRD        MRDBK 

NATOM    MRDBK      8

IATOM    MRDBK  C1     0
IATOM    MRDBK  C2     1
IATOM    MRDBK  O2     2
IATOM    MRDBK  CM     3
IATOM    MRDBK  C3     4
IATOM    MRDBK  C4     5
IATOM    MRDBK  O4     6
IATOM    MRDBK  C5     7

ATOMNAME MRDBK    0  C1 
ATOMNAME MRDBK    1  C2 
ATOMNAME MRDBK    2  O2 
ATOMNAME MRDBK    3  CM 
ATOMNAME MRDBK    4  C3 
ATOMNAME MRDBK    5  C4 
ATOMNAME MRDBK    6  O4 
ATOMNAME MRDBK    7  C5 

CONNECT  MRDBK  C1  ion        0    C2 
CONNECT  MRDBK  C2  ion        0    C1   0    O2   0    CM   0    C3 
CONNECT  MRDBK  O2  ion        0    C2 
CONNECT  MRDBK  CM  ion        0    C2 
CONNECT  MRDBK  C3  ion        0    C2   0    C4 
CONNECT  MRDBK  C4  ion        0    C3   0    O4   0    C5 
CONNECT  MRDBK  O4  ion        0    C4 
CONNECT  MRDBK  C5  ion        0    C4 

### This is a temporary parameter file made for residue IPA ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST IPA        IPABK 

NATOM    IPABK      4

IATOM    IPABK  C1     0
IATOM    IPABK  C2     1
IATOM    IPABK  C3     2
IATOM    IPABK  O2     3

ATOMNAME IPABK    0  C1 
ATOMNAME IPABK    1  C2 
ATOMNAME IPABK    2  C3 
ATOMNAME IPABK    3  O2 

CONNECT  IPABK  C1  ion        0    C2 
CONNECT  IPABK  C2  ion        0    C1   0    C3   0    O2 
CONNECT  IPABK  C3  ion        0    C2 
CONNECT  IPABK  O2  ion        0    C2 

