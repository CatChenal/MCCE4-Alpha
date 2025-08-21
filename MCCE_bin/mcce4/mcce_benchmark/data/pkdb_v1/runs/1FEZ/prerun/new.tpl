### This is a temporary parameter file made for residue _MG ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST _MG        _MGBK 

NATOM    _MGBK      1

IATOM    _MGBK MG      0

ATOMNAME _MGBK    0 MG  

CONNECT  _MGBK MG   ion      

### This is a temporary parameter file made for residue WO4 ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST WO4        WO4BK 

NATOM    WO4BK      5

IATOM    WO4BK  W      0
IATOM    WO4BK  O1     1
IATOM    WO4BK  O2     2
IATOM    WO4BK  O3     3
IATOM    WO4BK  O4     4

ATOMNAME WO4BK    0  W  
ATOMNAME WO4BK    1  O1 
ATOMNAME WO4BK    2  O2 
ATOMNAME WO4BK    3  O3 
ATOMNAME WO4BK    4  O4 

CONNECT  WO4BK  W   ion        0    O1   0    O2   0    O3   0    O4 
CONNECT  WO4BK  O1  ion        0    W  
CONNECT  WO4BK  O2  ion        0    W  
CONNECT  WO4BK  O3  ion        0    W  
CONNECT  WO4BK  O4  ion        0    W  

