### This is a temporary parameter file made for residue CSD ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST CSD        CSDBK 

NATOM    CSDBK      8

IATOM    CSDBK  N      0
IATOM    CSDBK  CA     1
IATOM    CSDBK  CB     2
IATOM    CSDBK  SG     3
IATOM    CSDBK  OD1    4
IATOM    CSDBK  OD2    5
IATOM    CSDBK  C      6
IATOM    CSDBK  O      7

ATOMNAME CSDBK    0  N  
ATOMNAME CSDBK    1  CA 
ATOMNAME CSDBK    2  CB 
ATOMNAME CSDBK    3  SG 
ATOMNAME CSDBK    4  OD1
ATOMNAME CSDBK    5  OD2
ATOMNAME CSDBK    6  C  
ATOMNAME CSDBK    7  O  

CONNECT  CSDBK  N   ion        0    CA 
CONNECT  CSDBK  CA  ion        0    N    0    CB   0    C  
CONNECT  CSDBK  CB  ion        0    CA   0    SG 
CONNECT  CSDBK  SG  ion        0    CB   0    OD1  0    OD2
CONNECT  CSDBK  OD1 ion        0    SG 
CONNECT  CSDBK  OD2 ion        0    SG 
CONNECT  CSDBK  C   ion        0    CA   0    O  
CONNECT  CSDBK  O   ion        0    C  

### This is a temporary parameter file made for residue EDO ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST EDO        EDOBK 

NATOM    EDOBK      4

IATOM    EDOBK  C1     0
IATOM    EDOBK  O1     1
IATOM    EDOBK  C2     2
IATOM    EDOBK  O2     3

ATOMNAME EDOBK    0  C1 
ATOMNAME EDOBK    1  O1 
ATOMNAME EDOBK    2  C2 
ATOMNAME EDOBK    3  O2 

CONNECT  EDOBK  C1  ion        0    O1   0    C2 
CONNECT  EDOBK  O1  ion        0    C1 
CONNECT  EDOBK  C2  ion        0    C1   0    O2 
CONNECT  EDOBK  O2  ion        0    C2 

### This is a temporary parameter file made for residue _MG ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST _MG        _MGBK 

NATOM    _MGBK      1

IATOM    _MGBK MG      0

ATOMNAME _MGBK    0 MG  

CONNECT  _MGBK MG   ion      

