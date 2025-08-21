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

