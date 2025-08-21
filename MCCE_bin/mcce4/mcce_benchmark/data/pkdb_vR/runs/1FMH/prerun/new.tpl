### This is a temporary parameter file made for residue ACE ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST ACE        ACEBK 

NATOM    ACEBK      6

IATOM    ACEBK  CH3    0
IATOM    ACEBK  CH3    1
IATOM    ACEBK  C      2
IATOM    ACEBK  O      3
IATOM    ACEBK  C      4
IATOM    ACEBK  O      5

ATOMNAME ACEBK    0  CH3
ATOMNAME ACEBK    1  CH3
ATOMNAME ACEBK    2  C  
ATOMNAME ACEBK    3  O  
ATOMNAME ACEBK    4  C  
ATOMNAME ACEBK    5  O  

CONNECT  ACEBK  CH3 ion        0    CH3  0    C    0    C  
CONNECT  ACEBK  CH3 ion        0    CH3  0    C    0    C  
CONNECT  ACEBK  C   ion        0    CH3  0    CH3  0    O    0    C    0    O  
CONNECT  ACEBK  O   ion        0    C    0    C    0    O  
CONNECT  ACEBK  C   ion        0    CH3  0    CH3  0    C    0    O    0    O  
CONNECT  ACEBK  O   ion        0    C    0    O    0    C  

### This is a temporary parameter file made for residue NH2 ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST NH2        NH2BK 

NATOM    NH2BK      1

IATOM    NH2BK  N      0

ATOMNAME NH2BK    0  N  

CONNECT  NH2BK  N   ion      

