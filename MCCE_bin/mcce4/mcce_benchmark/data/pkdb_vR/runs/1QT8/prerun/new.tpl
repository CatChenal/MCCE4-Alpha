### This is a temporary parameter file made for residue HED ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST HED        HEDBK 

NATOM    HEDBK      8

IATOM    HEDBK  C1     0
IATOM    HEDBK  O1     1
IATOM    HEDBK  C2     2
IATOM    HEDBK  S3     3
IATOM    HEDBK  S4     4
IATOM    HEDBK  C5     5
IATOM    HEDBK  C6     6
IATOM    HEDBK  O6     7

ATOMNAME HEDBK    0  C1 
ATOMNAME HEDBK    1  O1 
ATOMNAME HEDBK    2  C2 
ATOMNAME HEDBK    3  S3 
ATOMNAME HEDBK    4  S4 
ATOMNAME HEDBK    5  C5 
ATOMNAME HEDBK    6  C6 
ATOMNAME HEDBK    7  O6 

CONNECT  HEDBK  C1  ion        0    O1   0    C2 
CONNECT  HEDBK  O1  ion        0    C1 
CONNECT  HEDBK  C2  ion        0    C1   0    S3 
CONNECT  HEDBK  S3  ion        0    C2   0    S4 
CONNECT  HEDBK  S4  ion        0    S3   0    C5 
CONNECT  HEDBK  C5  ion        0    S4   0    C6 
CONNECT  HEDBK  C6  ion        0    C5   0    O6 
CONNECT  HEDBK  O6  ion        0    C6 

