### This is a temporary parameter file made for residue PCA ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST PCA        PCABK 

NATOM    PCABK      8

IATOM    PCABK  N_     0
IATOM    PCABK  CA     1
IATOM    PCABK  CB     2
IATOM    PCABK  CG     3
IATOM    PCABK  CD     4
IATOM    PCABK  OE     5
IATOM    PCABK  C_     6
IATOM    PCABK  O_     7

ATOMNAME PCABK    0  N_ 
ATOMNAME PCABK    1  CA 
ATOMNAME PCABK    2  CB 
ATOMNAME PCABK    3  CG 
ATOMNAME PCABK    4  CD 
ATOMNAME PCABK    5  OE 
ATOMNAME PCABK    6  C_ 
ATOMNAME PCABK    7  O_ 

CONNECT  PCABK  N_  ion        0    CA   0    CD 
CONNECT  PCABK  CA  ion        0    N_   0    CB   0    C_ 
CONNECT  PCABK  CB  ion        0    CA   0    CG 
CONNECT  PCABK  CG  ion        0    CB   0    CD 
CONNECT  PCABK  CD  ion        0    N_   0    CG   0    OE 
CONNECT  PCABK  OE  ion        0    CD 
CONNECT  PCABK  C_  ion        0    CA   0    O_ 
CONNECT  PCABK  O_  ion        0    C_ 

