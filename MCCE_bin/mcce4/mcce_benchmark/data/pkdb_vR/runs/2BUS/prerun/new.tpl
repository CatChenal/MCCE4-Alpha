### This is a temporary parameter file made for residue PCA ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST PCA        PCABK 

NATOM    PCABK      16

IATOM    PCABK  N_     0
IATOM    PCABK  CA     1
IATOM    PCABK  N_     2
IATOM    PCABK  CA     3
IATOM    PCABK  CB     4
IATOM    PCABK  CG     5
IATOM    PCABK  CD     6
IATOM    PCABK  OE     7
IATOM    PCABK  CB     8
IATOM    PCABK  CG     9
IATOM    PCABK  CD    10
IATOM    PCABK  OE    11
IATOM    PCABK  C_    12
IATOM    PCABK  O_    13
IATOM    PCABK  C_    14
IATOM    PCABK  O_    15

ATOMNAME PCABK    0  N_ 
ATOMNAME PCABK    1  CA 
ATOMNAME PCABK    2  N_ 
ATOMNAME PCABK    3  CA 
ATOMNAME PCABK    4  CB 
ATOMNAME PCABK    5  CG 
ATOMNAME PCABK    6  CD 
ATOMNAME PCABK    7  OE 
ATOMNAME PCABK    8  CB 
ATOMNAME PCABK    9  CG 
ATOMNAME PCABK   10  CD 
ATOMNAME PCABK   11  OE 
ATOMNAME PCABK   12  C_ 
ATOMNAME PCABK   13  O_ 
ATOMNAME PCABK   14  C_ 
ATOMNAME PCABK   15  O_ 

CONNECT  PCABK  N_  ion        0    CA   0    N_   0    CA   0    CD   0    CD 
CONNECT  PCABK  CA  ion        0    N_   0    N_   0    CA   0    CB   0    CB   0    C_   0    C_ 
CONNECT  PCABK  N_  ion        0    N_   0    CA   0    CA   0    CD   0    CD 
CONNECT  PCABK  CA  ion        0    N_   0    CA   0    N_   0    CB   0    CB   0    C_   0    C_ 
CONNECT  PCABK  CB  ion        0    CA   0    CA   0    CG   0    CB   0    CG 
CONNECT  PCABK  CG  ion        0    CB   0    CD   0    CB   0    CG   0    CD 
CONNECT  PCABK  CD  ion        0    N_   0    N_   0    CG   0    OE   0    CG   0    CD   0    OE 
CONNECT  PCABK  OE  ion        0    CD   0    CD   0    OE 
CONNECT  PCABK  CB  ion        0    CA   0    CA   0    CB   0    CG   0    CG 
CONNECT  PCABK  CG  ion        0    CB   0    CG   0    CD   0    CB   0    CD 
CONNECT  PCABK  CD  ion        0    N_   0    N_   0    CG   0    CD   0    OE   0    CG   0    OE 
CONNECT  PCABK  OE  ion        0    CD   0    OE   0    CD 
CONNECT  PCABK  C_  ion        0    CA   0    CA   0    O_   0    C_   0    O_ 
CONNECT  PCABK  O_  ion        0    C_   0    C_   0    O_ 
CONNECT  PCABK  C_  ion        0    CA   0    CA   0    C_   0    O_   0    O_ 
CONNECT  PCABK  O_  ion        0    C_   0    O_   0    C_ 

