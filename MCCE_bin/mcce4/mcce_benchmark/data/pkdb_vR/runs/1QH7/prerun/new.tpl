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

### This is a temporary parameter file made for residue XYP ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST XYP        XYPBK 

NATOM    XYPBK      10

IATOM    XYPBK  O1     0
IATOM    XYPBK  C1     1
IATOM    XYPBK  C2     2
IATOM    XYPBK  C3     3
IATOM    XYPBK  C4     4
IATOM    XYPBK  C5     5
IATOM    XYPBK  O2     6
IATOM    XYPBK  O3     7
IATOM    XYPBK  O4     8
IATOM    XYPBK  O5     9

ATOMNAME XYPBK    0  O1 
ATOMNAME XYPBK    1  C1 
ATOMNAME XYPBK    2  C2 
ATOMNAME XYPBK    3  C3 
ATOMNAME XYPBK    4  C4 
ATOMNAME XYPBK    5  C5 
ATOMNAME XYPBK    6  O2 
ATOMNAME XYPBK    7  O3 
ATOMNAME XYPBK    8  O4 
ATOMNAME XYPBK    9  O5 

CONNECT  XYPBK  O1  ion        0    C1 
CONNECT  XYPBK  C1  ion        0    O1   0    C2   0    O5 
CONNECT  XYPBK  C2  ion        0    C1   0    C3   0    O2 
CONNECT  XYPBK  C3  ion        0    C2   0    C4   0    O3 
CONNECT  XYPBK  C4  ion        0    C3   0    C5   0    O4 
CONNECT  XYPBK  C5  ion        0    C4   0    O5 
CONNECT  XYPBK  O2  ion        0    C2 
CONNECT  XYPBK  O3  ion        0    C3 
CONNECT  XYPBK  O4  ion        0    C4 
CONNECT  XYPBK  O5  ion        0    C1   0    C5 

