### This is a temporary parameter file made for residue GSH ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST GSH        GSHBK 

NATOM    GSHBK      20

IATOM    GSHBK  N1     0
IATOM    GSHBK  CA1    1
IATOM    GSHBK  C1     2
IATOM    GSHBK  O11    3
IATOM    GSHBK  O12    4
IATOM    GSHBK  CB1    5
IATOM    GSHBK  CG1    6
IATOM    GSHBK  CD1    7
IATOM    GSHBK  OE1    8
IATOM    GSHBK  N2     9
IATOM    GSHBK  CA2   10
IATOM    GSHBK  C2    11
IATOM    GSHBK  O2    12
IATOM    GSHBK  CB2   13
IATOM    GSHBK  SG2   14
IATOM    GSHBK  N3    15
IATOM    GSHBK  CA3   16
IATOM    GSHBK  C3    17
IATOM    GSHBK  O31   18
IATOM    GSHBK  O32   19

ATOMNAME GSHBK    0  N1 
ATOMNAME GSHBK    1  CA1
ATOMNAME GSHBK    2  C1 
ATOMNAME GSHBK    3  O11
ATOMNAME GSHBK    4  O12
ATOMNAME GSHBK    5  CB1
ATOMNAME GSHBK    6  CG1
ATOMNAME GSHBK    7  CD1
ATOMNAME GSHBK    8  OE1
ATOMNAME GSHBK    9  N2 
ATOMNAME GSHBK   10  CA2
ATOMNAME GSHBK   11  C2 
ATOMNAME GSHBK   12  O2 
ATOMNAME GSHBK   13  CB2
ATOMNAME GSHBK   14  SG2
ATOMNAME GSHBK   15  N3 
ATOMNAME GSHBK   16  CA3
ATOMNAME GSHBK   17  C3 
ATOMNAME GSHBK   18  O31
ATOMNAME GSHBK   19  O32

CONNECT  GSHBK  N1  ion        0    CA1
CONNECT  GSHBK  CA1 ion        0    N1   0    C1   0    CB1
CONNECT  GSHBK  C1  ion        0    CA1  0    O11  0    O12
CONNECT  GSHBK  O11 ion        0    C1 
CONNECT  GSHBK  O12 ion        0    C1 
CONNECT  GSHBK  CB1 ion        0    CA1  0    CG1
CONNECT  GSHBK  CG1 ion        0    CB1  0    CD1
CONNECT  GSHBK  CD1 ion        0    CG1  0    OE1  0    N2 
CONNECT  GSHBK  OE1 ion        0    CD1
CONNECT  GSHBK  N2  ion        0    CD1  0    CA2
CONNECT  GSHBK  CA2 ion        0    N2   0    C2   0    CB2
CONNECT  GSHBK  C2  ion        0    CA2  0    O2   0    N3 
CONNECT  GSHBK  O2  ion        0    C2 
CONNECT  GSHBK  CB2 ion        0    CA2  0    SG2
CONNECT  GSHBK  SG2 ion        0    CB2
CONNECT  GSHBK  N3  ion        0    C2   0    CA3
CONNECT  GSHBK  CA3 ion        0    N3   0    C3 
CONNECT  GSHBK  C3  ion        0    CA3  0    O31  0    O32
CONNECT  GSHBK  O31 ion        0    C3   0    O32
CONNECT  GSHBK  O32 ion        0    C3   0    O31

