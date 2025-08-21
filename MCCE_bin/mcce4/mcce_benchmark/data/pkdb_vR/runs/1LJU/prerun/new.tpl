### This is a temporary parameter file made for residue CSR ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST CSR        CSRBK 

NATOM    CSRBK      10

IATOM    CSRBK  N      0
IATOM    CSRBK  CA     1
IATOM    CSRBK  CB     2
IATOM    CSRBK  SG     3
IATOM    CSRBK AS      4
IATOM    CSRBK  O1     5
IATOM    CSRBK  O2     6
IATOM    CSRBK  O3     7
IATOM    CSRBK  C      8
IATOM    CSRBK  O      9

ATOMNAME CSRBK    0  N  
ATOMNAME CSRBK    1  CA 
ATOMNAME CSRBK    2  CB 
ATOMNAME CSRBK    3  SG 
ATOMNAME CSRBK    4 AS  
ATOMNAME CSRBK    5  O1 
ATOMNAME CSRBK    6  O2 
ATOMNAME CSRBK    7  O3 
ATOMNAME CSRBK    8  C  
ATOMNAME CSRBK    9  O  

CONNECT  CSRBK  N   ion        0    CA 
CONNECT  CSRBK  CA  ion        0    N    0    CB   0    C  
CONNECT  CSRBK  CB  ion        0    CA   0    SG 
CONNECT  CSRBK  SG  ion        0    CB 
CONNECT  CSRBK AS   ion        0    O1   0    O2   0    O3 
CONNECT  CSRBK  O1  ion        0   AS  
CONNECT  CSRBK  O2  ion        0   AS  
CONNECT  CSRBK  O3  ion        0   AS  
CONNECT  CSRBK  C   ion        0    CA   0    O  
CONNECT  CSRBK  O   ion        0    C  

### This is a temporary parameter file made for residue __K ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST __K        __KBK 

NATOM    __KBK      1

IATOM    __KBK  K      0

ATOMNAME __KBK    0  K  

CONNECT  __KBK  K   ion      

