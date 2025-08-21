### This is a temporary parameter file made for residue 2GP ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST 2GP        2GPBK 

NATOM    2GPBK      24

IATOM    2GPBK  P      0
IATOM    2GPBK  O1P    1
IATOM    2GPBK  O2P    2
IATOM    2GPBK  O3P    3
IATOM    2GPBK  O5'    4
IATOM    2GPBK  C5'    5
IATOM    2GPBK  C4'    6
IATOM    2GPBK  O4'    7
IATOM    2GPBK  C3'    8
IATOM    2GPBK  O3'    9
IATOM    2GPBK  C2'   10
IATOM    2GPBK  O2'   11
IATOM    2GPBK  C1'   12
IATOM    2GPBK  N9    13
IATOM    2GPBK  C8    14
IATOM    2GPBK  N7    15
IATOM    2GPBK  C5    16
IATOM    2GPBK  C6    17
IATOM    2GPBK  O6    18
IATOM    2GPBK  N1    19
IATOM    2GPBK  C2    20
IATOM    2GPBK  N2    21
IATOM    2GPBK  N3    22
IATOM    2GPBK  C4    23

ATOMNAME 2GPBK    0  P  
ATOMNAME 2GPBK    1  O1P
ATOMNAME 2GPBK    2  O2P
ATOMNAME 2GPBK    3  O3P
ATOMNAME 2GPBK    4  O5'
ATOMNAME 2GPBK    5  C5'
ATOMNAME 2GPBK    6  C4'
ATOMNAME 2GPBK    7  O4'
ATOMNAME 2GPBK    8  C3'
ATOMNAME 2GPBK    9  O3'
ATOMNAME 2GPBK   10  C2'
ATOMNAME 2GPBK   11  O2'
ATOMNAME 2GPBK   12  C1'
ATOMNAME 2GPBK   13  N9 
ATOMNAME 2GPBK   14  C8 
ATOMNAME 2GPBK   15  N7 
ATOMNAME 2GPBK   16  C5 
ATOMNAME 2GPBK   17  C6 
ATOMNAME 2GPBK   18  O6 
ATOMNAME 2GPBK   19  N1 
ATOMNAME 2GPBK   20  C2 
ATOMNAME 2GPBK   21  N2 
ATOMNAME 2GPBK   22  N3 
ATOMNAME 2GPBK   23  C4 

CONNECT  2GPBK  P   ion        0    O1P  0    O2P  0    O3P  0    O2'
CONNECT  2GPBK  O1P ion        0    P  
CONNECT  2GPBK  O2P ion        0    P  
CONNECT  2GPBK  O3P ion        0    P  
CONNECT  2GPBK  O5' ion        0    C5'
CONNECT  2GPBK  C5' ion        0    O5'  0    C4'
CONNECT  2GPBK  C4' ion        0    C5'  0    O4'  0    C3'
CONNECT  2GPBK  O4' ion        0    C4'  0    C1'
CONNECT  2GPBK  C3' ion        0    C4'  0    O3'  0    C2'
CONNECT  2GPBK  O3' ion        0    C3'
CONNECT  2GPBK  C2' ion        0    C3'  0    O2'  0    C1'
CONNECT  2GPBK  O2' ion        0    P    0    C2'
CONNECT  2GPBK  C1' ion        0    O4'  0    C2'  0    N9 
CONNECT  2GPBK  N9  ion        0    C1'  0    C8   0    C5   0    C4 
CONNECT  2GPBK  C8  ion        0    N9   0    N7   0    C5   0    C4 
CONNECT  2GPBK  N7  ion        0    C8   0    C5 
CONNECT  2GPBK  C5  ion        0    N9   0    C8   0    N7   0    C6   0    C4 
CONNECT  2GPBK  C6  ion        0    C5   0    O6   0    N1 
CONNECT  2GPBK  O6  ion        0    C6 
CONNECT  2GPBK  N1  ion        0    C6   0    C2 
CONNECT  2GPBK  C2  ion        0    N1   0    N2   0    N3 
CONNECT  2GPBK  N2  ion        0    C2 
CONNECT  2GPBK  N3  ion        0    C2   0    C4 
CONNECT  2GPBK  C4  ion        0    N9   0    C8   0    C5   0    N3 

### This is a temporary parameter file made for residue _CA ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST _CA        _CABK 

NATOM    _CABK      1

IATOM    _CABK CA      0

ATOMNAME _CABK    0 CA  

CONNECT  _CABK CA   ion      

