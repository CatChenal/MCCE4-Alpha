### This is a temporary parameter file made for residue _DG ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST _DG        _DGBK 

NATOM    _DGBK      20

IATOM    _DGBK  O5'    0
IATOM    _DGBK  C5'    1
IATOM    _DGBK  C4'    2
IATOM    _DGBK  O4'    3
IATOM    _DGBK  C3'    4
IATOM    _DGBK  O3'    5
IATOM    _DGBK  C2'    6
IATOM    _DGBK  C1'    7
IATOM    _DGBK  N9     8
IATOM    _DGBK  C8     9
IATOM    _DGBK  N7    10
IATOM    _DGBK  C5    11
IATOM    _DGBK  C6    12
IATOM    _DGBK  O6    13
IATOM    _DGBK  N1    14
IATOM    _DGBK  C2    15
IATOM    _DGBK  N2    16
IATOM    _DGBK  N3    17
IATOM    _DGBK  C4    18
IATOM    _DGBK HO5'   19

ATOMNAME _DGBK    0  O5'
ATOMNAME _DGBK    1  C5'
ATOMNAME _DGBK    2  C4'
ATOMNAME _DGBK    3  O4'
ATOMNAME _DGBK    4  C3'
ATOMNAME _DGBK    5  O3'
ATOMNAME _DGBK    6  C2'
ATOMNAME _DGBK    7  C1'
ATOMNAME _DGBK    8  N9 
ATOMNAME _DGBK    9  C8 
ATOMNAME _DGBK   10  N7 
ATOMNAME _DGBK   11  C5 
ATOMNAME _DGBK   12  C6 
ATOMNAME _DGBK   13  O6 
ATOMNAME _DGBK   14  N1 
ATOMNAME _DGBK   15  C2 
ATOMNAME _DGBK   16  N2 
ATOMNAME _DGBK   17  N3 
ATOMNAME _DGBK   18  C4 
ATOMNAME _DGBK   19 HO5'

CONNECT  _DGBK  O5' ion        0    C5'  0   HO5'
CONNECT  _DGBK  C5' ion        0    O5'  0    C4'  0   HO5'
CONNECT  _DGBK  C4' ion        0    C5'  0    O4'  0    C3'
CONNECT  _DGBK  O4' ion        0    C4'  0    C1'
CONNECT  _DGBK  C3' ion        0    C4'  0    O3'  0    C2'
CONNECT  _DGBK  O3' ion        0    C3'
CONNECT  _DGBK  C2' ion        0    C3'  0    C1'
CONNECT  _DGBK  C1' ion        0    O4'  0    C2'  0    N9 
CONNECT  _DGBK  N9  ion        0    C1'  0    C8   0    C5   0    C4 
CONNECT  _DGBK  C8  ion        0    N9   0    N7   0    C5 
CONNECT  _DGBK  N7  ion        0    C8   0    C5 
CONNECT  _DGBK  C5  ion        0    N9   0    C8   0    N7   0    C6   0    C4 
CONNECT  _DGBK  C6  ion        0    C5   0    O6   0    N1 
CONNECT  _DGBK  O6  ion        0    C6 
CONNECT  _DGBK  N1  ion        0    C6   0    C2 
CONNECT  _DGBK  C2  ion        0    N1   0    N2   0    N3 
CONNECT  _DGBK  N2  ion        0    C2 
CONNECT  _DGBK  N3  ion        0    C2   0    C4 
CONNECT  _DGBK  C4  ion        0    N9   0    C5   0    N3 
CONNECT  _DGBK HO5' ion        0    O5'  0    C5'

### This is a temporary parameter file made for residue _DC ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST _DC        _DCBK 

NATOM    _DCBK      19

IATOM    _DCBK  P      0
IATOM    _DCBK  OP1    1
IATOM    _DCBK  OP2    2
IATOM    _DCBK  O5'    3
IATOM    _DCBK  C5'    4
IATOM    _DCBK  C4'    5
IATOM    _DCBK  O4'    6
IATOM    _DCBK  C3'    7
IATOM    _DCBK  O3'    8
IATOM    _DCBK  C2'    9
IATOM    _DCBK  C1'   10
IATOM    _DCBK  N1    11
IATOM    _DCBK  C2    12
IATOM    _DCBK  O2    13
IATOM    _DCBK  N3    14
IATOM    _DCBK  C4    15
IATOM    _DCBK  N4    16
IATOM    _DCBK  C5    17
IATOM    _DCBK  C6    18

ATOMNAME _DCBK    0  P  
ATOMNAME _DCBK    1  OP1
ATOMNAME _DCBK    2  OP2
ATOMNAME _DCBK    3  O5'
ATOMNAME _DCBK    4  C5'
ATOMNAME _DCBK    5  C4'
ATOMNAME _DCBK    6  O4'
ATOMNAME _DCBK    7  C3'
ATOMNAME _DCBK    8  O3'
ATOMNAME _DCBK    9  C2'
ATOMNAME _DCBK   10  C1'
ATOMNAME _DCBK   11  N1 
ATOMNAME _DCBK   12  C2 
ATOMNAME _DCBK   13  O2 
ATOMNAME _DCBK   14  N3 
ATOMNAME _DCBK   15  C4 
ATOMNAME _DCBK   16  N4 
ATOMNAME _DCBK   17  C5 
ATOMNAME _DCBK   18  C6 

CONNECT  _DCBK  P   ion        0    OP1  0    OP2  0    O5'
CONNECT  _DCBK  OP1 ion        0    P  
CONNECT  _DCBK  OP2 ion        0    P  
CONNECT  _DCBK  O5' ion        0    P    0    C5'
CONNECT  _DCBK  C5' ion        0    O5'  0    C4'
CONNECT  _DCBK  C4' ion        0    C5'  0    O4'  0    C3'
CONNECT  _DCBK  O4' ion        0    C4'  0    C1'
CONNECT  _DCBK  C3' ion        0    C4'  0    O3'  0    C2'
CONNECT  _DCBK  O3' ion        0    C3'
CONNECT  _DCBK  C2' ion        0    C3'  0    C1'
CONNECT  _DCBK  C1' ion        0    O4'  0    C2'  0    N1 
CONNECT  _DCBK  N1  ion        0    C1'  0    C2   0    C6 
CONNECT  _DCBK  C2  ion        0    N1   0    O2   0    N3 
CONNECT  _DCBK  O2  ion        0    C2 
CONNECT  _DCBK  N3  ion        0    C2   0    C4 
CONNECT  _DCBK  C4  ion        0    N3   0    N4   0    C5 
CONNECT  _DCBK  N4  ion        0    C4 
CONNECT  _DCBK  C5  ion        0    C4   0    C6 
CONNECT  _DCBK  C6  ion        0    N1   0    C5 

### This is a temporary parameter file made for residue _DA ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST _DA        _DABK 

NATOM    _DABK      21

IATOM    _DABK  P      0
IATOM    _DABK  OP1    1
IATOM    _DABK  OP2    2
IATOM    _DABK  O5'    3
IATOM    _DABK  C5'    4
IATOM    _DABK  C4'    5
IATOM    _DABK  O4'    6
IATOM    _DABK  C3'    7
IATOM    _DABK  O3'    8
IATOM    _DABK  C2'    9
IATOM    _DABK  C1'   10
IATOM    _DABK  N9    11
IATOM    _DABK  C8    12
IATOM    _DABK  N7    13
IATOM    _DABK  C5    14
IATOM    _DABK  C6    15
IATOM    _DABK  N6    16
IATOM    _DABK  N1    17
IATOM    _DABK  C2    18
IATOM    _DABK  N3    19
IATOM    _DABK  C4    20

ATOMNAME _DABK    0  P  
ATOMNAME _DABK    1  OP1
ATOMNAME _DABK    2  OP2
ATOMNAME _DABK    3  O5'
ATOMNAME _DABK    4  C5'
ATOMNAME _DABK    5  C4'
ATOMNAME _DABK    6  O4'
ATOMNAME _DABK    7  C3'
ATOMNAME _DABK    8  O3'
ATOMNAME _DABK    9  C2'
ATOMNAME _DABK   10  C1'
ATOMNAME _DABK   11  N9 
ATOMNAME _DABK   12  C8 
ATOMNAME _DABK   13  N7 
ATOMNAME _DABK   14  C5 
ATOMNAME _DABK   15  C6 
ATOMNAME _DABK   16  N6 
ATOMNAME _DABK   17  N1 
ATOMNAME _DABK   18  C2 
ATOMNAME _DABK   19  N3 
ATOMNAME _DABK   20  C4 

CONNECT  _DABK  P   ion        0    OP1  0    OP2  0    O5'
CONNECT  _DABK  OP1 ion        0    P  
CONNECT  _DABK  OP2 ion        0    P  
CONNECT  _DABK  O5' ion        0    P    0    C5'
CONNECT  _DABK  C5' ion        0    O5'  0    C4'
CONNECT  _DABK  C4' ion        0    C5'  0    O4'  0    C3'
CONNECT  _DABK  O4' ion        0    C4'  0    C1'
CONNECT  _DABK  C3' ion        0    C4'  0    O3'  0    C2'
CONNECT  _DABK  O3' ion        0    C3'
CONNECT  _DABK  C2' ion        0    C3'  0    C1'
CONNECT  _DABK  C1' ion        0    O4'  0    C2'  0    N9 
CONNECT  _DABK  N9  ion        0    C1'  0    C8   0    C4 
CONNECT  _DABK  C8  ion        0    N9   0    N7   0    C5   0    C4 
CONNECT  _DABK  N7  ion        0    C8   0    C5 
CONNECT  _DABK  C5  ion        0    C8   0    N7   0    C6   0    C4 
CONNECT  _DABK  C6  ion        0    C5   0    N6   0    N1 
CONNECT  _DABK  N6  ion        0    C6 
CONNECT  _DABK  N1  ion        0    C6   0    C2 
CONNECT  _DABK  C2  ion        0    N1   0    N3 
CONNECT  _DABK  N3  ion        0    C2   0    C4 
CONNECT  _DABK  C4  ion        0    N9   0    C8   0    C5   0    N3 

### This is a temporary parameter file made for residue _DT ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST _DT        _DTBK 

NATOM    _DTBK      20

IATOM    _DTBK  P      0
IATOM    _DTBK  OP1    1
IATOM    _DTBK  OP2    2
IATOM    _DTBK  O5'    3
IATOM    _DTBK  C5'    4
IATOM    _DTBK  C4'    5
IATOM    _DTBK  O4'    6
IATOM    _DTBK  C3'    7
IATOM    _DTBK  O3'    8
IATOM    _DTBK  C2'    9
IATOM    _DTBK  C1'   10
IATOM    _DTBK  N1    11
IATOM    _DTBK  C2    12
IATOM    _DTBK  O2    13
IATOM    _DTBK  N3    14
IATOM    _DTBK  C4    15
IATOM    _DTBK  O4    16
IATOM    _DTBK  C5    17
IATOM    _DTBK  C7    18
IATOM    _DTBK  C6    19

ATOMNAME _DTBK    0  P  
ATOMNAME _DTBK    1  OP1
ATOMNAME _DTBK    2  OP2
ATOMNAME _DTBK    3  O5'
ATOMNAME _DTBK    4  C5'
ATOMNAME _DTBK    5  C4'
ATOMNAME _DTBK    6  O4'
ATOMNAME _DTBK    7  C3'
ATOMNAME _DTBK    8  O3'
ATOMNAME _DTBK    9  C2'
ATOMNAME _DTBK   10  C1'
ATOMNAME _DTBK   11  N1 
ATOMNAME _DTBK   12  C2 
ATOMNAME _DTBK   13  O2 
ATOMNAME _DTBK   14  N3 
ATOMNAME _DTBK   15  C4 
ATOMNAME _DTBK   16  O4 
ATOMNAME _DTBK   17  C5 
ATOMNAME _DTBK   18  C7 
ATOMNAME _DTBK   19  C6 

CONNECT  _DTBK  P   ion        0    OP1  0    OP2  0    O5'
CONNECT  _DTBK  OP1 ion        0    P  
CONNECT  _DTBK  OP2 ion        0    P  
CONNECT  _DTBK  O5' ion        0    P    0    C5'
CONNECT  _DTBK  C5' ion        0    O5'  0    C4'
CONNECT  _DTBK  C4' ion        0    C5'  0    O4'  0    C3'
CONNECT  _DTBK  O4' ion        0    C4'  0    C1'
CONNECT  _DTBK  C3' ion        0    C4'  0    O3'  0    C2'
CONNECT  _DTBK  O3' ion        0    C3'
CONNECT  _DTBK  C2' ion        0    C3'  0    C1'
CONNECT  _DTBK  C1' ion        0    O4'  0    C2'  0    N1 
CONNECT  _DTBK  N1  ion        0    C1'  0    C2   0    C6 
CONNECT  _DTBK  C2  ion        0    N1   0    O2   0    N3 
CONNECT  _DTBK  O2  ion        0    C2 
CONNECT  _DTBK  N3  ion        0    C2   0    C4 
CONNECT  _DTBK  C4  ion        0    N3   0    O4   0    C5 
CONNECT  _DTBK  O4  ion        0    C4 
CONNECT  _DTBK  C5  ion        0    C4   0    C7   0    C6 
CONNECT  _DTBK  C7  ion        0    C5 
CONNECT  _DTBK  C6  ion        0    N1   0    C5 

