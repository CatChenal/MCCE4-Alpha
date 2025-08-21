### This is a temporary parameter file made for residue _DC ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST _DC        _DCBK 

NATOM    _DCBK      16

IATOM    _DCBK  O5'    0
IATOM    _DCBK  C5'    1
IATOM    _DCBK  C4'    2
IATOM    _DCBK  O4'    3
IATOM    _DCBK  C3'    4
IATOM    _DCBK  O3'    5
IATOM    _DCBK  C2'    6
IATOM    _DCBK  C1'    7
IATOM    _DCBK  N1     8
IATOM    _DCBK  C2     9
IATOM    _DCBK  O2    10
IATOM    _DCBK  N3    11
IATOM    _DCBK  C4    12
IATOM    _DCBK  N4    13
IATOM    _DCBK  C5    14
IATOM    _DCBK  C6    15

ATOMNAME _DCBK    0  O5'
ATOMNAME _DCBK    1  C5'
ATOMNAME _DCBK    2  C4'
ATOMNAME _DCBK    3  O4'
ATOMNAME _DCBK    4  C3'
ATOMNAME _DCBK    5  O3'
ATOMNAME _DCBK    6  C2'
ATOMNAME _DCBK    7  C1'
ATOMNAME _DCBK    8  N1 
ATOMNAME _DCBK    9  C2 
ATOMNAME _DCBK   10  O2 
ATOMNAME _DCBK   11  N3 
ATOMNAME _DCBK   12  C4 
ATOMNAME _DCBK   13  N4 
ATOMNAME _DCBK   14  C5 
ATOMNAME _DCBK   15  C6 

CONNECT  _DCBK  O5' ion        0    C5'
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

### This is a temporary parameter file made for residue _DG ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST _DG        _DGBK 

NATOM    _DGBK      22

IATOM    _DGBK  P      0
IATOM    _DGBK  OP1    1
IATOM    _DGBK  OP2    2
IATOM    _DGBK  O5'    3
IATOM    _DGBK  C5'    4
IATOM    _DGBK  C4'    5
IATOM    _DGBK  O4'    6
IATOM    _DGBK  C3'    7
IATOM    _DGBK  O3'    8
IATOM    _DGBK  C2'    9
IATOM    _DGBK  C1'   10
IATOM    _DGBK  N9    11
IATOM    _DGBK  C8    12
IATOM    _DGBK  N7    13
IATOM    _DGBK  C5    14
IATOM    _DGBK  C6    15
IATOM    _DGBK  O6    16
IATOM    _DGBK  N1    17
IATOM    _DGBK  C2    18
IATOM    _DGBK  N2    19
IATOM    _DGBK  N3    20
IATOM    _DGBK  C4    21

ATOMNAME _DGBK    0  P  
ATOMNAME _DGBK    1  OP1
ATOMNAME _DGBK    2  OP2
ATOMNAME _DGBK    3  O5'
ATOMNAME _DGBK    4  C5'
ATOMNAME _DGBK    5  C4'
ATOMNAME _DGBK    6  O4'
ATOMNAME _DGBK    7  C3'
ATOMNAME _DGBK    8  O3'
ATOMNAME _DGBK    9  C2'
ATOMNAME _DGBK   10  C1'
ATOMNAME _DGBK   11  N9 
ATOMNAME _DGBK   12  C8 
ATOMNAME _DGBK   13  N7 
ATOMNAME _DGBK   14  C5 
ATOMNAME _DGBK   15  C6 
ATOMNAME _DGBK   16  O6 
ATOMNAME _DGBK   17  N1 
ATOMNAME _DGBK   18  C2 
ATOMNAME _DGBK   19  N2 
ATOMNAME _DGBK   20  N3 
ATOMNAME _DGBK   21  C4 

CONNECT  _DGBK  P   ion        0    OP1  0    OP2  0    O5'
CONNECT  _DGBK  OP1 ion        0    P  
CONNECT  _DGBK  OP2 ion        0    P  
CONNECT  _DGBK  O5' ion        0    P    0    C5'
CONNECT  _DGBK  C5' ion        0    O5'  0    C4'
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
CONNECT  _DABK  C2  ion        0    N1   0    N3   0    C4 
CONNECT  _DABK  N3  ion        0    C2   0    C4 
CONNECT  _DABK  C4  ion        0    N9   0    C8   0    C5   0    C2   0    N3 

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

### This is a temporary parameter file made for residue 2DT ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST 2DT        2DTBK 

NATOM    2DTBK      19

IATOM    2DTBK  P      0
IATOM    2DTBK  OP1    1
IATOM    2DTBK  OP2    2
IATOM    2DTBK  O5'    3
IATOM    2DTBK  N1     4
IATOM    2DTBK  C6     5
IATOM    2DTBK  C2     6
IATOM    2DTBK  O2     7
IATOM    2DTBK  N3     8
IATOM    2DTBK  C4     9
IATOM    2DTBK  O4    10
IATOM    2DTBK  C5    11
IATOM    2DTBK  C5M   12
IATOM    2DTBK  C2'   13
IATOM    2DTBK  C5'   14
IATOM    2DTBK  C4'   15
IATOM    2DTBK  O4'   16
IATOM    2DTBK  C1'   17
IATOM    2DTBK  C3'   18

ATOMNAME 2DTBK    0  P  
ATOMNAME 2DTBK    1  OP1
ATOMNAME 2DTBK    2  OP2
ATOMNAME 2DTBK    3  O5'
ATOMNAME 2DTBK    4  N1 
ATOMNAME 2DTBK    5  C6 
ATOMNAME 2DTBK    6  C2 
ATOMNAME 2DTBK    7  O2 
ATOMNAME 2DTBK    8  N3 
ATOMNAME 2DTBK    9  C4 
ATOMNAME 2DTBK   10  O4 
ATOMNAME 2DTBK   11  C5 
ATOMNAME 2DTBK   12  C5M
ATOMNAME 2DTBK   13  C2'
ATOMNAME 2DTBK   14  C5'
ATOMNAME 2DTBK   15  C4'
ATOMNAME 2DTBK   16  O4'
ATOMNAME 2DTBK   17  C1'
ATOMNAME 2DTBK   18  C3'

CONNECT  2DTBK  P   ion        0    OP1  0    OP2  0    O5'
CONNECT  2DTBK  OP1 ion        0    P  
CONNECT  2DTBK  OP2 ion        0    P  
CONNECT  2DTBK  O5' ion        0    P    0    C5'
CONNECT  2DTBK  N1  ion        0    C6   0    C2   0    C1'
CONNECT  2DTBK  C6  ion        0    N1   0    C5 
CONNECT  2DTBK  C2  ion        0    N1   0    O2   0    N3 
CONNECT  2DTBK  O2  ion        0    C2 
CONNECT  2DTBK  N3  ion        0    C2   0    C4 
CONNECT  2DTBK  C4  ion        0    N3   0    O4   0    C5 
CONNECT  2DTBK  O4  ion        0    C4 
CONNECT  2DTBK  C5  ion        0    C6   0    C4   0    C5M
CONNECT  2DTBK  C5M ion        0    C5 
CONNECT  2DTBK  C2' ion        0    C1'  0    C3'
CONNECT  2DTBK  C5' ion        0    O5'  0    C4'
CONNECT  2DTBK  C4' ion        0    C5'  0    O4'  0    C3'
CONNECT  2DTBK  O4' ion        0    C4'  0    C1'
CONNECT  2DTBK  C1' ion        0    N1   0    C2'  0    O4'
CONNECT  2DTBK  C3' ion        0    C2'  0    C4'

### This is a temporary parameter file made for residue EDO ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST EDO        EDOBK 

NATOM    EDOBK      4

IATOM    EDOBK  C1     0
IATOM    EDOBK  O1     1
IATOM    EDOBK  C2     2
IATOM    EDOBK  O2     3

ATOMNAME EDOBK    0  C1 
ATOMNAME EDOBK    1  O1 
ATOMNAME EDOBK    2  C2 
ATOMNAME EDOBK    3  O2 

CONNECT  EDOBK  C1  ion        0    O1   0    C2 
CONNECT  EDOBK  O1  ion        0    C1 
CONNECT  EDOBK  C2  ion        0    C1   0    O2 
CONNECT  EDOBK  O2  ion        0    C2 

### This is a temporary parameter file made for residue _MG ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST _MG        _MGBK 

NATOM    _MGBK      1

IATOM    _MGBK MG      0

ATOMNAME _MGBK    0 MG  

CONNECT  _MGBK MG   ion      

### This is a temporary parameter file made for residue D3T ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST D3T        D3TBK 

NATOM    D3TBK      28

IATOM    D3TBK  PA     0
IATOM    D3TBK  O1A    1
IATOM    D3TBK  O2A    2
IATOM    D3TBK  O5'    3
IATOM    D3TBK  C5'    4
IATOM    D3TBK  C4'    5
IATOM    D3TBK  O4'    6
IATOM    D3TBK  C1'    7
IATOM    D3TBK  N1     8
IATOM    D3TBK  C6     9
IATOM    D3TBK  C2    10
IATOM    D3TBK  O2    11
IATOM    D3TBK  N3    12
IATOM    D3TBK  C4    13
IATOM    D3TBK  O4    14
IATOM    D3TBK  C5    15
IATOM    D3TBK  C5M   16
IATOM    D3TBK  C2'   17
IATOM    D3TBK  C3'   18
IATOM    D3TBK  O3A   19
IATOM    D3TBK  PB    20
IATOM    D3TBK  O1B   21
IATOM    D3TBK  O2B   22
IATOM    D3TBK  O3B   23
IATOM    D3TBK  PG    24
IATOM    D3TBK  O1G   25
IATOM    D3TBK  O2G   26
IATOM    D3TBK  O3G   27

ATOMNAME D3TBK    0  PA 
ATOMNAME D3TBK    1  O1A
ATOMNAME D3TBK    2  O2A
ATOMNAME D3TBK    3  O5'
ATOMNAME D3TBK    4  C5'
ATOMNAME D3TBK    5  C4'
ATOMNAME D3TBK    6  O4'
ATOMNAME D3TBK    7  C1'
ATOMNAME D3TBK    8  N1 
ATOMNAME D3TBK    9  C6 
ATOMNAME D3TBK   10  C2 
ATOMNAME D3TBK   11  O2 
ATOMNAME D3TBK   12  N3 
ATOMNAME D3TBK   13  C4 
ATOMNAME D3TBK   14  O4 
ATOMNAME D3TBK   15  C5 
ATOMNAME D3TBK   16  C5M
ATOMNAME D3TBK   17  C2'
ATOMNAME D3TBK   18  C3'
ATOMNAME D3TBK   19  O3A
ATOMNAME D3TBK   20  PB 
ATOMNAME D3TBK   21  O1B
ATOMNAME D3TBK   22  O2B
ATOMNAME D3TBK   23  O3B
ATOMNAME D3TBK   24  PG 
ATOMNAME D3TBK   25  O1G
ATOMNAME D3TBK   26  O2G
ATOMNAME D3TBK   27  O3G

CONNECT  D3TBK  PA  ion        0    O1A  0    O2A  0    O5'  0    O3A
CONNECT  D3TBK  O1A ion        0    PA 
CONNECT  D3TBK  O2A ion        0    PA 
CONNECT  D3TBK  O5' ion        0    PA   0    C5'
CONNECT  D3TBK  C5' ion        0    O5'  0    C4'
CONNECT  D3TBK  C4' ion        0    C5'  0    O4'  0    C3'
CONNECT  D3TBK  O4' ion        0    C4'  0    C1'
CONNECT  D3TBK  C1' ion        0    O4'  0    N1   0    C2'
CONNECT  D3TBK  N1  ion        0    C1'  0    C6   0    C2 
CONNECT  D3TBK  C6  ion        0    N1   0    C5 
CONNECT  D3TBK  C2  ion        0    N1   0    O2   0    N3 
CONNECT  D3TBK  O2  ion        0    C2 
CONNECT  D3TBK  N3  ion        0    C2   0    C4 
CONNECT  D3TBK  C4  ion        0    N3   0    O4   0    C5 
CONNECT  D3TBK  O4  ion        0    C4 
CONNECT  D3TBK  C5  ion        0    C6   0    C4   0    C5M
CONNECT  D3TBK  C5M ion        0    C5 
CONNECT  D3TBK  C2' ion        0    C1'  0    C3'
CONNECT  D3TBK  C3' ion        0    C4'  0    C2'
CONNECT  D3TBK  O3A ion        0    PA   0    PB 
CONNECT  D3TBK  PB  ion        0    O3A  0    O1B  0    O2B  0    O3B
CONNECT  D3TBK  O1B ion        0    PB 
CONNECT  D3TBK  O2B ion        0    PB 
CONNECT  D3TBK  O3B ion        0    PB   0    PG 
CONNECT  D3TBK  PG  ion        0    O3B  0    O1G  0    O2G  0    O3G
CONNECT  D3TBK  O1G ion        0    PG 
CONNECT  D3TBK  O2G ion        0    PG 
CONNECT  D3TBK  O3G ion        0    PG 

