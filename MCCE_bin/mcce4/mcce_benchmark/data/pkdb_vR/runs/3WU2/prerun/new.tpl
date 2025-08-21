### This is a temporary parameter file made for residue HSK ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST HSK        HSKBK 

NATOM    HSKBK      10

IATOM    HSKBK  CA     0
IATOM    HSKBK  N      1
IATOM    HSKBK  CB     2
IATOM    HSKBK  CG     3
IATOM    HSKBK  CD2    4
IATOM    HSKBK  NE2    5
IATOM    HSKBK  CE1    6
IATOM    HSKBK  ND1    7
IATOM    HSKBK  C      8
IATOM    HSKBK  OXT    9

ATOMNAME HSKBK    0  CA 
ATOMNAME HSKBK    1  N  
ATOMNAME HSKBK    2  CB 
ATOMNAME HSKBK    3  CG 
ATOMNAME HSKBK    4  CD2
ATOMNAME HSKBK    5  NE2
ATOMNAME HSKBK    6  CE1
ATOMNAME HSKBK    7  ND1
ATOMNAME HSKBK    8  C  
ATOMNAME HSKBK    9  OXT

CONNECT  HSKBK  CA  ion        0    N    0    CB   0    C  
CONNECT  HSKBK  N   ion        0    CA 
CONNECT  HSKBK  CB  ion        0    CA   0    CG 
CONNECT  HSKBK  CG  ion        0    CB   0    CD2  0    NE2  0    CE1  0    ND1
CONNECT  HSKBK  CD2 ion        0    CG   0    NE2
CONNECT  HSKBK  NE2 ion        0    CG   0    CD2  0    CE1  0    ND1
CONNECT  HSKBK  CE1 ion        0    CG   0    NE2  0    ND1
CONNECT  HSKBK  ND1 ion        0    CG   0    NE2  0    CE1
CONNECT  HSKBK  C   ion        0    CA   0    OXT
CONNECT  HSKBK  OXT ion        0    C  

### This is a temporary parameter file made for residue FME ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST FME        FMEBK 

NATOM    FMEBK      10

IATOM    FMEBK  N      0
IATOM    FMEBK  CA     1
IATOM    FMEBK  CN     2
IATOM    FMEBK  O1     3
IATOM    FMEBK  CB     4
IATOM    FMEBK  CG     5
IATOM    FMEBK  SD     6
IATOM    FMEBK  CE     7
IATOM    FMEBK  C      8
IATOM    FMEBK  O      9

ATOMNAME FMEBK    0  N  
ATOMNAME FMEBK    1  CA 
ATOMNAME FMEBK    2  CN 
ATOMNAME FMEBK    3  O1 
ATOMNAME FMEBK    4  CB 
ATOMNAME FMEBK    5  CG 
ATOMNAME FMEBK    6  SD 
ATOMNAME FMEBK    7  CE 
ATOMNAME FMEBK    8  C  
ATOMNAME FMEBK    9  O  

CONNECT  FMEBK  N   ion        0    CA   0    CN 
CONNECT  FMEBK  CA  ion        0    N    0    CB   0    C  
CONNECT  FMEBK  CN  ion        0    N    0    O1 
CONNECT  FMEBK  O1  ion        0    CN 
CONNECT  FMEBK  CB  ion        0    CA   0    CG 
CONNECT  FMEBK  CG  ion        0    CB   0    SD 
CONNECT  FMEBK  SD  ion        0    CG   0    CE 
CONNECT  FMEBK  CE  ion        0    SD 
CONNECT  FMEBK  C   ion        0    CA   0    O  
CONNECT  FMEBK  O   ion        0    C  

### This is a temporary parameter file made for residue OEC ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST OEC        OECBK 

NATOM    OECBK      10

IATOM    OECBK  O1     0
IATOM    OECBK CA1     1
IATOM    OECBK MN1     2
IATOM    OECBK  O2     3
IATOM    OECBK MN2     4
IATOM    OECBK  O3     5
IATOM    OECBK MN3     6
IATOM    OECBK  O4     7
IATOM    OECBK MN4     8
IATOM    OECBK  O5     9

ATOMNAME OECBK    0  O1 
ATOMNAME OECBK    1 CA1 
ATOMNAME OECBK    2 MN1 
ATOMNAME OECBK    3  O2 
ATOMNAME OECBK    4 MN2 
ATOMNAME OECBK    5  O3 
ATOMNAME OECBK    6 MN3 
ATOMNAME OECBK    7  O4 
ATOMNAME OECBK    8 MN4 
ATOMNAME OECBK    9  O5 

CONNECT  OECBK  O1  ion        0   MN1   0   MN2 
CONNECT  OECBK CA1  ion      
CONNECT  OECBK MN1  ion        0    O1   0    O3 
CONNECT  OECBK  O2  ion        0   MN2   0   MN3 
CONNECT  OECBK MN2  ion        0    O1   0    O2   0    O3 
CONNECT  OECBK  O3  ion        0   MN1   0   MN2 
CONNECT  OECBK MN3  ion        0    O2   0    O4 
CONNECT  OECBK  O4  ion        0   MN3   0   MN4 
CONNECT  OECBK MN4  ion        0    O4 
CONNECT  OECBK  O5  ion      

### This is a temporary parameter file made for residue CT1 ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST CT1        CT1BK 

NATOM    CT1BK      5

IATOM    CT1BK  C1     0
IATOM    CT1BK  C2     1
IATOM    CT1BK  C3     2
IATOM    CT1BK  C4     3
IATOM    CT1BK  C5     4

ATOMNAME CT1BK    0  C1 
ATOMNAME CT1BK    1  C2 
ATOMNAME CT1BK    2  C3 
ATOMNAME CT1BK    3  C4 
ATOMNAME CT1BK    4  C5 

CONNECT  CT1BK  C1  ion        0    C2 
CONNECT  CT1BK  C2  ion        0    C1   0    C3 
CONNECT  CT1BK  C3  ion        0    C2   0    C4   0    C5 
CONNECT  CT1BK  C4  ion        0    C3 
CONNECT  CT1BK  C5  ion        0    C3 

### This is a temporary parameter file made for residue CT2 ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST CT2        CT2BK 

NATOM    CT2BK      5

IATOM    CT2BK  C6     0
IATOM    CT2BK  C7     1
IATOM    CT2BK  C8     2
IATOM    CT2BK  C9     3
IATOM    CT2BK  C10    4

ATOMNAME CT2BK    0  C6 
ATOMNAME CT2BK    1  C7 
ATOMNAME CT2BK    2  C8 
ATOMNAME CT2BK    3  C9 
ATOMNAME CT2BK    4  C10

CONNECT  CT2BK  C6  ion        0    C7 
CONNECT  CT2BK  C7  ion        0    C6   0    C8 
CONNECT  CT2BK  C8  ion        0    C7   0    C9   0    C10
CONNECT  CT2BK  C9  ion        0    C8 
CONNECT  CT2BK  C10 ion        0    C8 

### This is a temporary parameter file made for residue CT3 ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST CT3        CT3BK 

NATOM    CT3BK      5

IATOM    CT3BK  C11    0
IATOM    CT3BK  C12    1
IATOM    CT3BK  C13    2
IATOM    CT3BK  C14    3
IATOM    CT3BK  C15    4

ATOMNAME CT3BK    0  C11
ATOMNAME CT3BK    1  C12
ATOMNAME CT3BK    2  C13
ATOMNAME CT3BK    3  C14
ATOMNAME CT3BK    4  C15

CONNECT  CT3BK  C11 ion        0    C12
CONNECT  CT3BK  C12 ion        0    C11  0    C13
CONNECT  CT3BK  C13 ion        0    C12  0    C14  0    C15
CONNECT  CT3BK  C14 ion        0    C13
CONNECT  CT3BK  C15 ion        0    C13

### This is a temporary parameter file made for residue CT4 ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST CT4        CT4BK 

NATOM    CT4BK      5

IATOM    CT4BK  C16    0
IATOM    CT4BK  C17    1
IATOM    CT4BK  C18    2
IATOM    CT4BK  C19    3
IATOM    CT4BK  C20    4

ATOMNAME CT4BK    0  C16
ATOMNAME CT4BK    1  C17
ATOMNAME CT4BK    2  C18
ATOMNAME CT4BK    3  C19
ATOMNAME CT4BK    4  C20

CONNECT  CT4BK  C16 ion        0    C17
CONNECT  CT4BK  C17 ion        0    C16  0    C18
CONNECT  CT4BK  C18 ion        0    C17  0    C19  0    C20
CONNECT  CT4BK  C19 ion        0    C18
CONNECT  CT4BK  C20 ion        0    C18

### This is a temporary parameter file made for residue SQD ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST SQD        SQDBK 

NATOM    SQDBK      26

IATOM    SQDBK  O6     0
IATOM    SQDBK  C44    1
IATOM    SQDBK  C45    2
IATOM    SQDBK  C46    3
IATOM    SQDBK  O47    4
IATOM    SQDBK  C7     5
IATOM    SQDBK  O49    6
IATOM    SQDBK  C8     7
IATOM    SQDBK  O48    8
IATOM    SQDBK  C23    9
IATOM    SQDBK  O10   10
IATOM    SQDBK  C24   11
IATOM    SQDBK  C1    12
IATOM    SQDBK  C2    13
IATOM    SQDBK  O2    14
IATOM    SQDBK  C3    15
IATOM    SQDBK  O3    16
IATOM    SQDBK  C4    17
IATOM    SQDBK  O4    18
IATOM    SQDBK  C5    19
IATOM    SQDBK  C6    20
IATOM    SQDBK  O5    21
IATOM    SQDBK  S     22
IATOM    SQDBK  O7    23
IATOM    SQDBK  O8    24
IATOM    SQDBK  O9    25

ATOMNAME SQDBK    0  O6 
ATOMNAME SQDBK    1  C44
ATOMNAME SQDBK    2  C45
ATOMNAME SQDBK    3  C46
ATOMNAME SQDBK    4  O47
ATOMNAME SQDBK    5  C7 
ATOMNAME SQDBK    6  O49
ATOMNAME SQDBK    7  C8 
ATOMNAME SQDBK    8  O48
ATOMNAME SQDBK    9  C23
ATOMNAME SQDBK   10  O10
ATOMNAME SQDBK   11  C24
ATOMNAME SQDBK   12  C1 
ATOMNAME SQDBK   13  C2 
ATOMNAME SQDBK   14  O2 
ATOMNAME SQDBK   15  C3 
ATOMNAME SQDBK   16  O3 
ATOMNAME SQDBK   17  C4 
ATOMNAME SQDBK   18  O4 
ATOMNAME SQDBK   19  C5 
ATOMNAME SQDBK   20  C6 
ATOMNAME SQDBK   21  O5 
ATOMNAME SQDBK   22  S  
ATOMNAME SQDBK   23  O7 
ATOMNAME SQDBK   24  O8 
ATOMNAME SQDBK   25  O9 

CONNECT  SQDBK  O6  ion        0    C44  0    C1 
CONNECT  SQDBK  C44 ion        0    O6   0    C45
CONNECT  SQDBK  C45 ion        0    C44  0    C46  0    O47
CONNECT  SQDBK  C46 ion        0    C45  0    O48
CONNECT  SQDBK  O47 ion        0    C45  0    C7 
CONNECT  SQDBK  C7  ion        0    O47  0    O49  0    C8 
CONNECT  SQDBK  O49 ion        0    C7 
CONNECT  SQDBK  C8  ion        0    C7 
CONNECT  SQDBK  O48 ion        0    C46  0    C23
CONNECT  SQDBK  C23 ion        0    O48  0    O10  0    C24
CONNECT  SQDBK  O10 ion        0    C23
CONNECT  SQDBK  C24 ion        0    C23
CONNECT  SQDBK  C1  ion        0    O6   0    C2   0    O5 
CONNECT  SQDBK  C2  ion        0    C1   0    O2   0    C3 
CONNECT  SQDBK  O2  ion        0    C2 
CONNECT  SQDBK  C3  ion        0    C2   0    O3   0    C4 
CONNECT  SQDBK  O3  ion        0    C3 
CONNECT  SQDBK  C4  ion        0    C3   0    O4   0    C5 
CONNECT  SQDBK  O4  ion        0    C4 
CONNECT  SQDBK  C5  ion        0    C4   0    C6   0    O5 
CONNECT  SQDBK  C6  ion        0    C5   0    S  
CONNECT  SQDBK  O5  ion        0    C1   0    C5 
CONNECT  SQDBK  S   ion        0    C6   0    O7   0    O8   0    O9 
CONNECT  SQDBK  O7  ion        0    S  
CONNECT  SQDBK  O8  ion        0    S  
CONNECT  SQDBK  O9  ion        0    S  

### This is a temporary parameter file made for residue sqd ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST sqd        sqdBK 

NATOM    sqdBK      28

IATOM    sqdBK  C9     0
IATOM    sqdBK  C10    1
IATOM    sqdBK  C11    2
IATOM    sqdBK  C12    3
IATOM    sqdBK  C13    4
IATOM    sqdBK  C14    5
IATOM    sqdBK  C15    6
IATOM    sqdBK  C16    7
IATOM    sqdBK  C17    8
IATOM    sqdBK  C18    9
IATOM    sqdBK  C19   10
IATOM    sqdBK  C20   11
IATOM    sqdBK  C21   12
IATOM    sqdBK  C22   13
IATOM    sqdBK  C25   14
IATOM    sqdBK  C26   15
IATOM    sqdBK  C27   16
IATOM    sqdBK  C28   17
IATOM    sqdBK  C29   18
IATOM    sqdBK  C30   19
IATOM    sqdBK  C31   20
IATOM    sqdBK  C32   21
IATOM    sqdBK  C33   22
IATOM    sqdBK  C34   23
IATOM    sqdBK  C35   24
IATOM    sqdBK  C36   25
IATOM    sqdBK  C37   26
IATOM    sqdBK  C38   27

ATOMNAME sqdBK    0  C9 
ATOMNAME sqdBK    1  C10
ATOMNAME sqdBK    2  C11
ATOMNAME sqdBK    3  C12
ATOMNAME sqdBK    4  C13
ATOMNAME sqdBK    5  C14
ATOMNAME sqdBK    6  C15
ATOMNAME sqdBK    7  C16
ATOMNAME sqdBK    8  C17
ATOMNAME sqdBK    9  C18
ATOMNAME sqdBK   10  C19
ATOMNAME sqdBK   11  C20
ATOMNAME sqdBK   12  C21
ATOMNAME sqdBK   13  C22
ATOMNAME sqdBK   14  C25
ATOMNAME sqdBK   15  C26
ATOMNAME sqdBK   16  C27
ATOMNAME sqdBK   17  C28
ATOMNAME sqdBK   18  C29
ATOMNAME sqdBK   19  C30
ATOMNAME sqdBK   20  C31
ATOMNAME sqdBK   21  C32
ATOMNAME sqdBK   22  C33
ATOMNAME sqdBK   23  C34
ATOMNAME sqdBK   24  C35
ATOMNAME sqdBK   25  C36
ATOMNAME sqdBK   26  C37
ATOMNAME sqdBK   27  C38

CONNECT  sqdBK  C9  ion        0    C10
CONNECT  sqdBK  C10 ion        0    C9   0    C11
CONNECT  sqdBK  C11 ion        0    C10  0    C12
CONNECT  sqdBK  C12 ion        0    C11  0    C13
CONNECT  sqdBK  C13 ion        0    C12  0    C14
CONNECT  sqdBK  C14 ion        0    C13  0    C15
CONNECT  sqdBK  C15 ion        0    C14  0    C16
CONNECT  sqdBK  C16 ion        0    C15  0    C17
CONNECT  sqdBK  C17 ion        0    C16  0    C18
CONNECT  sqdBK  C18 ion        0    C17  0    C19
CONNECT  sqdBK  C19 ion        0    C18  0    C20
CONNECT  sqdBK  C20 ion        0    C19  0    C21
CONNECT  sqdBK  C21 ion        0    C20  0    C22
CONNECT  sqdBK  C22 ion        0    C21
CONNECT  sqdBK  C25 ion        0    C26
CONNECT  sqdBK  C26 ion        0    C25  0    C27
CONNECT  sqdBK  C27 ion        0    C26  0    C28
CONNECT  sqdBK  C28 ion        0    C27  0    C29
CONNECT  sqdBK  C29 ion        0    C28  0    C30
CONNECT  sqdBK  C30 ion        0    C29  0    C31
CONNECT  sqdBK  C31 ion        0    C30  0    C32
CONNECT  sqdBK  C32 ion        0    C31  0    C33
CONNECT  sqdBK  C33 ion        0    C32  0    C34
CONNECT  sqdBK  C34 ion        0    C33  0    C35
CONNECT  sqdBK  C35 ion        0    C34  0    C36
CONNECT  sqdBK  C36 ion        0    C35  0    C37
CONNECT  sqdBK  C37 ion        0    C36  0    C38
CONNECT  sqdBK  C38 ion        0    C37

### This is a temporary parameter file made for residue LMG ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST LMG        LMGBK 

NATOM    LMGBK      21

IATOM    LMGBK  C1     0
IATOM    LMGBK  O1     1
IATOM    LMGBK  C2     2
IATOM    LMGBK  O2     3
IATOM    LMGBK  C3     4
IATOM    LMGBK  O3     5
IATOM    LMGBK  C4     6
IATOM    LMGBK  O4     7
IATOM    LMGBK  C5     8
IATOM    LMGBK  O5     9
IATOM    LMGBK  C6    10
IATOM    LMGBK  O6    11
IATOM    LMGBK  C7    12
IATOM    LMGBK  C8    13
IATOM    LMGBK  C9    14
IATOM    LMGBK  O7    15
IATOM    LMGBK  C10   16
IATOM    LMGBK  O9    17
IATOM    LMGBK  O8    18
IATOM    LMGBK  C28   19
IATOM    LMGBK  O10   20

ATOMNAME LMGBK    0  C1 
ATOMNAME LMGBK    1  O1 
ATOMNAME LMGBK    2  C2 
ATOMNAME LMGBK    3  O2 
ATOMNAME LMGBK    4  C3 
ATOMNAME LMGBK    5  O3 
ATOMNAME LMGBK    6  C4 
ATOMNAME LMGBK    7  O4 
ATOMNAME LMGBK    8  C5 
ATOMNAME LMGBK    9  O5 
ATOMNAME LMGBK   10  C6 
ATOMNAME LMGBK   11  O6 
ATOMNAME LMGBK   12  C7 
ATOMNAME LMGBK   13  C8 
ATOMNAME LMGBK   14  C9 
ATOMNAME LMGBK   15  O7 
ATOMNAME LMGBK   16  C10
ATOMNAME LMGBK   17  O9 
ATOMNAME LMGBK   18  O8 
ATOMNAME LMGBK   19  C28
ATOMNAME LMGBK   20  O10

CONNECT  LMGBK  C1  ion        0    O1   0    C2   0    O6 
CONNECT  LMGBK  O1  ion        0    C1   0    C7 
CONNECT  LMGBK  C2  ion        0    C1   0    O2   0    C3 
CONNECT  LMGBK  O2  ion        0    C2 
CONNECT  LMGBK  C3  ion        0    C2   0    O3   0    C4 
CONNECT  LMGBK  O3  ion        0    C3 
CONNECT  LMGBK  C4  ion        0    C3   0    O4   0    C5 
CONNECT  LMGBK  O4  ion        0    C4 
CONNECT  LMGBK  C5  ion        0    C4   0    C6   0    O6 
CONNECT  LMGBK  O5  ion        0    C6 
CONNECT  LMGBK  C6  ion        0    C5   0    O5 
CONNECT  LMGBK  O6  ion        0    C1   0    C5 
CONNECT  LMGBK  C7  ion        0    O1   0    C8 
CONNECT  LMGBK  C8  ion        0    C7   0    C9   0    O7 
CONNECT  LMGBK  C9  ion        0    C8   0    O8 
CONNECT  LMGBK  O7  ion        0    C8   0    C10
CONNECT  LMGBK  C10 ion        0    O7   0    O9 
CONNECT  LMGBK  O9  ion        0    C10
CONNECT  LMGBK  O8  ion        0    C9   0    C28
CONNECT  LMGBK  C28 ion        0    O8   0    O10
CONNECT  LMGBK  O10 ion        0    C28

### This is a temporary parameter file made for residue lmg ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST lmg        lmgBK 

NATOM    lmgBK      30

IATOM    lmgBK  C11    0
IATOM    lmgBK  C12    1
IATOM    lmgBK  C13    2
IATOM    lmgBK  C14    3
IATOM    lmgBK  C15    4
IATOM    lmgBK  C16    5
IATOM    lmgBK  C17    6
IATOM    lmgBK  C18    7
IATOM    lmgBK  C19    8
IATOM    lmgBK  C20    9
IATOM    lmgBK  C21   10
IATOM    lmgBK  C22   11
IATOM    lmgBK  C23   12
IATOM    lmgBK  C24   13
IATOM    lmgBK  C25   14
IATOM    lmgBK  C29   15
IATOM    lmgBK  C30   16
IATOM    lmgBK  C31   17
IATOM    lmgBK  C32   18
IATOM    lmgBK  C33   19
IATOM    lmgBK  C34   20
IATOM    lmgBK  C35   21
IATOM    lmgBK  C36   22
IATOM    lmgBK  C37   23
IATOM    lmgBK  C38   24
IATOM    lmgBK  C39   25
IATOM    lmgBK  C40   26
IATOM    lmgBK  C41   27
IATOM    lmgBK  C42   28
IATOM    lmgBK  C43   29

ATOMNAME lmgBK    0  C11
ATOMNAME lmgBK    1  C12
ATOMNAME lmgBK    2  C13
ATOMNAME lmgBK    3  C14
ATOMNAME lmgBK    4  C15
ATOMNAME lmgBK    5  C16
ATOMNAME lmgBK    6  C17
ATOMNAME lmgBK    7  C18
ATOMNAME lmgBK    8  C19
ATOMNAME lmgBK    9  C20
ATOMNAME lmgBK   10  C21
ATOMNAME lmgBK   11  C22
ATOMNAME lmgBK   12  C23
ATOMNAME lmgBK   13  C24
ATOMNAME lmgBK   14  C25
ATOMNAME lmgBK   15  C29
ATOMNAME lmgBK   16  C30
ATOMNAME lmgBK   17  C31
ATOMNAME lmgBK   18  C32
ATOMNAME lmgBK   19  C33
ATOMNAME lmgBK   20  C34
ATOMNAME lmgBK   21  C35
ATOMNAME lmgBK   22  C36
ATOMNAME lmgBK   23  C37
ATOMNAME lmgBK   24  C38
ATOMNAME lmgBK   25  C39
ATOMNAME lmgBK   26  C40
ATOMNAME lmgBK   27  C41
ATOMNAME lmgBK   28  C42
ATOMNAME lmgBK   29  C43

CONNECT  lmgBK  C11 ion        0    C12
CONNECT  lmgBK  C12 ion        0    C11  0    C13
CONNECT  lmgBK  C13 ion        0    C12  0    C14
CONNECT  lmgBK  C14 ion        0    C13  0    C15
CONNECT  lmgBK  C15 ion        0    C14  0    C16
CONNECT  lmgBK  C16 ion        0    C15  0    C17
CONNECT  lmgBK  C17 ion        0    C16  0    C18
CONNECT  lmgBK  C18 ion        0    C17  0    C19
CONNECT  lmgBK  C19 ion        0    C18  0    C20
CONNECT  lmgBK  C20 ion        0    C19  0    C21
CONNECT  lmgBK  C21 ion        0    C20  0    C22
CONNECT  lmgBK  C22 ion        0    C21  0    C23
CONNECT  lmgBK  C23 ion        0    C22  0    C24
CONNECT  lmgBK  C24 ion        0    C23  0    C25
CONNECT  lmgBK  C25 ion        0    C24
CONNECT  lmgBK  C29 ion        0    C30
CONNECT  lmgBK  C30 ion        0    C29  0    C31
CONNECT  lmgBK  C31 ion        0    C30  0    C32
CONNECT  lmgBK  C32 ion        0    C31  0    C33
CONNECT  lmgBK  C33 ion        0    C32  0    C34
CONNECT  lmgBK  C34 ion        0    C33  0    C35
CONNECT  lmgBK  C35 ion        0    C34  0    C36
CONNECT  lmgBK  C36 ion        0    C35  0    C37
CONNECT  lmgBK  C37 ion        0    C36  0    C38
CONNECT  lmgBK  C38 ion        0    C37  0    C39
CONNECT  lmgBK  C39 ion        0    C38  0    C40
CONNECT  lmgBK  C40 ion        0    C39  0    C41
CONNECT  lmgBK  C41 ion        0    C40  0    C42
CONNECT  lmgBK  C42 ion        0    C41  0    C43
CONNECT  lmgBK  C43 ion        0    C42

### This is a temporary parameter file made for residue UNL ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST UNL        UNLBK 

NATOM    UNLBK      36

IATOM    UNLBK  O9     0
IATOM    UNLBK  C7     1
IATOM    UNLBK  C8     2
IATOM    UNLBK  C9     3
IATOM    UNLBK  C10    4
IATOM    UNLBK  C11    5
IATOM    UNLBK  C12    6
IATOM    UNLBK  C13    7
IATOM    UNLBK  C14    8
IATOM    UNLBK  C15    9
IATOM    UNLBK  C16   10
IATOM    UNLBK  C17   11
IATOM    UNLBK  C18   12
IATOM    UNLBK  O7    13
IATOM    UNLBK  C5    14
IATOM    UNLBK  C4    15
IATOM    UNLBK  O6    16
IATOM    UNLBK  C6    17
IATOM    UNLBK  O8    18
IATOM    UNLBK  C23   19
IATOM    UNLBK  O10   20
IATOM    UNLBK  C24   21
IATOM    UNLBK  C25   22
IATOM    UNLBK  C26   23
IATOM    UNLBK  C27   24
IATOM    UNLBK  C28   25
IATOM    UNLBK  C29   26
IATOM    UNLBK  C30   27
IATOM    UNLBK  C31   28
IATOM    UNLBK  C32   29
IATOM    UNLBK  C33   30
IATOM    UNLBK  C34   31
IATOM    UNLBK  C35   32
IATOM    UNLBK  C36   33
IATOM    UNLBK  C37   34
IATOM    UNLBK  C38   35

ATOMNAME UNLBK    0  O9 
ATOMNAME UNLBK    1  C7 
ATOMNAME UNLBK    2  C8 
ATOMNAME UNLBK    3  C9 
ATOMNAME UNLBK    4  C10
ATOMNAME UNLBK    5  C11
ATOMNAME UNLBK    6  C12
ATOMNAME UNLBK    7  C13
ATOMNAME UNLBK    8  C14
ATOMNAME UNLBK    9  C15
ATOMNAME UNLBK   10  C16
ATOMNAME UNLBK   11  C17
ATOMNAME UNLBK   12  C18
ATOMNAME UNLBK   13  O7 
ATOMNAME UNLBK   14  C5 
ATOMNAME UNLBK   15  C4 
ATOMNAME UNLBK   16  O6 
ATOMNAME UNLBK   17  C6 
ATOMNAME UNLBK   18  O8 
ATOMNAME UNLBK   19  C23
ATOMNAME UNLBK   20  O10
ATOMNAME UNLBK   21  C24
ATOMNAME UNLBK   22  C25
ATOMNAME UNLBK   23  C26
ATOMNAME UNLBK   24  C27
ATOMNAME UNLBK   25  C28
ATOMNAME UNLBK   26  C29
ATOMNAME UNLBK   27  C30
ATOMNAME UNLBK   28  C31
ATOMNAME UNLBK   29  C32
ATOMNAME UNLBK   30  C33
ATOMNAME UNLBK   31  C34
ATOMNAME UNLBK   32  C35
ATOMNAME UNLBK   33  C36
ATOMNAME UNLBK   34  C37
ATOMNAME UNLBK   35  C38

CONNECT  UNLBK  O9  ion        0    C7 
CONNECT  UNLBK  C7  ion        0    O9   0    C8   0    O7 
CONNECT  UNLBK  C8  ion        0    C7   0    C9 
CONNECT  UNLBK  C9  ion        0    C8   0    C10
CONNECT  UNLBK  C10 ion        0    C9   0    C11
CONNECT  UNLBK  C11 ion        0    C10  0    C12
CONNECT  UNLBK  C12 ion        0    C11  0    C13
CONNECT  UNLBK  C13 ion        0    C12  0    C14
CONNECT  UNLBK  C14 ion        0    C13  0    C15
CONNECT  UNLBK  C15 ion        0    C14  0    C16
CONNECT  UNLBK  C16 ion        0    C15  0    C17
CONNECT  UNLBK  C17 ion        0    C16  0    C18
CONNECT  UNLBK  C18 ion        0    C17
CONNECT  UNLBK  O7  ion        0    C7   0    C5 
CONNECT  UNLBK  C5  ion        0    O7   0    C4   0    C6 
CONNECT  UNLBK  C4  ion        0    C5   0    O6 
CONNECT  UNLBK  O6  ion        0    C4 
CONNECT  UNLBK  C6  ion        0    C5   0    O8 
CONNECT  UNLBK  O8  ion        0    C6   0    C23
CONNECT  UNLBK  C23 ion        0    O8   0    O10  0    C24
CONNECT  UNLBK  O10 ion        0    C23
CONNECT  UNLBK  C24 ion        0    C23  0    C25
CONNECT  UNLBK  C25 ion        0    C24  0    C26
CONNECT  UNLBK  C26 ion        0    C25  0    C27
CONNECT  UNLBK  C27 ion        0    C26  0    C28
CONNECT  UNLBK  C28 ion        0    C27  0    C29
CONNECT  UNLBK  C29 ion        0    C28  0    C30
CONNECT  UNLBK  C30 ion        0    C29  0    C31
CONNECT  UNLBK  C31 ion        0    C30  0    C32
CONNECT  UNLBK  C32 ion        0    C31  0    C33
CONNECT  UNLBK  C33 ion        0    C32  0    C34
CONNECT  UNLBK  C34 ion        0    C33  0    C35
CONNECT  UNLBK  C35 ion        0    C34  0    C36
CONNECT  UNLBK  C36 ion        0    C35  0    C37
CONNECT  UNLBK  C37 ion        0    C36  0    C38
CONNECT  UNLBK  C38 ion        0    C37

### This is a temporary parameter file made for residue LMT ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST LMT        LMTBK 

NATOM    LMTBK      24

IATOM    LMTBK  C1B    0
IATOM    LMTBK  C2B    1
IATOM    LMTBK  C3B    2
IATOM    LMTBK  C4B    3
IATOM    LMTBK  C5B    4
IATOM    LMTBK  C6B    5
IATOM    LMTBK  O1B    6
IATOM    LMTBK  O2B    7
IATOM    LMTBK  O3B    8
IATOM    LMTBK  O4'    9
IATOM    LMTBK  O5B   10
IATOM    LMTBK  O6B   11
IATOM    LMTBK  C1'   12
IATOM    LMTBK  C2'   13
IATOM    LMTBK  C3'   14
IATOM    LMTBK  C4'   15
IATOM    LMTBK  C5'   16
IATOM    LMTBK  C6'   17
IATOM    LMTBK  O1'   18
IATOM    LMTBK  O2'   19
IATOM    LMTBK  O3'   20
IATOM    LMTBK  O5'   21
IATOM    LMTBK  O6'   22
IATOM    LMTBK  C1    23

ATOMNAME LMTBK    0  C1B
ATOMNAME LMTBK    1  C2B
ATOMNAME LMTBK    2  C3B
ATOMNAME LMTBK    3  C4B
ATOMNAME LMTBK    4  C5B
ATOMNAME LMTBK    5  C6B
ATOMNAME LMTBK    6  O1B
ATOMNAME LMTBK    7  O2B
ATOMNAME LMTBK    8  O3B
ATOMNAME LMTBK    9  O4'
ATOMNAME LMTBK   10  O5B
ATOMNAME LMTBK   11  O6B
ATOMNAME LMTBK   12  C1'
ATOMNAME LMTBK   13  C2'
ATOMNAME LMTBK   14  C3'
ATOMNAME LMTBK   15  C4'
ATOMNAME LMTBK   16  C5'
ATOMNAME LMTBK   17  C6'
ATOMNAME LMTBK   18  O1'
ATOMNAME LMTBK   19  O2'
ATOMNAME LMTBK   20  O3'
ATOMNAME LMTBK   21  O5'
ATOMNAME LMTBK   22  O6'
ATOMNAME LMTBK   23  C1 

CONNECT  LMTBK  C1B ion        0    C2B  0    O1B  0    O5B
CONNECT  LMTBK  C2B ion        0    C1B  0    C3B  0    O2B
CONNECT  LMTBK  C3B ion        0    C2B  0    C4B  0    O3B
CONNECT  LMTBK  C4B ion        0    C3B  0    C5B  0    O4'
CONNECT  LMTBK  C5B ion        0    C4B  0    C6B  0    O5B
CONNECT  LMTBK  C6B ion        0    C5B  0    O6B
CONNECT  LMTBK  O1B ion        0    C1B  0    C4'
CONNECT  LMTBK  O2B ion        0    C2B
CONNECT  LMTBK  O3B ion        0    C3B
CONNECT  LMTBK  O4' ion        0    C4B
CONNECT  LMTBK  O5B ion        0    C1B  0    C5B
CONNECT  LMTBK  O6B ion        0    C6B
CONNECT  LMTBK  C1' ion        0    C2'  0    O1'  0    O5'
CONNECT  LMTBK  C2' ion        0    C1'  0    C3'  0    O2'
CONNECT  LMTBK  C3' ion        0    C2'  0    C4'  0    O3'
CONNECT  LMTBK  C4' ion        0    O1B  0    C3'  0    C5'
CONNECT  LMTBK  C5' ion        0    C4'  0    C6'  0    O5'
CONNECT  LMTBK  C6' ion        0    C5'  0    O6'
CONNECT  LMTBK  O1' ion        0    C1'  0    C1 
CONNECT  LMTBK  O2' ion        0    C2'
CONNECT  LMTBK  O3' ion        0    C3'
CONNECT  LMTBK  O5' ion        0    C1'  0    C5'
CONNECT  LMTBK  O6' ion        0    C6'
CONNECT  LMTBK  C1  ion        0    O1'

### This is a temporary parameter file made for residue lmt ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST lmt        lmtBK 

NATOM    lmtBK      11

IATOM    lmtBK  C2     0
IATOM    lmtBK  C3     1
IATOM    lmtBK  C4     2
IATOM    lmtBK  C5     3
IATOM    lmtBK  C6     4
IATOM    lmtBK  C7     5
IATOM    lmtBK  C8     6
IATOM    lmtBK  C9     7
IATOM    lmtBK  C10    8
IATOM    lmtBK  C11    9
IATOM    lmtBK  C12   10

ATOMNAME lmtBK    0  C2 
ATOMNAME lmtBK    1  C3 
ATOMNAME lmtBK    2  C4 
ATOMNAME lmtBK    3  C5 
ATOMNAME lmtBK    4  C6 
ATOMNAME lmtBK    5  C7 
ATOMNAME lmtBK    6  C8 
ATOMNAME lmtBK    7  C9 
ATOMNAME lmtBK    8  C10
ATOMNAME lmtBK    9  C11
ATOMNAME lmtBK   10  C12

CONNECT  lmtBK  C2  ion        0    C3 
CONNECT  lmtBK  C3  ion        0    C2   0    C4 
CONNECT  lmtBK  C4  ion        0    C3   0    C5 
CONNECT  lmtBK  C5  ion        0    C4   0    C6 
CONNECT  lmtBK  C6  ion        0    C5   0    C7 
CONNECT  lmtBK  C7  ion        0    C6   0    C8 
CONNECT  lmtBK  C8  ion        0    C7   0    C9 
CONNECT  lmtBK  C9  ion        0    C8   0    C10
CONNECT  lmtBK  C10 ion        0    C9   0    C11
CONNECT  lmtBK  C11 ion        0    C10  0    C12
CONNECT  lmtBK  C12 ion        0    C11

### This is a temporary parameter file made for residue GOL ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST GOL        GOLBK 

NATOM    GOLBK      6

IATOM    GOLBK  C1     0
IATOM    GOLBK  O1     1
IATOM    GOLBK  C2     2
IATOM    GOLBK  O2     3
IATOM    GOLBK  C3     4
IATOM    GOLBK  O3     5

ATOMNAME GOLBK    0  C1 
ATOMNAME GOLBK    1  O1 
ATOMNAME GOLBK    2  C2 
ATOMNAME GOLBK    3  O2 
ATOMNAME GOLBK    4  C3 
ATOMNAME GOLBK    5  O3 

CONNECT  GOLBK  C1  ion        0    O1   0    C2 
CONNECT  GOLBK  O1  ion        0    C1 
CONNECT  GOLBK  C2  ion        0    C1   0    O2   0    C3 
CONNECT  GOLBK  O2  ion        0    C2 
CONNECT  GOLBK  C3  ion        0    C2   0    O3 
CONNECT  GOLBK  O3  ion        0    C3 

### This is a temporary parameter file made for residue HTG ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST HTG        HTGBK 

NATOM    HTGBK      19

IATOM    HTGBK  C1     0
IATOM    HTGBK  S1     1
IATOM    HTGBK  C2     2
IATOM    HTGBK  O2     3
IATOM    HTGBK  C3     4
IATOM    HTGBK  O3     5
IATOM    HTGBK  C4     6
IATOM    HTGBK  O4     7
IATOM    HTGBK  C5     8
IATOM    HTGBK  O5     9
IATOM    HTGBK  C6    10
IATOM    HTGBK  O6    11
IATOM    HTGBK  C1'   12
IATOM    HTGBK  C2'   13
IATOM    HTGBK  C3'   14
IATOM    HTGBK  C4'   15
IATOM    HTGBK  C5'   16
IATOM    HTGBK  C6'   17
IATOM    HTGBK  C7'   18

ATOMNAME HTGBK    0  C1 
ATOMNAME HTGBK    1  S1 
ATOMNAME HTGBK    2  C2 
ATOMNAME HTGBK    3  O2 
ATOMNAME HTGBK    4  C3 
ATOMNAME HTGBK    5  O3 
ATOMNAME HTGBK    6  C4 
ATOMNAME HTGBK    7  O4 
ATOMNAME HTGBK    8  C5 
ATOMNAME HTGBK    9  O5 
ATOMNAME HTGBK   10  C6 
ATOMNAME HTGBK   11  O6 
ATOMNAME HTGBK   12  C1'
ATOMNAME HTGBK   13  C2'
ATOMNAME HTGBK   14  C3'
ATOMNAME HTGBK   15  C4'
ATOMNAME HTGBK   16  C5'
ATOMNAME HTGBK   17  C6'
ATOMNAME HTGBK   18  C7'

CONNECT  HTGBK  C1  ion        0    S1   0    C2   0    O5 
CONNECT  HTGBK  S1  ion        0    C1   0    C1'
CONNECT  HTGBK  C2  ion        0    C1   0    O2   0    C3 
CONNECT  HTGBK  O2  ion        0    C2 
CONNECT  HTGBK  C3  ion        0    C2   0    O3   0    C4 
CONNECT  HTGBK  O3  ion        0    C3 
CONNECT  HTGBK  C4  ion        0    C3   0    O4   0    C5 
CONNECT  HTGBK  O4  ion        0    C4 
CONNECT  HTGBK  C5  ion        0    C4   0    O5   0    C6 
CONNECT  HTGBK  O5  ion        0    C1   0    C5 
CONNECT  HTGBK  C6  ion        0    C5   0    O6 
CONNECT  HTGBK  O6  ion        0    C6 
CONNECT  HTGBK  C1' ion        0    S1   0    C2'
CONNECT  HTGBK  C2' ion        0    C1'  0    C3'
CONNECT  HTGBK  C3' ion        0    C2'  0    C4'
CONNECT  HTGBK  C4' ion        0    C3'  0    C5'
CONNECT  HTGBK  C5' ion        0    C4'  0    C6'
CONNECT  HTGBK  C6' ion        0    C5'  0    C7'
CONNECT  HTGBK  C7' ion        0    C6'

### This is a temporary parameter file made for residue DGA ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST DGA        DGABK 

NATOM    DGABK      21

IATOM    DGABK  C1A    0
IATOM    DGABK  C2A    1
IATOM    DGABK  O1A    2
IATOM    DGABK  C1B    3
IATOM    DGABK  C2B    4
IATOM    DGABK  O1B    5
IATOM    DGABK  O1G    6
IATOM    DGABK  C1G    7
IATOM    DGABK  C2G    8
IATOM    DGABK  O2G    9
IATOM    DGABK  C3G   10
IATOM    DGABK  O3G   11
IATOM    DGABK  C1D   12
IATOM    DGABK  C2D   13
IATOM    DGABK  O2D   14
IATOM    DGABK  C3D   15
IATOM    DGABK  O3D   16
IATOM    DGABK  C4D   17
IATOM    DGABK  O4D   18
IATOM    DGABK  C5D   19
IATOM    DGABK  O6D   20

ATOMNAME DGABK    0  C1A
ATOMNAME DGABK    1  C2A
ATOMNAME DGABK    2  O1A
ATOMNAME DGABK    3  C1B
ATOMNAME DGABK    4  C2B
ATOMNAME DGABK    5  O1B
ATOMNAME DGABK    6  O1G
ATOMNAME DGABK    7  C1G
ATOMNAME DGABK    8  C2G
ATOMNAME DGABK    9  O2G
ATOMNAME DGABK   10  C3G
ATOMNAME DGABK   11  O3G
ATOMNAME DGABK   12  C1D
ATOMNAME DGABK   13  C2D
ATOMNAME DGABK   14  O2D
ATOMNAME DGABK   15  C3D
ATOMNAME DGABK   16  O3D
ATOMNAME DGABK   17  C4D
ATOMNAME DGABK   18  O4D
ATOMNAME DGABK   19  C5D
ATOMNAME DGABK   20  O6D

CONNECT  DGABK  C1A ion        0    C2A  0    O1A  0    O1G
CONNECT  DGABK  C2A ion        0    C1A
CONNECT  DGABK  O1A ion        0    C1A
CONNECT  DGABK  C1B ion        0    C2B  0    O1B  0    O2G
CONNECT  DGABK  C2B ion        0    C1B
CONNECT  DGABK  O1B ion        0    C1B
CONNECT  DGABK  O1G ion        0    C1A  0    C1G
CONNECT  DGABK  C1G ion        0    O1G  0    C2G
CONNECT  DGABK  C2G ion        0    C1G  0    O2G  0    C3G
CONNECT  DGABK  O2G ion        0    C1B  0    C2G
CONNECT  DGABK  C3G ion        0    C2G  0    O3G
CONNECT  DGABK  O3G ion        0    C3G  0    C1D
CONNECT  DGABK  C1D ion        0    O3G  0    C2D  0    O6D
CONNECT  DGABK  C2D ion        0    C1D  0    O2D  0    C3D
CONNECT  DGABK  O2D ion        0    C2D
CONNECT  DGABK  C3D ion        0    C2D  0    O3D  0    C4D
CONNECT  DGABK  O3D ion        0    C3D
CONNECT  DGABK  C4D ion        0    C3D  0    O4D  0    C5D
CONNECT  DGABK  O4D ion        0    C4D
CONNECT  DGABK  C5D ion        0    C4D  0    O6D
CONNECT  DGABK  O6D ion        0    C1D  0    C5D

### This is a temporary parameter file made for residue dgd ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST dgd        dgdBK 

NATOM    dgdBK      28

IATOM    dgdBK  C3A    0
IATOM    dgdBK  C4A    1
IATOM    dgdBK  C5A    2
IATOM    dgdBK  C6A    3
IATOM    dgdBK  C7A    4
IATOM    dgdBK  C8A    5
IATOM    dgdBK  C9A    6
IATOM    dgdBK  CAA    7
IATOM    dgdBK  CBA    8
IATOM    dgdBK  CCA    9
IATOM    dgdBK  CDA   10
IATOM    dgdBK  CEA   11
IATOM    dgdBK  CFA   12
IATOM    dgdBK  CGA   13
IATOM    dgdBK  C3B   14
IATOM    dgdBK  C4B   15
IATOM    dgdBK  C5B   16
IATOM    dgdBK  C6B   17
IATOM    dgdBK  C7B   18
IATOM    dgdBK  C8B   19
IATOM    dgdBK  C9B   20
IATOM    dgdBK  CAB   21
IATOM    dgdBK  CBB   22
IATOM    dgdBK  CCB   23
IATOM    dgdBK  CDB   24
IATOM    dgdBK  CEB   25
IATOM    dgdBK  CFB   26
IATOM    dgdBK  CGB   27

ATOMNAME dgdBK    0  C3A
ATOMNAME dgdBK    1  C4A
ATOMNAME dgdBK    2  C5A
ATOMNAME dgdBK    3  C6A
ATOMNAME dgdBK    4  C7A
ATOMNAME dgdBK    5  C8A
ATOMNAME dgdBK    6  C9A
ATOMNAME dgdBK    7  CAA
ATOMNAME dgdBK    8  CBA
ATOMNAME dgdBK    9  CCA
ATOMNAME dgdBK   10  CDA
ATOMNAME dgdBK   11  CEA
ATOMNAME dgdBK   12  CFA
ATOMNAME dgdBK   13  CGA
ATOMNAME dgdBK   14  C3B
ATOMNAME dgdBK   15  C4B
ATOMNAME dgdBK   16  C5B
ATOMNAME dgdBK   17  C6B
ATOMNAME dgdBK   18  C7B
ATOMNAME dgdBK   19  C8B
ATOMNAME dgdBK   20  C9B
ATOMNAME dgdBK   21  CAB
ATOMNAME dgdBK   22  CBB
ATOMNAME dgdBK   23  CCB
ATOMNAME dgdBK   24  CDB
ATOMNAME dgdBK   25  CEB
ATOMNAME dgdBK   26  CFB
ATOMNAME dgdBK   27  CGB

CONNECT  dgdBK  C3A ion        0    C4A
CONNECT  dgdBK  C4A ion        0    C3A  0    C5A
CONNECT  dgdBK  C5A ion        0    C4A  0    C6A
CONNECT  dgdBK  C6A ion        0    C5A  0    C7A
CONNECT  dgdBK  C7A ion        0    C6A  0    C8A
CONNECT  dgdBK  C8A ion        0    C7A  0    C9A
CONNECT  dgdBK  C9A ion        0    C8A  0    CAA
CONNECT  dgdBK  CAA ion        0    C9A  0    CBA
CONNECT  dgdBK  CBA ion        0    CAA  0    CCA
CONNECT  dgdBK  CCA ion        0    CBA  0    CDA
CONNECT  dgdBK  CDA ion        0    CCA  0    CEA
CONNECT  dgdBK  CEA ion        0    CDA  0    CFA
CONNECT  dgdBK  CFA ion        0    CEA  0    CGA
CONNECT  dgdBK  CGA ion        0    CFA
CONNECT  dgdBK  C3B ion        0    C4B
CONNECT  dgdBK  C4B ion        0    C3B  0    C5B
CONNECT  dgdBK  C5B ion        0    C4B  0    C6B
CONNECT  dgdBK  C6B ion        0    C5B  0    C7B
CONNECT  dgdBK  C7B ion        0    C6B  0    C8B
CONNECT  dgdBK  C8B ion        0    C7B  0    C9B
CONNECT  dgdBK  C9B ion        0    C8B  0    CAB
CONNECT  dgdBK  CAB ion        0    C9B  0    CBB
CONNECT  dgdBK  CBB ion        0    CAB  0    CCB
CONNECT  dgdBK  CCB ion        0    CBB  0    CDB
CONNECT  dgdBK  CDB ion        0    CCB  0    CEB
CONNECT  dgdBK  CEB ion        0    CDB  0    CFB
CONNECT  dgdBK  CFB ion        0    CEB  0    CGB
CONNECT  dgdBK  CGB ion        0    CFB

### This is a temporary parameter file made for residue DGB ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST DGB        DGBBK 

NATOM    DGBBK      13

IATOM    DGBBK  O5D    0
IATOM    DGBBK  C6D    1
IATOM    DGBBK  C1E    2
IATOM    DGBBK  C2E    3
IATOM    DGBBK  O2E    4
IATOM    DGBBK  C3E    5
IATOM    DGBBK  O3E    6
IATOM    DGBBK  C4E    7
IATOM    DGBBK  O4E    8
IATOM    DGBBK  C5E    9
IATOM    DGBBK  O6E   10
IATOM    DGBBK  C6E   11
IATOM    DGBBK  O5E   12

ATOMNAME DGBBK    0  O5D
ATOMNAME DGBBK    1  C6D
ATOMNAME DGBBK    2  C1E
ATOMNAME DGBBK    3  C2E
ATOMNAME DGBBK    4  O2E
ATOMNAME DGBBK    5  C3E
ATOMNAME DGBBK    6  O3E
ATOMNAME DGBBK    7  C4E
ATOMNAME DGBBK    8  O4E
ATOMNAME DGBBK    9  C5E
ATOMNAME DGBBK   10  O6E
ATOMNAME DGBBK   11  C6E
ATOMNAME DGBBK   12  O5E

CONNECT  DGBBK  O5D ion        0    C6D  0    C1E
CONNECT  DGBBK  C6D ion        0    O5D
CONNECT  DGBBK  C1E ion        0    O5D  0    C2E  0    O6E
CONNECT  DGBBK  C2E ion        0    C1E  0    O2E  0    C3E
CONNECT  DGBBK  O2E ion        0    C2E
CONNECT  DGBBK  C3E ion        0    C2E  0    O3E  0    C4E
CONNECT  DGBBK  O3E ion        0    C3E
CONNECT  DGBBK  C4E ion        0    C3E  0    O4E  0    C5E
CONNECT  DGBBK  O4E ion        0    C4E
CONNECT  DGBBK  C5E ion        0    C4E  0    O6E  0    C6E
CONNECT  DGBBK  O6E ion        0    C1E  0    C5E
CONNECT  DGBBK  C6E ion        0    C5E  0    O5E
CONNECT  DGBBK  O5E ion        0    C6E

### This is a temporary parameter file made for residue LHG ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST LHG        LHGBK 

NATOM    LHGBK      21

IATOM    LHGBK  O1     0
IATOM    LHGBK  C1     1
IATOM    LHGBK  C2     2
IATOM    LHGBK  O2     3
IATOM    LHGBK  C3     4
IATOM    LHGBK  O3     5
IATOM    LHGBK  P      6
IATOM    LHGBK  O4     7
IATOM    LHGBK  O5     8
IATOM    LHGBK  O6     9
IATOM    LHGBK  C4    10
IATOM    LHGBK  C5    11
IATOM    LHGBK  C6    12
IATOM    LHGBK  O7    13
IATOM    LHGBK  C7    14
IATOM    LHGBK  O9    15
IATOM    LHGBK  C8    16
IATOM    LHGBK  O8    17
IATOM    LHGBK  C23   18
IATOM    LHGBK  O10   19
IATOM    LHGBK  C24   20

ATOMNAME LHGBK    0  O1 
ATOMNAME LHGBK    1  C1 
ATOMNAME LHGBK    2  C2 
ATOMNAME LHGBK    3  O2 
ATOMNAME LHGBK    4  C3 
ATOMNAME LHGBK    5  O3 
ATOMNAME LHGBK    6  P  
ATOMNAME LHGBK    7  O4 
ATOMNAME LHGBK    8  O5 
ATOMNAME LHGBK    9  O6 
ATOMNAME LHGBK   10  C4 
ATOMNAME LHGBK   11  C5 
ATOMNAME LHGBK   12  C6 
ATOMNAME LHGBK   13  O7 
ATOMNAME LHGBK   14  C7 
ATOMNAME LHGBK   15  O9 
ATOMNAME LHGBK   16  C8 
ATOMNAME LHGBK   17  O8 
ATOMNAME LHGBK   18  C23
ATOMNAME LHGBK   19  O10
ATOMNAME LHGBK   20  C24

CONNECT  LHGBK  O1  ion        0    C1 
CONNECT  LHGBK  C1  ion        0    O1   0    C2 
CONNECT  LHGBK  C2  ion        0    C1   0    O2   0    C3 
CONNECT  LHGBK  O2  ion        0    C2 
CONNECT  LHGBK  C3  ion        0    C2   0    O3 
CONNECT  LHGBK  O3  ion        0    C3   0    P  
CONNECT  LHGBK  P   ion        0    O3   0    O4   0    O5   0    O6 
CONNECT  LHGBK  O4  ion        0    P  
CONNECT  LHGBK  O5  ion        0    P  
CONNECT  LHGBK  O6  ion        0    P    0    C4 
CONNECT  LHGBK  C4  ion        0    O6   0    C5 
CONNECT  LHGBK  C5  ion        0    C4   0    C6   0    O7 
CONNECT  LHGBK  C6  ion        0    C5   0    O8 
CONNECT  LHGBK  O7  ion        0    C5   0    C7 
CONNECT  LHGBK  C7  ion        0    O7   0    O9   0    C8 
CONNECT  LHGBK  O9  ion        0    C7 
CONNECT  LHGBK  C8  ion        0    C7 
CONNECT  LHGBK  O8  ion        0    C6   0    C23  0    O10
CONNECT  LHGBK  C23 ion        0    O8   0    O10  0    C24
CONNECT  LHGBK  O10 ion        0    O8   0    C23
CONNECT  LHGBK  C24 ion        0    C23

### This is a temporary parameter file made for residue lhg ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST lhg        lhgBK 

NATOM    lhgBK      28

IATOM    lhgBK  C9     0
IATOM    lhgBK  C10    1
IATOM    lhgBK  C11    2
IATOM    lhgBK  C12    3
IATOM    lhgBK  C13    4
IATOM    lhgBK  C14    5
IATOM    lhgBK  C15    6
IATOM    lhgBK  C16    7
IATOM    lhgBK  C17    8
IATOM    lhgBK  C18    9
IATOM    lhgBK  C19   10
IATOM    lhgBK  C20   11
IATOM    lhgBK  C21   12
IATOM    lhgBK  C22   13
IATOM    lhgBK  C25   14
IATOM    lhgBK  C26   15
IATOM    lhgBK  C27   16
IATOM    lhgBK  C28   17
IATOM    lhgBK  C29   18
IATOM    lhgBK  C30   19
IATOM    lhgBK  C31   20
IATOM    lhgBK  C32   21
IATOM    lhgBK  C33   22
IATOM    lhgBK  C34   23
IATOM    lhgBK  C35   24
IATOM    lhgBK  C36   25
IATOM    lhgBK  C37   26
IATOM    lhgBK  C38   27

ATOMNAME lhgBK    0  C9 
ATOMNAME lhgBK    1  C10
ATOMNAME lhgBK    2  C11
ATOMNAME lhgBK    3  C12
ATOMNAME lhgBK    4  C13
ATOMNAME lhgBK    5  C14
ATOMNAME lhgBK    6  C15
ATOMNAME lhgBK    7  C16
ATOMNAME lhgBK    8  C17
ATOMNAME lhgBK    9  C18
ATOMNAME lhgBK   10  C19
ATOMNAME lhgBK   11  C20
ATOMNAME lhgBK   12  C21
ATOMNAME lhgBK   13  C22
ATOMNAME lhgBK   14  C25
ATOMNAME lhgBK   15  C26
ATOMNAME lhgBK   16  C27
ATOMNAME lhgBK   17  C28
ATOMNAME lhgBK   18  C29
ATOMNAME lhgBK   19  C30
ATOMNAME lhgBK   20  C31
ATOMNAME lhgBK   21  C32
ATOMNAME lhgBK   22  C33
ATOMNAME lhgBK   23  C34
ATOMNAME lhgBK   24  C35
ATOMNAME lhgBK   25  C36
ATOMNAME lhgBK   26  C37
ATOMNAME lhgBK   27  C38

CONNECT  lhgBK  C9  ion        0    C10
CONNECT  lhgBK  C10 ion        0    C9   0    C11
CONNECT  lhgBK  C11 ion        0    C10  0    C12
CONNECT  lhgBK  C12 ion        0    C11  0    C13
CONNECT  lhgBK  C13 ion        0    C12  0    C14
CONNECT  lhgBK  C14 ion        0    C13  0    C15
CONNECT  lhgBK  C15 ion        0    C14  0    C16
CONNECT  lhgBK  C16 ion        0    C15  0    C17
CONNECT  lhgBK  C17 ion        0    C16  0    C18
CONNECT  lhgBK  C18 ion        0    C17  0    C19
CONNECT  lhgBK  C19 ion        0    C18  0    C20
CONNECT  lhgBK  C20 ion        0    C19  0    C21
CONNECT  lhgBK  C21 ion        0    C20  0    C22
CONNECT  lhgBK  C22 ion        0    C21
CONNECT  lhgBK  C25 ion        0    C26
CONNECT  lhgBK  C26 ion        0    C25  0    C27
CONNECT  lhgBK  C27 ion        0    C26  0    C28
CONNECT  lhgBK  C28 ion        0    C27  0    C29
CONNECT  lhgBK  C29 ion        0    C28  0    C30
CONNECT  lhgBK  C30 ion        0    C29  0    C31
CONNECT  lhgBK  C31 ion        0    C30  0    C32
CONNECT  lhgBK  C32 ion        0    C31  0    C33
CONNECT  lhgBK  C33 ion        0    C32  0    C34
CONNECT  lhgBK  C34 ion        0    C33  0    C35
CONNECT  lhgBK  C35 ion        0    C34  0    C36
CONNECT  lhgBK  C36 ion        0    C35  0    C37
CONNECT  lhgBK  C37 ion        0    C36  0    C38
CONNECT  lhgBK  C38 ion        0    C37

### This is a temporary parameter file made for residue RRX ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST RRX        RRXBK 

NATOM    RRXBK      41

IATOM    RRXBK  C1     0
IATOM    RRXBK  C2     1
IATOM    RRXBK  O2     2
IATOM    RRXBK  C3     3
IATOM    RRXBK  C40    4
IATOM    RRXBK  C30    5
IATOM    RRXBK  C39    6
IATOM    RRXBK  C29    7
IATOM    RRXBK  C28    8
IATOM    RRXBK  C27    9
IATOM    RRXBK  C26   10
IATOM    RRXBK  C38   11
IATOM    RRXBK  C25   12
IATOM    RRXBK  C24   13
IATOM    RRXBK  C23   14
IATOM    RRXBK  C22   15
IATOM    RRXBK  C37   16
IATOM    RRXBK  C21   17
IATOM    RRXBK  C20   18
IATOM    RRXBK  C19   19
IATOM    RRXBK  C18   20
IATOM    RRXBK  C36   21
IATOM    RRXBK  C17   22
IATOM    RRXBK  C16   23
IATOM    RRXBK  C15   24
IATOM    RRXBK  C14   25
IATOM    RRXBK  C13   26
IATOM    RRXBK  C35   27
IATOM    RRXBK  C12   28
IATOM    RRXBK  C11   29
IATOM    RRXBK  C10   30
IATOM    RRXBK  C9    31
IATOM    RRXBK  C34   32
IATOM    RRXBK  C8    33
IATOM    RRXBK  C7    34
IATOM    RRXBK  C6    35
IATOM    RRXBK  C32   36
IATOM    RRXBK  C31   37
IATOM    RRXBK  C5    38
IATOM    RRXBK  C33   39
IATOM    RRXBK  C4    40

ATOMNAME RRXBK    0  C1 
ATOMNAME RRXBK    1  C2 
ATOMNAME RRXBK    2  O2 
ATOMNAME RRXBK    3  C3 
ATOMNAME RRXBK    4  C40
ATOMNAME RRXBK    5  C30
ATOMNAME RRXBK    6  C39
ATOMNAME RRXBK    7  C29
ATOMNAME RRXBK    8  C28
ATOMNAME RRXBK    9  C27
ATOMNAME RRXBK   10  C26
ATOMNAME RRXBK   11  C38
ATOMNAME RRXBK   12  C25
ATOMNAME RRXBK   13  C24
ATOMNAME RRXBK   14  C23
ATOMNAME RRXBK   15  C22
ATOMNAME RRXBK   16  C37
ATOMNAME RRXBK   17  C21
ATOMNAME RRXBK   18  C20
ATOMNAME RRXBK   19  C19
ATOMNAME RRXBK   20  C18
ATOMNAME RRXBK   21  C36
ATOMNAME RRXBK   22  C17
ATOMNAME RRXBK   23  C16
ATOMNAME RRXBK   24  C15
ATOMNAME RRXBK   25  C14
ATOMNAME RRXBK   26  C13
ATOMNAME RRXBK   27  C35
ATOMNAME RRXBK   28  C12
ATOMNAME RRXBK   29  C11
ATOMNAME RRXBK   30  C10
ATOMNAME RRXBK   31  C9 
ATOMNAME RRXBK   32  C34
ATOMNAME RRXBK   33  C8 
ATOMNAME RRXBK   34  C7 
ATOMNAME RRXBK   35  C6 
ATOMNAME RRXBK   36  C32
ATOMNAME RRXBK   37  C31
ATOMNAME RRXBK   38  C5 
ATOMNAME RRXBK   39  C33
ATOMNAME RRXBK   40  C4 

CONNECT  RRXBK  C1  ion        0    C2   0    C6   0    C32  0    C31
CONNECT  RRXBK  C2  ion        0    C1   0    C3 
CONNECT  RRXBK  O2  ion        0    C28
CONNECT  RRXBK  C3  ion        0    C2   0    C4 
CONNECT  RRXBK  C40 ion        0    C30
CONNECT  RRXBK  C30 ion        0    C40  0    C39  0    C29  0    C25
CONNECT  RRXBK  C39 ion        0    C30
CONNECT  RRXBK  C29 ion        0    C30  0    C28
CONNECT  RRXBK  C28 ion        0    O2   0    C29  0    C27
CONNECT  RRXBK  C27 ion        0    C28  0    C26
CONNECT  RRXBK  C26 ion        0    C27  0    C38  0    C25
CONNECT  RRXBK  C38 ion        0    C26
CONNECT  RRXBK  C25 ion        0    C30  0    C26  0    C24
CONNECT  RRXBK  C24 ion        0    C25  0    C23
CONNECT  RRXBK  C23 ion        0    C24  0    C22
CONNECT  RRXBK  C22 ion        0    C23  0    C37  0    C21
CONNECT  RRXBK  C37 ion        0    C22
CONNECT  RRXBK  C21 ion        0    C22  0    C20
CONNECT  RRXBK  C20 ion        0    C21  0    C19
CONNECT  RRXBK  C19 ion        0    C20  0    C18
CONNECT  RRXBK  C18 ion        0    C19  0    C36  0    C17
CONNECT  RRXBK  C36 ion        0    C18
CONNECT  RRXBK  C17 ion        0    C18  0    C16
CONNECT  RRXBK  C16 ion        0    C17  0    C15
CONNECT  RRXBK  C15 ion        0    C16  0    C14
CONNECT  RRXBK  C14 ion        0    C15  0    C13
CONNECT  RRXBK  C13 ion        0    C14  0    C35  0    C12
CONNECT  RRXBK  C35 ion        0    C13
CONNECT  RRXBK  C12 ion        0    C13  0    C11
CONNECT  RRXBK  C11 ion        0    C12  0    C10
CONNECT  RRXBK  C10 ion        0    C11  0    C9 
CONNECT  RRXBK  C9  ion        0    C10  0    C34  0    C8 
CONNECT  RRXBK  C34 ion        0    C9 
CONNECT  RRXBK  C8  ion        0    C9   0    C7 
CONNECT  RRXBK  C7  ion        0    C8   0    C6 
CONNECT  RRXBK  C6  ion        0    C1   0    C7   0    C5 
CONNECT  RRXBK  C32 ion        0    C1 
CONNECT  RRXBK  C31 ion        0    C1 
CONNECT  RRXBK  C5  ion        0    C6   0    C33  0    C4 
CONNECT  RRXBK  C33 ion        0    C5 
CONNECT  RRXBK  C4  ion        0    C3   0    C5 

### This is a temporary parameter file made for residue _MG ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST _MG        _MGBK 

NATOM    _MGBK      1

IATOM    _MGBK MG      0

ATOMNAME _MGBK    0 MG  

CONNECT  _MGBK MG   ion      

### This is a temporary parameter file made for residue SO4 ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST SO4        SO4BK 

NATOM    SO4BK      5

IATOM    SO4BK  S      0
IATOM    SO4BK  O1     1
IATOM    SO4BK  O2     2
IATOM    SO4BK  O3     3
IATOM    SO4BK  O4     4

ATOMNAME SO4BK    0  S  
ATOMNAME SO4BK    1  O1 
ATOMNAME SO4BK    2  O2 
ATOMNAME SO4BK    3  O3 
ATOMNAME SO4BK    4  O4 

CONNECT  SO4BK  S   ion        0    O1   0    O2   0    O3   0    O4 
CONNECT  SO4BK  O1  ion        0    S  
CONNECT  SO4BK  O2  ion        0    S  
CONNECT  SO4BK  O3  ion        0    S  
CONNECT  SO4BK  O4  ion        0    S  

