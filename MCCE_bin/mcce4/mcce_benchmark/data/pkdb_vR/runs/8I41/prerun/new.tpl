### This is a temporary parameter file made for residue PLC ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST PLC        PLCBK 

NATOM    PLCBK      36

IATOM    PLCBK  C1     0
IATOM    PLCBK  C2     1
IATOM    PLCBK  C3     2
IATOM    PLCBK  C'     3
IATOM    PLCBK  C1'    4
IATOM    PLCBK  C2'    5
IATOM    PLCBK  C3'    6
IATOM    PLCBK  C4'    7
IATOM    PLCBK  C5'    8
IATOM    PLCBK  C6'    9
IATOM    PLCBK  C7'   10
IATOM    PLCBK  C8'   11
IATOM    PLCBK  C9'   12
IATOM    PLCBK  CA'   13
IATOM    PLCBK  CB'   14
IATOM    PLCBK  CB    15
IATOM    PLCBK  C1B   16
IATOM    PLCBK  C2B   17
IATOM    PLCBK  C3B   18
IATOM    PLCBK  C4B   19
IATOM    PLCBK  C5B   20
IATOM    PLCBK  C6B   21
IATOM    PLCBK  C7B   22
IATOM    PLCBK  C8B   23
IATOM    PLCBK  C9B   24
IATOM    PLCBK  CAA   25
IATOM    PLCBK  CBA   26
IATOM    PLCBK  O'    27
IATOM    PLCBK  OB    28
IATOM    PLCBK  O2    29
IATOM    PLCBK  O3    30
IATOM    PLCBK  O1P   31
IATOM    PLCBK  O2P   32
IATOM    PLCBK  O3P   33
IATOM    PLCBK  O4P   34
IATOM    PLCBK  P     35

ATOMNAME PLCBK    0  C1 
ATOMNAME PLCBK    1  C2 
ATOMNAME PLCBK    2  C3 
ATOMNAME PLCBK    3  C' 
ATOMNAME PLCBK    4  C1'
ATOMNAME PLCBK    5  C2'
ATOMNAME PLCBK    6  C3'
ATOMNAME PLCBK    7  C4'
ATOMNAME PLCBK    8  C5'
ATOMNAME PLCBK    9  C6'
ATOMNAME PLCBK   10  C7'
ATOMNAME PLCBK   11  C8'
ATOMNAME PLCBK   12  C9'
ATOMNAME PLCBK   13  CA'
ATOMNAME PLCBK   14  CB'
ATOMNAME PLCBK   15  CB 
ATOMNAME PLCBK   16  C1B
ATOMNAME PLCBK   17  C2B
ATOMNAME PLCBK   18  C3B
ATOMNAME PLCBK   19  C4B
ATOMNAME PLCBK   20  C5B
ATOMNAME PLCBK   21  C6B
ATOMNAME PLCBK   22  C7B
ATOMNAME PLCBK   23  C8B
ATOMNAME PLCBK   24  C9B
ATOMNAME PLCBK   25  CAA
ATOMNAME PLCBK   26  CBA
ATOMNAME PLCBK   27  O' 
ATOMNAME PLCBK   28  OB 
ATOMNAME PLCBK   29  O2 
ATOMNAME PLCBK   30  O3 
ATOMNAME PLCBK   31  O1P
ATOMNAME PLCBK   32  O2P
ATOMNAME PLCBK   33  O3P
ATOMNAME PLCBK   34  O4P
ATOMNAME PLCBK   35  P  

CONNECT  PLCBK  C1  ion        0    C2   0    O3P
CONNECT  PLCBK  C2  ion        0    C1   0    C3   0    O2 
CONNECT  PLCBK  C3  ion        0    C2   0    O3 
CONNECT  PLCBK  C'  ion        0    C1'  0    O'   0    O2 
CONNECT  PLCBK  C1' ion        0    C'   0    C2'
CONNECT  PLCBK  C2' ion        0    C1'  0    C3'
CONNECT  PLCBK  C3' ion        0    C2'  0    C4'
CONNECT  PLCBK  C4' ion        0    C3'  0    C5'
CONNECT  PLCBK  C5' ion        0    C4'  0    C6'
CONNECT  PLCBK  C6' ion        0    C5'  0    C7'
CONNECT  PLCBK  C7' ion        0    C6'  0    C8'
CONNECT  PLCBK  C8' ion        0    C7'  0    C9'
CONNECT  PLCBK  C9' ion        0    C8'  0    CA'
CONNECT  PLCBK  CA' ion        0    C9'  0    CB'
CONNECT  PLCBK  CB' ion        0    CA'
CONNECT  PLCBK  CB  ion        0    C1B  0    OB   0    O3 
CONNECT  PLCBK  C1B ion        0    CB   0    C2B
CONNECT  PLCBK  C2B ion        0    C1B  0    C3B
CONNECT  PLCBK  C3B ion        0    C2B  0    C4B
CONNECT  PLCBK  C4B ion        0    C3B  0    C5B
CONNECT  PLCBK  C5B ion        0    C4B  0    C6B
CONNECT  PLCBK  C6B ion        0    C5B  0    C7B
CONNECT  PLCBK  C7B ion        0    C6B  0    C8B
CONNECT  PLCBK  C8B ion        0    C7B  0    C9B
CONNECT  PLCBK  C9B ion        0    C8B  0    CAA
CONNECT  PLCBK  CAA ion        0    C9B  0    CBA
CONNECT  PLCBK  CBA ion        0    CAA
CONNECT  PLCBK  O'  ion        0    C' 
CONNECT  PLCBK  OB  ion        0    CB 
CONNECT  PLCBK  O2  ion        0    C2   0    C' 
CONNECT  PLCBK  O3  ion        0    C3   0    CB 
CONNECT  PLCBK  O1P ion        0    P  
CONNECT  PLCBK  O2P ion        0    P  
CONNECT  PLCBK  O3P ion        0    C1   0    P  
CONNECT  PLCBK  O4P ion        0    P  
CONNECT  PLCBK  P   ion        0    O1P  0    O2P  0    O3P  0    O4P

