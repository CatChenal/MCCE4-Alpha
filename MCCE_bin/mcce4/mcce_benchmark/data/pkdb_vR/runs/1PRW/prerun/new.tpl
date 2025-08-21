### This is a temporary parameter file made for residue ACE ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST ACE        ACEBK 

NATOM    ACEBK      6

IATOM    ACEBK  CH3    0
IATOM    ACEBK  CH3    1
IATOM    ACEBK  C      2
IATOM    ACEBK  O      3
IATOM    ACEBK  C      4
IATOM    ACEBK  O      5

ATOMNAME ACEBK    0  CH3
ATOMNAME ACEBK    1  CH3
ATOMNAME ACEBK    2  C  
ATOMNAME ACEBK    3  O  
ATOMNAME ACEBK    4  C  
ATOMNAME ACEBK    5  O  

CONNECT  ACEBK  CH3 ion        0    CH3  0    C    0    C  
CONNECT  ACEBK  CH3 ion        0    CH3  0    C    0    C  
CONNECT  ACEBK  C   ion        0    CH3  0    CH3  0    O    0    C    0    O  
CONNECT  ACEBK  O   ion        0    C    0    C    0    O  
CONNECT  ACEBK  C   ion        0    CH3  0    CH3  0    C    0    O    0    O  
CONNECT  ACEBK  O   ion        0    C    0    O    0    C  

### This is a temporary parameter file made for residue M3L ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST M3L        M3LBK 

NATOM    M3LBK      12

IATOM    M3LBK  N      0
IATOM    M3LBK  CA     1
IATOM    M3LBK  CB     2
IATOM    M3LBK  CG     3
IATOM    M3LBK  CD     4
IATOM    M3LBK  CE     5
IATOM    M3LBK  NZ     6
IATOM    M3LBK  CM1    7
IATOM    M3LBK  CM2    8
IATOM    M3LBK  CM3    9
IATOM    M3LBK  C     10
IATOM    M3LBK  O     11

ATOMNAME M3LBK    0  N  
ATOMNAME M3LBK    1  CA 
ATOMNAME M3LBK    2  CB 
ATOMNAME M3LBK    3  CG 
ATOMNAME M3LBK    4  CD 
ATOMNAME M3LBK    5  CE 
ATOMNAME M3LBK    6  NZ 
ATOMNAME M3LBK    7  CM1
ATOMNAME M3LBK    8  CM2
ATOMNAME M3LBK    9  CM3
ATOMNAME M3LBK   10  C  
ATOMNAME M3LBK   11  O  

CONNECT  M3LBK  N   ion        0    CA 
CONNECT  M3LBK  CA  ion        0    N    0    CB   0    C  
CONNECT  M3LBK  CB  ion        0    CA   0    CG 
CONNECT  M3LBK  CG  ion        0    CB   0    CD 
CONNECT  M3LBK  CD  ion        0    CG   0    CE 
CONNECT  M3LBK  CE  ion        0    CD   0    NZ 
CONNECT  M3LBK  NZ  ion        0    CE   0    CM1  0    CM2  0    CM3
CONNECT  M3LBK  CM1 ion        0    NZ 
CONNECT  M3LBK  CM2 ion        0    NZ 
CONNECT  M3LBK  CM3 ion        0    NZ 
CONNECT  M3LBK  C   ion        0    CA   0    O  
CONNECT  M3LBK  O   ion        0    C  

