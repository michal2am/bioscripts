# PySCeS input file
# Stochastic Simulation Algorithm input format
# JWM model GABAAR model

R1:
    R > AR
    kon2*R
R2:
    AR > R
    koff*AR
R3:
    AR > A2R
    kon*AR
R4:
    A2R > AR
    koff2*A2R
R5:
    A2R > A2D
    d*A2R
R6:
    A2D > A2R
    r*A2D
R7:
    A2R > A2O
    b*A2R
R8:
    A2O > A2R
    a*A2O

# InitPar
kon = 99.0
kon2 = 192.0
koff = 330
koff2 = 660
d = 5.27
r = 0.07
a = 1.03
b = 7.67

# InitVar
R = 100
AR = 0
A2R = 0
A2D = 0
A2O = 0
