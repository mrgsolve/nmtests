$PROB RUN# 103

$INPUT C ID TIME EVID AMT CMT SS II ADDL RATE LAGT MODE DUR2 RAT2 BIOAV TVCL DV

$DATA ../../data/103.csv IGNORE=C

$SUBROUTINES ADVAN2 TRANS2

$PK

CL=TVCL*EXP(ETA(1))

TVV2=THETA(2)
V=TVV2*EXP(ETA(2))

TVKA=THETA(3)
KA=TVKA*EXP(ETA(3))

ALAG2 = LAGT
F2 = BIOAV

IF(MODE.EQ.1) R2 = RAT2
IF(MODE.EQ.2) D2 = DUR2

$ERROR
IPRED=A(2)/(V/1000)
Y=IPRED*EXP(ERR(1))

CP = IPRED

$THETA
(5,   FIX) ;; CL
(30,  FIX) ;; V
(1.5, FIX) ;; KA

$OMEGA
0.0 FIX
0.0 FIX
0.0 FIX

$SIGMA
0.00 FIX

$TABLE FILE=TAB TIME EVID CP IPRED PRED DV NOPRINT ONEHEADER NOAPPEND

$SIMULATION (2674474) ONLYSIMULATION

