
;$PARAM TVKG =0.6, TVKS0=0.4, TVGAMMA=0.8, TVBASE =70, DOSE=10
;
;$OMEGA 0.005 0.03 0.003 0.01
;
;$CMT RESP
;
;$MAIN
;double KG = TVKG * exp(ETA(1));
;double KS0 = TVKS0 * exp(ETA(2));
;double GAMMA = TVGAMMA * exp(ETA(3));
;double BASE = TVBASE*exp(ETA(4));
;RESP_0 = BASE;

;$ODE
;double KS = KS0 *exp( -GAMMA * SOLVERTIME);
;dxdt_RESP = KG * RESP - KS*log(DOSE)*RESP;
;$CAPTURE
;RESP_0 BASE KS KS0 KG GAMMA ETA(4)

$PROB RUN# c001

$INPUT C ID TIME DOSE CMT EVID MDV DV

$DATA data/claret001.csv IGNORE=C

$SUBROUTINES ADVAN13 TOL=8

$MODEL NCOMP=1

$PK
KG = 0.6;
KS0 = 0.4;
GAMMA=0.8;
BASE=70;
A_0(1) = BASE;

$DES
KS = KS0 *EXP( -GAMMA * T);
DADT(1) = KG * A(1) - KS*LOG(DOSE)*A(1)

$ERROR
IPRED=A(1);
Y=IPRED+ERR(1);

$THETA
0, FIX

$OMEGA
0.0 FIX
0.0 FIX
0.0 FIX

$SIGMA
0.00 FIX

$TABLE FILE=TAB ID TIME Y NOPRINT ONEHEADER NOAPPEND

$SIMULATION (2674474) ONLYSIMULATION
