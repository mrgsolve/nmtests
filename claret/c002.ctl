
;$PRED
;double KG = TVKG * exp(ETA(1));
;double KS0 = TVKS0 * exp(ETA(2));
;double GAMMA = TVGAMMA * exp(ETA(3));
;double BASE=TVBASE * exp(ETA(4));
;double RESP_0=BASE;

;double KS = KS0*exp(-GAMMA * TIME);
;double RESP = RESP_0 *exp((KG*TIME)-(KS0*log(DOSE)/GAMMA)*(1-exp(-GAMMA*TIME)));


$PROB RUN# c001

$INPUT C ID TIME DOSE CMT EVID MDV DV

$DATA data/claret001.csv IGNORE=C

$PRED
KG = 0.6;
KS0 = 0.4;
GAMMA=0.8;
BASE=70;
RESP_0 = BASE;

KS = KS0 *EXP( -GAMMA * TIME);
RESP = RESP_0 * EXP((KG*TIME) - (KS0 * LOG(DOSE)/GAMMA) * (1-EXP(-GAMMA*TIME)))

IPRED=RESP;
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

