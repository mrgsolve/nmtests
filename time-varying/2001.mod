$SET req = ""
$PARAM TVCL = 1.1, V = 20, KA = 1.5
LAGT = 0, MODE = 0, DUR2 = 2, RAT2 = 10, BIOAV = 1, 
WT = 70
$PKMODEL cmt = "GUT CENT", depot = TRUE
$MAIN
double CL = TVCL*pow(WT/70,0.75); 
$TABLE
capture DV = (CENT/(V/1000));
capture CP = DV;
$CAPTURE LAGT MODE DUR2 RAT2 BIOAV WT CL
