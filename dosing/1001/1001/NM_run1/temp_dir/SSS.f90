!*********************************COPYRIGHT******************************************
!                                                                                   !
!       THE NONMEM SYSTEM MAY BE DISTRIBUTED ONLY BY ICON DEVELOPMENT               !
!       SOLUTIONS.                                                                  !
!                                                                                   !
!       COPYRIGHT BY ICON DEVELOPMENT SOLUTIONS                                     !
!       2009-2020 ALL RIGHTS RESERVED.                                              !
!                                                                                   !
!       DO NOT ATTEMPT TO MODIFY CODE WITHOUT FIRST CONSULTING WITH                 !
!       ICON DEVELOPMENT SOLUTIONS.                                                 !
!                                                                                   !
!************************************************************************************
!
!-----------------------------HISTORY------------------------------------------------
! VERSION     : NONMEM VII
! AUTHOR      : ALISON J. BOECKMANN
! CREATED ON  : JAN/1984
! LANGUAGE    : FORTRAN 90/95
! LAST UPDATE : JUL/2008 - COMMON BLOCKS REPLACED WITH MODULES
!               NOV/2008 - INTRODUCED HEADER INFORMATIONS AND RESTRUCTURED AS PER
!                          THE NONMEM STANDARDS
!               APR/2009 - IQUIT ADDED TO HANDLE ERROR MESSAGES
!               OCT/2009 - BUG FIX WHEN USING THE FEATURE "SS=2".  
!               FEB/2010 - CHANGED SIZES TO PRSIZES
!               AUG/2012 - INTEGRATED NONMEM7.3ALPHA6.3 CHANGES
!
!----------------------------- SSS.F90 ----------------------------------------------
!
! SUBROUTINE SSS
!
! DESCRIPTION : Supervisor of steady state
!
! ARGUMENTS   : NONE
!
! CALLED BY   : PREDI - Initialization-Finalization routine of PRED
!               PRED  - This is PREDPP main program. Provides prediction, partial 
!                       derivatives of the statistical model with respect to ETA
!                       and EPSILON random variables and stores them in the G
!                       and H arguments of the PRED routine.
!
! CALLS       : SS - ??
!
! ALGORITHM   : - If ICALLD=1, execute initialization call for DES at ICALL=1
!               - If (SSID /= 0) Normal entry; Else
!                 - CALL SS; RETURN
!               - Normal entry
!                 - When I_SS feature is active, use nothing from current event record
!                 - Initialize all SS variables to 0
!                 - Get multiple infusion
!                   - Duration is modelled
!                   - Get multiple serial infusions
!                   - Rate is modelled
!                 - Else, get multiple bolus doses.
!                 - Get constant infusion
!                   - Rate is modelled
!                 - Execute test for ADVAN6, 8, and 10 with SS. Otherwise, COMPAC=.TRUE.
!                   is default
!                 - Start of ON/OFF feature:
!                   - General model, some compartment(s) OFF. 
!                   - Compress state vector, status vector, and set up mapping
!                   - When I_SS feature is active, use nothing from current event record
!                   - Redo the mapping for compact arrays
!                   - End of ON/OFF compression
!                 - If mapped:
!                   - General model, some compartment(s) OFF
!                   - Uncompress state vector, status vector, and restore 1-1 mapping
!                   - End of ON/OFF feature
!                 - Set up infusions
!                 - Set up MI+1 infusions with constant values
!                 - Add variable information, shortest infusion first
!                 - Compute for remaining infusions: shortest to longest
!               - RETURN
!
! MODULES USED: PRSIZES,NMPRD_INT,PRCM_INT,PRCOM_INT,PRMOD_INT,PROCM_INT,PRCM_REAL,
!               PRCOM_REAL,PROCM_REAL,NMPRD_CHAR,PRCM_LOG,PRCOM_LOG
!
! CONTAINS    : NONE
!
! LOCAL'S     : AMT,DEL1,G2,GG,I,IF,II,IX,J,JJ,K,KD,KP,KR,MI,MIC,MM,RATE
!
!---------------------------- END OF HEADER -----------------------------------------
!
      SUBROUTINE SSS
!
      USE PRSIZES,    ONLY: ISIZE,DPSIZE,MAXIC
! INTEGER
      USE NMPRD_INT,  ONLY: IERPRD,IQUIT
      USE PRCM_INT,   ONLY: AA,IAI,IAJ,IAK,IPI,IPJ,IPK,ITI,ITK,LNCM1,LNCM2,MAP,PP,TT,&
                            MAPINV,NBRON,NIAI,NIAI1,NIPI,NIPI1,NIT,NIT1,XAA,XIAI,XNC,&
                            XIAJ,XIAK,XIPI,XIPJ,XIPK,XITI,XITK,XLNCM1,XLNCM2,XNCM1,  &
                            XNIAI,XNIAI1,XNIPI,XNIPI1,XNIT,XNIT1,XPP,XTT
      USE PRCOM_INT,  ONLY: IBACK,IBF,ID,IHEAD,INEXT,IRGG,IRR,JAMT,JDELTA,JRATE,JSS, &
                            BETA,IP,IPOOL,MCOMP,NC,NCM1,NPETAS,SSC,SSID,SV
      USE PRMOD_INT,  ONLY: ISSNOW,ICALLD
      USE PROCM_INT,  ONLY: NEVENT,IDXETA
! REAL
      USE PRCM_REAL,  ONLY: SSAE2,SSRE2,SDELE2,RHOE2
      USE PRCOM_REAL, ONLY: D2DTE,DT,DTE,I2AEA,IA,IAA,IAEA,IDA,IDEA,IRA,IREA,SAMT,   &
                            SDEL1,I2DEA,I2REA,R,RHO,RHOE,SDEL,SDELE,SSA,SSAE,SSR,    &
                            SSRE,XD,XR,G3,ZERO
      USE PROCM_REAL, ONLY: AMNT,DAETA,D2AETA,EVTREC
! CHARACTER
      USE NMPRD_CHAR, ONLY: ETEXT
! LOGICAL
      USE PRCM_LOG,   ONLY: DOFINL,GENMOD,MAPPED,COMPAC
      USE PRCOM_LOG,  ONLY: NEWWAY,NOETAS,SECOND
!
      IMPLICIT NONE
!
      SAVE
!
!------------------------------------------------------------------------------------
!     COMMON /NMPRD1/ IERPRD,NETEXT
!     COMMON /NMPRD2/ ETEXT(3)
!     INTEGER IERPRD,NETEXT
!     CHARACTER*132 ETEXT
!     COMMON /PRCOM0/ NP,NBP,YFORM
!     COMMON /PRCOM0/ MAXKF,IFORM
!     COMMON /PRCOM0/ IDC,IDO,MAXIC,ISV,IINST,ITURN
!     COMMON /PRCOM0/ JTIME,JCONT,JEVENT,JAMT,JRATE,JSS,JDELTA
!     COMMON /PRCOM0/ JCOMPT,JCOMPF,JERROR,SSC,KREC,JMORE,JDUM
!     COMMON /PROCM7/ EVTREC(5,PD+1),N
!     DOUBLE PRECISION EVTREC
!     COMMON /PRCOM1/ NOETAS,SECOND
!     COMMON /PRCOM2/ IBF,IRR,IS,ID,ITSC,IFR,ILAG
!     COMMON /PRCOM3/ ITRANS,IRGG,IREV,NPETAS,NPEPS
!     COMMON /PRCOM4/ G3,HH,DELTA,DT,DTE
!     COMMON /PRCOM4/ YMULT,ZERO,ONE,XR,XD,TSTART,DTSTAR
!     COMMON /PRCOM4/ DDELTA,D2DELT,ADTE,D2ADTE
!     COMMON /PRCOM5/ ISPEC,DCTR,BETA,DD
!     COMMON /PRCOM5/ IP,IPOOL,IHEAD,INEXT,IBACK,SV
!     COMMON /PRCOM6/ IA,IAA,IAEA,IRA,IREA,IDA,IDEA,R,RE
!     COMMON /PRCOM6/ RHO,RHOE,SDEL,SDELE,SSA,SSAE,SSR,SSRE
!     COMMON /PRCOM6/ SAMT,SDEL1
!     COMMON /PRCOM6/ I2AEA,I2REA,I2DEA,R2E,D2DTE,D2TSTA
!     COMMON /PRCM6A/ SSAE2,SSRE2,SDELE2,RHOE2
!     COMMON /PRCOM7/ ADVID,SSID
!     COMMON /PRCOMN/ LOGUNT,NC
!     DOUBLE PRECISION DELTA,G3,HH
!     DOUBLE PRECISION SDELE,RHOE,SSAE,SSRE,YMULT,ZERO,XR,XD
!     DOUBLE PRECISION SSAE2,SSRE2,SDELE2,RHOE2
!     DOUBLE PRECISION ONE,TSTART,DTSTAR(PE)
!     DOUBLE PRECISION DDELTA(PE),D2DELT(PE,PE),ADTE(PE),D2ADTE(PE,PE)
!     DOUBLE PRECISION IA(90),IAA(90),IAEA(90,PE),IRA(90),IDA(90)
!     DOUBLE PRECISION IREA(90,PE),IDEA(90,PE),R(PC),RE(PC,PE)
!     DOUBLE PRECISION I2REA(90,PE,PE),I2DEA(90,PE,PE),I2AEA(90,PE,PE)
!     DOUBLE PRECISION R2E(PC,PE,PE),D2DTE(PE,PE)
!     DOUBLE PRECISION D2TSTA(PE,PE)
!     DOUBLE PRECISION DT,DTE(PE),RHO,SDEL,SSA,SSR
!     DOUBLE PRECISION SAMT,SDEL1
!     DIMENSION SDELE(PE),RHOE(PE),SSAE(PE),SSRE(PE)
!     DIMENSION SSAE2(PE,PE),SSRE2(PE,PE),SDELE2(PE,PE),RHOE2(PE,PE)
!     DIMENSION G3(PG+1,PE+1,PE+1),HH(PE,PE)
!     INTEGER LOGUNT,IRGG,IREV,ITRANS
!     INTEGER NPETAS,NPEPS,N
!     INTEGER JCONT,JTIME,JEVENT,JAMT,JRATE,JSS,JDELTA
!     INTEGER JCOMPT,JCOMPF,JERROR
!     INTEGER NC,IDC,IDO,NP,NBP,SSC,KREC,JMORE,JDUM
!     INTEGER ISV(PC),IBF(PC),IRR(PC),SV(PC)
!     INTEGER IINST(PC),ITURN(PC),ITSC,IFR,ILAG(PC),IS(PC),ID(PC)
!     INTEGER ADVID,SSID,MAXKF,IFORM(PG+1),YFORM,MAXIC
!     INTEGER BETA(90),IPOOL(90),IP,IHEAD,INEXT(90),IBACK(90)
!     INTEGER ISPEC,DD(90),DCTR
!     LOGICAL NOETAS,SECOND
!     COMMON /PROCM4/ A,DAETA,D2AETA
!     DOUBLE PRECISION A,DAETA,D2AETA
!     DIMENSION A(PC),DAETA(PC,PE),D2AETA(PC,PE,PE)
!     COMMON /PRCMX1/ GENMOD,MAPPED,COMPAC
!     LOGICAL GENMOD,MAPPED,COMPAC
!     COMMON /PRCMX2/ XNC,XNCM1,MAP,MAPINV,NBRON,XNBRON
!     INTEGER XNC,XNCM1,MAP(PC),MAPINV(PC),NBRON,XNBRON
!     COMMON /PRCOME/ NRD,MCOMP,NCM1,IH
!     INTEGER NRD,MCOMP,NCM1,IH
!     COMMON /PRCMLS/ XLNCM1,XLNCM2,LNCM1,LNCM2
!     INTEGER XLNCM1,XLNCM2,LNCM1,LNCM2
! BIOAVAILABILITY FLAG
!     COMMON /PRCOMU/ NEWWAY
!     LOGICAL NEWWAY
!     COMMON /PRDDE1/ICALLD,IDEFD(2),IDEFA(2)
!     INTEGER ICALLD,IDEFD,IDEFA
! COUNTERS FOR INDICES FOR THE COMPACT ARRAYS (X VERSIONS ARE UNMAPPED)
!     COMMON /PRCMDE/ NIAI,NIPI,NIT,XNIAI,XNIPI,XNIT
!     INTEGER NIAI,NIPI,NIT,XNIAI,XNIPI,XNIT
!     COMMON /PRCMDE/ IAI(PIR),IAJ(PIR),IAK(PIR)
!     COMMON /PRCMDE/ IPI(PIR),IPJ(PIR),IPK(PIR)
!     COMMON /PRCMDE/ ITI(PIR),ITK(PIR)
!     COMMON /PRCMDE/ XIAI(PIR),XIAJ(PIR),XIAK(PIR)
!     COMMON /PRCMDE/ XIPI(PIR),XIPJ(PIR),XIPK(PIR)
!     COMMON /PRCMDE/ XITI(PIR),XITK(PIR)
!     COMMON /PRCMDE/ IAC(PIR),IPC(PIR),ITC(PIR)
!     INTEGER IAI,IAJ,IAK,IPI,IPJ,IPK,ITI,ITK
!     INTEGER XIAI,XIAJ,XIAK,XIPI,XIPJ,XIPK,XITI,XITK
!     INTEGER IAC,IPC,ITC
!     COMMON /PRCMDE/ NIAI1,NIPI1,NIT1,XNIAI1,XNIPI1,XNIT1
!     COMMON /PRCMDE/ AA,PP,TT,XAA,XPP,XTT
!     INTEGER NIAI1,NIPI1,NIT1,XNIAI1,XNIPI1,XNIT1
!     INTEGER AA(PIR),PP(PIR),TT(PIR),XAA(PIR),XPP(PIR),XTT(PIR)
! FOR LVOUT FEATURE
!     COMMON /PROCM5/ NACTIV,M(0:PE)
!     INTEGER NACTIV,M
!     COMMON /PRCM03/ DIDCAA,DIDDES,DIDAES,DOFINL
!     LOGICAL DIDCAA,DIDDES,DIDAES,DOFINL
!     COMMON /PRDPK4/ I_SS,ISSNOW,ISSMOD
!     INTEGER I_SS,ISSNOW,ISSMOD
!------------------------------------------------------------------------
!
! Local Variables           
!
      INTEGER(KIND=ISIZE) :: I,IF,II,IX,J,JJ,K,KD,KP,KR,MI,MM 
!     
      REAL(KIND=DPSIZE)   :: AMT,DEL1,G2,GG,RATE 
!     
      CHARACTER(LEN=4)    :: MIC
!     
      G2(I,J,K)=G3(I,IDXETA(J)+1,IDXETA(K)+1) ! Statement function
      GG(I,J)=G3(I,IDXETA(J-1)+1,1)
!
! Steady state variables:
! From PRED: 
! SSC      = Compartment number for dose
! DT       = Capital Delta (Dosing interval)
! SDEL     = Small Delta (Duration of infusion)
! SSA      = Dose amount
! SSR      = Rate of constant infusion
! RHO      = Rate of multiple infusion
! RHOE     = Partial derivative of RHO wrt ETA(K)
! DTE(K)   = Partial derivative of DT wrt ETA(K)
! SDELE(K) = Partial derivative of SDEL wrt ETA(K)
! SSAE(K)  = Partial derivative of SSA wrt ETA(K)
! SSRE(K)  = Partial derivative of SSR wrt ETA(K)
!
! For implementation of general ADVANS' compartment ON/OFF feature
! GENMOD true for general ADVAN (set by PREDI)
! XNC, XNCM1, XNBRON are original values of NC, NCM1, NBRON (set by PREDI)
! NBRON tells how many compartments are now on (not counting output)
! MAPPED is true when a MAPPING in effect (set by SADVAN, SSS)
! MAP, MAPINV: MAPS 'REAL' Compartment Nos. To 'REDUCED SET' and V.V.
!
! LNCM1 = Number of DES-TYPE compartments (XLNCM1 is original No.)
! LNCM2 = Number of AES-TYPE compartments (XLNCM2 is original No.)
!
! Indices for the compact arrays (X versions are unmapped)
! IAC maps the contents of the compact array DAC
! IPC maps the contents of the compact array DPC
! ITC maps the contents of the compact array DTC
!
! ICALLD=1 Initialization call for DES at ICALL=1    
! Normal entry
      IF (ICALLD == 1 .OR. SSID == 0 .OR. DOFINL) THEN
        CALL SS; GO TO 999
      END IF
!
! Added 5/2008; When I_SS feature is active, use nothing from current event record
      IF (ISSNOW > 0) THEN
        DEL1=ZERO
        RATE=ZERO
        AMT=ZERO
      ELSE
        DEL1=EVTREC(NEVENT,JDELTA)
        RATE=EVTREC(NEVENT,JRATE)
        AMT= EVTREC(NEVENT,JAMT)
      END IF
!      
! Initialize all SS variables to 0
      SSA=ZERO; SDEL=ZERO;  DT=ZERO
      SSR=ZERO; SDEL1=ZERO; RHO=ZERO
!      
      IF (.NOT. NOETAS) THEN 
        SSAE(1:NPETAS)=ZERO
        SSRE(1:NPETAS)=ZERO
        DTE(1:NPETAS)=ZERO
        SDELE(1:NPETAS)=ZERO
        RHOE(1:NPETAS)=ZERO
        IF (SECOND) THEN
          DO K=1,NPETAS
            DO J=K,NPETAS
              SSAE2(J,K)=ZERO
              SSRE2(J,K)=ZERO
              D2DTE(J,K)=ZERO
              SDELE2(J,K)=ZERO
              RHOE2(J,K)=ZERO
            END DO
          END DO
        END IF
      END IF
!      
      IF (RATE /= ZERO) THEN
        IF (DEL1 /= ZERO) THEN      ! Multiple infusion
          KP=IBF(SSC)
          IF (GG(KP,1) < ZERO .AND. (RATE < ZERO .OR. NEWWAY)) THEN
            ETEXT(2)='PK PARAMETER FOR BIOAVAILABILITY FRACTION IS NEGATIVE'
            IERPRD=1; GO TO 999
          END IF
          SAMT=AMT*GG(KP,1)
          IF (RATE /= XR) THEN 
            IF (RATE == XD) THEN    ! Duration is modelled
              KD=ID(SSC)
              SDEL1=GG(KD,1)
              IF (SDEL1 <= ZERO)  THEN
                ETEXT(2)='PK PARAMETER FOR DURATION IS NON-POSITIVE'
                IERPRD=1; GO TO 999
              END IF
              RHO=SAMT/SDEL1
            ELSE                    ! Multiple serial infusions
              RHO=RATE
              IF (NEWWAY) THEN
                SDEL1=SAMT/RHO
              ELSE
                SDEL1=AMT/RHO
              END IF
            END IF 
          ELSE                      ! Rate is modelled
            KR=IRR(SSC)
            RHO=GG(KR,1)
            IF (RHO <= ZERO) THEN
              ETEXT(2)='PK PARAMETER FOR RATE IS NON-POSITIVE'
              IERPRD=1; GO TO 999
            END IF
            SDEL1=SAMT/RHO
          END IF
          MI=INT(SDEL1/DEL1)
! Fixes a bug seen on the cray, in SP version. 11/95
! It can happen that A=N*B, where A, B, and N are integers,
! But INT(A/B)=N-1  instead of INT(A/B)=N.
! This causes incorrect amounts to be calculated in SS
!          
          IF (SDEL1-MI*DEL1 == DEL1) MI=MI+1
!          
          SDEL=SDEL1-MI*DEL1
          DT=DEL1-SDEL
          SSR=MI*RHO
!          
          IF (NOETAS .OR. (RATE > 0 .AND. (.NOT. NEWWAY .OR. KP == IRGG))) GO TO 295
!          
          DO K=1,NPETAS
            IF (RATE == XR) THEN  
              RHOE(K)=GG(KR,K+1)
              SDELE(K)=(AMT*GG(KP,K+1)-SDEL1*RHOE(K))/RHO
            ELSE
              IF (RATE /= XD) THEN
                SDELE(K)=AMT*GG(KP,K+1)/RHO
              ELSE
                SDELE(K)=GG(KD,K+1)
                RHOE(K)=(AMT*GG(KP,K+1)-RHO*SDELE(K))/SDEL1
              END IF
            END IF
            DTE(K)=-SDELE(K)
            SSRE(K)=MI*RHOE(K)
          END DO
!          
          IF (SECOND) THEN
            DO K=1,NPETAS
              DO J=K,NPETAS
                IF (RATE /= XR) THEN
                  IF (RATE == XD) THEN   
                    SDELE2(J,K)=G2(KD,J,K)
                    RHOE2(J,K)=(AMT*G2(KP,J,K)-RHO*SDELE2(J,K)-SDELE(K)*RHOE(J)-RHOE(K)*SDELE(J))/SDEL1
                  ELSE 
                    SDELE2(J,K)=AMT*G2(KP,J,K)/RHO
                  END IF  
                ELSE
                  RHOE2(J,K)=G2(KR,J,K)
                  SDELE2(J,K)=(AMT*G2(KP,J,K)-SDEL1*RHOE2(J,K)-SDELE(J)*RHOE(K)-SDELE(K)*RHOE(J))/RHO
                END IF 
                D2DTE(J,K)=-SDELE2(J,K)
                SSRE2(J,K)=MI*RHOE2(J,K)
              END DO
            END DO
          END IF
          GO TO 295  
        END IF 
      ELSE     ! Multiple bolus doses
        KP=IBF(SSC)
        IF (GG(KP,1) < ZERO) THEN
          ETEXT(2)='PK PARAMETER FOR BIOAVAILABILITY FRACTION IS NEGATIVE'
          IERPRD=1; GO TO 999
        END IF
        DT=DEL1
        SSA=AMT*GG(KP,1)
        IF (.NOT. NOETAS) THEN 
          DO K=1,NPETAS
            SSAE(K)=AMT*GG(KP,K+1)
          END DO  
          IF (SECOND) THEN
            DO K=1,NPETAS
              DO J=K,NPETAS
                SSAE2(J,K)=AMT*G2(KP,J,K)
              END DO
            END DO
          END IF
        END IF
        GO TO 295      
      END IF
!      
! Constant infusion
      IF (RATE /= XR) THEN
        SSR=RATE
      ELSE          ! Rate is modelled
        KR=IRR(SSC)
        SSR=GG(KR,1)
        IF (SSR < ZERO) THEN
          ETEXT(2)='PK PARAMETER FOR RATE IS NEGATIVE'
          IERPRD=1; GO TO 999
        END IF
!        
        IF (.NOT. NOETAS) THEN
          DO K=1,NPETAS
            SSRE(K)=GG(KR,K+1)
          END DO  
          IF (SECOND) THEN
            DO K=1,NPETAS
              DO J=K,NPETAS
                SSRE2(J,K)=G2(KR,J,K)
              END DO
            END DO
          END IF
        END IF
      END IF
!
  295 CONTINUE
! This test for ADVAN6, 8, and 10 with SS. otherwise, COMPAC=.TRUE. is default
      IF (SECOND  .AND. .NOT. COMPAC) THEN
        ETEXT(2)='COMPACT DES ARRAYS REQUIRED WITH LAPLACIAN AND ANALYTICAL SECOND DERIVATIVES, OR CHOOSE NUMERICAL.'
        IERPRD=2; GO TO 999
      END IF 
!  
! ON/OFF feature
      IF (GENMOD .AND. NBRON /= NCM1) THEN
! General model, some compartment(S) OFF
! Compress state vector, status vector, and set up mapping
        LNCM1=XLNCM1
        JJ=0
        DO II=1,NCM1
          IF (SV(II) == 0) THEN
            MAP(II)=0
            IF (II <= XLNCM1) LNCM1=LNCM1-1
          ELSE
            JJ=JJ+1
            MAP(II)=JJ
            MAPINV(JJ)=II
! Changed 5/2008
! When I_SS feature is active, use nothing from current event record
!           IF (II /= JJ .AND. (EVTREC(N,JSS) == 2 .OR. ISSNOW == 2)) THEN
!            IF (II /= JJ .AND. ((ISSNOW == 0 .AND. EVTREC(NEVENT,JSS) == 2) .OR. ISSNOW == 2)) THEN
            IF (II /= JJ .AND. ((ISSNOW <= 0 .AND. EVTREC(NEVENT,JSS) == 2) .OR. ISSNOW == 2)) THEN !7.1.1

              AMNT(JJ)=AMNT(II)
              IF (.NOT. NOETAS) THEN
                DAETA(JJ,1:NPETAS)=DAETA(II,1:NPETAS)
                IF (SECOND) THEN
                  DO K=1,NPETAS
                    DO J=K,NPETAS
                      D2AETA(JJ,J,K)=D2AETA(II,J,K)
                    END DO
                  END DO
                END IF
              END IF
            END IF
            SV(JJ)=SV(II)
          END IF
        END DO
        NC=NBRON+1
        NCM1=NBRON
        LNCM2=NCM1-LNCM1
        MAPPED=.TRUE.
! Redo the mapping for compact arrays
        IF (COMPAC) THEN
          NIAI=0; NIAI1=0; NIT1=0
          NIPI=0; NIPI1=0; NIT=0
          DO I=1,XNIAI
            IF (MAP(XIAI(I)) == 0 .OR. MAP(XIAJ(I)) == 0) CYCLE
            IF (XIAK(I) /= 0 .AND. XIAK(I) <= MCOMP+1) THEN
              IF (MAP(XIAK(I)) == 0) CYCLE
            END IF
            NIAI=NIAI+1
            AA(NIAI)=XAA(I)
            IAI(NIAI)=MAP(XIAI(I))
            IAJ(NIAI)=MAP(XIAJ(I))
            IAK(NIAI)=XIAK(I)
            IF (IAK(NIAI) == 0) THEN
              NIAI1=NIAI1+1; CYCLE
            END IF
            IF (IAK(NIAI) <= MCOMP+1) IAK(NIAI)=MAP(XIAK(I))
          END DO
          DO I=1,XNIPI
            IF (MAP(XIPI(I)) == 0) CYCLE
            NIPI=NIPI+1
            PP(NIPI)=XPP(I)
            IPI(NIPI)=MAP(XIPI(I))
            IPJ(NIPI)=XIPJ(I)
            IPK(NIPI)=XIPK(I)
            IF (IPK(NIPI) == 0) NIPI1=NIPI1+1
          END DO
          DO I=1,XNIT
            IF (MAP(XITI(I)) == 0) CYCLE
            IF (XITK(I) /= 0 .AND. XITK(I) <= MCOMP+1) THEN
              IF (MAP(XITK(I)) == 0) CYCLE
            END IF
            NIT=NIT+1
            TT(NIT)=XTT(I)
            ITI(NIT)=MAP(XITI(I))
            ITK(NIT)=XITK(I)
            IF (ITK(NIT) == 0) THEN
              NIT1=NIT1+1; CYCLE
            END IF
            IF (ITK(NIT) <= MCOMP+1) ITK(NIT)=MAP(XITK(I))
          END DO
        END IF
        SSC=MAP(SSC)
      END IF        ! End ON/OFF compression
!
      CALL SS
      IF (IQUIT == 1) GO TO 999
!      
      IF (MAPPED) THEN  
! General model, some compartment(S) OFF
        SSC=MAPINV(SSC)
! Uncompress state vector, status vector, and restore 1-1 mapping
        DO II=XNCM1,1,-1
          JJ=MAP(II)
          IF (JJ == 0 .AND. II <= NBRON) THEN
            AMNT(II)=ZERO
            IF (.NOT.NOETAS) THEN
              DAETA(II,1:NPETAS)=ZERO
              IF (SECOND) THEN
                DO K=1,NPETAS
                  DO J=K,NPETAS
                    D2AETA(II,J,K)=ZERO
                  END DO
                END DO
              END IF
            END IF
            SV(II)=0
          ELSE IF (JJ /= 0 .AND. JJ /= II) THEN
            AMNT(II)=AMNT(JJ)
            IF (.NOT. NOETAS) THEN
              DAETA(II,1:NPETAS)=DAETA(JJ,1:NPETAS)
              IF (SECOND) THEN
                DO K=1,NPETAS
                  DO J=K,NPETAS
                    D2AETA(II,J,K)=D2AETA(JJ,J,K)
                  END DO
                END DO
              END IF
            END IF
            SV(II)=1
          END IF
          MAP(II)=II
          MAPINV(II)=II
        END DO
        NC=XNC
        NCM1=XNCM1
        LNCM1=XLNCM1
        LNCM2=XLNCM2
        MAPPED=.FALSE.
        IF (COMPAC) THEN
          NIAI=XNIAI; NIAI1=XNIAI1; NIT=XNIT
          NIPI=XNIPI; NIPI1=XNIPI1; NIT1=XNIT1
!                   
          AA(1:NIAI) =XAA(1:NIAI)
          IAI(1:NIAI)=XIAI(1:NIAI)
          IAJ(1:NIAI)=XIAJ(1:NIAI)
          IAK(1:NIAI)=XIAK(1:NIAI)
!          
          PP(1:NIPI) =XPP(1:NIPI)
          IPI(1:NIPI)=XIPI(1:NIPI)
          IPJ(1:NIPI)=XIPJ(1:NIPI)
          IPK(1:NIPI)=XIPK(1:NIPI)
!          
          TT(1:XNIT) =XTT(1:XNIT)
          ITI(1:XNIT)=XITI(1:XNIT)
          ITK(1:XNIT)=XITK(1:XNIT)
        END IF
      END IF    ! End ON/OFF feature
!   
      IF (IERPRD > 0 .OR. RHO == ZERO) GO TO 999
!      
! Set up infusions
      IF (IP+MI+1 > MAXIC) THEN
        WRITE (MIC,'(I4)') MAXIC
        ETEXT(2)='NUMBER OF ACTIVE INFUSIONS'//                           &
                 ' IMPLIED BY STEADY STATE DOSE RECORD'// ' EXCEEDS THE'
        ETEXT(3)='MAXIMUM ('//MIC//'). PERHAPS THE DURATION GREATLY'//    &
                 ' EXCEEDS THE INTERDOSE INTERVAL II'
        IERPRD=1; GO TO 999
      END IF
! Set up MI+1 infusions with constant values
      MM=MI+1
      DO J=1,MM
        IP=IP+1
        IF=IPOOL(IP)
        BETA(IF)=SSC
        INEXT(IF)=IHEAD
        IF (IHEAD /= 0) IBACK(IHEAD)=IF
        IHEAD=IF
        IA(IF)=RATE
        IF (RATE == XR) CYCLE
        IRA(IF)=RHO
        IF (NOETAS) CYCLE
! This works even if (RATE > 0 && NOT NEWWAY) because RHOE, SDELE are zero
        IREA(IF,1:NPETAS)=RHOE(1:NPETAS)
        IDEA(IF,1:NPETAS)=SDELE(1:NPETAS)
        IF (SECOND) THEN
          DO K=1,NPETAS
            DO JJ=K,NPETAS
              I2REA(IF,JJ,K)=RHOE2(JJ,K)
              I2DEA(IF,JJ,K)=SDELE2(JJ,K)
            END DO
          END DO
        END IF
      END DO
!    
! ADD variable information, shortest infusion first
      IF (RATE /= XR) THEN  
        IDA(IF)=SDEL
      ELSE     
        IAA(IF)=RHO*SDEL
        IF (.NOT. NOETAS) THEN
          IAEA(IF,1:NPETAS)=RHOE(1:NPETAS)*SDEL+RHO*SDELE(1:NPETAS) 
          IF (SECOND) THEN
            DO K=1,NPETAS
              DO J=K,NPETAS
                I2AEA(IF,J,K)=RHOE(K)*SDELE(J)+RHOE2(J,K)*SDEL+RHO*SDELE2(J,K)+RHOE(J)*SDELE(K)
              END DO
            END DO
          END IF
        END IF
      END IF
!        
! Remaining infusions shortest to longest
      IF (MI /= 0) THEN  
        DO J=1,MI
          IX=IF
          IF=INEXT(IX)
          IF (RATE == XR) THEN 
            IAA(IF)=IAA(IX)+DEL1*RHO
            IF (NOETAS) CYCLE
            IAEA(IF,1:NPETAS)=IAEA(IX,1:NPETAS)+DEL1*RHOE(1:NPETAS)
            IF (SECOND) THEN
              DO K=1,NPETAS
                DO JJ=K,NPETAS
                  I2AEA(IF,JJ,K)=I2AEA(IX,JJ,K)+DEL1*RHOE2(JJ,K)
                END DO
              END DO
            END IF
          ELSE
            IDA(IF)=IDA(IX)+DEL1
          END IF  
        END DO
      END IF
!
  999 RETURN
!
      END SUBROUTINE SSS
!