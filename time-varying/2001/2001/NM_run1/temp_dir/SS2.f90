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
!               FEB/2010 - CHANGED SIZES TO PRSIZES
!               FEB/2011 - INTEGRATED 7.2BETA5.8B MODIFICATIONS
!
!----------------------------- SS2.f90 ----------------------------------------------
!
! SUBROUTINE SS
!
! DESCRIPTION : One compartment model with absorption and with urine compartment
!
! ARGUMENTS   : NONE
!
! CALLED BY   : SSS - Supervisor of steady state
!
! CALLS       : NONE
!
! ALGORITHM   : - If SSID /= 0, then Normal entry; Else
!                 - Set SSID=2 and RETURN
!               - Normal entry
!                 - Central compartment
!                 - Check for errors in basic PK parameters
!               - Dose is to drug depot
!                 - Check for errors in basic PK parameters
!                 - Compute multiple infusions
!                 - Else, compute constant infusion
!                 - Compute multiple bolus doses
!
! MODULES USED: PRSIZES,NMPRD_INT,PRCOM_INT,PRMOD_INT,PROCM_INT,PRCM_REAL,PRCOM_REAL,
!               PROCM_REAL,NMPRD_CHAR,PRCOM_LOG 
!
! CONTAINS    : NONE
!
! LOCAL'S     : AA,AAE,AAE2,DEN,DENE,DENE2,EDEL,EDELE,EDELE2,EX,EXD,EXDE,EXDE2,EXE,
!               EXE2,G2,GG,I,J,K,KB,KBE,KBE2,L,OMEDEL,OMEX,PE,PR,PRE,PRE2,Q1,Q1E,Q1E2,
!               TEMP,TEMPE,TEMPE2,U,UE,UE2,V,VE,VE2,W,WE,WE2,X,XE,XE2
!
!---------------------------- END OF HEADER -----------------------------------------
!
      SUBROUTINE SS
!      
      USE PRSIZES,    ONLY: ISIZE,DPSIZE,PE   
! INTEGER
      USE NMPRD_INT,  ONLY: IERPRD
      USE PRCOM_INT,  ONLY: IDO,IDC,ITSC,NPETAS,SSID,SSC
      USE PRMOD_INT,  ONLY: IKE,IKA
      USE PROCM_INT,  ONLY: IDXETA
! REAL       
      USE PRCM_REAL,  ONLY: SSAE2,SSRE2,SDELE2,RHOE2
      USE PRCOM_REAL, ONLY: D2DTE,DT,DTE,RHO,RHOE,SDEL,SDELE,SSA,SSAE,SSR,SSRE,G3, &
                            ZERO,ONE
      USE PROCM_REAL, ONLY: AMNT,DAETA,D2AETA
! CHARACTER
      USE NMPRD_CHAR, ONLY: ETEXT
! LOGICAL       
      USE PRCOM_LOG,  ONLY: NOETAS,SECOND     
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
!     COMMON /PRCOM1/ NOETAS,SECOND
!     COMMON /PRCOM2/ IBF,IRR,IS,ID,ITSC,IFR,ILAG
!     COMMON /PRCOM3/ ITRANS,IRGG,IREV,NPETAS,NPEPS
!     COMMON /PRCOM4/ G3,HH,DELTA,DT,DTE
!     COMMON /PRCOM4/ YMULT,ZERO,ONE,XR,XD,TSTART,DTSTAR
!     COMMON /PRCOM4/ DDELTA,D2DELT,ADTE,D2ADTE
!     COMMON /PRCOM6/ IA,IAA,IAEA,IRA,IREA,IDA,IDEA,R,RE
!     COMMON /PRCOM6/ RHO,RHOE,SDEL,SDELE,SSA,SSAE,SSR,SSRE
!     COMMON /PRCOM6/ SAMT,SDEL1
!     COMMON /PRCOM6/ I2AEA,I2REA,I2DEA,R2E,D2DTE,D2TSTA
!     COMMON /PRCM6A/ SSAE2,SSRE2,SDELE2,RHOE2
!     COMMON /PRCOM7/ ADVID,SSID
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
!     INTEGER IRGG,IREV,ITRANS,NPETAS,NPEPS
!     INTEGER JCONT,JTIME,JEVENT,JAMT,JRATE,JSS,JDELTA
!     INTEGER JCOMPT,JCOMPF,JERROR
!     INTEGER IDC,IDO,NP,NBP,SSC,KREC,JMORE,JDUM
!     INTEGER ISV(PC),IBF(PC),IRR(PC)
!     INTEGER IINST(PC),ITURN(PC),ITSC,IFR,ILAG(PC),IS(PC),ID(PC)
!     INTEGER ADVID,SSID,MAXKF,IFORM(PG+1),YFORM,MAXIC
!     LOGICAL NOETAS,SECOND
!     COMMON /PROCM4/ A,DAETA,D2AETA
!     DOUBLE PRECISION A,DAETA,D2AETA
!     DIMENSION A(PC),DAETA(PC,PE),D2AETA(PC,PE,PE)
!     COMMON /PRMON3/ IOUT,IV,IKE,IKA
!     INTEGER IOUT,IKE,IKA,IV
!------------------------------------------------------------------------------------
!
! Local Variables  
!
      INTEGER(KIND=ISIZE) :: I,J,K,L
!     
      REAL(KIND=DPSIZE)   :: AA(2),AAE(2,PE),AAE2(2,PE,PE),DEN,DENE(PE),DENE2,EDEL(2), &
                             EDELE(2,PE),EDELE2(2),EX(2),EXD(2),EXDE(2,PE),EXDE2(2),   &  
                             EXE(2,PE),EXE2(2),G2,GG,KB,KBE(PE),KBE2(PE,PE),OMEDEL(2), &
                             OMEX(2),PR(2),PRE(2,PE),PRE2(2),Q1(2),Q1E(2,PE),Q1E2(2),  &
                             TEMP,TEMPE(PE),TEMPE2,U(2),UE(2,PE),UE2(2,PE,PE),V(2),    &
                             VE(2,PE),VE2(2),W(2),WE(2,PE),WE2(2),X(2),XE(2,PE),XE2(2)
!    
      G2(I,J,K)=G3(I,IDXETA(J)+1,IDXETA(K)+1) ! Statement function
      GG(I,J)=G3(I,IDXETA(J-1)+1,1)
!
! Steady state variables:
! From PRED: 
! SSC     = Compartment number for dose
! DT      = Captial delta (dosing interval)
! SDEL    = Small delta (duration of infusion)
! SSA     = Dose amount
! SSR     = Rate of constant infusion
! RHO     = Rate of multiple infusion
! RHOE    = Partial derivative of RHO wrt ETA(K)
! DTE(K)  = Partial derivative of DT wrt ETA(K)
! SDELE(K)= Partial derivative of SDEL wrt ETA(K)
! SSAE(K) = Partial derivative of SSA wrt ETA(K)
! SSRE(K) = Partial derivative of SSR wrt ETA(K)
!
      IF (SSID == 0) THEN
        SSID=2; GO TO 999  
      END IF
!      
! Normal entry; Central compartment
      AA(2)=GG(ITSC,1)*GG(IKE,1)
!
      IF (AA(2) <= ZERO) THEN   ! Check for errors in basic PK parameters
        IF (GG(ITSC,1) <= ZERO) THEN
          ETEXT(2)='PK PARAMETER FOR TIME SCALE IS NON-POSITIVE'
        ELSE
          ETEXT(2)='PK PARAMETER FOR K IS NON-POSITIVE'
        END IF
        IERPRD=1; GO TO 999  
      END IF
!          
      IF (.NOT. NOETAS) THEN
        DO K=1,NPETAS
          AAE(2,K)=GG(ITSC,K+1)*GG(IKE,1)+GG(ITSC,1)*GG(IKE,K+1)
        END DO  
        IF (SECOND) THEN
          DO K=1,NPETAS
            DO J=K,NPETAS
              AAE2(2,J,K) = GG(ITSC,J+1)*GG(IKE,K+1)+GG(ITSC,1)*G2(IKE,J,K)  &
                           +GG(ITSC,K+1)*GG(IKE,J+1)+G2(ITSC,J,K)*GG(IKE,1)
            END DO
          END DO
        END IF     
      END IF
!  
      I=2
      IF (SSC /= IDO) THEN      ! Dose is to drug depot
        I=1
        IF (GG(IKA,1)-GG(IKE,1) == ZERO) THEN
          ETEXT(2)='PK PARAMETERS FOR KA AND K ARE EQUAL'
          IERPRD=1; GO TO 999  
        END IF
        KB=GG(IKA,1)/(GG(IKA,1)-GG(IKE,1))
        AA(1)=GG(ITSC,1)*GG(IKA,1)
!
        IF (AA(1) <= ZERO) THEN ! Check for errors in basic PK parameters
          IF (GG(ITSC,1) <= ZERO) THEN
            ETEXT(2)='PK PARAMETER FOR TIME SCALE IS NON-POSITIVE'
          ELSE
            ETEXT(2)='PK PARAMETER FOR KA IS NON-POSITIVE'
          END IF
          IERPRD=1; GO TO 999 
        END IF
        IF (AA(1) == AA(2)) THEN
          ETEXT(2)='PK PARAMETERS FOR KA AND K ARE EQUAL'
          IERPRD=1; GO TO 999 
        END IF
!       
        IF (.NOT. NOETAS) THEN
          DO K=1,NPETAS
            AAE(1,K) = GG(ITSC,K+1)*GG(IKA,1)+GG(ITSC,1)*GG(IKA,K+1)
            KBE(K)   = (GG(IKA,K+1)-KB*(GG(IKA,K+1)-GG(IKE,K+1)))/(GG(IKA,1)-GG(IKE,1))
          END DO
!          
          IF (SECOND) THEN
            DO K=1,NPETAS
              DO J=K,NPETAS
                AAE2(1,J,K) = GG(ITSC,J+1)*GG(IKA,K+1)+GG(ITSC,1)*G2(IKA,J,K)+GG(ITSC,K+1)  & 
                             *GG(IKA,J+1)+G2(ITSC,J,K)*GG(IKA,1)
                KBE2(J,K)   = (G2(IKA,J,K)-KBE(J)*(GG(IKA,K+1)-GG(IKE,K+1))-KB*(G2(IKA,J,K) & 
                             -G2(IKE,J,K))-KBE(K)*(GG(IKA,J+1)-GG(IKE,J+1)))/(GG(IKA,1)-GG(IKE,1))
              END DO
            END DO
          END IF
        END IF
      END IF
!        
      IF (DT == ZERO .OR. RHO /= ZERO) THEN       
        DO J=I,2
          U(J)=SSR/AA(J)
          IF (NOETAS) CYCLE
          DO K=1,NPETAS
            UE(J,K)=(SSRE(K)-U(J)*AAE(J,K))/AA(J)
          END DO
          IF (SECOND) THEN
            DO K=1,NPETAS
              DO L=K,NPETAS
                UE2(J,L,K)=(SSRE2(L,K)-UE(J,L)*AAE(J,K)-U(J)*AAE2(J,L,K)-UE(J,K)*AAE(J,L))/AA(J)
              END DO
            END DO
          END IF
        END DO
!        
        IF (DT /= ZERO) THEN    ! Multiple infusions
          DO J=I,2
            EX(J)=EXP(-AA(J)*SDEL)
            OMEX(J)=ONE-EX(J)
            EXD(J)=EXP(-AA(J)*DT)
            EDEL(J)=EX(J)*EXD(J)
            OMEDEL(J)=ONE-EDEL(J)
            IF (OMEDEL(J) == ZERO) THEN
              ETEXT(1)='ELIMINATION FROM SOME COMPARTMENT IS NEGLIGIBLE.'
              ETEXT(2)='STEADY-STATE IS NOT ACHIEVABLE.'
              IERPRD=1; GO TO 999  
            END IF
            V(J) = RHO/AA(J)
            Q1(J)= OMEX(J)/OMEDEL(J)
            PR(J)= V(J)*Q1(J)
            PR(J)= PR(J)*EXD(J)
          END DO
          AMNT(SSC)=AMNT(SSC)+U(I)+PR(I)
!         
          IF (SSC /= IDO) THEN
            W(1:2)=U(1:2)*KB  
            X(1:2)=V(1:2)*KB
            TEMP = (W(1)+X(1)*EXD(1)*Q1(1))*(EDEL(2)-EDEL(1))+U(2)+W(1)*EDEL(1)-W(2) & 
                  *EDEL(2)+V(2)*EXD(2)-X(2)*EDEL(2)+X(1)*(EXD(2)-EXD(1)+EDEL(1))
            DEN  = TEMP/OMEDEL(2)
            AMNT(IDO) = AMNT(IDO)+DEN
          END IF
!          
          IF (NOETAS) GO TO 999
!          
          DO K=1,NPETAS
            DO J=I,2
              EXE(J,K)  = -EX(J)*(AA(J)*SDELE(K)+AAE(J,K)*SDEL)
              EXDE(J,K) = -EXD(J)*(AA(J)*DTE(K)+AAE(J,K)*DT)
              EDELE(J,K)= EX(J)*EXDE(J,K)+EXE(J,K)*EXD(J)
              VE(J,K)   = (RHOE(K)-V(J)*AAE(J,K))/AA(J)
              Q1E(J,K)  = (-EXE(J,K)+Q1(J)*EDELE(J,K))/OMEDEL(J)
              PRE(J,K)  = VE(J,K)*Q1(J)+V(J)*Q1E(J,K)
              PRE(J,K)  = PRE(J,K)*EXD(J)+V(J)*Q1(J)*EXDE(J,K)
            END DO
            DAETA(SSC,K)=DAETA(SSC,K)+UE(I,K)+PRE(I,K)
            IF (SSC == IDO) CYCLE
            DO J=1,2
              WE(J,K)=U(J)*KBE(K)+UE(J,K)*KB
              XE(J,K)=V(J)*KBE(K)+VE(J,K)*KB
            END DO
            TEMPE(K) = (WE(1,K)+XE(1,K)*EXD(1)*Q1(1)+X(1)*EXDE(1,K)*Q1(1)+X(1)*EXD(1)     & 
                      *Q1E(1,K))*(EDEL(2)-EDEL(1))
            TEMPE(K) = TEMPE(K)+(W(1)+X(1)*EXD(1)*Q1(1))*(EDELE(2,K)-EDELE(1,K))+UE(2,K)  & 
                      +WE(1,K)*EDEL(1)+W(1)*EDELE(1,K)-WE(2,K)*EDEL(2)-W(2)*EDELE(2,K)    &
                      +VE(2,K)*EXD(2)+V(2)*EXDE(2,K)-XE(2,K)*EDEL(2)-X(2)*EDELE(2,K)      &
                      +XE(1,K)*(EXD(2)-EXD(1)+EDEL(1))+X(1)*(EXDE(2,K)-EXDE(1,K)+EDELE(1,K))
            DENE(K) = (TEMPE(K)+DEN*EDELE(2,K))/OMEDEL(2)
            DAETA(IDO,K)=DAETA(IDO,K)+DENE(K)
          END DO
!                    
          IF (SECOND) THEN
            DO K=1,NPETAS
              DO L=K,NPETAS
                DO J=I,2
                  EXE2(J)  = -EXE(J,L)*(AA(J)*SDELE(K)+AAE(J,K)*SDEL)-EX(J)*(AAE(J,L)*SDELE(K)   & 
                             +AAE2(J,L,K)*SDEL+AA(J)*SDELE2(L,K)+AAE(J,K)*SDELE(L))
                  EXDE2(J) = -EXDE(J,L)*(AA(J)*DTE(K)+AAE(J,K)*DT)-EXD(J)*(AAE(J,L)*DTE(K)       & 
                             +AAE2(J,L,K)*DT)-EXD(J)*(AA(J)*D2DTE(L,K)+AAE(J,K)*DTE(L))
                  EDELE2(J)= EXE(J,L)*EXDE(J,K)+EXE2(J)*EXD(J)+EX(J)*EXDE2(J)+EXE(J,K)*EXDE(J,L)
                  VE2(J)   = (RHOE2(L,K)-VE(J,L)*AAE(J,K)-V(J)*AAE2(J,L,K)-VE(J,K)*AAE(J,L))/AA(J)
                  Q1E2(J)  = (-EXE2(J)+Q1E(J,L)*EDELE(J,K)+Q1(J)*EDELE2(J)                       &
                             +Q1E(J,K)*EDELE(J,L))/OMEDEL(J)
                  PRE2(J)  = (VE2(J)*Q1(J)+VE(J,L)*Q1E(J,K)+VE(J,K)*Q1E(J,L)                     &
                             +V(J)*Q1E2(J))*EXD(J)+(VE(J,K)*Q1(J)+V(J)*Q1E(J,K))*EXDE(J,L)       &  
                             +VE(J,L)*Q1(J)*EXDE(J,K)+V(J)*Q1E(J,L)*EXDE(J,K)+V(J)*Q1(J)*EXDE2(J)
                END DO
                D2AETA(SSC,L,K) = D2AETA(SSC,L,K)+UE2(I,L,K)+PRE2(I)
                IF (SSC == IDO) CYCLE
                DO J=1,2
                  WE2(J) = UE(J,L)*KBE(K)+UE2(J,L,K)*KB+U(J)*KBE2(L,K)+UE(J,K)*KBE(L)
                  XE2(J) = VE(J,L)*KBE(K)+VE2(J)*KB+V(J)*KBE2(L,K)+VE(J,K)*KBE(L)
                END DO
                TEMPE2 = (WE2(1)+XE2(1)*EXD(1)*Q1(1)+XE(1,K)*EXDE(1,L)*Q1(1)                     &
                        +XE(1,K)*EXD(1)*Q1E(1,L)+XE(1,L)*EXDE(1,K)*Q1(1)+X(1)*EXDE2(1)*Q1(1)     &
                        +X(1)*EXDE(1,K)*Q1E(1,L)+XE(1,L)*EXD(1)*Q1E(1,K)+X(1)*EXDE(1,L)*Q1E(1,K) &
                        +X(1)*EXD(1)*Q1E2(1))*(EDEL(2)-EDEL(1))+(WE(1,K)+XE(1,K)*EXD(1)*Q1(1)    &
                        +X(1)*EXDE(1,K)*Q1(1)+X(1)*EXD(1)*Q1E(1,K))*(EDELE(2,L)-EDELE(1,L))
                TEMPE2 = TEMPE2+(WE(1,L)+XE(1,L)*EXD(1)*Q1(1)+X(1)*EXDE(1,L)*Q1(1)               &
                        +X(1)*EXD(1)*Q1E(1,L))*(EDELE(2,K)-EDELE(1,K))+(W(1)                     &
                        +X(1)*EXD(1)*Q1(1))*(EDELE2(2)-EDELE2(1))+UE2(2,L,K)                     &
                        +WE2(1)*EDEL(1)+WE(1,L)*EDELE(1,K)+WE(1,K)*EDELE(1,L)                    &
                        +W(1)*EDELE2(1)-WE2(2)*EDEL(2)-WE(2,L)*EDELE(2,K)-WE(2,K)*EDELE(2,L)     &
                        -W(2)*EDELE2(2)+VE2(2)*EXD(2)+VE(2,L)*EXDE(2,K)+VE(2,K)*EXDE(2,L)        &
                        +V(2)*EXDE2(2)-XE2(2)*EDEL(2)-XE(2,L)*EDELE(2,K)-XE(2,K)*EDELE(2,L)      &
                        -X(2)*EDELE2(2)+XE2(1)*(EXD(2)-EXD(1)+EDEL(1))+XE(1,K)*(EXDE(2,L)        &
                        -EXDE(1,L)+EDELE(1,L))+XE(1,L)*(EXDE(2,K)-EXDE(1,K)+EDELE(1,K))          &
                        +X(1)*(EXDE2(2)-EXDE2(1)+EDELE2(1))
                D2AETA(IDO,L,K) = D2AETA(IDO,L,K)+(TEMPE2+DENE(L)*EDELE(2,K)+DEN*EDELE2(2)       &
                                 +DENE(K)*EDELE(2,L))/OMEDEL(2)
              END DO
            END DO
          END IF
        ELSE      ! Constant infusion
          AMNT(SSC)=AMNT(SSC)+U(I)       
          IF (.NOT. NOETAS) THEN
            FORALL (K=1:NPETAS) DAETA(SSC,K)=DAETA(SSC,K)+UE(I,K)
            IF (SECOND) THEN
              DO K=1,NPETAS
                DO J=K,NPETAS
                  D2AETA(SSC,J,K)=D2AETA(SSC,J,K)+UE2(I,J,K)
                END DO
              END DO
            END IF
          END IF
          IF (SSC /= IDO) THEN
            AMNT(IDO)=AMNT(IDO)+U(2)
            IF (.NOT. NOETAS) THEN
              FORALL (K=1:NPETAS) DAETA(IDO,K)=DAETA(IDO,K)+UE(2,K)
              IF (SECOND) THEN
                DO K=1,NPETAS
                  DO J=K,NPETAS
                    D2AETA(IDO,J,K)=D2AETA(IDO,J,K)+UE2(2,J,K)
                  END DO
                END DO
              END IF
            END IF
          END IF  
        END IF  
        GO TO 999
      END IF  
!        
! Multiple bolus doses
      DO J=I,2
        EX(J)=EXP(-AA(J)*DT)
        OMEX(J)=ONE-EX(J)
        IF (OMEX(J) == ZERO) THEN
          ETEXT(1)='ELIMINATION FROM SOME COMPARTMENT IS NEGLIGIBLE.'
          ETEXT(2)='STEADY-STATE IS NOT ACHIEVABLE.'
          IERPRD=1; GO TO 999  
        END IF
        Q1(J)=SSA/OMEX(J)
      END DO
      AMNT(SSC)=AMNT(SSC)+Q1(I)
      IF (SSC == IDC) THEN
        DEN=ONE/(OMEX(1)*OMEX(2))
        AMNT(IDO)=AMNT(IDO)+KB*SSA*(EX(2)-EX(1))*DEN
      END IF
!      
      IF (.NOT. NOETAS) THEN
        DO K=1,NPETAS
          DO J=I,2
            EXE(J,K) = -EX(J)*(AA(J)*DTE(K)+AAE(J,K)*DT)
            Q1E(J,K) = (SSAE(K)+Q1(J)*EXE(J,K))/OMEX(J)
          END DO
          DAETA(SSC,K)=DAETA(SSC,K)+Q1E(I,K)
          IF (SSC == IDC) THEN
            DENE(K) = (EXE(1,K)*OMEX(2)+OMEX(1)*EXE(2,K))*DEN*DEN
            DAETA(IDO,K) = DAETA(IDO,K)+(KBE(K)*SSA*(EX(2)-EX(1))            &
                          +KB*SSAE(K)*(EX(2)-EX(1))                          &  
                          +KB*SSA*(EXE(2,K)-EXE(1,K)))*DEN                   &
                          +KB*SSA*(EX(2)-EX(1))*DENE(K)
          END IF
        END DO
!        
        IF (SECOND) THEN
          DO K=1,NPETAS
            DO L=K,NPETAS
              DO J=I,2
                EXE2(J) = -EXE(J,L)*(AA(J)*DTE(K)                            &
                          +AAE(J,K)*DT)-EX(J)*(AAE(J,L)*DTE(K)               &
                          +AAE2(J,L,K)*DT+AA(J)*D2DTE(L,K)                   &
                          +AAE(J,K)*DTE(L))
                Q1E2(J) = (SSAE2(L,K)+Q1E(J,L)*EXE(J,K)+Q1(J)*EXE2(J)        &
                         +Q1E(J,K)*EXE(J,L))/OMEX(J)
              END DO
              D2AETA(SSC,L,K) = D2AETA(SSC,L,K)+Q1E2(I)
              IF (SSC == IDC) THEN
                DENE2 = (EXE(1,K)*OMEX(2)+OMEX(1)*EXE(2,K))*2*DEN*DENE(L)    &
                       +(EXE2(1)*OMEX(2)-EXE(1,L)*EXE(2,K)-EXE(1,K)*EXE(2,L) &
                       +OMEX(1)*EXE2(2))*DEN*DEN
                D2AETA(IDO,L,K) = D2AETA(IDO,L,K)                            &
                                 +(KBE2(L,K)*SSA*(EX(2)-EX(1))               &
                                 +KBE(K)*SSAE(L)*(EX(2)-EX(1))               &
                                 +KBE(K)*SSA*(EXE(2,L)-EXE(1,L))             &
                                 +KBE(L)*SSAE(K)*(EX(2)-EX(1))               &
                                 +KB*SSAE2(L,K)*(EX(2)-EX(1))                &
                                 +KB*SSAE(K)*(EXE(2,L)-EXE(1,L))             &
                                 +KBE(L)*SSA*(EXE(2,K)-EXE(1,K))             &
                                 +KB*SSAE(L)*(EXE(2,K)-EXE(1,K))             &
                                 +KB*SSA*(EXE2(2)-EXE2(1)))*DEN              &
                                 +(KBE(K)*SSA*(EX(2)-EX(1))                  &
                                 +KB*SSAE(K)*(EX(2)-EX(1))                   &
                                 +KB*SSA*(EXE(2,K)-EXE(1,K)))*DENE(L)        &
                                 +KB*SSA*(EX(2)-EX(1))*DENE2                 &
                                 +(KBE(L)*SSA*(EX(2)-EX(1))                  &
                                 +KB*SSAE(L)*(EX(2)-EX(1))                   &
                                 +KB*SSA*(EXE(2,L)-EXE(1,L)))*DENE(K)
              END IF
            END DO
          END DO
        END IF
      END IF      
!     
  999 RETURN 
!                        
      END SUBROUTINE SS
