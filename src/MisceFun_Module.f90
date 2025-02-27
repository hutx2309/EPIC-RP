MODULE MisceFun_Module
USE PARM
IMPLICIT NONE

CONTAINS
! ------------------------------------1. Check Files ---------------------------------------
SUBROUTINE OPENV(NUM,FNAM,ID,NMSO)
      !     EPIC1102      
      !     VERIFIES THE EXISTENCE OF A FILE BEFORE OPENING IT
      IMPLICIT NONE
      ! local variables  
      INTEGER:: JRT, ID, NMSO, NUM 
      CHARACTER(LEN = 80)::FNAM,ADIR
      CHARACTER(LEN = 160)::FNM
      LOGICAL::XMIS

      JRT=0
      SELECT CASE(ID)
            CASE(0) 
                  ADIR='G:\EEPIC_PilotStudy\LM\'
            CASE(2)
                  ADIR='C:\SITE\'
            CASE(3)
                  ADIR='C:\SOIL\'
            CASE(4)
                  ADIR='C:\OPSC\'
            CASE DEFAULT
                  FNM=FNAM
             
      END SELECT
      
      FNM=ADJUSTR(ADIR)//TRIM(ADJUSTL(FNAM))
      FNM=TRIM(ADJUSTL(FNM)) 
      
      INQUIRE(FILE=FNM,EXIST=XMIS)

      IF(XMIS)THEN
            OPEN(NUM,FILE=FNM)
      ELSE
            IF(runMode==0)THEN
                  WRITE(*,'(/A/)')'File '//TRIM(FNM)//' IS MISSING.'
                  ERROR STOP
            ELSE
                  WRITE(NMSO,'(A,A8,1X,A96,A)')' !!!!! ',siteName,TRIM(FNM),&
                        ' IS MISSING.'
                  ERROR STOP
            END IF
      END IF 
END SUBROUTINE OPENV    
! ------------------ 1.1 Split real numbers ----------------------------
REAL FUNCTION split(X)
      REAL:: X
      split = AINT(X) ! AINT get the integer part of a real number      
      X = X - split
      
END FUNCTION split
    
! ------------------- 1.2 Split integers --------------------------------
    
SUBROUTINE split_int(X1, X2)
      IMPLICIT NONE
      INTEGER:: X1, X2
      REAL:: XX
      XX=FLOAT(X1)
      XX=XX*.1
      X2=INT(XX)
      X1=X1-X2*10

END SUBROUTINE split_int

! ----------------- 2. S Curve coefficients ---------------------------
SUBROUTINE Coeff_Scurve(X1,X2,X3,X4)
      !     EPIC1102
      !     THIS SUBPROGRAM COMPUTES S CURVE PARMS GIVEN 2 (X,Y) POINTS.
      IMPLICIT NONE
      REAL:: X1, X2, X3, X4, XX
      
      XX=LOG(X3/X1-X3)
      X2=(XX-LOG(X4/X2-X4))/(X4-X3)
      X1=XX+X3*X2
END SUBROUTINE Coeff_Scurve

! ----------------- 3. Leap year check -------------------------------------
SUBROUTINE Leap_Yr_Check(year, leap, leapYearFlag)

      IMPLICIT NONE
      ! global variables
      INTEGER:: year, leap, leapYearFlag
      ! leapYearFlag = 0 : leap year is considered 
      ! for common year, leap = 1; for leap year, leap = 0
      leap = 1    
      IF(MOD(year,100)/=0)THEN
          IF(MOD(year,4)/=0)RETURN           ! common year
      ELSE
          IF(MOD(year,400)/=0)RETURN         ! common year 
      END IF
      leap = leapYearFlag
 
END SUBROUTINE Leap_Yr_Check

! -------- 4.1 Random number generator ------------------
REAL FUNCTION Generate_Random(J)
      IMPLICIT NONE 
      INTEGER:: J,K
      K=IX(J)/127773
      IX(J)=16807*(IX(J)-K*127773)-K*2836
      IF(IX(J)<0)IX(J)=IX(J)+2147483647
      Generate_Random=IX(J)*4.656612875D-10
END FUNCTION Generate_Random

! ----------------------- 4.2 Generate numbers from triangular distribution -------------------    
REAL FUNCTION Triangle_Num(BLM,QMN,UpperLimit,KK)
      !     EPIC1102
      !     THIS SUBPROGRAM GENERATES NUMBERS FROM A TRIANGULAR DISTRIBUTION
      !     GIVEN X AXIS POINTS AT START & END AND PEAK Y VALUE.
      IMPLICIT NONE
    
      REAL:: BLM, QMN, U3, RN, Y, UpperLimit, B1, B2, X1, AMN
      INTEGER:: KK
      
      U3=QMN-BLM
      RN=Generate_Random(randID(KK))
      Y=2./(UpperLimit-BLM)
      B2=UpperLimit-QMN
      B1=RN/Y
      X1=Y*U3/2.
      IF(RN>X1)THEN
            Triangle_Num=UpperLimit-SQRT(B2*B2-2.*B2*(B1-.5*U3))
      ELSE
            Triangle_Num=SQRT(2.*B1*U3)+BLM
      END IF    
      IF(KK/=7.AND.KK/=4)RETURN
      AMN=(UpperLimit+QMN+BLM)/3.
      Triangle_Num=Triangle_Num*QMN/AMN
      IF(Triangle_Num>=1.)Triangle_Num=.99
END FUNCTION Triangle_Num

! ------------------- 5. Shuffle Data randomly -------------------------
SUBROUTINE Shuffle
     !     EPIC1102
     !     THIS SUBPROGRAM SHUFFLES DATA RANDOMLY. (BRATLEY,FOX,SCHRAGE,P.34)
      IMPLICIT NONE
      ! local variables
      INTEGER:: I, II, K
      REAL:: RN
      DO I=20,2,-1
          II=randID(I)
          RN=Generate_Random(21)
          K=INT(I*RN+1)
          randID(I)=randID(K)
          randID(K)=II
      END DO
END SUBROUTINE Shuffle

! ------------- 6.1 Sort Intergers into ascending order ---------------------
SUBROUTINE Sort_IntNum(NZ, NY, M) 
      !     EPIC1102
      !     THIS SUBPROGRAM SORTS INTEGERS INTO ASCENDING ORDER USING
      !     RIPPLE SORT
      IMPLICIT NONE
	! local variables potential problems with local variables (SAVE attribute)
      INTEGER:: M, NB, K, I, J, MK, K1, N1, NY(M), NZ(M)
 
      NB=M-1
      J=M
      DO I=1,NB
          J=J-1
          MK=0
          DO K=1,J
              K1=K+1
              IF(NZ(NY(K))<=NZ(NY(K1)))CYCLE
              N1=NY(K1)
              NY(K1)=NY(K)
              NY(K)=N1
              MK=1
          END DO
          IF(MK==0)EXIT
      END DO
END SUBROUTINE Sort_IntNum

! ------------------ 6.2 Sort real numbers into ascending order ---------------------
SUBROUTINE Sort_RealNum(X, NXX, M)
      !     THIS SUBPROGRAM SORTS REAL NUMBERS INTO ASCENDING ORDER USING
      !     RIPPLE SORT      potential problems with local variables (SAVE attribute)
      INTEGER:: M, NXX(200), NB, I, J, K, KP1, MK, N1
      REAL:: X(200)
      NB=M-1
      J=M
      DO I=1,NB
          J=J-1
          MK=0
          DO K=1,J
              KP1=K+1
              IF(X(NXX(K))>=X(NXX(KP1)))CYCLE
              N1=NXX(KP1)
              NXX(KP1)=NXX(K)
              NXX(K)=N1
              MK=1
          END DO
          IF(MK==0)EXIT
      END DO
END SUBROUTINE Sort_RealNum

! ---------------- 7. Acout: output variables --------------------------------------- 
SUBROUTINE ACOUT(XX, QQ, GKG)
	!     EPIC1102
      REAL:: GKG , QQ
      REAL(KIND = 8):: XX,X1,X2,X3,XI, N2, N1
      INTEGER:: I
      IF(XX<1.E-10)RETURN
      X1=.1*GKG*XX/(QQ+1.E-10)
      IF(X1<1000.)THEN
            DO I=1,4
                  N2=X1
                  IF(N2>0)EXIT
                  X1=X1*1.D3
            END DO
            I=MIN(I,4)
            XI=I
      ELSE
            XI=0.
            X1=.001*X1
      END IF
      N2=XX+.5
      X2=N2
      N1=X1+.5
      X1=N1
      X3=1.D-3*X1+1.D-4*XI
      XX=X2+X3
END SUBROUTINE ACOUT

! ------------------- 8. Adjust SCS runoff curve numbers ------------------------
SUBROUTINE Adjust_CN(CNII, X1)
      !     EPIC1102
      !     THIS SUBPROGRAM ADJUSTS THE 2 CONDITION SCS RUNOFF CURVE
      !     NUMBER FOR WATERSHED SLOPE AND COMPUTES CN1 AND CN3.
      IMPLICIT NONE
 
      ! local variables
      INTEGER:: J,L1, N1  
      REAL:: CNII, C2, X1, S2, SUM, TOT, RTO, Z1, Z2, adjCN 
     
      C2=100.-CNII          ! Eq. 33a         CNII = CN_2s
      CN3=CNII*EXP(.006729*C2)     
      adjCN=(1.-2.*EXP(-13.86*uplandSteep))/3.  ! Eq 2.3 in EPIC theoretical doc              
      X1=(CN3-CNII)*adjCN+CNII
      C2=100.-X1             ! Eq. 33a
      CN1=MAX(.4*X1,X1-20.*C2/(C2+EXP(2.533-.0636*C2))) ! Eq. 32
      SMX=254.*(100./CN1-1.)   ! SMX == s in Eq. 30 in APEX-doc Runoff volume the SCS curve number method
      CN3=X1*EXP(.006729*C2)   ! X1 == CN_2s in Eq. 33 in APEX-doc
      S3=254.*(100./CN3-1.)
      S2=254.*(100./X1-1.)
      SUM=0.
      TOT=0.
      DO J=1,Actual_SoilLayers
          ISL=Layer_ID(J)
          IF(Z(ISL)>1.)GO TO 2
          SUM=SUM+fieldCapacity(ISL)-Wilt_Point(ISL)
          TOT=TOT+Porosity(ISL)-Wilt_Point(ISL)
      END DO
      GO TO 3
    2 L1=Layer_ID(J-1)
      RTO=(1.-Z(L1))/(Z(ISL)-Z(L1))
      SUM=SUM+RTO*(fieldCapacity(ISL)-Wilt_Point(ISL))
      TOT=TOT+RTO*(Porosity(ISL)-Wilt_Point(ISL))
    3 N1=INT(100.+S_Curve(30,2)*(TOT/SUM-1.)+.5)
      S_Curve(4,1)=1.-S2/SMX+S_Curve(30,1)
      S_Curve(4,2)=1.-S3/SMX+N1
      Z1=split(S_Curve(4,1))
      Z2=split(S_Curve(4,2))
      CALL Coeff_Scurve(S_Curve(4,1),S_Curve(4,2),Z1,Z2) 
END SUBROUTINE Adjust_CN 

!--------------------------- 9. Command Line ----------------
SUBROUTINE Get_CommendLine(CM)
      !     THIS COMMAND LINE ROUTINE FOR MICROSOFT/UNIX/SALFORD
      IMPLICIT NONE
      INTEGER::NNARG, I  ! , NCHAR  not needed in gfortran
      CHARACTER (LEN = 80), DIMENSION(9):: CM

      CM=' '

      ! windows system
      ! NNARG=NARGS()
      ! linux system
      NNARG = IARGC()

      IF(NNARG>9)NNARG=9

      DO I=1,NNARG
          ! intel Fortran
          ! CALL GETARG(I-1,CM(I),NCHAR)
          ! gfortran
           CALL GETARG(I-1, CM(I)) 
      END DO
END SUBROUTINE Get_CommendLine

! ---------------------10. Cal Day of Year ----------------------
! given the month and days in the month cal day of year 
INTEGER FUNCTION Cal_DOY(Mon, I, leapY)    
      IMPLICIT NONE 
      INTEGER, INTENT(IN):: Mon, I, leapY
      
      Cal_DOY=NC(Mon)+I
      IF(Mon>2) Cal_DOY=Cal_DOY-leapY
END FUNCTION Cal_DOY

! ----------------------11. Cal Month --------------------------   
! Give the day of year, cal the month     
SUBROUTINE Cal_Mon(doy, MOX)  
      IMPLICIT NONE
      ! local variables 
      INTEGER:: M1, NDA, MOX, doy
     
      IF(doy>NC(2))THEN
            DO MOX = 2,12
                  M1= MOX + 1
                  NDA=NC(M1)-leapYr
                  IF(doy<=NDA) RETURN 
            END DO
      END IF
      MOX = 1
END SUBROUTINE Cal_Mon 
 
! ----------------------- 12. Cal Day of Month -----------------    
!     THIS SUBPROGRAM COMPUTES THE DAY OF THE MONTH, GIVEN THE MONTH AND
!     THE DAY OF THE YEAR.
INTEGER FUNCTION  Cal_DOM(doy, monx, leapYear) 
      IMPLICIT NONE
      INTEGER:: doy, monx, leapYear

      Cal_DOM = doy -NC(monx)
        
      IF(monx>2) Cal_DOM= Cal_DOM + leapYear
END FUNCTION Cal_DOM
! -------------------- 3. Heat Accumulation ---------------------
REAL FUNCTION Cal_Heat(J,K,BASE,NHS)
      ! EPIC1102
      ! THIS SUBPROGRAM ACCUMULATES HEAT UNITS FOR USE IN CPTHU.
      IMPLICIT NONE
      ! local variables
      INTEGER:: J, K, NHS 
      REAL:: BASE, TA, TGX
      
      Cal_Heat=0.
      MO=1
      DO DayOfYear=J,K
            CALL Cal_Mon(DayOfYear, MO)
            IF(JDHU<=366)THEN
                  CALL Cal_DL_MaxRd 
                  IF(DayLength<winterDormancy.AND.NHS==0)CYCLE
            END IF
            TA=Inted_DL_Rd(monAveT)
            TGX=TA-BASE
            IF(TGX>0.)Cal_Heat=Cal_Heat+TGX
      END DO
END FUNCTION Cal_Heat
! ------------------1. calcualte daily maximum solar ratidation ------------------
SUBROUTINE Cal_DL_MaxRd 
      !     EPIC1102
      !     THIS SUBPROGRAM COMPUTES DAY LENGTH & MAX SOLAR RADIATION AT THE
      !     EARTHS SURFACE.
      IMPLICIT NONE
 
      ! local variables
      INTEGER:: XI
      REAL:: SD, CH, H
      
      XI=DayOfYear
      SD=.4102*SIN((XI-80.25)/PIT)   ! 0.4102 radians = 23.5 degree   SD: the sun's declination angle   see EPIC-doc 2.195 P47
      CH=-YTN*TAN(SD)
      IF(CH>=1.)THEN
            H=0.
      ELSE
            IF(CH<=-1.)THEN
                  H=3.1416
            ELSE
                  H=ACOS(CH)
            END IF
      END IF
      DD=1.+.0335*SIN((XI+88.2)/PIT)      
      DayLength=7.72*H
      Delta_Daylength = DayLength-HR0      ! cal the daylenght change between day i and i-1
      HR0=DayLength
      maxRad=30.*DD*(H*YLTS*SIN(SD)+YLTC*COS(SD)*SIN(H))
     
END SUBROUTINE Cal_DL_MaxRd 

! ---------------- 2. Intepolate Montly Rd to Daily Rd -----------------------
REAL FUNCTION Inted_DL_Rd(DL_Rd)
   !  EPIC1102
   !  THIS SUBPROGRAM INTERPOLATES MONTHLY VALUES OF DAYLIGHT HOURS AND
   !  MAXIMUM SOLAR RADIATION TO PROVIDE DAILY VALUES.
      IMPLICIT NONE
    
      ! local variables  potential problems with local variables (SAVE attribute)
      INTEGER:: M1, N2, N1
      REAL:: DL_Rd(12), RTO, XX, X1, X2    
 
      M1=MO+1
      N2=NC(M1)
      IF(MO==2)N2=N2-leapYr
      N1=NC(MO)
      X1=DayOfYear-N1
      X2=N2-N1
      RTO=X1/X2
      IF(M1==13)M1=1
      XX=DL_Rd(M1)-DL_Rd(MO)
      !DRV=XX/X2    ! not used
      Inted_DL_Rd=XX*RTO + DL_Rd(MO)
  
END FUNCTION Inted_DL_Rd
! ------------------- 13. Integrate modified expoential equations -------------------
SUBROUTINE Integrate_ExpotentialEqs(WW,SUM)
      !     EPIC1102
      !     THIS SUBPROGRAM INTEGRATES THE MODIFIED EXPONENTIAL EQ.
      IMPLICIT NONE
      ! local variables
      REAL:: X1, X2, XY, DX, SUM, Y1, Y2, WW  
      X1=1.
      DX=.1
      SUM=0.
      Y1=0.
      !     WRITE(KW(1),3)X1,DX,Y1
      DO WHILE(DX>1.E-4)
            XY=0.
            DO WHILE(XY<.1)
                  X2=X1-DX
                  IF(X2<=0.)EXIT
                  Y2=(-LOG(X2))**WW
                  XY=(Y2+Y1)*DX
                  SUM=SUM+XY
                  X1=X2
                  Y1=Y2
            END DO
            DX=DX*.5
      END DO
      SUM=SUM*.5
END SUBROUTINE Integrate_ExpotentialEqs
 
END MODULE MisceFun_Module