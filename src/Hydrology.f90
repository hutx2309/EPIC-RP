MODULE Hydrology
USE PARM
USE MisceFun_Module
 
IMPLICIT NONE
CONTAINS
    
! 1. ------------------- Snow melting model -----------------------------------------------
                                
SUBROUTINE Snow_Melt                            
      !     EPIC1102
      !     THIS SUBPROGRAM PREDICTS DAILY SNOW MELT WHEN THE AVERAGE IrrWater TEMPERATURE EXCEEDS 0 DEGREES C.
      !     Ref : APEX-doc Eq 109 a-d    
      
      IMPLICIT NONE
 
      ! local variables       
      REAL:: X1, X2, F, SNPKT
 
      X2=MIN(surfTem,soilTem(Layer_ID(2)))
      X1=SQRT(TMX*SRAD)
      SNPKT=.3333*(2.*X2+TX)
      F= snowPackAge /(snowPackAge+EXP(S_Curve(16,1)-S_Curve(16,2)*snowPackAge))  ! Snow pack age in days
      snowMelting=MAX(0.,X1*(1.52+.54*F*SNPKT))            ! SNPKT: snow pack temperature
      snowMelting=MIN(snowMelting,snowWater)                  ! snowMelting: mm/d
      snowWater=snowWater-snowMelting                       ! snowWater: snow present in mm of water 
      IF(peakRainRate<1.E-5)peakRainRate=.042*snowMelting
      RainfallIntesity_5hr=.042
      Rainfall_Duration=24.
END SUBROUTINE Snow_Melt
                           
!  2. ---------------   ESTIMATES THE USLE RAINFALL ENERGY FACTOR, GIVEN DAILY RAINFALL.------
SUBROUTINE RainEnergy_Factor(IVR)
      ! Ref:  APEX-doc Eqs 116 -127 
      IMPLICIT NONE
 
      ! local variables
      REAL:: ALMN = 0.02083, AJP, X1 
      INTEGER:: IVR
    
      IF(IVR==0)THEN
            AJP=1.-EXP(-125./(Rainfall+5.))
            alpha05L=Triangle_Num(ALMN,halfHrRainInfo(NWI,MO),AJP,4)
            X1=-2.*LOG(1.-alpha05L)
            RainfallIntesity_5hr=2.*Rainfall*alpha05L+.001
            IF(peakRainRate<1.E-5)peakRainRate=X1*Rainfall+.001
            Rainfall_Duration=MIN(24.,4.605/X1)
      ELSE
            peakRainRate=Rainfall/24.     ! Rainfall: mm/d
            Rainfall_Duration=24.
            RainfallIntesity_5hr=peakRainRate
            X1=RainfallIntesity_5hr
      END IF    
      rainEnergyFactor=MAX(0.,Rainfall*(12.1+8.9*(LOG10(peakRainRate)-.4343))*RainfallIntesity_5hr/1000.)  ! APEX-doc  Eq. 119
      RainfallIntesity_5hr=X1                ! Peak rainfall rate (mm/h)  RainfallIntesity_5hr  maximum 0.5-h rainfall intensity
END SUBROUTINE RainEnergy_Factor
    
! 3. ---------- SOLVES THE GREEN & AMPT INFILTRATION EQ ASSUMING----------------------------
      !     F1 IS INCREMENTED BY TOTAL RAIN DURING DT
SUBROUTINE solve_infi(A, PT, RX)
      ! Ref: APEX-doc Eq. 40    !SATK: soil saturated conductivity in mm/h
      IMPLICIT NONE 
      ! local variables
      REAL:: F1, PT, X1, ZI, RX, A, Q1
    
      F1=PT-Runoff                      ! accumulated infitration in mm
      X1=satuateCondv0
      IF(soilTem(Layer_ID(2))<-1.)X1=.01*X1
      ZI=X1*(SCN/(F1+1.))               ! Eq. 40  
      IF(RX>ZI)THEN
            Q1=A*(RX-ZI)/RX
      ELSE
            Q1=0.
      END IF
      Runoff=Runoff+Q1
END SUBROUTINE solve_infi
    
! 4. ------------------ expotential rainfall ----------------------------------------
SUBROUTINE Runoff_Expo 
      ! Ref: APEX-doc Eqs 42-50
      IMPLICIT NONE
   
      ! local varialbes
      REAL:: DRFV = 2.0, PT, QMN, BLM, R1, RTP, XK1, XK2, X1, XKP1, XKP2, RX, A
 
      PT=0.
      maxRH=.95
      QMN=.25
      BLM=.05
      R1=Triangle_Num(BLM,QMN,maxRH,8)
      RTP=R1*Rainfall
      XK1=R1/4.605
      XK2=XK1*(1.-R1)/R1
      Rainfall_Duration=Rainfall/(peakRainRate*(XK1+XK2))
      X1=peakRainRate*Rainfall_Duration
      XKP1=XK1*X1
      XKP2=XK2*X1
      RX=0.
      PT=0.
      WRITE(KW(4),1)IYR,MO,DayOfMon,SCN,XK1,XK2,Rainfall,RTP,peakRainRate,Rainfall_Duration,R1
      DO WHILE(PT<RTP)
            PT=PT+DRFV
            RX=peakRainRate*(1.-(RTP-PT)/XKP1)
            CALL solve_infi(DRFV,PT,RX )
      END DO
      A=RTP-PT+DRFV
      PT=RTP
      RX=peakRainRate
      CALL solve_infi(A,PT,RX )
      DO WHILE(PT<Rainfall)
            PT=PT+DRFV
            RX=peakRainRate*(1.-(PT-RTP)/XKP2)
            CALL solve_infi(DRFV,PT,RX )
      END DO
      RX=0.
      A=Rainfall-PT
      PT=Rainfall
      CALL solve_infi(A,PT,RX )
 
    1 FORMAT(1X,3I4,9F10.2)
 
END SUBROUTINE Runoff_Expo 
    
! 5. ---------------  uniformal rainfall ----------------------------------------
SUBROUTINE Runoff_Unif 
      !     EPIC1102
      !     THIS SUBPROGRAM DISTRIBUTES DAILY RAINFALL UNIFORMLY &
      !     FURNISHES THE GREEN & AMPT SUBPROGRAM RAIN INCREMENTS OF EQUAL
      !     VOLUME = DRFV
      IMPLICIT NONE
    
      ! local variables
      REAL:: DRFV = 5.0, PT, A   
    
      PT=0.
      WRITE(KW(1),1)SCN,Rainfall,peakRainRate,Rainfall_Duration
      DO WHILE(PT<Rainfall)
            PT=PT+DRFV
            CALL solve_infi(DRFV, PT, peakRainRate )    ! returned variable
            WRITE(KW(1),1)PT,Runoff
      END DO
      A=Rainfall-PT+DRFV
      PT=Rainfall
      CALL solve_infi(A, PT, peakRainRate )
      WRITE(KW(1),1)PT,Runoff
      WRITE(KW(1),13)IYR,MO,DayOfMon,SCN,Rainfall,PT,Runoff,peakRainRate,Rainfall_Duration
       
    1 FORMAT(1X,10E13.5)
   13 FORMAT(1X,3I4,10F10.2)
END SUBROUTINE Runoff_Unif
    
! 6. ------------------ Furrow storage of water ------------------------------
SUBROUTINE Furrow_Store 
      ! Ref: APEX-doc 318 -319
      IMPLICIT NONE   
 
      ! local variables
      REAL:: DH, H2, X1, TW, BW , DI, D2, D3, TW2, TW3, A2, A3,XX, ZZ
    
      DH=dikeHeight*.001                 !  DH: dike height in m    DHT: dike height
      H2=2.*DH              
      X1 = 0.001 * DKHL
      TW=Ridge_IntervalX-X1   
      BW=MAX(TW-4.*X1,.1*TW)
      DI=dikeInterval-X1
      D2=DH*(1.-2.*uplandSteep)     ! Eq. 319 
      D3=DH-uplandSteep*(DI-H2)
      X1=(TW-BW)/DH
      TW2=BW+D2*X1
      TW3=BW+D3*X1
      A2=.5*D2*(TW2+BW)
      A3=.5*D3*(TW3+BW)
      XX=DH/uplandSteep
      ZZ=DI-H2
      X1=furrowDikeSafeFactor/(Ridge_IntervalX*dikeInterval)
      IF(XX>ZZ)THEN
            Dike_Volum=X1*(A2*DH+.5*(A2+A3)*(DI-4.*DH)+A3*D3)    ! the APEX model has divided by Ridge_IntervalX*dikeInterval
      ELSE
            Dike_Volum=X1*A2*(DH+.5*(XX-H2))
      END IF
 
END SUBROUTINE Furrow_Store
    
REAL FUNCTION HQP(X1,ITP,INT0)
      IMPLICIT NONE
  
      ! local variables
      REAL::X1
      INTEGER:: ITP,INT0
    
      HQP=coeffTR55(1,INT0,ITP)+X1*(coeffTR55(2,INT0,ITP)+X1*(coeffTR55(3,INT0,ITP)+X1*&
        &(coeffTR55(4,INT0,ITP)+X1*(coeffTR55(5,INT0,ITP)+X1*(coeffTR55(6,INT0,ITP)+X1*&
        &(coeffTR55(7,INT0,ITP)+X1*coeffTR55(8,INT0,ITP)))))))
    
END FUNCTION HQP
 
!7. --------------------HTR55 Method to cal runoff---------------------------------------
    
SUBROUTINE TR55
      IMPLICIT NONE
      ! questionable XOX and INT should be real or integer?
      ! local variable    
      REAL, DIMENSION(18):: PIAF = [0.0,0.1,0.2,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0]
      REAL:: RTO, X1, Y1, Y, Y2     
      INTEGER:: XOX, INT0, INT1
      
      IF(peakRunoffRate>.35)THEN
            XOX=INT((peakRunoffRate-.35)/.05+5.)
      ELSE
            XOX=INT(peakRunoffRate/.1+1.)
      END IF
      INT0=XOX
      INT1=INT0+1
      RTO=(peakRunoffRate-PIAF(INT0))/(PIAF(INT1)-PIAF(INT0))
      X1=LOG(TC)
      Y1=HQP(X1,peakRateMethod ,INT0)
      IF(INT0<17)THEN
            Y2=HQP(X1,peakRateMethod ,INT1)
            Y=Y1+(Y2-Y1)*RTO
            Y=EXP(Y)
      ELSE
            Y=EXP(Y1)*(1.0-RTO)
      END IF
      peakRunoffRate=Y*Rainfall
END SUBROUTINE TR55
 
!8 . ---------------------------Predict daily runoff and peak runoff------------------------------------------
SUBROUTINE Predict_Runoff(IVR)  
      !Ref: APEX-doc Eqs 29- 76  RUNOFF VOLUME
      IMPLICIT NONE
                               
      ! local variables
      INTEGER:: IVR,II, JJ, J, L, L1, I
      REAL:: ADD, SUM, XX, ZZ, RTO, BLM, BB, TOT,X1, X2, X4, ALTC  ! where is XHSM allocated?
      ! -------- if 1------------
      IF(IVR==0)THEN
            SUM=0.
            ADD=0.
            IF(landUseCode==35)THEN   ! Pavement and urban area
                  SCN=25400./CN0-254.    ! Eq.30
                  CN=CN0
            ELSE
                  II=CN_Method+1           ! methods  of calculate CN
                  SELECT CASE(II)
                        CASE(1)
                              XX=0.
                              DO JJ=1,Actual_SoilLayers
                                    J=Layer_ID(JJ)
                                    IF(Z(J)>1.)EXIT
                                    ZZ=(Z(J)-XX)/Z(J)
                                    SUM=SUM+(soilWater(J)-Wilt_Point(J))*ZZ/(fieldCapacity(J)-Wilt_Point(J))
                                    ADD=ADD+ZZ
                                    XX=Z(J)
                              END DO
                              IF(JJ<=Actual_SoilLayers)THEN
                                    ZZ=1.-XX
                                    SUM=SUM+(soilWater(J)-Wilt_Point(J))*ZZ/(fieldCapacity(J)-Wilt_Point(J))
                                    ADD=ADD+ZZ
                              END IF
                              SUM=SUM/ADD
                              IF(SUM>0.)THEN
                                    SUM=SUM*100.
                                    SCN=SMX*(1.-SUM/(SUM+EXP(S_Curve(4,1)-S_Curve(4,2)*SUM)))
                              ELSE
                                    SCN=SMX*(1.-SUM)**2
                              END IF
                              SCN=modelPara(81)*SCN
                        CASE(2)
                              DO JJ=1,Actual_SoilLayers
                                    L=Layer_ID(JJ)
                                    IF(Z(L)>1.)EXIT
                                    SUM=SUM+soilWater(L)-Wilt_Point(L)
                                    ADD=ADD+fieldCapacity(L)-Wilt_Point(L)
                                    L1=L
                              END DO
                              IF(JJ<=Actual_SoilLayers)THEN
                                    RTO=(1.-Z(L1))/(Z(L)-Z(L1))
                                    SUM=SUM+(soilWater(L)-Wilt_Point(L))*RTO
                                    ADD=ADD+(fieldCapacity(L)-Wilt_Point(L))*RTO
                              END IF
                              SUM=SUM/ADD
                              IF(SUM>0.)THEN
                                    SUM=SUM*100.
                                    SCN=SMX*(1.-SUM/(SUM+EXP(S_Curve(4,1)-S_Curve(4,2)*SUM)))
                              ELSE
                                    SCN=SMX*(1.-SUM)**2
                              END IF
                        CASE(3)
                              DO JJ=1,Actual_SoilLayers
                                    ISL=Layer_ID(JJ)
                                    SUM=SUM+soilWater(ISL)-Wilt_Point(ISL)
                                    ADD=ADD+fieldCapacity(ISL)-Wilt_Point(ISL)
                                    IF(Z(ISL)>1.)EXIT
                              END DO
                              SUM=SUM/ADD
                              IF(SUM<0.)THEN
                                    SCN=SMX*(1.-SUM)**2
                              ELSE
                                    RTO=MIN(.98,SUM)
                                    SCN=SMX*(1.-RTO)
                              END IF
                        CASE(4)
                              SCN=25400./CN0-254.
                              CN=CN0
                        CASE(5)
                              SCN=MAX(3.,SCI)
                  END SELECT
                  IF(II/=4)THEN
                        SCN=(SCN-SMX)*EXP(modelPara(75)*(1.-cropResidu(LD1)))+SMX  
                        IF(soilTem(Layer_ID(2))<-1.)SCN=SCN*modelPara(65)
                        IF(ICV>0)SCN=SCN*SQRT(1.-FCV) 
                        CN=25400./(SCN+254.)
                        IF(curveNumFlag==0)THEN
                              maxRH=MIN(99.5,CN+5.)
                              BLM=MAX(1.,CN-5.)
                              CN=Triangle_Num(BLM,CN,maxRH,8)
                        END IF
                        SCN=25400./CN-254.
                        IF(SCN<3.)THEN
                              SCN=3. 
                              CN=25400./(SCN+254.)
                        END IF
                  END IF
            END IF    
            BB=.2*SCN
            TOT=100.
            DO I=1,9
                  TOT=TOT-5.
                  IF(CN>TOT)EXIT
            END DO
            initialNO3(I)=initialNO3(I)+1.
            RTO=MIN(1.,SCN/SMX)
            crackFlow= modelPara(17)*Rainfall*RTO
            Rainfall=Rainfall-crackFlow
            IF(totRunoff>0.)THEN
                  SELECT CASE(runoff_Method) ! runoff estimation methods 
                        CASE(1)
                              X1=totRunoff-BB
                              IF(X1>0.)THEN
                                    Runoff=X1*X1/(totRunoff+.8*SCN)
                              ELSE
                                    Runoff=0.
                              END IF
                        CASE(2,3)
                              CALL Runoff_Expo 
                        CASE(4)
                              CALL Runoff_Unif 
                  END SELECT
            END IF    
      ELSE
            Runoff=IrrRunoffRatio*totRunoff    
      END IF ! ---------- end of if 1 -----------
     ! -------------if 2 --------------
      IF(furrowFlag>0.AND.Runoff>0.)THEN
            IF(dikeHeight>.01.AND.Ridge_IntervalX>.01)THEN
                  CALL Furrow_Store 
                  X1=MAX(0.,soilWater(LD1)-Porosity(LD1))
                  IF(NOP>0) WRITE(KW(1),17)IYR,MO,DayOfMon,dikeHeight,RidgeHeight2,Runoff,Dike_Volum,X1,XHSM
                  IF(Runoff<=Dike_Volum-X1)THEN
                        Runoff=0.
                  ELSE
                        dikeHeight=0.
                        IF(ABS(IDRL) < 1.0E-8 .AND. CPHT(Crop_Num) < 1.)dikeHeight=DKHL    
                  END IF
            END IF
      END IF  ! ------------- end of if 2 -----------------
      ! --------------- if 3 Peak runoff rate-------------
      IF(peakRateMethod ==0)THEN   ! /////// Modified rational eqs method ////////
            X2=Runoff/Rainfall_Duration       ! qcl = X2: average flow rate in m3/s 
            IF(X2>1.)THEN
                  X2=X2**.25
            ELSE
                  X2=1.
            END IF   
            X4=MIN(uplandSlopLen/360.,TCS/X2)
            TC=X4+TCC/X2                      ! TC: time of concentration in h   
            ALTC=1.-EXP(-TC*RainfallIntesity_5hr)
            X1=ALTC*Runoff
            peakRunoffRate=X1/TC
            Flow_Rate=X1/X4
      ELSE            !///////// TR55 method ///////////
            TC=TCC+TCS/SQRT(totRunoff)                   ! totRunoff: tot runoff volume 
            peakRunoffRate=MIN(.95,BB/totRunoff)
            CALL TR55 
      END IF ! ---------- end of if 3 ----------
      
      IF(KFL(4)>0)WRITE(KW(4),27)IYR,MO,DayOfMon,CN,totRunoff,Runoff,TC,peakRunoffRate,Rainfall_Duration,ALTC,alpha05L
   27 FORMAT(1X,3I4,9F10.2)
   17 FORMAT(1X,3I4,12X,'dikeHeight=',F5.0,'MM',2X,'RidgeHeight2=',F5.0,'MM',2X,'Q=',F5.1,'MM',2X,&
            'Dike_Volum=',F5.1,'MM',2X,'soilWater=',F5.1,'MM',2X,'PotentHU_frac=',F6.2)
END SUBROUTINE Predict_Runoff
  
! 9. - ---------------------Water Table Dynamics ----------------------
SUBROUTINE WaterTable_Dyn 
      IMPLICIT NONE
      !     EPIC1102
      !     THIS SUBPROGRAM SIMULATES WATER TABLE DYNAMICS AS A FUNCTION OF RAIN AND EVAP.
 
      ! local variables
      INTEGER:: NN, K
      REAL:: RTO, SUM, TOT, XX
 
      RTO=100.*groundWaterStor/groundWaterMaxStor
      waterTableHigt=waterTableMaxDep+(waterTableMinDep-waterTableMaxDep)*RTO/(RTO+EXP(S_Curve(19,1)-S_Curve(19,2)*RTO))
      SUM=0.
      TOT=0.
      IF(waterTableHigt<=Z(Layer_ID(Actual_SoilLayers)))THEN
            XX=0.
            NN=0
            DO K=1,Actual_SoilLayers
                  ISL=Layer_ID(K)
                  SUM=SUM+soilWater(ISL)
                  IF(waterTableHigt<=Z(ISL))THEN
                        IF(NN>0)THEN
                              soilWater(ISL)=Porosity(ISL)
                        ELSE
                              NN=1
                              IF(soilWater(ISL)>fieldCapacity(ISL))soilWater(ISL)=fieldCapacity(ISL)
                              soilWater(ISL)=(soilWater(ISL)*(waterTableHigt-XX)+Porosity(ISL)*(Z(ISL)-waterTableHigt))/&
                              &(Z(ISL)-XX)
                        END IF
                  END IF
                  TOT=TOT+soilWater(ISL)
                  XX=Z(ISL)
            END DO
      END IF
      XX=TOT-SUM
      Inflow(MO)=Inflow(MO)+XX      ! inflow to root zone from water table
      SMM(20,MO)=SMM(20,MO)+XX      ! SMM(20,MO), VAR(20) :  inflow
      VAR(20)=XX
END SUBROUTINE WaterTable_Dyn
              
!10.  ----------------------LAGOON --------------------------
SUBROUTINE Lagoon( JRT)   
      IMPLICIT NONE
      !     EPIC1102
      !     THIS SUBPROGRAM ESTIMATES INFLOW, STORAGE, & IRRIGATION FROM
      !     LAGOONS. RUNOFF IS ESTIMATED WITH 90 CN.
      ! local varialbes
      INTEGER:: JRT
      REAL:: X1, X2, DP, TW, SA, EV, XX, Runoff_toLagoon
 
      X2=10.*DALG
      DP=.1677*lagoonVol**.3333
      TW=18.*DP
      SA=.0001*TW*TW
      EV=6.*EO*SA
      lagoonVol=lagoonVol-EV+lagoonInput
      EV=EV/X2
      SMM(21,MO)=SMM(21,MO)+EV                    ! SMM(21,MO), VAR(21) : lagoon evperation
      VAR(21)=EV
      XX=lagoonInput/X2
      SMM(22,MO)=SMM(22,MO)+XX                    ! SMM(22,MO), VAR(22) : water wash to lagoon 
      VAR(22)=XX
      X1=totRunoff-5.64
      IF(X1>0.)THEN
            Runoff_toLagoon=X1*X1/(totRunoff+22.6)
            Runoff_toLagoon=10.*(Runoff_toLagoon*(DALG-SA)+totRunoff*SA)
            lagoonVol=lagoonVol+Runoff_toLagoon
            Runoff_toLagoon=Runoff_toLagoon/X2
            SMM(23,MO)=SMM(23,MO)+Runoff_toLagoon     ! SMM(23,MO), VAR(23) : runoff to lagoon              (mm)
            VAR(23)=Runoff_toLagoon
      END IF
      IF(lagoonVol<= lagoonVol0)THEN
            JRT=1                                     ! JRT ==1: water from lagoon was used for irrigation
            RETURN
      END IF
      IF(lagoonVol>lagoonMaxVol)THEN
            X1=(lagoonVol-lagoonMaxVol)/X2
            WRITE(KW(1),4)IYR,MO,DayOfMon,X1
            SMM(24,MO)=SMM(24,MO)+X1                    ! SMM(24,MO), VAR(24) : lagoon overflow             (mm)
            VAR(24)=X1
            lagoonVol=lagoonMaxVol
      END IF
      IF(soilWaterRZ>=potentialAvailWater)THEN
            JRT=1
            RETURN
      END IF
      XX=10.*areaWshed*irrFromLagoon
      X1=XX/X2
      SMM(25,MO)=SMM(25,MO)+X1                        ! SMM(25,MO), VAR(25) : Irrigation from a lagoon    (mm)
      VAR(25)=X1
      lagoonVol=MAX(1.E-5,lagoonVol-XX)
      JRT=0
 
    4 FORMAT(T10,'***** LAGOON OVERFLOWED ',3I4,F4.0,' MM')
END SUBROUTINE Lagoon

! ------------------------------ 10.2 Check Lagoon water balance --------------
SUBROUTINE Check_LagoonWater(Q,EV,O,RG,VLGE,VLGB,WW)
      !     EPIC1102
      !     THIS SUBPROGRAM CHECKS THE LAGOON WATER BALANCE AT THE END OF A
      !     SIMULATION.
      IMPLICIT NONE
      ! local variables
      REAL:: DF, Q,EV,O,RG,VLGE,VLGB,WW, PER
      
      WRITE(KW(1),'(T10,A)')'LAGOON WATER BALANCE'
      DF=VLGB+Q-EV-O-VLGE-RG+WW
      PER=200.*DF/(VLGB+VLGE)
      WRITE(KW(1),3)DF,VLGB,Q,EV,O,VLGE,RG,WW
      VLGB=VLGE
  
    3 FORMAT(8E16.6)
END SUBROUTINE Check_LagoonWater
	
! ----------------------------- 10.3 Lagoon manure balance ---------------------
SUBROUTINE Check_LagoonMaure(WTI,WTO,WTB,WTE)
      !     EPIC1102
      !     THIS SUBPROGRAM THE LAGOON MANURE BALANCE
      IMPLICIT NONE
      ! local variables
      REAL:: DF, WTI,WTO,WTB,WTE, PER
 
      DF=WTB+WTI-WTO-WTE
      PER=200.*DF/(WTB+WTE)
      WRITE(KW(1),2)
      WRITE(KW(1),1)PER,DF,WTB,WTI,WTO,WTE     
 
    1 FORMAT(5X,'PER =',E13.6,2X,'DF  =',E13.6,2X,'WTB =',E13.6,2X,&
     &'WTI =',E13.6,2X,'WTO =',E13.6,2X,'WTE =',E13.6)
    2 FORMAT(T10,'LAGOON MANURE BALANCE')
 
END SUBROUTINE Check_LagoonMaure

! ----------------------------- EP --------------------------------------------
REAL FUNCTION  Sat_Vapor(TK)
      IMPLICIT NONE
      !     THIS SUBPROGRAM COMPUTES THE SATURATION VAPOR PRESSURE GIVEN 
      !     TEMPERATURE.
      REAL:: TK
      Sat_Vapor=.1*EXP(54.879-5.029*LOG(TK)-6790.5/TK)  ! Eq. 93i in APEX-doc
      
END FUNCTION Sat_Vapor
 
SUBROUTINE Cal_EVP 
      IMPLICIT NONE
      !     EPIC1102
      !     THIS SUBPROGRAM ESTIMATES DAILY EVAPOTRANSPIRATION.  THERE ARE FOUR OPTIONS FOR COMPUTING POTENTIAL EVAP
      !     (PENMAN-MONTEITH, PENMAN, PRIESTLEY-TAYLOR, & HARGREAVES)
      ! APEX-doc Pg.28    
      ! (ICV may cause problems, FCV is factor of vegetative cover? or machinary efficiency?)
   
      ! local variables
      INTEGER:: NN,K, K1, J
      REAL:: SUM, potenPlantEP,  CHMX, X1, F, EAJ, TK, RSO, EA, ED, RALB1, RAMM, DLT, XX, XL,TK4, RBO, RTO, &     
           RN, X2, FWV, X3, Air_Density, UZZ, ZZ, Z0, ZD, RV, FVPD, G1,RC,XZ, XY, Z1, TOT
      SUM=.01
      potenPlantEP=0.
      NN=NBC(IRO)      ! NN: number of crops?
      CHMX=0.
      DO K=1,NN
            K1=JE(K)
            IF(K1>MNC)CYCLE
            SUM=SUM+Current_LAI(K1)
            IF(CPHT(K1)>CHMX)CHMX=CPHT(K1)
            EP(K1)=0.      ! Potential evaporation for each crop      mm/d
      END DO
      X1=MAX(.4*SUM,modelPara(41)*(groundCover+.1))  ! groundCover: all above ground plant material  (t/ha)
      EAJ=EXP(-X1)                                     ! Eq. 100: EAJ: soil cover index
      IF(snowWater>5.)THEN   ! if snow is 5 mm or more, set albedo 0.6 and EAJ 0.5 for estimate EO and snow is evaporated at that rate.
            Albedo=.6
            EAJ=.5
      ELSE
            IF(ICV>0)THEN  ! ICV >0 : no crop  
                  X1=1.-EAJ
                  Albedo=(.23*X1+PALB*FCV)/(FCV+X1)
                  EAJ=MIN(EAJ,1.-FCV) 
            ELSE           ! ICV =0:  crop are growing
                  Albedo=.23*(1.-EAJ)+soilAlbedo*EAJ
            END IF
      END IF          
      TK=TX+273.
      RSO=maxRad
      XL=2.501-2.2E-3*TX   ! XL == HV: Latent heat of vaporization in MJ/kg
      EA=.1*EXP(54.879-5.029*LOG(TK)-6790.5/TK) ! Saturation vapor pressure (kPa)
      ED=EA*RHD                                 ! Vapor pressure
      VPD=EA-ED            
      SMM(9,MO)=SMM(9,MO)+VPD                   ! SMM(9,MO), VAR(9): VPD
      VAR(9)=VPD
      RALB1=SRAD*(1.-Albedo)                    ! Eq. 93a
      DLT=EA*(6790.5/TK-5.029)/TK               ! Slope of the saturation vapor pressure curve in kPa/C
      XX=DLT+psychrometerConst                                ! GMA:  psychrometer constant in kPa/C    
      TK4=TK**4
      RBO=(.34-.14*SQRT(ED))*4.9E-9*TK4         ! RBO:  net outgoing long wave radiation in MJ/m2/d
      RTO=MIN(.99,SRAD/(RSO+.1))         
      RN=RALB1-RBO*(.9*RTO+.1)                  ! RN: net radiation in MJ/m2/d  
      X2=RN*DLT                          
      SELECT CASE(ET_Method)
            CASE(5)
                  ! BAIER-ROBERTSON PET METHOD
                  EO=MAX(0.,.288*TMX-.144*TMN+.139*RSO-4.931)
            CASE(4)
                  ! HARGREAVES PET METHOD
                  RAMM=RSO/XL
                  EO=MAX(0.,modelPara(38)*RAMM*(TX+17.8)*(TMX-TMN)**HGX)    ! Potential Evaportation mm/d  
            CASE(3)
                  ! PRIESTLEY-TAYLOR PET METHOD
                  RAMM=RALB1/XL
                  EO=1.28*RAMM*DLT/XX
            CASE(2)
                  ! PENMAN PET METHOD
                  FWV=2.7+1.63*U10
                  X3=psychrometerConst*FWV*VPD
                  X1=X2/XL+X3
                  EO=X1/XX
            CASE(1)
                  ! PENMAN-MONTEITH PET METHOD
                  Air_Density=.01276*barometricPressure/(1.+.00367*TX)
                  IF(IGO==0)GO TO 12    ! growing crop?
                  IF(CHMX<8.)THEN
                        UZZ=U10           ! daily mean wind speed
                        ZZ=10.
                  ELSE
                        ZZ=CHMX+2.
                        UZZ=U10*LOG(ZZ/.0005)/9.9035
                  END IF
                  X1=LOG10(CHMX+.01)
                  Z0=10.**(.997*X1-.883)     ! surface roughtness
                  ZD=10.**(.979*X1-.154)     ! displacement height of the crop  
                  RV=6.25*(LOG((ZZ-ZD)/Z0))**2/UZZ  ! aerodynamic resistance for heat and vapor transfer  s/m
                  X3=VPD-VPDThreshold(Crop_Num)
                  IF(X3>0.)THEN
                        FVPD=MAX(.1,1.-VPDPara2(Crop_Num)*X3)
                  ELSE
                        FVPD=1.
                  END IF
                  G1=maxStomaCond(Crop_Num)*FVPD ! maxStomaCond: crops leaf resistance  s/m
                  RC=modelPara(1)/((SUM+.01)*G1*EXP(.00155*(330.-CO2)))
                  potenPlantEP=modelPara(74)*(X2+86.66*Air_Density*VPD/RV)/(XL*(DLT+psychrometerConst*(1.+RC/RV)))  ! Potential plant evaporation in mm/d, this is what we finally want 
              12  RV=350./U10
                  EO=modelPara(74)*(X2+86.66*Air_Density*VPD/RV)/(XL*XX)
                  IF(potenPlantEP>EO)EO=potenPlantEP
            CASE DEFAULT              
                  ! HARGREAVES PET METHOD
                  RAMM=RSO/XL
                  EO=MAX(0.,modelPara(38)*RAMM*(TX+17.8)*(TMX-TMN)**HGX)
                  EO=MIN(9.,EO)    
      END SELECT
 
      IF(ET_Method>1) potenPlantEP=MIN(EO, EO*SUM/3.)
      
      IF(IGO>0)THEN      ! growing crops
            XX=potenPlantEP/SUM
            DO K=1,NN
                  K1=JE(K)
                  IF(K1>MNC)CYCLE
                  EP(K1)=Current_LAI(K1)*XX
            END DO
      END IF
      
      soilEvp=EO*EAJ                       ! soilEvp: potential soil evaporation in APEX-doc P31
      soilRad=RALB1
      soilEvp=MIN(soilEvp,soilEvp*EO/(soilEvp+potenPlantEP+1.E-10))  !Eq. 105 in APEX-doc  P31
      
      IF(snowWater>=soilEvp)GO TO 21
      XX=soilEvp-snowWater
      soilEvp=snowWater
      snowWater=0.
      TOT=0.
      DO J=1,Actual_SoilLayers
            ISL=Layer_ID(J)
            RTO=1000.*Z(ISL)     ! m  to mm
            SUM=XX*RTO/(RTO+EXP(S_Curve(2,1)-S_Curve(2,2)*RTO)) !S_curve(2) : governing soil evaporation as a function of soil depth 
            XZ=fieldCapacity(ISL)-Wilt_Point(ISL)
            IF(soilWater(ISL)<fieldCapacity(ISL))THEN
                  F=EXP(modelPara(12)*(soilWater(ISL)-fieldCapacity(ISL))/XZ)
            ELSE
                  F=1.
            END IF
            ZZ=SUM-modelPara(61)*TOT-(1.-modelPara(61))*soilEvp  ! para(61) : weighting factor 
            soilElement(ISL)=ZZ*F                                ! soilElement is used to store soil evaporation at each layer since here
            XY=modelPara(5)*Wilt_Point(ISL)                      ! para(5): soil water lower limit of water content
            IF(Z(ISL)>.2)GO TO 20
            IF(soilWater(ISL)-soilElement(ISL)<XY)soilElement(ISL)=soilWater(ISL)-XY-1.E-5
            soilEvp=soilEvp+soilElement(ISL)
            soilWater(ISL)=soilWater(ISL)-soilElement(ISL)
            TOT=SUM
      END DO
      J=Actual_SoilLayers
   20 NEV=J
      Z1=Z(Layer_ID(J-1))
      RTO=(.2-Z1)/(Z(ISL)-Z1)
      X1=RTO*soilWater(ISL)
      X2=RTO*XY
      IF(X1-soilElement(ISL)<X2)soilElement(ISL)=X1-X2
      soilEvp=soilEvp+soilElement(ISL)
      soilWater(ISL)=MAX(1.E-5,soilWater(ISL)-soilElement(ISL))
      GO TO 23
   21 snowWater=snowWater-soilEvp
      NEV=1
   23 XX=MAX(0.,EO-soilEvp)
      IF(potenPlantEP>XX)THEN
            X1=XX/potenPlantEP
            potenPlantEP=XX
            DO K=1,NN
                  K1=JE(K)
                  IF(K1>MNC)CYCLE
                  EP(K1)=EP(K1)*X1
            END DO
      END IF
      VAR(12)=potenPlantEP                   ! SMM(12, MO), VAR(12) : Potential plant evaporation in mm/d
      SMM(12,MO)=SMM(12,MO)+potenPlantEP 
 
END SUBROUTINE Cal_EVP
 
! ----------------------------- check soil water balance ------------------------------------
SUBROUTINE Print_WaterVar(P,Q,ET,SubSurfFlow_lat,O,RG,SWW,QINT,RGDL) 
      ! EPIC1102
      ! THIS SUBPROGRAM CHECKS THE SOIL WATER BALANCE AT THE END OF A SIMULATION.
      IMPLICIT NONE
      ! arguments list
      REAL:: SubSurfFlow_lat
      ! local variables
      REAL:: SWW, P, Q, ET, O,  RG, QINT, RGDL, DF, PER 
      
      WRITE(KW(1),'(T10,A)')'SOIL WATER BALANCE'
      DF=SWW+P-Q-ET-O-totSoilWater-SubSurfFlow_lat+RG-snowWater+QINT-RGDL
      PER=100.*DF/totSoilWater
      WRITE(KW(1),1)PER,DF,SWW,P,Q,ET,O,SubSurfFlow_lat,RG,RGDL,snowWater,QINT,totSoilWater
      SWW=totSoilWater+snowWater
         
    1 FORMAT(5X,'PER =',E13.6,2X,'DF  =',E13.6,2X,'BSW =',E13.6,2X,&
            'PCP =',E13.6,2X,'Q   =',E13.6,2X,'ET  =',E13.6/5X,'percolateFlow=',E13.6,&
            2X,'subsurfLaterFlow =',E13.6,2X,'IRG =',E13.6,2X,'IRDL=',E13.6,2X,'snowWater =',&
            E13.6,2X,'Inflow =',E13.6,2X/5X,'FSW =',E13.6)
END SUBROUTINE Print_WaterVar
 
END MODULE Hydrology