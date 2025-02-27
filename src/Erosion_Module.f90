MODULE Erosion_Module
USE PARM
USE MisceFun_Module
USE Soil_Module
IMPLICIT NONE
   CONTAINS 
    
! ************************************I. water erosion *****************************************   
    !-------------------1 RUSLE C factor --------------------
REAL FUNCTION RUSLE_C() 
      IMPLICIT NONE
      !EPCM1102.FOR
      !     THIS SUBPROGRAM ESTIMATES THE RUSLE C FACTOR DAILY.
      ! Ref: APEX-doc  114a-d
      ! questionable: should NN be integer or real?
      ! local variables
      INTEGER:: I, L, K, KK, NN
      REAL:: SUM, TOT, RTO, PLU, FRUF, FRSD , X1, FBIO, FPPL
  
      NN=NBC(IRO)
      SUM=0.
      TOT=0.
      DO I=2,Actual_SoilLayers
            L=Layer_ID(I)
            IF(Z(L)>plawDepth)GO TO 3
            IF(I>1)TOT=TOT+cropResidu(L)
            DO K=1,NN
                  IF(JE(K)>MNC)CYCLE
                  SUM=SUM+RWT(L,JE(K))
            END DO
      END DO
      GO TO 4
    3 KK=Layer_ID(I-1)
      RTO=(plawDepth-Z(KK))/(Z(L)-Z(KK))
      TOT=TOT+RTO*cropResidu(L)
      DO K=1,NN
            IF(JE(K)>MNC)CYCLE
            SUM=SUM+RTO*RWT(L,JE(K))
      END DO
    4 SUM=SUM/plawDepth
      TOT=TOT/(plawDepth-Z(LD1))
      PLU=.951*RCF*EXP(-.0451*SUM-.00943*TOT/SQRT(RCF))            
      FRUF=MIN(1.,EXP(-.026*(RandRough-6.1)))                 !Eq. 114c
      IF(abvGroundResidue<15.)THEN
            FRSD=EXP(-modelPara(23)*abvGroundResidue)   ! C residule effect in RUSLE eqs
      ELSE
            FRSD=.0001
      END IF
      X1=MAX(GroundCover_Frac,GroundCover_StandLiveBio)    
      FBIO=1.-X1*EXP(-modelPara(26)*CPHT(Crop_Num))          ! ground biomass factor
      FPPL=.9*(1.-coverFactorPlantPopu)+.1
      RUSLE_C=MAX(1.E-10,FRSD*FBIO*FRUF*FPPL)                  ! CropManage_Factor: crop management factor 
END FUNCTION RUSLE_C
                    
! 2. ---------------------------USLE C factor -------------------------                    
REAL FUNCTION USLE_C ()
      !     THIS SUBPROGRAM ESTIMATES THE USLE C FACTOR BASED ON PLANT POP & BIOMASS & RESIDUE COVER
      IMPLICIT NONE
  
      ! local variables
      REAL:: X1
    
      X1=MIN(10.,modelPara(60)*groundCover)
      USLE_C=(.8*EXP(-X1)+minCFactor)*(.9*(1.-coverFactorPlantPopu)+.1)
      USLE_C=USLE_C*EXP(-.05*rockFrac(LD1))
     
END FUNCTION USLE_C
                    
! 3. ---------------------- soil loss by water erosion -----------                     
Subroutine Water_Erosion 
      IMPLICIT NONE
      !     EPIC1102
      !     THIS SUBPROGRAM PREDICTS DAILY SOIL LOSS CAUSED BY WATER EROSION
      !     AND ESTIMATES THE NUTRIENT ENRICHMENT RATIO.
 
      ! local  variables
      INTEGER:: I
      REAL, DIMENSION(5):: PSZ = [200.,10.,2.,30.,500.]
      REAL:: CVX, F, XX, YLM, X1, B,B1, B2, RXM, RLFX, YI, T2, QQ, CY, DR1,SUM, &
      YSW, Wshed_Area2 
 
      CropManage_Factor = RUSLE_C( )  
    
      coverFactor = USLE_C( )

      CVX=coverFactor
      IF(erosionC_factor==0)CVX=CropManage_Factor
      F=1.
      XX=soilTem(Layer_ID(2))     
      IF(XX<=0.)THEN                  ! The effects of frozen soils. XX is temperature of the second layer.
            XX=273.+XX                  ! turn temperature to K       
            F=XX/(XX+EXP(S_Curve(18,1)-S_Curve(18,2)*XX))  ! Express soil temperature effect on erosion of frozen soils
      END IF
      XX=CVX*USL*F
      YLM=Runoff
      IF(Rainfall>12.7)THEN
            Sediment(3)=MIN(YLM,rainEnergyFactor*XX*1.292)
            SMM(29,MO)=SMM(29,MO)+coverFactor*rainEnergyFactor ! SMM(29, MO): 
            SMM(28,MO)=SMM(28,MO)+rainEnergyFactor              ! SMM(28, MO): rainfall energy factor
            X1=CropManage_Factor*rainEnergyFactor
            SMM(37,MO)=SMM(37,MO)+X1                               ! SMM(37, MO): 
            Sediment(7)=MIN(YLM,X1*RSLK)
            IF(Sediment(7)<1.E-10)GO TO 2
            B=BETA
            RXM=B/(1.+B)                                           ! APEX-doc Eq. 113f
            RLFX=slopeLength**RXM
            SMM(90,MO)=SMM(90,MO)+RLFX*RSF*rainEnergyFactor      ! SMM(90, MO) : SL*EI
            SMM(91,MO)=SMM(91,MO)+REK*rainEnergyFactor
            YI=MIN(YLM,.5*X1*RSK)
 
            DO I=1,3                                                                            
                  PSZ(I)=.411*PSZ(I)*PSZ(I)     ! PSZ: particle size                                           
            END DO       
            SUM=PSZ(1)*sandFrac(LD1)     
            SUM=SUM+PSZ(2)*siltFrac(LD1)
            SUM=SUM+PSZ(3)*clayFrac(LD1)
            SUM=.01*modelPara(71)*SUM/(Flow_Rate+1.E-5)
            T2=modelPara(72)*Flow_Rate*uplandSteep
            IF(T2>YI)THEN
                  Sediment(8)=1.5*X1*RLFX*RSK
                  RUSM(1)=RUSM(1)+Sediment(8)
            ELSE
                  Sediment(8)=MAX(0.,YI+SUM*(T2-YI))
                  RUSM(2)=RUSM(2)+Sediment(8)
                  RUSM(3)=RUSM(3)+YI
            END IF
            Sediment(8)=MIN(YLM,Sediment(8))
      END IF
    2 IF(Runoff<1.)RETURN
      YSW=.79*WSX**.009               ! Eq. 111e in APEX doc without *Q*q      MUSS, in EYSED: QQ = Q*q                                     
      Wshed_Area2=1.586*WSX**.12      ! Eq. 111c without *Q*q:             MUSLE, in EYSED: QQ = Q*q                  
      peakRainRate=peakRainRate*(Runoff/totRunoff)**.1     ! totRunoff : rainfall, snowmetling, runoff
      !Sediment(2): Sediment by Onstad-Foster model 
      Sediment(2)=MIN(YLM,(.646*rainEnergyFactor+.45*Runoff*peakRunoffRate**.3333)*XX) 
      QQ=Runoff*peakRunoffRate                         ! APEX-doc Eq. 111e 
      ! Sediment(4): Sediment by MUSS model
      Sediment(4)=MIN(YLM,YSW*QQ**.65*XX)   
      ! Sediment(1): Sediment by MUST model
      Sediment(1)=MIN(YLM,2.5*SQRT(QQ)*XX)        
      CX(MO)=CX(MO)+1.
      totHfHourRain(MO)=totHfHourRain(MO)+alpha05L
      ! Sediment(5): Sediment by MUSLE model
      Sediment(5)=MIN(YLM,Wshed_Area2*QQ**.56*XX)
      DR=SQRT(peakRunoffRate/peakRainRate)
      CY=.1*Sediment(4)/Runoff+1.E-10
      IF(Enrich_Method>0)THEN
            enrichRatio=.78*CY**(-.2468)
      ELSE
            DR1=1./DR
            B2=-LOG10(DR1)/2.699
            B1=1./.1**B2
            enrichRatio=MAX(1.,B1*(CY+1.E-4)**B2)
            enrichRatio=MIN(enrichRatio,3.)
      END IF
      SMM(58,MO)=SMM(58,MO)+enrichRatio                       ! SMM(58,MO): Enrichment ratio
   
END SUBROUTINE Water_Erosion
    
! *************************** II . Wind Erosion *************************************
   !  2.1 -----------------Soil erodibility by wind-------------------------
REAL FUNCTION Wind_Erodibility(sand,silt, clay, CaCo3f)
      !     THIS SUBPROGRAM ESTIMATES THE SOIL ERODIBILITY FACTOR FOR THE WIND
      !     EROSION EQ.
      ! questionable: Caco3, will wind_erodibility can not be referenced?
      IMPLICIT NONE
      REAL:: sand,silt, clay, CaCo3f
 
      IF(sand>85.+.5*clay)THEN
            Wind_Erodibility=1.
            RETURN
      END IF
      IF(sand>70.+clay)THEN
            Wind_Erodibility=.43
            RETURN
      END IF
      IF(silt>80..AND.clay<12.)THEN
            Wind_Erodibility=.12
            RETURN
      END IF
      IF(CaCo3f >0.)THEN
            IF(sand <45..OR.clay <20..OR.silt >28.)THEN
                  Wind_Erodibility=.28
                  RETURN
            ELSE
                  Wind_Erodibility=.18
                  RETURN
            END IF
      END IF  
      IF(clay <7.)THEN
            IF(silt <50.)THEN
                  Wind_Erodibility=.28
                  RETURN
            ELSE
                  Wind_Erodibility=.18
                  RETURN
            END IF
      END IF        
      IF(clay <20.)THEN
            IF(sand >52.)THEN
                  Wind_Erodibility=.28
                  RETURN
            ELSE
                  Wind_Erodibility=.18
                  RETURN
            END IF
      END IF        
      IF(clay <27.)THEN
            IF(silt <28.)THEN
                  Wind_Erodibility=.18
                  RETURN
            ELSE
                  Wind_Erodibility=.16
                  RETURN
            END IF
      END IF        
      IF(clay <35..AND.sand <20.)THEN
            Wind_Erodibility=.12
            RETURN
      END IF
      IF(clay <35.)THEN
            IF(sand <45.)THEN        
                  Wind_Erodibility=.16
                  RETURN
            ELSE
                  Wind_Erodibility=.18
                  RETURN
            END IF
      END IF        
      IF(sand >45.)THEN
            Wind_Erodibility=.18
            RETURN
      ELSE
            Wind_Erodibility=.28
            RETURN
      END IF
 
END FUNCTION Wind_Erodibility
    
! ------------- 2.2 Potential wind erosion rate ----------------------
    
REAL FUNCTION PWindEro_Rate(Y1, ALG, USTW)
      ! Ref: APEX-doc Eq 137  
      IMPLICIT NONE 
      ! local variable
      REAL:: ALG, DU10, USTR, X1, YWR, Y1, USTW    
      DU10=U10*Y1
      IF(DU10>U10MX(MO))U10MX(MO)=DU10
      USTR=.0408*DU10
      X1=USTR*USTR-USTW
      IF(X1<0.)THEN
            PWindEro_Rate=0.
      ELSE
            YWR=.255*X1**1.5          ! erosion rate  (kg/m/s) 
            PWindEro_Rate=YWR*ALG     ! APEX does not multiply ALG
      END IF
 
END FUNCTION PWindEro_Rate
 
! -------------- 2.3 Potential Daily Wind Erosion --------------------
    
SUBROUTINE PWind_Erosion(ALG, USTW, USTRT )
      IMPLICIT NONE   
      ! local variables
      REAL:: ALG, USTRT, USTW, WW, SUM, DU10,Y1, XX, X1, DX, Z1, XY, X2, Y2, Z2, XZ  
  
      WW=MIN(.5,Triangle_Num(.1,.35,.6,9)*UAVM(MO)/U10)
      SUM=.6424*WW**(-.1508)*EXP(.4336*WW)
      DU10=USTRT/.0408
      Y1=DU10/U10
      XX=LOG10(Y1)/WW
      IF(XX>1.3)THEN
            windErodSoil=0.
            RETURN
      END IF
      IF(XX<-3.)GO TO 1
      XX=10**XX
      GO TO 6    
    1 X1=1.
      GO TO 4 
    6 X1=EXP(-XX)
    4 DX=.1
      windErodSoil=0.
      Z1=0.
      DO WHILE(DX>1.E-4)                          ! integration of wind erosion rate
            XY=0.
            DO WHILE(XY<.1)
                  X2=X1-DX
                  IF(X2<=0.)EXIT
                  Y2=(-LOG(X2))**WW/SUM
                  Z2=PWindEro_Rate(Y2, ALG, USTW)      ! wind erosion rate     kg/m/s
                  XY=(Y1+Y2)*DX
                  XZ=(Z2+Z1)*DX
                  windErodSoil=windErodSoil+XZ
                  X1=X2
                  Y1=Y2
                  Z1=Z2
            END DO
            DX=DX*.5
      END DO
      windErodSoil=.5*windErodSoil
END SUBROUTINE PWind_Erosion
    
! ---------------------2.4 Daily Soil Loss by Wind Erosion ----------------------
SUBROUTINE SoiLoss_WindEro(JRT,  WindDir )  
      IMPLICIT NONE  
      !     EPIC1102
      !     THIS SUBPROGRAM ESTIMATES DAILY SOIL LOSS CAUSED BY WIND EROSION,
      !     GIVEN THE AVERAGE WIND SPEED AND DIRECTION.
      !Ref: APEX-doc 135 - 143
      INTEGER:: JRT 
 
      ! local varialbes
      REAL:: ALG, WindDir, BT,RIF,RFB, RFC,RRF, X1, RTO, ROKF, USTW, USTT, USTRT, WK1,VF  
    
      JRT=0
      IF(U10<modelPara(67))THEN  ! wind erosion threshold wind speed(normal: 6.0, range: 4.0-10.0)
            JRT=1
            RETURN
      END IF
      ! --------------- unsheltered distance factor ----------------
      BT=PI2/4.+WindDir-fieldAngle                     ! Eq. 143b: angle of the wind relative to ridges   
      ALG=fieldLength*fieldWidth/(fieldLength*ABS(COS(BT))+fieldWidth*ABS(SIN(BT))) ! Eq. 143a
      ALG=1.-EXP(-ALG/.07)                      ! Eq. 143
      ! ---------------- roughtness factor  ------------------------
      RRF=11.9*(1.-EXP(-(RandRough/9.8)**1.3))  ! clod roghtness factor  RandRough: random roughtness (mm)
      RIF=ABS(SIN(BT))*(1.27*RidgeHeight2**.52) ! ridge roughtness factor, but !!! APEX does not use sin(BT)
      RFB=MAX(.1,RRF+RIF)
      RFC=.77*1.002**RidgeHeight2               ! APEX-doc Eq. 140
      Roughtness_Factor=1.
      X1=(10./RFB)**RFC
      IF(X1<10.)Roughtness_Factor=1.-EXP(-X1)  
      ! --------------- vegetative cove factor ----------------------
      X1=vegCoverFactor                   
      X1=MIN(10.,X1+windEroCrop(3,JD)*cropResidu(LD1))   ! Eq. 142a  Vegetative equivalent
      VF=1.-X1/(X1+EXP(S_Curve(13,1)-S_Curve(13,2)*X1))                ! Eq. 142:  wind erosion cover factor
      ! --------------- ST/WP ---------------------------------------      
      RTO=soilWater(LD1)/Wilt_Point(LD1)
      
      SMM(68,MO)=SMM(68,MO)+soilWater(LD1)             ! SMM(68, MO), VAR(MO): soil water content of top layer (10mm thick)
      VAR(68)=soilWater(LD1)
      NWDA=NWDA+1
      ! --------------- integrate wind erosion rate -----------------
      USTRT=.0161*SQRT(soilDiameter)                   ! Eq. 136b in APEX doc                                         
      USTT=USTRT*USTRT   
      USTW=USTT+.5*RTO*RTO
      ROKF=EXP(-.047*rockFrac(LD1))                    ! Coarse fraction factor
      CALL PWind_Erosion(ALG, USTW, USTRT )
      ! ----------------- soil erodibility factor by wind -----------
      WK1=Wind_Erod*EXP(-.001*TLMF)
      SMM(38,MO)=SMM(38,MO)+WK1                         ! SMM(38, MO), VAR(MO): Wind erosion soil erodibility factor
      VAR(38)=WK1
      ! ----------------- final equation of soil erosion by wind ------------
      windErodSoil=8640.*windErodSoil*Roughtness_Factor*VF*ROKF*WK1   ! APEX-doc Eq 135   kg/m   
      TLMF=TLMF+windErodSoil
END SUBROUTINE SoiLoss_WindEro          
    
! ------------------------------ 3. Thickness of soil removed by each erosion -------------------------------
! ---------------------------- 3.1 distribute variables between layers ------------------------
REAL FUNCTION DistributeVar(X1,X2,Z1,Z2,ZZ)
      IMPLICIT NONE
      REAL:: X1, X2, Z1, Z2, ZZ

      DistributeVar = (Z1*X1 + Z2*X2)/ZZ

END FUNCTION DistributeVar

! ---------------------------------- Thickness change --------------------------
SUBROUTINE Thickness_Ero(JRT)
      !     EPIC1102
      !     THIS SUBPROGRAM CALCULATES THE THICKNESS OF SOIL REMOVED BY EACH
      !     EROSION EVENT, MOVES THE TOP LAYER INTO THE SECOND LAYER BY A DIST
      !     EQUAL THE ERODED THICKNESS, AND ADJUSTS THE TOP LAYER SOIL
      !     PROPERTIES BY INTERPOLATION.  WHEN THE SECOND LAYER IS REDUCED TO
      !     TO A THICKNESS OF 10 MM, IT IS PLACED INTO THE THIRD LAYER.
      IMPLICIT NONE
      ! local variables
      INTEGER:: J, LD2, MXZ, JRT, MX1, L1
      REAL:: TK, W1, W2, Z1, Z2, VV, ZZ, ZMX, XX, RTO, RX, X3, XY  
    
      TK=.1*totSoilErosion/bulkDensity(LD1)
      THK=THK+TK
      TK=MIN(TK*.001,Z(LD1))
      J=Layer_ID(Actual_SoilLayers)
      W1=Z(LD1)-TK
      W2=TK
      WT(LD1)=WT(LD1)-totSoilErosion
      IF(Z(Layer_ID(2))-Z(LD1)-TK>.01)GO TO 10
      ! REMOVE LAYER 2 AND PLACE SMALL REMAINING CONTENTS IN LAYER 3
      ISL=Layer_ID(3)
      LD2=Layer_ID(2)
      Z1=Z(LD2)-Z(LD1)
      Z2=Z(ISL)-Z(LD2)
      VV=Z1+Z2
      ZZ=.01/VV
      sorbRatioP(ISL)=DistributeVar(sorbRatioP(LD2),sorbRatioP(ISL),Z1,Z2,VV)
      mineralBulkDensity(ISL)=DistributeVar(mineralBulkDensity(LD2),mineralBulkDensity(ISL),Z1,Z2,VV)
      clayFrac(ISL)=DistributeVar(clayFrac(LD2),clayFrac(ISL),Z1,Z2,VV)
      siltFrac(ISL)=DistributeVar(siltFrac(LD2),siltFrac(ISL),Z1,Z2,VV)
      sandFrac(ISL)=100.-clayFrac(ISL)-siltFrac(ISL)
      rockFrac(ISL)=DistributeVar(rockFrac(LD2),rockFrac(ISL),Z1,Z2,VV)
      satuateCondv(ISL)=DistributeVar(satuateCondv(LD2),satuateCondv(ISL),Z1,Z2,VV)
      PH(ISL)=DistributeVar(PH(LD2),PH(ISL),Z1,Z2,VV)
      WT(ISL)=WT(ISL)+WT(LD2)
      bulkDensity(ISL)=.0001*WT(ISL)/VV
      NO3_N_Soil(ISL)=NO3_N_Soil(ISL)+NO3_N_Soil(LD2)
      NPassiveHumus(ISL)=NPassiveHumus(ISL)+NPassiveHumus(LD2)
      NSlowHumus(ISL)=NSlowHumus(ISL)+NSlowHumus(LD2)
      NBiomass(ISL)=NBiomass(ISL)+NBiomass(LD2)
      NStructLitt(ISL)=NStructLitt(ISL)+NStructLitt(LD2)
      NMetabLitt(ISL)=NMetabLitt(ISL)+NMetabLitt(LD2)
      CPassiveHumus(ISL)=CPassiveHumus(ISL)+CPassiveHumus(LD2)
      CSlowHumus(ISL)=CSlowHumus(ISL)+CSlowHumus(LD2)
      CBiomass(ISL)=CBiomass(ISL)+CBiomass(LD2)
      structLitt(ISL)=structLitt(ISL)+structLitt(LD2)
      metabLitt(ISL)=metabLitt(ISL)+metabLitt(LD2)
      lgStructLitt(ISL)=lgStructLitt(ISL)+lgStructLitt(LD2)
      CStructLitt(ISL)=CStructLitt(ISL)+CStructLitt(LD2)
      CMetabLitt(ISL)=CMetabLitt(ISL)+CMetabLitt(LD2)
      CLgStructLitt(ISL)=CLgStructLitt(ISL)+CLgStructLitt(LD2)
      NLgStructLitt(ISL)=CStructLitt(ISL)-CLgStructLitt(ISL)
      SOP(ISL)=SOP(ISL)+SOP(LD2)
      cropResidu(ISL)=.001*(structLitt(ISL)+metabLitt(ISL))
      labileP(ISL)=labileP(ISL)+labileP(LD2)
      activeMineralP(ISL)=activeMineralP(ISL)+activeMineralP(LD2)
      FreshOrgP_Residu(ISL)=FreshOrgP_Residu(ISL)+FreshOrgP_Residu(LD2)
      stableMineralP(ISL)=stableMineralP(ISL)+stableMineralP(LD2)
      Wilt_Point(ISL)=Wilt_Point(ISL)+Wilt_Point(LD2)
      fieldCapacity(ISL)=fieldCapacity(ISL)+fieldCapacity(LD2)
      Porosity(ISL)=Porosity(ISL)+Porosity(LD2)
      CALL WiltPoint_Sat(ISL)
      soilWater(ISL)=soilWater(ISL)+soilWater(LD2)
      LORG(LD1)=LORG(LD2)
      ! SPLIT LAYER NEAREST SURFACE WITH THICKNESS > 0.15 M IN HALF
      ZMX=0.
      L1=LD1
      MXZ=LD2
      DO J=3,Actual_SoilLayers
            ISL=Layer_ID(J)
            ZZ=Z(ISL)-Z(L1)
            IF(ZZ>.15)THEN
                  MXZ=J
                  ZMX=ZZ
                  GO TO 7
            END IF
            L1=ISL
            IF(ZZ<=ZMX+.01)CYCLE
            ZMX=ZZ
            MXZ=J
      END DO
      ISL=Layer_ID(MXZ)
      L1=Layer_ID(MXZ-1)
      IF(ZMX>minThickMaxLayer)GO TO 7
      Actual_SoilLayers=Actual_SoilLayers-1
      IF(Actual_SoilLayers<=2)THEN
            JRT=1
            RETURN
      END IF
      DO J=2,Actual_SoilLayers
            Layer_ID(J)=Layer_ID(J+1)
      END DO
      GO TO 10
    7 MX1=MXZ-1
      IF(MX1>2)THEN
            DO J=2,MX1
                  Layer_ID(J)=Layer_ID(J+1)
            END DO
      END IF
      Layer_ID(MX1)=LD2
      LORG(LD2)=LORG(Layer_ID(MXZ))
      CALL Soil_TopLayerMove(ISL,L1,LD2,1,.5)
      ! ADJUST LAYER 1 BY INTERPOLATING BETWEEN LAYER 2 USING ERODED THICK
   10 ISL=Layer_ID(2)
      sorbRatioP(LD1)=DistributeVar(sorbRatioP(LD1),sorbRatioP(ISL),W1,W2,Z(LD1))
      mineralBulkDensity(LD1)=DistributeVar(mineralBulkDensity(LD1),mineralBulkDensity(ISL),W1,W2,Z(LD1))
      clayFrac(LD1)=DistributeVar(clayFrac(LD1),clayFrac(ISL),W1,W2,Z(LD1))
      siltFrac(LD1)=DistributeVar(siltFrac(LD1),siltFrac(ISL),W1,W2,Z(LD1))
      sandFrac(LD1)=100.-clayFrac(LD1)-siltFrac(LD1)
      bulkDensity(LD1)=DistributeVar(bulkDensity(LD1),bulkDensity(ISL),W1,W2,Z(LD1))
      satuateCondv(LD1)=DistributeVar(satuateCondv(LD1),satuateCondv(ISL),W1,W2,Z(LD1))
      ! LateralFlow_Cond(LD1)=satuateCondv(LD1)*uplandSteep
      PH(LD1)=DistributeVar(PH(LD1),PH(ISL),W1,W2,Z(LD1))
      XX=Z(ISL)-Z(LD1)
      RootZone_Depth=XX-TK
      RTO=RootZone_Depth/XX
      RX=1.
      IF(rockFrac(LD1)>0.)THEN
            X3=rockFrac(LD1)
            rockFrac(LD1)=DistributeVar(rockFrac(LD1),rockFrac(ISL),W1,W2,Z(LD1))
            RX=(100.-rockFrac(LD1))/(100.-X3)
      END IF
      W1=W1/Z(LD1)
      W2=W2/XX
      Wilt_Point(LD1)=DistributeVar(Wilt_Point(LD1),Wilt_Point(ISL),W1,W2,RX)
      fieldCapacity(LD1)=DistributeVar(fieldCapacity(LD1),fieldCapacity(ISL),W1,W2,RX)
      Porosity(ISL)=Porosity(ISL)*RTO
      Wilt_Point(ISL)=Wilt_Point(ISL)*RTO
      fieldCapacity(ISL)=fieldCapacity(ISL)*RTO
      CALL WiltPoint_Sat(ISL )
      WT(ISL)=WT(ISL)*RTO
      Porosity(LD1)=DistributeVar(Porosity(LD1),Porosity(ISL),W1,W2,RX)
      CALL WiltPoint_Sat(LD1)
      NO3_N_Soil(LD1)=NO3_N_Soil(LD1)+Soil_Move(NO3_N_Soil(ISL),W2)
      NPassiveHumus(LD1)=NPassiveHumus(LD1)+Soil_Move(NPassiveHumus(ISL),W2)
      NSlowHumus(LD1)=NSlowHumus(LD1)+Soil_Move(NSlowHumus(ISL),W2)
      NBiomass(LD1)=NBiomass(LD1)+Soil_Move(NBiomass(ISL),W2)
      NStructLitt(LD1)=NStructLitt(LD1)+Soil_Move(NStructLitt(ISL),W2)
      NMetabLitt(LD1)=NMetabLitt(LD1)+Soil_Move(NMetabLitt(ISL),W2)
      CPassiveHumus(LD1)=CPassiveHumus(LD1)+Soil_Move(CPassiveHumus(ISL),W2)
      CSlowHumus(LD1)=CSlowHumus(LD1)+Soil_Move(CSlowHumus(ISL),W2)
      CBiomass(LD1)=CBiomass(LD1)+Soil_Move(CBiomass(ISL),W2)
      structLitt(LD1)=structLitt(LD1)+Soil_Move(structLitt(ISL),W2)
      metabLitt(LD1)=metabLitt(LD1)+Soil_Move(metabLitt(ISL),W2)
      lgStructLitt(LD1)=lgStructLitt(LD1)+Soil_Move(lgStructLitt(ISL),W2)
      CStructLitt(LD1)=CStructLitt(LD1)+Soil_Move(CStructLitt(ISL),W2)
      CMetabLitt(LD1)=CMetabLitt(LD1)+Soil_Move(CMetabLitt(ISL),W2)
      CLgStructLitt(LD1)=CLgStructLitt(LD1)+Soil_Move(CLgStructLitt(ISL),W2)
      NLgStructLitt(LD1)=CStructLitt(LD1)-CLgStructLitt(LD1)
      SOP(LD1)=SOP(LD1)+Soil_Move(SOP(ISL),W2)
      cropResidu(LD1)=.001*(structLitt(LD1)+metabLitt(LD1))
      labileP(LD1)=labileP(LD1)+Soil_Move(labileP(ISL),W2)
      activeMineralP(LD1)=activeMineralP(LD1)+Soil_Move(activeMineralP(ISL),W2)
      FreshOrgP_Residu(LD1)=FreshOrgP_Residu(LD1)+Soil_Move(FreshOrgP_Residu(ISL),W2)
      stableMineralP(LD1)=stableMineralP(LD1)+Soil_Move(stableMineralP(ISL),W2)
      solubleK(LD1)=solubleK(LD1)+Soil_Move(solubleK(ISL),W2)
      exchangeK(LD1)=exchangeK(LD1)+Soil_Move(exchangeK(ISL),W2)
      fixK(LD1)=fixK(LD1)+Soil_Move(fixK(ISL),W2)
      soilWater(LD1)=soilWater(LD1)+Soil_Move(soilWater(ISL),W2)
      XY=.01*(1.-.01*rockFrac(LD1))
      XX=XY
      WT(LD1)=bulkDensity(LD1)*100.
      aveBulkDensity=WT(LD1)
      DO J=2,Actual_SoilLayers
            ISL=Layer_ID(J)
            Z(ISL)=Z(ISL)-TK
            aveBulkDensity=aveBulkDensity+WT(ISL)
      END DO
      XX=Z(Layer_ID(Actual_SoilLayers))
      IF(XX<=minThickProfile)THEN
            JRT=1
            RETURN
      END IF
      IF(BIG>=XX)THEN
            BIG=XX
            plawDepth=XX
      END IF
      aveBulkDensity=aveBulkDensity*1.E-4/Z(Layer_ID(Actual_SoilLayers))
      layersEqualThick=INT(Z(Layer_ID(Actual_SoilLayers))/layerThickness+1.)
END SUBROUTINE Thickness_Ero 

END MODULE Erosion_Module