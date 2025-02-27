MODULE Crop_Module
USE PARM
USE Weather_Generator
USE Hydrology
USE Management 
IMPLICIT NONE

    CONTAINS

!------------------------- updata crop parms --------------------
SUBROUTINE Update_CropParms(LRG)
      IMPLICIT NONE
      ! arguments list 
      INTEGER, INTENT(INOUT):: LRG
      ! local variables
      REAL, DIMENSION(3,3):: AKX(3, 3) = &
      RESHAPE([1.0307,-3.0921,2.0614,-.040816,4.1224,-4.0816,.010101,-1.0303,2.0202], [3, 3]) 
      REAL, DIMENSION(3):: XTP
      INTEGER:: I, J, K
      REAL:: X1, X2

      cropLoop: DO I=1,LC   ! number of crops                                                                 
            NHU(I)=IHU(I)                             ! What is IHU, NHU?  
          
            IF(NHU(I)>LRG)LRG=NHU(I)                  ! What is LRG?  
          
            IF(maxRootDep(I)>RZ)RZ=maxRootDep(I)    ! maxRootDep: maximum root depth in m for crop i  
          
            IF(cropCode(I)==plantCategoryCode(7).OR.cropCode(I)==plantCategoryCode(8).OR.cropCode(I)==&
                        plantCategoryCode(10)) XMTU(I)=(1.-EXP(-turnoutFracCottn(I)/XMTU(I)))/144.  

            ligninFrac(3,I)=ligninFrac(2,I)  
            ligninFrac(1,I)=ligninFrac(1,I)/ligninFrac(2,I) 
            ligninFrac(2,I)=.99     
          
            CALL Coeff_Scurve(ligninFrac(1,I),ligninFrac(2,I),.5,1.)
          
            IF(NP_uptakeCurve==0)THEN
                  CALL Cal_NP_Biomass(uptakeParaN(1,I),uptakeParaN(2,I),uptakeParaN(3,I),uptakeParaN(4,I)) 
                  CALL Cal_NP_Biomass(uptakeParaP(1,I),uptakeParaP(2,I),uptakeParaP(3,I),uptakeParaP(4,I))
            ELSE           
                  uptakeParaN(4,I)=uptakeParaN(1,I)  
                  X1=uptakeParaN(1,I)-uptakeParaN(3,I)  
                  uptakeParaN(1,I)=1.-(uptakeParaN(2,I)-uptakeParaN(3,I))/X1  
                  uptakeParaN(2,I)=1.-.00001/X1  
                  CALL Coeff_Scurve(uptakeParaN(1,I),uptakeParaN(2,I),.5,1.)  
                  uptakeParaP(4,I)=uptakeParaP(1,I)  
                  X1=uptakeParaP(1,I)-uptakeParaP(3,I)  
                  uptakeParaP(1,I)=1.-(uptakeParaP(2,I)-uptakeParaP(3,I))/X1  
                  uptakeParaP(2,I)=1.-.00001/X1   
                  CALL Coeff_Scurve(uptakeParaP(1,I),uptakeParaP(2,I),.5,1.)   
            END IF 
          
            DO K=1,3  
                  XTP(K)=0.  
                  DO J=1,3  
                        XTP(K)=XTP(K)+uptakeParaK(J,I)*AKX(K,J)        ! AKX are coefficients  
                  END DO                                                                            
            END DO   
          
            uptakeParaK(1,I)=XTP(1)    
            uptakeParaK(2,I)=XTP(2)   
            uptakeParaK(4,I)=XTP(3)  
            IHU(I)=1                                                                       
            X1=split(pointsLAIDevp(1,I))*.01                             ! LAI development curve  
            X2=split(pointsLAIDevp(2,I))*.01                                                        
            CALL Coeff_Scurve(pointsLAIDevp(1,I),pointsLAIDevp(2,I),X1,X2)                                          
            X1=split(pointsFrostDamg(1,I))                               ! Frost damage curve  
            X2=split(pointsFrostDamg(2,I))                                                            
            CALL Coeff_Scurve(pointsFrostDamg(1,I),pointsFrostDamg(2,I),X1,X2)  
            CO2EffOnBio(1,I)=RUE(I)*.01  
            X2=split(CO2EffOnBio(2,I))                                   ! CO2 effects on biomass curve  
            CALL Coeff_Scurve(CO2EffOnBio(1,I),CO2EffOnBio(2,I),330.,X2)                                        
            X2=split(VPDPara2(I))                                                              
            VPDPara2(I)=(1.-VPDPara2(I))/(X2-VPDThreshold(I))                                              
            UNA(I)=modelPara(39)*uptakeParaN(3,I)*RUE(I)*modelPara(28)   !    
            ULYN(I)=UNA(I)                                                                      
            BLYN(I)=0.
      END DO   cropLoop

END SUBROUTINE Update_CropParms

! --------------------- 6. Relationships between NP and RUE Accumulation ----------
SUBROUTINE Cal_NP_Biomass(P0,P5,P1,A)
      !     THIS SUBPROGRAM COMPUTES PARAMETERS OF AN EQUATION DESCRIBING THE
      !     N AND P RELATIONS TO BIOMASS ACCUMULATION.
      IMPLICIT NONE
      ! local variables
      INTEGER:: I
      REAL:: A, A5, EA, EA1, EG, EG1, P0, P0G, PEG, P1, P01, P5, X1, PG5, FU, DFDA
 
      A=5.0
      DO I=1,10                  !repeat the calculation 10 times, and I is not involved
            A5=A*.5
            EA=EXP(A)
            EA1=EA-1.
            EG=EXP(-A5)
            P0G=P0*EG
            EG1=EXP(A5)
            PEG=P1*(EA-EG1)
            P01=P0*(1.-EG)
            X1=PEG-P01
            PG5=.5*P0G
            FU=X1/EA1+P0G-P5
            IF(ABS(FU)<1.E-7)GO TO 3
            DFDA=(EA1*(P1*(EA-.5*EG1)-PG5)-EA*X1)/(EA1*EA1)-PG5
            A=A-FU/DFDA
      END DO
      WRITE(KW(1),4)A,FU
    3 P5=(P1*EA-P0)/EA1
      P0=P0-P5
 
    4 FORMAT(//T10,'Relation_NP_Biomass DID NOT CONVERGE',2E16.6)
 
END SUBROUTINE Cal_NP_Biomass          
! ------------------------- 1. Cal LAI --------------------------------
SUBROUTINE LAI_HU_RootDep(JRT )
      !     EPIC1102
      !     THIS SUBPROGRAM CALCUALTES LEAF AREA INDEX, HEAT UNITS, ROOT DEPTH
      !     AND TEMPERATURE STRESS FOR THE.
      ! LAI simulation is different from SWAT, seems like APEX
      ! Plant height is the same as SWAT
      IMPLICIT NONE
      ! arguments list
      INTEGER, INTENT(OUT):: JRT                        
      ! local variables
      INTEGER:: I, K1, N1
      REAL, SAVE:: SLA0(12)         ! The SAVE attribute is very important
      REAL:: X1, X2, X3, X4, XX, XY, XZ, XW, XUN, XUP, CPR, CNR, CKR, XPHU, TGX, F, FF, F1, F2, F3, &
       SLAX, SUM,TOT, ADD, RTO, DLAX, W1, FHR    

      ! --------------------- HU ------------------------ 
                      
      JRT=0
      X1=totCropBio(Crop_Num)+1.E-10
      CPR=actualCropP(Crop_Num)/X1            ! UP1: the actual P content of the crop in kg/ha 
      CNR=actualCropN(Crop_Num)/X1            ! UN1: the actual N content of the crop in kg/ha
      CKR=MAX(1.E-5,actualCropK(Crop_Num)/X1) ! UK1
      AJWA(Crop_Num) =1.
      XPHU=potentialHeatUnit(Crop_Num,IHU(Crop_Num)) ! Why PHU is so large an array? (30,150)  30: for 30 plants; 150 for what?
      X4=optimalTem(Crop_Num)-plantMinT(Crop_Num)
      TGX=TX-plantMinT(Crop_Num)
      HU(Crop_Num)=HU(Crop_Num)+MAX(0.,TGX)    ! HU: number of heat units accumulated during a day in oC
      IF(DayOfYear==JDHU.AND.cropCode(Crop_Num)/=plantCategoryCode(7).AND.cropCode(Crop_Num)/=plantCategoryCode(8)&
       .AND.cropCode(Crop_Num)/=plantCategoryCode(10).AND.cropCode(Crop_Num)/=plantCategoryCode(4))THEN
            HU(Crop_Num)=XPHU*modelPara(19)    ! 19: fraction of maturity at spring growth initiation allows fall growing crops to reset heat unit index
            pestGrowIndx=MIN(0.,pestGrowIndx) ! PSTS: accumulated daily pest index
            GrowingSeason_Days=0                ! IPST: Growing season lengths in day                 
      END IF
      HUI(Crop_Num)=HU(Crop_Num)/XPHU         ! HUI: heat unit index is daily HU values devided by potential heat units of a crop
      IF(HU(Crop_Num)>XPHU)THEN
            WaterFrac_Yield=MAX(waterFracYield(Crop_Num),WaterFrac_Yield-EO*.002)
            IF(cropCode(Crop_Num)==plantCategoryCode(3).OR.cropCode(Crop_Num)==plantCategoryCode(6))THEN
                  ! 3 - prennial legume ; 6 - perennial
                  HU(Crop_Num)=0.
                  JRT=2
            ELSE
                  IF(cropCode(Crop_Num)==plantCategoryCode(7).OR.cropCode(Crop_Num)==plantCategoryCode(8).OR. &
                              cropCode(Crop_Num)==plantCategoryCode(10))THEN
                        ! 7 - Evergreen tree; 8 - Deciduous tree; 10 - N-fixing tree
                        JRT=0
                        HU(Crop_Num)=.9*HU(Crop_Num)
                  ELSE
                        JRT=1
                  END IF
            END IF
            RETURN
      END IF      
      ! ----------------------- LAI simulation ----------------------
      
      F2=HUI(Crop_Num)/(HUI(Crop_Num)+EXP(pointsLAIDevp(1,Crop_Num)-pointsLAIDevp(2,Crop_Num)*HUI(Crop_Num))) ! F2: = HUF in APEX-doc Eqs 276a
	  ! F2 is also the fr_LAImax in SWAT-doc  Eq 5:2.1.10
      IF(cropCode(Crop_Num)==plantCategoryCode(7).OR.cropCode(Crop_Num)==plantCategoryCode(8).OR.&
         cropCode(Crop_Num)==plantCategoryCode(10))THEN
            ! 7 - Evergreen tree; 8 - Deciduous tree; 10 - N-fixing tree
            X1=HSM/yearTotHeat 
            F1=X1/(X1+EXP(pointsLAIDevp(1,Crop_Num)-pointsLAIDevp(2,Crop_Num)*X1))
            F=F1
            XLAI(Crop_Num)=MAX(.1,DMLX(Crop_Num)* F2)
            F3=SQRT(HUI(Crop_Num)+1.E-10)
      ELSE
            F=F2 
            F3=SQRT(F2+1.E-10)
      END IF
	  
      FF=F-WLV(Crop_Num)    ! WLV: LAI fraction of last day
      XX=FF*XLAI(Crop_Num)  ! XLAI: maximum LAI for the plant
      X2=1.
      SLAX=0.
      X3=(Current_LAI(Crop_Num)+.001)*CPHT(Crop_Num)
      IF(IGO==1)THEN
            SUM=X3
            TOT=abvGroundBiom(Crop_Num)
            ADD=Current_LAI(Crop_Num)
      ELSE    
            SUM=0.
            TOT=0.
            ADD=0.
            DO I=1,IGO
                  K1=JE(I)    
                  IF(K1>MNC)CYCLE
                  IF(Current_LAI(K1)>SLAX)SLAX=Current_LAI(K1)
                  SUM=SUM+Current_LAI(K1)*CPHT(Crop_Num)
                  TOT=TOT+abvGroundBiom(K1)
                  ADD=ADD+Current_LAI(K1)
            END DO
            IF(SLAX>2.)X2=X3/SUM
      END IF
      IF(XX>0.)THEN
            X1=XX*X2*(1.+Delta_Daylength)**modelPara(70)  ! Ref: EPIC-doc Eq 2.193 Pg52, this equation is for biomass but here used for LAI?
            IF(cropCode(Crop_Num)/=plantCategoryCode(7).AND.cropCode(Crop_Num)/=plantCategoryCode(8).AND.&
                 cropCode(Crop_Num)/=plantCategoryCode(10)) X1=X1*SQRT(Min_Stress_Factor(Crop_Num))*SHRL 
                  ! 7, 8 , 10 are trees
                  ! REG : plant growth factor = 1 - max(stress factors)   Ref: SWAT-doc 2009 Eq 5:3.2.3
                     
            Current_LAI(Crop_Num)=MIN(XLAI(Crop_Num),Current_LAI(Crop_Num)+X1)    ! X1: delta_LAI  Current_LAI: real LAI
          
      END IF
      WLV(Crop_Num)=F
      RTO=TGX/X4
      IF(TGX>0..AND.RTO<2.)THEN 
            Min_Stress_Factor(Crop_Num)=SIN(1.5707*RTO)    ! same equation as RootGrowth_Factor 1.5707= pi/2
      ELSE
            Min_Stress_Factor(Crop_Num)=0.
      END IF
      
      ! --------------------- crop height ---------------------
      
      IF(Current_LAI(Crop_Num)<.001)THEN
            Current_LAI(Crop_Num)=.001
      ELSE         
            FF=MAX(0.,F3-WCHT(Crop_Num))
            CPHT(Crop_Num)=MIN(maxCropHeight(Crop_Num),CPHT(Crop_Num)+FF*maxCropHeight(Crop_Num))
            ! LfWidth was added by TXH: max leaf width is 0.1 m
            LfWidth(Crop_Num)=MIN(0.1,LfWidth(Crop_Num)+FF*0.1)
          
            IF(HUI(Crop_Num)>XDLAI(Crop_Num).AND.DayLength>winterDormancy.AND.XDLA0(Crop_Num)>0.)THEN
                  XX=(1.-HUI(Crop_Num))/XDLA0(Crop_Num)    ! ref: APEX-doc Eq 277
                  IF(XX>1.E-5)THEN
                        XX=LOG10(XX)
                  ELSE
                        XX=-5.
                  END IF
                  IF(cropCode(Crop_Num)/=plantCategoryCode(7).AND.cropCode(Crop_Num)/=plantCategoryCode(8).AND.&
                             cropCode(Crop_Num)/=plantCategoryCode(10))THEN
                        RTO=factorLAIDecline(Crop_Num)*XX
                        IF(RTO<-10.)RTO=-10.
                        X1=SLA0(Crop_Num)*10.**RTO
                        IF(X1<Current_LAI(Crop_Num))THEN
                              DLAX=Current_LAI(Crop_Num)-X1
                              Current_LAI(Crop_Num)=X1
                              Current_LAI(Crop_Num)=MAX(.001,X1)
                              X2=DLAX*abvGroundBiom(Crop_Num)*leafWeightFrac(Crop_Num)
                              XZ=MIN(actualCropN(Crop_Num),1000.*fracYieldN(Crop_Num)*X2)
                              XY=MIN(actualCropP(Crop_Num),1000.*fracYieldP(Crop_Num)*X2)
                              W1=1000.*fracYieldK(Crop_Num)*X2
                      
                              IF(X2>0.) CALL CN_TopSoil(X2,XZ,0 )
                      
                              FreshOrgP_Residu(LD1)=FreshOrgP_Residu(LD1)+XY
                              actualCropN(Crop_Num)=actualCropN(Crop_Num)-XZ
                              actualCropP(Crop_Num)=actualCropP(Crop_Num)-XY
                        END IF    
                  END IF
                  RTO=factorBioDecline(Crop_Num)*XX
                  IF(RTO<-10.)RTO=-10.
                  AJWA(Crop_Num)=10.**RTO
            ELSE
                  XDLA0(Crop_Num)=1.-XDLAI(Crop_Num)
                  SLA0(Crop_Num)=Current_LAI(Crop_Num)
            END IF              
      END IF
      WCHT(Crop_Num)=F3
      XX=MAX(CPHT(Crop_Num),RD(Crop_Num),2.5*maxRootDep(Crop_Num)*HUI(Crop_Num))
      RD(Crop_Num)=MIN(maxRootDep(Crop_Num),Z(Layer_ID(Actual_SoilLayers)),XX)
      GroundCover_Frac=ADD/(ADD+EXP(S_Curve(23,1)-S_Curve(23,2)*ADD)) ! FGC: fraction plant ground cover
      GroundCover_StandLiveBio=TOT/(TOT+EXP(S_Curve(26,1)-S_Curve(26,2)*TOT))   ! FGSL: plant ground cover as a function of standing live biomass
      CLG=ligninFrac(3,Crop_Num)*HUI(Crop_Num)/(HUI(Crop_Num)+EXP(ligninFrac(1,Crop_Num)-ligninFrac(2,Crop_Num)*&
      HUI(Crop_Num)))
      SHRL=1.
      FHR=0.
      IF(DayLength+1.E-5>winterDormancy)THEN
            SRA=SRA+SRAD
            GSVP=GSVP+VPD
      ELSE
            SHRL=0.
            FHR=1.-DayLength/winterDormancy
      END IF
      !IF(cropCode(Crop_Num)/=plantCategoryCode(7))THEN
      IF(TMN>-1.)THEN             ! ??? questionable no end if
            IF(SHRL>0.)GO TO 15
                  F=FHR
            ELSE
                  XX=ABS(TMN)
                  F=XX/(XX+EXP(pointsFrostDamg(1,Crop_Num)-pointsFrostDamg(2,Crop_Num)*XX))
                  F=MAX(F,FHR)
            END IF
            IF(abvGroundBiom(Crop_Num)>0..AND.abvGroundBiom(Crop_Num)<totCropBio(Crop_Num).AND. &
                        cropCode(Crop_Num)/=plantCategoryCode(7).AND.cropCode(Crop_Num)/=plantCategoryCode(8).AND.&
                        cropCode(Crop_Num)/=plantCategoryCode(10))THEN
                  abvGroundBiom(Crop_Num)=MAX(0.,totCropBio(Crop_Num)-totRootWeight(Crop_Num))
                  XX=F*abvGroundBiom(Crop_Num)
                  abvGroundBiom(Crop_Num)=abvGroundBiom(Crop_Num)-XX
                  totCropBio(Crop_Num)=totRootWeight(Crop_Num)+abvGroundBiom(Crop_Num)
                  standCropResi(Crop_Num)=standCropResi(Crop_Num)+XX
                  standDeadLignin=standDeadLignin+CLG*XX
                  XY=XX*CNR
                  XZ=XX*CPR
                  XW=XX*CKR
                  XUN=actualCropN(Crop_Num)
                  IF(XUN-XY<.01) XY=XUN-.01
                  XUP=actualCropP(Crop_Num)
                  IF(XUP-XZ<.01) XZ=XUP-.01
                  standDeadResiK=standDeadResiK+XW
                  standDeadResiN(Crop_Num)=standDeadResiN(Crop_Num)+XY
                  standDeadResiP=standDeadResiP+XZ
                  !actualCropK(Crop_Num)=actualCropK(Crop_Num)-XW
                  actualCropN(Crop_Num)=XUN-XY
                  actualCropP(Crop_Num)=XUP-XZ
            END IF
            IF(cropCode(Crop_Num)==plantCategoryCode(7)) F=F*F
            Current_LAI(Crop_Num)=MAX(.001,Current_LAI(Crop_Num)*(1.-F))
            !END IF
         15 IF(Min_Stress_Factor(Crop_Num)>0.) RETURN
            N1=NCP(Crop_Num)
            Growing_Stress(5,N1,Crop_Num)=Growing_Stress(5,N1,Crop_Num)+1.
            SFMO(5,Crop_Num)=SFMO(5,Crop_Num)+1.
            growStressFactor(5,Crop_Num)=1.
            BLYN(Crop_Num)=BLYN(Crop_Num)+1.
            SMMC(17,Crop_Num,MO)=SMMC(17,Crop_Num,MO)+1.
            cropVar(17,Crop_Num)=0.
            JRT=1
      RETURN
 
END SUBROUTINE LAI_HU_RootDep

! ---------------------------------- 2. Water and nutrient usage -------------------------------------------
! ------------------- 2.2 determine minimum stress fractor for root and total biomass ------------
SUBROUTINE Determine_MinStress(I, J, F, RootGrowth_Factor, X, JRT)
      !     EPIC1102
      !     DETERMINES MINIMUM STRESS FACTORS FOR ROOT AND TOTAL BIOMASS
      !     GROWTH.
      IMPLICIT NONE
      ! local variables
      INTEGER:: JRT, J, I
      REAL:: F, RootGrowth_Factor, X
      JRT=0
      IF(F>=RootGrowth_Factor)RETURN
      RootGrowth_Factor=MAX(0.,F)
      J=I
      IF(RootGrowth_Factor<=X)JRT=1
      RETURN
END SUBROUTINE

! ----------------------------------2.3 Root Growth Stress ----------------------------------
SUBROUTINE Root_Stress(JRT, RootGrowth_Factor)
      !     EPIC1102
      !     THIS SUBPROGRAM CALCULATES ROOT GROWTH STRESSES CAUSED BY
      !     TEMPERATURE,  ALUMINUM TOXICITY, AND SOIL STRENGTH AND DETERMINES
      !     THE ACTIVE CONSTRAINT ON ROOT GROWTH (THE MINIMUM STRESS FACTOR).
      IMPLICIT NONE
 
      ! local variables 
      INTEGER:: II, JRT
      REAL:: RootGrowth_Factor, XX, A0, F

      RootGrowth_Factor=1.
      IF(modelPara(2)>1.99)RETURN
      ! root growing stress by Temperature
      II=3
      XX=soilTem(ISL)/optimalTem(Crop_Num)
      IF(XX<=0.)THEN
            RootGrowth_Factor=0.
            GO TO 4
      END IF
      IF(XX<1.)RootGrowth_Factor=SIN(1.5708*XX)
      ! root growing stress by Al concentration
      A0=10.+(cropTolerantAl(Crop_Num)-1.)*20.                   ! Ref: APEX-doc Eq. 301a
      IF(ALS(ISL)>A0)THEN
            F=(100.-ALS(ISL))/(100.-A0)
            CALL Determine_MinStress(2,II,F,RootGrowth_Factor,.1,JRT)
            IF(JRT>0)GO TO 4
      END IF
      CALL BulkDens_Effect(postTillBulkDensity(ISL),modelPara(2),F,ISL,3 )
      XX=rockFrac(ISL)
      F=F*(1.-XX/(XX+EXP(S_Curve(1,1)-S_Curve(1,2)*XX)))  ! effect of course fragment content on root growth restriction
      CALL Determine_MinStress(1,II,F,RootGrowth_Factor,.1,JRT)
    4 STDA(II,Crop_Num)=STDA(II,Crop_Num)+(1.-RootGrowth_Factor)/Actual_SoilLayers
      IF(ISL==Layer_ID(Layer_RD_Reach))RGSM=RootGrowth_Factor
      RETURN

END SUBROUTINE Root_Stress
! ---------------------------------- 2.4 Cal water usage ---------------------------------------
SUBROUTINE Cal_WaterUsage(CPWU, RootGrowth_Factor, JRT)
      !     EPIC1102
      !     THIS SUBPROGRAM DISTRIBUTES PLANT EVAPORATION THROUGH THE ROOT
      !     ZONE AND CALCULATES ACTUAL PLANT WATER USE BASED ON SOIL WATER
      !     AVAILABILITY.
      !  Ref: APEX-doc  Eqs. 283
      IMPLICIT NONE
 					 
      ! local variables
      INTEGER:: JRT
      REAL:: BLM, CPWU, SUM, TOS, X1, XX, F, RootGrowth_Factor

      BLM=Wilt_Point(ISL)
      IF(Z(ISL)<=.5)BLM=modelPara(5)*Wilt_Point(ISL) ! para5: soil water lower limit of water content in top 0.5 m soil depth 
      IF(ISL/=LD1)THEN
            CALL Root_Stress( JRT, RootGrowth_Factor)
		 
            CPWU=CPWU*RootGrowth_Factor   ! product of all RGF above depth Z
      END IF
      SUM=EP(Crop_Num)*(1.-EXP(-UB1*RootZone_Depth/RD(Crop_Num)))/UOB ! APEX-doc Eq.283a : potential water use rate in mm/d at depth Z
      TOS=36.*elecCondv(ISL)
      XX=LOG10(Wilt_Point(ISL))
      X1=3.1761-1.6576*(LOG10(soilWater(ISL))-XX)/(LOG10(fieldCapacity(ISL))-XX)
      IF(X1<4.)THEN
            WTN=MAX(5.,10.**X1)
            XX=TOS+WTN
            IF(XX<5000.)THEN
                  F=1.-XX/(XX+EXP(S_Curve(21,1)-S_Curve(21,2)*XX)) ! XX: gravimetric + osmotic tension 
                  wiltPoint(ISL)=MIN(SUM-CPWU*plantEP-(1.-CPWU)*Poten_WaterUse_Rate,soilWater(ISL)-BLM)*&
                  RootGrowth_Factor*F
                  wiltPoint(ISL)=MAX(0.,wiltPoint(ISL))    ! mm/d
            END IF
      ELSE
            WTN=10000.
      END IF    
      Poten_WaterUse_Rate=SUM
END SUBROUTINE Cal_WaterUsage

! ---------------------------------- 2. Water and nutrient usage -----------------------------------
SUBROUTINE Water_Nutrient(JRT )
      !     EPIC1102
      !     THIS SUBPROGRAM IS THE MASTER WATER AND NUTRIENT USE SUBPROGRAM.
      !     CALLS HSWU AND NUPPO FOR EACH SOIL LAYER.
      IMPLICIT NONE
 
      ! local variables
      INTEGER:: IAR, J, JRT
      REAL:: dep, TOT, CPWU, RootGrowth_Factor, X1, X2, X3, X4, RTO, F  

      Layer_RD_Reach=1
      IAR=0
      Poten_WaterUse_Rate=0.      ! mm/d
      dep=0.
      TOT=0.
      CPWU=1.
      RootGrowth_Factor=1.
      DO J=1,Actual_SoilLayers
            ISL=Layer_ID(J)
            IF(Z(ISL)<1.)THEN
                  SAT=SAT+soilWater(ISL)
                  TOT=TOT+Porosity(ISL)
                  dep=Z(ISL)
            ELSE
                  IF(IAR==0)THEN
                        IAR=1
                        X3=1.- dep
                        X4=Z(ISL)- dep
                        RTO=X3/X4
                        IF(waterTableHigt<=Z(ISL))THEN
                              X1=Porosity(ISL)*(Z(ISL)-waterTableHigt)/X4
                              X2=soilWater(ISL)-X1
                              IF(waterTableHigt>1.)THEN
                                    SAT=SAT+X2*X3/(waterTableHigt- dep)
                              ELSE
                                    SAT=SAT+X2+Porosity(ISL)*(1.-waterTableHigt)/X4
                              END IF
                        ELSE
                              SAT=SAT+RTO*soilWater(ISL)
                        END IF    
                        TOT=TOT+RTO*Porosity(ISL)
                  END IF    
            END IF
            IF(Layer_RD_Reach>1)CYCLE
            IF(RD(Crop_Num)>Z(ISL))THEN
                  RootZone_Depth=Z(ISL)
            ELSE
                  RootZone_Depth=RD(Crop_Num)
                  Layer_RD_Reach=MAX(Layer_RD_Reach,J)
            END IF
            CALL Cal_WaterUsage(CPWU, RootGrowth_Factor, JRT )
          
            plantEP = plantEP + wiltPoint(ISL)
      END DO
      IF(Layer_RD_Reach==0)Layer_RD_Reach=Actual_SoilLayers
      ! Aeration stress
      IF(soilWaterRZ>potentialAvailWater)THEN
            RTO=MIN(1.,SAT/TOT)
            F=100.*(RTO-areationThreshold(Crop_Num))/(1.0001-areationThreshold(Crop_Num))
            IF(F>0.)THEN
                  SAT=1.-F/(F+EXP(S_Curve(7,1)-S_Curve(7,2)*F))
            ELSE
                  SAT=1.
            END IF    
      END IF
END SUBROUTINE Water_Nutrient

! --------------------------------- 3. Predicts Potential Growth and roots -------------------------

SUBROUTINE Poten_Growth0 
      !     EPIC1102
      !     THIS SUBPROGRAM PREDICTS DAILY POTENTIAL GROWTH OF TOTAL PLANT
      !     BIOMASS AND ROOTS.
      IMPLICIT NONE
 
      ! local variables
      INTEGER:: I, J
      REAL:: XX, X1, X2, X3, X4, SUM, XLA1, RTO, SLA1, TLI, XXC, PAR_RUE, VPDX
  
      XX=.685-.209*ROSP
      IF(IGO>1)THEN
            SUM=0.
            DO I=1,LC
                  J=KG(I)
                  IF(J==0)CYCLE
                  SUM=SUM+PPL0(J)
            END DO    
            XLA1=maxLAI(Crop_Num)*SUM/(SUM+EXP(coeffPopuCurve(1,Crop_Num)-coeffPopuCurve(2,Crop_Num)*SUM))
            RTO=XLA1/XLAI(Crop_Num)
            SLA1=RTO*Current_LAI(Crop_Num)
            TLI=1.-EXP(-XX*SLA1)
            XXC=-LOG(1.-TLI/REAL(IGO))/Current_LAI(Crop_Num)
      ELSE    
            XXC=XX    
      END IF    
      IF(PAR_Method==0)THEN
            X1=1.-EXP(-XXC*Current_LAI(Crop_Num))
      ELSE
            X1=EVI
      END IF
      ! ////////////////////////////////////////////////////////////////////////////
      !   my model works like this part
      PAR_RUE=.0005*SRAD*X1       ! PAR: intercepted photosynthetic active radiation in MJ/m2/d 0.0005 = 0.5*0.001
      X2=PAR_RUE*AJWA(Crop_Num)   !  
      IF(biomsConver_Method==0)THEN
            X1=MAX(VPD-1.,-.5)
            XX=RUE(Crop_Num)-VPDPara(Crop_Num)*X1   ! vpd is 0 for the ERL paper
            potentialBioIncrese(Crop_Num)=XX*X2    ! t/ha/day
            IF(potentialBioIncrese(Crop_Num)<1.E-10)potentialBioIncrese(Crop_Num)=0.
          !//////////////////////////////////////////////////////////////////////////////
      ELSE
            X2=X2*RUE(Crop_Num)
            VPDX=.67*(Sat_Vapor(TMX+273)-RHD*Sat_Vapor(TX+273.))
            X3=.01*WUE(Crop_Num)*VPDX**(-.5)
            X4= plantEP * X3
            potentialBioIncrese(Crop_Num)=MIN(X2,X4) 
      END IF
      totCropBio(Crop_Num)=totCropBio(Crop_Num)+potentialBioIncrese(Crop_Num)  
      ! WRITE(KW(51),"(1X, I4, 1X, I4, 3F8.3 )") IYR, Date, X2, potentialBioIncrese(Crop_Num), totCropBio(Crop_Num)
 
END SUBROUTINE Poten_Growth0

! ---------------------------------- 4. Nutrition Demand -----------------------------------------
! -------------------------- 4.1 N demand -----------------------------------------
SUBROUTINE N_Demand 
      !     EPIC1102
      !     THIS SUBPROGRAM CALCULATES THE DAILY N DEMAND FOR OPTIMAL PLANT GROWTH.
      ! Ref: SWAT-doc 2009 Section 5:2.3.1  N Uptake
      IMPLICIT NONE

      ! local variables
      REAL:: CNT
      
      IF(NP_uptakeCurve==0)THEN
            CNT=uptakeParaN(2,Crop_Num)+uptakeParaN(1,Crop_Num)*EXP(-uptakeParaN(4,Crop_Num)*HUI(Crop_Num))
      ELSE
            CNT=(uptakeParaN(4,Crop_Num)-uptakeParaN(3,Crop_Num))*(1.-HUI(Crop_Num)/(HUI(Crop_Num)+EXP(uptakeParaN(1,Crop_Num)-&
            uptakeParaN(2,Crop_Num)*HUI(Crop_Num))))+uptakeParaN(3,Crop_Num)
      END IF
      Optimal_N(Crop_Num)=CNT*totCropBio(Crop_Num)*1000.
      !UNMX(Crop_Num)=MAX(Optimal_N(Crop_Num),UNMX(Crop_Num))
      potentialDemandN=MIN(4000.*uptakeParaN(3,Crop_Num)*potentialBioIncrese(Crop_Num),Optimal_N(Crop_Num)-actualCropN(Crop_Num))
      potentialDemandN=MAX(0.,potentialDemandN)
      !potentialDemandN=Optimal_N(Crop_Num)-actualCropN(Crop_Num)
END SUBROUTINE N_Demand

! -------------------------------4.2 P demand --------------------------
SUBROUTINE P_Demand 
      !     EPIC1102     
      !     THIS SUBPROGRAM CALCULATES THE DAILY P DEMAND FOR OPTIMAL PLANT GROWTH.
      ! Ref: SWAT-doc  Eqs 5:2.3.19
      IMPLICIT NONE
 
      ! local variables
      REAL:: CPT

      IF(NP_uptakeCurve==0)THEN
            CPT=uptakeParaP(2,Crop_Num)+uptakeParaP(1,Crop_Num)*EXP(-uptakeParaP(4,Crop_Num)*HUI(Crop_Num))
      ELSE
            CPT=(uptakeParaP(4,Crop_Num)-uptakeParaP(3,Crop_Num))*(1.-HUI(Crop_Num)/(HUI(Crop_Num)+EXP(uptakeParaP(1,Crop_Num)-&
            &uptakeParaP(2,Crop_Num)*HUI(Crop_Num))))+uptakeParaP(3,Crop_Num)
      END IF
      Optimal_P(Crop_Num)=CPT*totCropBio(Crop_Num)*1000.    ! kg P/ha   totCropBio is t/ha
      UPMX(Crop_Num)=MAX(Optimal_P(Crop_Num),UPMX(Crop_Num))
      potentialDemandP=MIN(4000.*uptakeParaP(3,Crop_Num)*potentialBioIncrese(Crop_Num),Optimal_P(Crop_Num)-actualCropP(Crop_Num))
      potentialDemandP=MAX(0.,potentialDemandP)         ! Why not multiply 1.5?  SWAT-doc 2009 Eq. 5: 2.3.23
END SUBROUTINE P_Demand

! 4.3 --------------------------4.3 K demand ---------------------------------
SUBROUTINE K_Demand 
      !     EPIC1102
      !     THIS SUBPROGRAM CALCULATES THE DAILY K DEMAND FOR OPTIMAL PLANT
      !     GROWTH.
      IMPLICIT NONE
 
      ! local variables
      REAL:: CKT
 
      CKT=uptakeParaK(1,Crop_Num)+HUI(Crop_Num)*(uptakeParaK(2,Crop_Num)+HUI(Crop_Num)*uptakeParaK(4,Crop_Num))
      UK2(Crop_Num)=CKT*totCropBio(Crop_Num)*1000.
      IF(UK2(Crop_Num)<actualCropK(Crop_Num))UK2(Crop_Num)=actualCropK(Crop_Num)
      PlantK_Demand=MIN(4000.*uptakeParaK(3,Crop_Num)*potentialBioIncrese(Crop_Num),UK2(Crop_Num)-actualCropK(Crop_Num))
END SUBROUTINE K_Demand

! -----------------------------4.4 P supply ------------------------------------
SUBROUTINE NPK_Supply 
      !     EPIC1102
      !     THIS SUBPROGRAM CALCULATES THE DAILY POTENTIAL SOIL SUPPLY OF P
      !     FOR EACH LAYER.
      ! Ref: APEX-doc Eqs 289 - 291
      IMPLICIT NONE
 				   
      ! local variables
      INTEGER::J
      REAL:: XX,Amount_NO3N, F
 
      IF(totRootWeight(Crop_Num)>0.)THEN
            XX=1.5*potentialDemandP/totRootWeight(Crop_Num)  ! totRootWeight: root weight in layer 1 in kg/ha  
      ELSE
            XX=1.
      END IF
      DO J=1,Layer_RD_Reach
            ISL=Layer_ID(J)
            Amount_NO3N=modelPara(91)*(NO3_N_Soil(ISL)-.001*modelPara(27)*WT(ISL))+(1.-modelPara(91))*NH3_Weight(ISL) ! amount NO3-N in layer ISL (kg/ha)
            soilSupplyRateN(ISL)=MAX(0., Amount_NO3N*wiltPoint(ISL)/(soilWater(ISL)+.001))
            SUN=SUN+soilSupplyRateN(ISL)
            soilSupplyRateK(ISL)=solubleK(ISL)*wiltPoint(ISL)/(soilWater(ISL)+.001)
            SUK=SUK+soilSupplyRateK(ISL)
            F=1000.*labileP(ISL)/WT(ISL) !F == CSP, concentration of labile P in soil layer 1 in g/t 
            IF(F>30.)THEN
                  F=1.
            ELSE
                  F=F/(F+EXP(8.0065-.3604*F))   ! Eq. 289a in APEX-doc P63  
            END IF
            soilSupplyRateP(ISL)=XX*F*RWT(ISL,Crop_Num)   ! RWT: root weight   
            IF(soilSupplyRateP(ISL)>=labileP(ISL))soilSupplyRateP(ISL)=.9*labileP(ISL)
            SUP=SUP+soilSupplyRateP(ISL)        ! kg/ha/d
      END DO
END SUBROUTINE NPK_Supply

! ----------------------------- 5 Overall Stress -----------------------------
! -----------------------------5.1 N fix ---------------------------
SUBROUTINE N_Fix 
      !     EPIC1102
      !     THIS SUBPROGRAM ESTIMATES N FIXATION FOR LEGUMES.
      ! Ref: APEX-doc  Eqs 288 (a-i)
      ! Ref: SWAT-doc  Eqs 5:2.3.9
      IMPLICIT NONE
 
      ! local variables
      INTEGER:: J, L1
      REAL:: SUM, TOT, ADD, RTO, X1, FXN, FXW, FXG, FXS, FXP, FIXR

      IF(HUI(Crop_Num)<.15.OR.HUI(Crop_Num)>.75)GO TO 8
      SUM=0.
      TOT=0.
      ADD=0.
      DO J=1,Actual_SoilLayers
            ISL=Layer_ID(J)
            IF(Z(ISL)>.3)GO TO 3
            SUM=SUM+soilWater(ISL)-Wilt_Point(ISL)
            TOT=TOT+fieldCapacity(ISL)-Wilt_Point(ISL)
      END DO
      GO TO 4
    3 L1=Layer_ID(J-1)
      RTO=(.3-Z(L1))/(Z(ISL)-Z(L1))
      SUM=SUM+(soilWater(ISL)-Wilt_Point(ISL))*RTO
      TOT=TOT+(fieldCapacity(ISL)-Wilt_Point(ISL))*RTO
    4 X1=SUM/TOT
      IF(X1<=.25)GO TO 8
      DO J=1,Actual_SoilLayers
            ISL=Layer_ID(J)
            IF(Z(ISL)>RD(Crop_Num))GO TO 6
            ADD=ADD+NO3_N_Soil(ISL)   ! ADD is TNO3 in Eq 288h  APEX-doc
      END DO
      GO TO 7
    6 L1=Layer_ID(J-1)
      RTO=(RD(Crop_Num)-Z(L1))/(Z(ISL)-Z(L1))
      ADD=ADD+NO3_N_Soil(ISL)*RTO
    7 FXN=1.5-.005*ADD/RD(Crop_Num)
      IF(FXN>0.)THEN
            FXW=1.333*X1-.333
            FXG=(HUI(Crop_Num)-.1)*5.
            FXS=4.-5.*HUI(Crop_Num)
            FXP=MIN(FXG,FXS,1.)      ! FXP: growth stage (0-1.)
            FIXR=MIN(FXW,FXN,1.)*FXP ! FXW: soil water factor (0-1.)  FXN: soil nitrate factor (0-1.)
            FixedN_Final=FIXR*potentialDemandN
      END IF
    8 FixedN_Final=MAX(0.,modelPara(7)*FixedN_Final+(1.-modelPara(7))*potentialDemandN)  !  Ref: Eq 288 in APEX-doc   kg/ha/d
      IF(FixedN_Final>modelPara(68))FixedN_Final=modelPara(68)
      SMM(50,MO)=SMM(50,MO)+FixedN_Final
      dailyFixedN=dailyFixedN+FixedN_Final
      potentialDemandN=potentialDemandN-FixedN_Final
END SUBROUTINE N_Fix

! ----------------------------- 5.2 Actual N uptake ------------------
SUBROUTINE N_Uptake(UU,XNO3,XNH3,DMD,AJF,SUPL)
      !     EPIC1102
      !     THIS SUBPROGRAM COMPUTES ACTUAL N PLANT UPTAKE FROM EACH
      !     LAYER (UPTAKE = MINIMUM OF PLANT DEMAND AND SOIL SUPPLY)
      IMPLICIT NONE
 
      ! local variables potential problems with local variables (SAVE attribute UU)
      INTEGER:: J, K
      REAL:: UU(15), XNO3(15), XNH3(15), XX, X1, X2, X21, X3, DMD, SUPL, AD1, AD2, AJF, RTO
 
      X2=DMD/(SUPL+1.E-20)
      AD1=0.
      IF(X2<1..AND.DMD>0.)THEN
            DO J=1,Layer_RD_Reach
                  K=Layer_ID(J)
                  UU(K)=UU(K)*X2
                  AD1=AD1+UU(K)
            END DO
            SUPL=AD1
      ELSE
            AD1=0.
            X2=AJF*(DMD-SUPL)
            X21=X2
            DO J=1,Layer_RD_Reach
                  K=Layer_ID(J)
                  XX=UU(K)+X2
                  X3=.001*modelPara(27)*WT(K)
                  X1=modelPara(91)*(XNO3(K)-X3)+(1.-modelPara(91))*XNH3(K)
                  IF(XX<X1)EXIT
                  IF(X1>0.)THEN
                        X2=X2-X1+UU(K)
                        UU(K)=X1
                        AD1=AD1+UU(K)
                  ELSE
                        UU(K)=0.
                  END IF
            END DO
            IF(J>Layer_RD_Reach)THEN
                  SUPL=AD1
            ELSE
                  UU(K)=XX
                  SUPL=SUPL+X21
            END IF
      END IF
      AD1=0.
      AD2=0.
      DO J=1,Layer_RD_Reach
            K=Layer_ID(J)
            X3=.001*modelPara(27)*WT(K)
            XX=(1.-modelPara(91))*XNH3(K)
            RTO=XX/(modelPara(91)*(XNO3(K)-X3)+XX)
            X1=RTO*UU(K)
            AD1=AD1+X1
            X2=UU(K)-X1
            AD2=AD2+X2
            XNO3(K)=XNO3(K)-X2
            XNH3(K)=XNH3(K)-X1
      END DO
      SMM(108,MO)=SMM(108,MO)+AD2
      VAR(108)=AD2
      SMM(109,MO)=SMM(109,MO)+AD1
      VAR(109)=AD1    
END SUBROUTINE N_Uptake

! ------------------------------ 5.3 ACTUAL P uptake ------------------
SUBROUTINE P_Uptake(UU,AN,DMD,AJF,SUPL)
      ! THIS SUBPROGRAM COMPUTES ACTUAL P PLANT UPTAKE FROM EACH
      ! LAYER (UPTAKE = MINIMUM OF PLANT DEMAND AND SOIL SUPPLY).
      IMPLICIT NONE
 
      ! local variables potential problems with local variables (SAVE attribute UU)
      INTEGER:: J, K
      REAL:: UU(15),AN(15), SUM, X1, X2, X21, XX, AJF, DMD, SUPL 
      
      SUM=0.
      X2=AJF*(DMD-SUPL)
      X21=X2
      DO J=1,Layer_RD_Reach
            K=Layer_ID(J)
            XX=UU(K)+X2
            X1=AN(K)-.001*modelPara(27)*WT(K)
            IF(XX<X1)GO TO 6
            IF(X1>0.)THEN
                  X2=X2-X1+UU(K)
                  UU(K)=X1
                  SUM=SUM+UU(K)
            ELSE
                  UU(K)=0.
            END IF
      END DO
      SUPL=SUM
      RETURN
    6 UU(K)=XX
      SUPL=SUPL+X21
END SUBROUTINE P_Uptake
! ------------------------------ 5.4 N or P Stress ------------------
REAL FUNCTION StressFactor(U1,U2)
      ! EPIC1102
      ! THIS SUBPROGRAM CALCULATES THE PLANT STRESS FACTOR CAUSED BY LIMITED SUPPLY OF N OR P.
      !  Ref: SWAT-doc 2009  Eq. 5.3.1.6
      IMPLICIT NONE
      ! local variables
      REAL:: U1, U2      ! Kg (N/P)/ha
 
      StressFactor=200.*(U1/U2-.5)   
      IF(StressFactor>0.)THEN
            StressFactor=StressFactor/(StressFactor+EXP(S_Curve(8,1)-S_Curve(8,2)*StressFactor))  ! Ref: APEX-doc Eq 296  but the coefficients are very different! Watch out
      ELSE
            StressFactor=0.
      END IF
END FUNCTION StressFactor

SUBROUTINE Overall_Stress(JRT)
      ! EPIC1102
      ! THIS SUBPROGRAM ESTIMATES PLANT STRESS FACTORS CAUSED BY LIMITED
      ! N, P, AIR, AND WATER AND DETERMINES THE ACTIVE CONSTRAINT
      ! (MINIMUM STRESS FACTOR--N, P, WATER, OR TEMPERATURE).  CALLS
      ! NFIX AND Apply_Fert (AUTOMATIC FERTILIZER OPTION).
      IMPLICIT NONE
 
      ! local variables
      INTEGER:: I, J1, J3, N1, JRT
      INTEGER, DIMENSION(7):: JFS = [(i, i = 13, 19)]
      REAL:: X1, X2, X3, XX, ZZ   

      J3=5
      cropVar(17,Crop_Num)=Min_Stress_Factor(Crop_Num)
      ! Water stress
      IF(EP(Crop_Num)>0.)THEN
            IF(soilWaterRZ>0.)THEN
                  waterStress=100.*soilWaterRZ/potentialAvailWater
                  waterStress=waterStress/(waterStress+EXP(S_Curve(11,1)-S_Curve(11,2)*waterStress))
            ELSE
                  waterStress=0.
            END IF
            waterStress=(1.-modelPara(35))*waterStress+modelPara(35)*plantEP/(EP(Crop_Num)+1.E-10)
      END IF
      cropVar(13,Crop_Num)=waterStress
      FixedN_Final=0.
      IF(cropCode(Crop_Num)==plantCategoryCode(1).OR.cropCode(Crop_Num)==plantCategoryCode(2).OR.&
      cropCode(Crop_Num)==plantCategoryCode(3).OR.cropCode(Crop_Num)==plantCategoryCode(10))&
      CALL N_Fix 
      
      CALL N_Uptake(soilSupplyRateN,NO3_N_Soil,NH3_Weight,potentialDemandN,1.,SUN)
					
      X1=SUN/(potentialDemandN+1.E-10)
      potentialDemandN=SUN+FixedN_Final
      actualCropN(Crop_Num)=actualCropN(Crop_Num)+potentialDemandN
      IF(potentialDemandP>SUP)CALL P_Uptake(soilSupplyRateP,labileP,potentialDemandP,1.,SUP)
      
      X2=MIN(1.,SUP/(potentialDemandP+1.E-10))
      potentialDemandP=SUP
      actualCropP(Crop_Num)=actualCropP(Crop_Num)+potentialDemandP
      !CALL NAJN(soilSupplyRateK,solubleK,PlantK_Demand,SUK,1.,0)
      X3=SUK/(PlantK_Demand+1.E-10)
      PlantK_Demand=SUK
      !actualCropK(Crop_Num)=actualCropK(Crop_Num)+PlantK_Demand
      N_StressFactor = StressFactor(actualCropN(Crop_Num),Optimal_N(Crop_Num)) ! UN2: optimal N content for crop in kg/ha
      N_StressFactor=MAX(X1,N_StressFactor)
      cropVar(14,Crop_Num)=N_StressFactor
      P_StressFactor = StressFactor(actualCropP(Crop_Num),Optimal_P(Crop_Num)) ! UP2: optimal P content for crop in kg/ha
      P_StressFactor=MAX(X2,P_StressFactor)
      cropVar(15,Crop_Num)=P_StressFactor
      ! CALL NUTS(actualCropK(Crop_Num),UK2(Crop_Num),SK)
      ! SK=MAX(X3,SK)
      ! cropVar(16,Crop_Num)=SK
      cropVar(18,Crop_Num)=SAT
      X1=.15625*Tot_Salt_RootZone/Tot_Water_RootZone
      XX=X1-salinityThreshold(2,Crop_Num)
      IF(XX>0.)THEN
            SSLT=MAX(0.,1.-salinityThreshold(1,Crop_Num)*XX)
      ELSE
            SSLT=1.
      END IF
      cropVar(19,Crop_Num)=SSLT
      DO
            CALL Determine_MinStress(6,J3,SAT,Min_Stress_Factor(Crop_Num),0.,JRT)
            IF(JRT>0)EXIT
            IF(biomsConver_Method==0)CALL Determine_MinStress(1,J3,waterStress,Min_Stress_Factor(Crop_Num),0.,JRT)
            IF(JRT>0)EXIT
            CALL Determine_MinStress(3,J3,P_StressFactor,Min_Stress_Factor(Crop_Num),0.,JRT)
            IF(JRT>0)EXIT
            ! CALL CFRG(4,J3,SK,Min_Stress_Factor(Crop_Num),0.,JRT)
            ! IF(JRT>0)EXIT
            CALL Determine_MinStress(7,J3,SSLT,Min_Stress_Factor(Crop_Num),0.,JRT)
            IF(JRT>0)EXIT
            ZZ=Min_Stress_Factor(Crop_Num)
            CALL Determine_MinStress(2,J3,N_StressFactor,ZZ,Min_Stress_Factor(Crop_Num),JRT)
            IF(JRT==0)EXIT
            Min_Stress_Factor(Crop_Num)=ZZ
            IF(cropCode(Crop_Num)==plantCategoryCode(1).OR.cropCode(Crop_Num)==plantCategoryCode(2).OR. &
                  cropCode(Crop_Num)==plantCategoryCode(3).OR.cropCode(Crop_Num)==plantCategoryCode(10))EXIT
            IF(NStressTrigger>N_StressFactor.AND.NFA>=fertApplyMinInterval) CALL Apply_Fert(3, IAUF)
            EXIT          ! BFT: N stress factor to trigger automatic fertilizer
      END DO    
      XX=1.-Min_Stress_Factor(Crop_Num)
      N1=NCP(Crop_Num)
      SFMO(J3,Crop_Num)=SFMO(J3,Crop_Num)+XX
      Growing_Stress(J3,N1,Crop_Num)=Growing_Stress(J3,N1,Crop_Num)+XX
      growStressFactor(J3,Crop_Num)=XX
      BLYN(Crop_Num)=BLYN(Crop_Num)+1.
      J1=JFS(J3)
      SMMC(J1,Crop_Num,MO)=SMMC(J1,Crop_Num,MO)+XX
      cropVar(J1,Crop_Num)=Min_Stress_Factor(Crop_Num)
      SMMC(20,Crop_Num,MO)=SMMC(20,Crop_Num,MO)+XX
      cropVar(20,Crop_Num)=Min_Stress_Factor(Crop_Num)
END SUBROUTINE Overall_Stress

! ---------------------------------------- 6. RUE Accumulation ----------------------------------------------
SUBROUTINE Actual_Growth 
      ! EPIC1102
      ! THIS SUBPROGRAM CALCULATES THE DAILY INCREASE IN PLANT BIOMASS,
      ! ROOT WEIGHT, AND YIELD BY ADJUSTING THE POTENTIAL VALUES WITH THE
      ! ACTIVE STRESS CONSTRAINT.
      ! Ref: SWAT-doc 2009 Eqs 5:2.4.1
      IMPLICIT NONE
      ! local variables   
      INTEGER:: I, J
      REAL:: XTP(15), ADD, X1, X2, YX, XX, XY, XZ, XW, ZZ, RWL, RGD, RF, F, FF, Leaf_Fall, DRW, SUM, &
       W1, DMD, UTO, SPL
 
      XX=Min_Stress_Factor(Crop_Num)*SHRL
      RWL=totRootWeight(Crop_Num)
      RGD=potentialBioIncrese(Crop_Num)*XX
      DRWX=0.
      X1=100.*HUI(Crop_Num)
      AJHI(Crop_Num)=HI(Crop_Num)*X1/(X1+EXP(S_Curve(3,1)-S_Curve(3,2)*X1))
      XX=100.*sumActuralEP(Crop_Num)/(sumPotentialEP(Crop_Num)+1.E-5) ! sumActuralEP: sum of actual EP on a given day; sumPotentialEP: sum of potential EP 
      F=XX/(XX+EXP(S_Curve(10,1)-S_Curve(10,2)*XX))      ! SWAT-doc 2009  Eq 5: 3.3.1
      XX=MAX(AJHI(Crop_Num)-lowerLimitHI(Crop_Num),0.)
      HIX=F*XX+lowerLimitHI(Crop_Num)                        ! actual harvest index
      YLDX=HIX*abvGroundBiom(Crop_Num)             ! crop yield   kg/ha  
      YX=totCropBio(Crop_Num)-potentialBioIncrese(Crop_Num)
      XX=MAX(1.E-5,YX+RGD)
      X2=RGD*420.
      SM99=SM99+X2
      VAR(99)=VAR(99)+X2
      SMM(99,MO)=SMM(99,MO)+X2
      RF=rootPartition(1,Crop_Num)*(1.-HUI(Crop_Num))+rootPartition(2,Crop_Num)*HUI(Crop_Num)
      IF(RF<.2)THEN
            RF=.2
      ELSE
            RF=MIN(RF,.99)
      END IF
      totRootWeight(Crop_Num)=modelPara(92)*RWL+(1.-modelPara(92))*RF*XX        
      DRW=totRootWeight(Crop_Num)-RWL
      totCropBio(Crop_Num)=XX
      totCropBioX(Crop_Num)=MAX(totCropBioX(Crop_Num),XX)              
      abvGroundBiom(Crop_Num)=MAX(0.,totCropBio(Crop_Num)-totRootWeight(Crop_Num))
      XX=0.
      SUM=0.
      X2=2.*RZ
      DO I=1,Layer_RD_Reach
            ISL=Layer_ID(I)
            XTP(ISL)=(Z(ISL)-XX)*EXP(-modelPara(56)*Z(ISL)/RZ)
            SUM=SUM+XTP(ISL)
            XX=Z(ISL)
      END DO
      IF(cropCode(Crop_Num)==plantCategoryCode(7).OR.cropCode(Crop_Num)==plantCategoryCode(8).OR. &
           cropCode(Crop_Num)==plantCategoryCode(10))THEN
            X1=100.*HSM/yearTotHeat    ! for trees
            IF(X1>75.)THEN
                  AboveGround_Biomass0(Crop_Num)=MAX(AboveGround_Biomass0(Crop_Num),abvGroundBiom(Crop_Num))
                  F=X1/(X1+EXP(S_Curve(29,1)-S_Curve(29,2)*X1))
                  FF=F-FLF0(Crop_Num)
                  Leaf_Fall=AboveGround_Biomass0(Crop_Num)*turnoutFracCottn(Crop_Num)*FF
                  !ADFL=ADFL+Leaf_Fall
                  FLF0(Crop_Num)=F
                  SMM(88,MO)=SMM(88,MO)+Leaf_Fall
	            !Current_LAI(Crop_Num)=MAX(.001,Current_LAI(Crop_Num)-Leaf_Fall)
                  totCropBio(Crop_Num)=totCropBio(Crop_Num)-Leaf_Fall
                  abvGroundBiom(Crop_Num)=abvGroundBiom(Crop_Num)-Leaf_Fall
                  X1=Leaf_Fall*1000.
                  XZ=fracYieldN(Crop_Num)*X1
                  XY=fracYieldP(Crop_Num)*X1
                  W1=fracYieldK(Crop_Num)*X1
                  CALL CN_TopSoil(Leaf_Fall,XZ,0 )
                  FreshOrgP_Residu(LD1)=FreshOrgP_Residu(LD1)+XY
                  actualCropN(Crop_Num)=MAX(1.E-5,actualCropN(Crop_Num)-XZ)
                  actualCropP(Crop_Num)=actualCropP(Crop_Num)-XY
                  !WRITE(KW(1),'(1X,A,3I4,5F10.4)')'~~~~~',IYR,MO,DayOfMon,F,FF,AboveGround_Biomass0(Crop_Num),Leaf_Fall,ADFL
            END IF
      END IF
      IF(cropCode(Crop_Num)==plantCategoryCode(3).OR.cropCode(Crop_Num)==plantCategoryCode(6))THEN  ! prennial
            ZZ=.01*(HUI(Crop_Num)+.01)**10*abvGroundBiom(Crop_Num)
            abvGroundBiom(Crop_Num)=abvGroundBiom(Crop_Num)-ZZ
            totCropBio(Crop_Num)=totCropBio(Crop_Num)-ZZ
            standCropResi(Crop_Num)=standCropResi(Crop_Num)+ZZ
            standDeadLignin=standDeadLignin+CLG*ZZ
            XY=ZZ*uptakeParaN(3,Crop_Num)
            XZ=ZZ*uptakeParaP(3,Crop_Num)
            XW=ZZ*uptakeParaK(3,Crop_Num)
            standDeadResiK=standDeadResiK+XW
            standDeadResiN(Crop_Num)=standDeadResiN(Crop_Num)+XY
            standDeadResiP=standDeadResiP+XZ
            !actualCropK(Crop_Num)=actualCropK(Crop_Num)-XW
            actualCropN(Crop_Num)=MAX(1.E-5,actualCropN(Crop_Num)-XY)
            actualCropP(Crop_Num)=actualCropP(Crop_Num)-XZ
            IF(HUI(Crop_Num)>.6.AND.abvGroundBiom(Crop_Num)<.1)HU(Crop_Num)=0.
      END IF
      ADD=0.
      DMD=0.
      DO J=1,Layer_RD_Reach
            ISL=Layer_ID(J)
            !NO3_N_Soil(ISL)=MAX(1.E-5,NO3_N_Soil(ISL)-soilSupplyRateN(ISL))
            IF(DRW<0.)THEN
                  UTO=RWT(ISL,Crop_Num)/(RWL+1.E-10)
                  DMD=-DRW*UTO+DMD
                  SPL=RWT(ISL,Crop_Num)
                  IF(DMD>=SPL)THEN
                        X1=SPL
                        DMD=DMD-SPL
                  ELSE
                        X1=DMD
                        DMD=0.                      
                  END IF
                  X2=MIN(actualCropN(Crop_Num),1000.*uptakeParaN(3,Crop_Num)*X1)
                  actualCropN(Crop_Num)=actualCropN(Crop_Num)-X2
                  X1=MIN(X1,totCropBio(Crop_Num))
                  CALL CN_TopSoil(X1, X2, 1 )
                  totCropBio(Crop_Num)=totCropBio(Crop_Num)-X1
                  X1=-X1
            ELSE
                  UTO=modelPara(55)*wiltPoint(ISL)/(plantEP+1.E-20)+(1.-modelPara(55))*XTP(ISL)/SUM
                  X1=DRW*UTO
            END IF
            soilWater(ISL)=soilWater(ISL)-wiltPoint(ISL)
            labileP(ISL)=labileP(ISL)-soilSupplyRateP(ISL)
            solubleK(ISL)=solubleK(ISL)-soilSupplyRateK(ISL)
            DRWX(ISL)=X1
            RWT(ISL,Crop_Num)=RWT(ISL,Crop_Num)+DRWX(ISL)
            ADD=ADD+RWT(ISL,Crop_Num)
      END DO
      totRootWeight(Crop_Num)=ADD
      abvGroundBiom(Crop_Num)=totCropBio(Crop_Num)-ADD
      IF(totRootWeight(Crop_Num)<=RWX(Crop_Num))RETURN
      RWX(Crop_Num)=totRootWeight(Crop_Num)
      M21=MO
      K21=DayOfMon
      DO I=1,Layer_RD_Reach
            ISL=Layer_ID(I)
            RWTX(ISL,Crop_Num)=RWT(ISL,Crop_Num)
      END DO
END SUBROUTINE Actual_Growth

END MODULE Crop_Module