MODULE ReadVarTable_Module
USE PARM
USE MisceFun_Module 
IMPLICIT NONE
CONTAINS
    
! -------------------1. CN Table ---------------------------------
SUBROUTINE CN2_Table
      IMPLICIT NONE
      !  EPIC1102
      !     THIS SUBPROGRAM CONTAINS THE SCS HYDROLOGIC SOIL GROUP-CURVE
      !     NUMBER TABLE.    
      REAL, DIMENSION(35,4):: CNX(35,4) = RESHAPE([77.,72.,67.,70.,65.,66.,62.,65.,63.,63.,61.,61.,59.,66.,&
      58.,64.,55.,63.,51.,68.,49.,39.,47.,25.,6.,30.,45.,36.,25.,59.,72.,&
      74.,43.,39.,98.,86.,81.,78.,79.,75.,74.,71.,76.,75.,74.,73.,72.,&
      70.,77.,72.,75.,69.,73.,67.,79.,69.,61.,67.,59.,35.,58.,66.,60.,&
      55.,74.,82.,84.,65.,61.,98.,91.,88.,85.,84.,82.,80.,78.,84.,83.,&
      82.,81.,79.,78.,85.,81.,83.,78.,80.,76.,86.,79.,74.,81.,75.,70.,&
      71.,77.,73.,70.,82.,87.,90.,77.,74.,98.,94.,91.,89.,88.,86.,82.,&
      81.,88.,87.,85.,84.,82.,81.,89.,85.,85.,83.,83.,80.,89.,84.,80.,&
      88.,83.,79.,78.,83.,79.,77.,86.,89.,92.,82.,80.,98.], [35,4])
      CN2=CNX(landUseCode,hydroGroup)
      RETURN
END SUBROUTINE
    
! ---------------------- 2. Fertalizer Table --------------------------
SUBROUTINE Fert_Table(L)
      IMPLICIT NONE
      !local variables potential problems with local variables (SAVE attribute XTP)
      REAL:: XTP(11) 
      INTEGER:: I, L, K
     
      IF(NDF>0)THEN
            DO L=1,NDF
                  IF(KDF(L)==JX(7))RETURN
            END DO
      END IF
      NDF=NDF+1
      KDF(NDF)=JX(7)
      KDF1(JX(7))=NDF
      READ(KR(9),'()')
      READ(KR(9),'()')
      DO
            READ(KR(9),395,IOSTAT=NFL)I,fertName(NDF),(XTP(K),K=2,11)
            IF(NFL/=0)THEN
                  IF(runMode==0)THEN
                        WRITE(*,*)'FERT NO = ',JX(7),' NOT IN FERT LIST FILE'
                        ERROR STOP
                  ELSE
                        WRITE(KW(MSO),'(A,A8,A,I4,A)')'!!!!! ',siteName,'FERT NO = ',&
                              JX(7),' NOT IN FERT LIST FILE'
                  END IF
                  STOP
            END IF
            IF(I==JX(7))EXIT
      END DO
      FN(NDF)=XTP(2)
      FP(NDF)=XTP(3)
      FK(NDF)=XTP(4)
      FNO(NDF)=XTP(5)
      FPO(NDF)=XTP(6)
      FNH3(NDF)=XTP(7)
      FOC(NDF)=XTP(8)
      FSLT(NDF)=XTP(9)
      FCST(NDF)=XTP(10)
      FCEM(NDF)=XTP(11)
      REWIND KR(9)
 
  395 FORMAT(1X,I4,1X,A8,11F8.0)
END SUBROUTINE Fert_Table
! ------------------------------ 3. Pesticide Table ---------------------------
SUBROUTINE Pest_Table
    
      IMPLICIT NONE
      ! arguments list
 
      ! local variables 
      REAL:: XTP(8)
      INTEGER:: L, J1
      ! ////////////////////////////////////////////////////////////

      IF(NDP>0)THEN
            DO L=1,NDP
                  IF(KDP(L)==JX(7))RETURN
            END DO
      END IF
      NDP=NDP+1
      KDP(NDP)=JX(7)
      KDP1(JX(7))=NDP
 
      J1=-1
      DO WHILE(J1/=JX(7))
            READ(KR(8),1,IOSTAT=NFL)J1,pestName(NDP),(XTP(L),L=2,8)
            IF(NFL/=0)THEN
                  IF(runMode==0)THEN
                        WRITE(*,*)'PEST NO = ',JX(7),' NOT IN PEST LIST FILE'
                        ERROR STOP
                  ELSE
                        WRITE(KW(MSO),'(A,A8,A,I4,A)')'!!!!! ',siteName,' PEST NO = ',&
                              JX(7),' NOT IN PEST LIST FILE'
                  END IF
                  STOP
            END IF
      END DO
      pestSolubity(NDP)=XTP(2)
      pestHfLifeSoil(NDP)=XTP(3)
      pestHfLifePlant(NDP)=XTP(4)
      WashOff_Frac(NDP)=XTP(5)
      Coeff_OrgCarbon(NDP)=XTP(6)
      Pest_Cost(NDP)=XTP(7)
      Pest_C_emission(NDP)=XTP(8)
      REWIND KR(8)
 
    1 FORMAT(I5,1X,A16,7E16.6)    
END SUBROUTINE Pest_Table

! ----------------------- 4. Tillage Table --------------------------------
SUBROUTINE Tillage_Table
    ! EPIC1102
    ! THIS SUBPROGRAM READS EQUIPMENT TABLE TO DETERMINE PARAMETERS OF INPUT OPERATIONS 
    ! AND COMPUTES OPERATION COSTS(EQUIPMENT + TRACTOR)
      IMPLICIT NONE
      ! argument list
 
      ! local variables
      CHARACTER(LEN =4), DIMENSION(5):: PCAT(5)= ['POWE','NON ','SELF','IRRI','CUST']
      CHARACTER(LEN = 4):: PCD
      CHARACTER(LEN = 8):: TILX
      INTEGER:: L, NN, JJ, J2, J, I, K
      REAL:: XTP(13), FLAB, XLP, HRY, PRIC, HRL, PWR, WDTE, WDT, SPD, SPDE, RC1, &
            RC2, XLB, FCM, RFV1, RFV2, X1, RTI, YR, AMF, SALV, CSTA, CSTI, CSTO, &
            CSTW, CAPM, TOCS, TAR, CSFU 
       ! management practice database 
 
      IF(NDT>0)THEN
            DO L=1,MXT                   ! initialized in AINLZ, what is it?
                  IF(NBE(L)/=JX(4).OR.NBT(L)/=JX(5))CYCLE
                  NDT=L
                  RETURN
            END DO
      END IF
      NDT=MXT+1
      MXT=NDT
      NBE(NDT)=JX(4)  ! ???
      NBT(NDT)=JX(5)  ! Tractor ID
      NN=JX(4)        ! Equipment ID
      COTL(NDT)=0.    ! ???
      COOP(NDT)=0.    ! ???
      FLAB=1.        
      JJ=JX(4)        ! 
      DO J=1,2                     
            READ(KR(3),25)TILX              ! Read tillage table
            READ(KR(3),25)TILX
            ! READ EQUIPMENT DATA TABLE
            ! 17  RR   = RANDOM SURFACE ROUGHNESS CREATED BY TILLAGE OPERATION (MM)
            J2=-1
            searchFile: DO WHILE(J2/=JJ)
                              READ(KR(3),18,IOSTAT=NFL)J2,TILX,PCD, &  ! PCD  = POWER CODE
                                    PRIC, &       !  1  PRIC = PURCHASE PRICE($)--EXCEPTION CUSTOM = COST($/ha)
                                    XLP,  &       !  2  XLP  = INITIAL LIST PRICE IN CURRENT($)
                                    HRY,  &       !  3  HRY  = ANNUAL USE(H)
                                    HRL,  &       !  4  HRL  = LIFE OF EQUIP(h)
                                    PWR,  &       !  5  PWR  = POWER OF UNIT(kw)
                                    WDT,  &
                                    SPD,  &       
                                    RC1,  &       !  8  RC1  = REPAIR COST COEF 1
                                    RC2,  &       !  9  RC2  = REPAIR COST COEF 2
                                    XLB,  &       ! 10  XLB  = LUBRICANT FACTOR
                                    FCM,  &       ! 11  FCM  = FUEL CONSUMPTION MULTIPLIER
                                    RFV1, &       ! 12  RFV1 = REMAINING FARM VALUE PARM 1
                                    RFV2, &       ! 13  RFV2 = REMAINING FARM VALUE PARM 2
                                    X1,   &
                                    RTI,  &       ! 15  RTI  = ANNUAL REAL INTEREST RATE($/$)
                                    (XTP(L),L=1,13)
                              IF(NFL/=0)THEN
                                    IF(runMode==0)THEN
                                          WRITE(*,*)'TILLAGE NO = ',JJ,' NOT IN TILL LIST FILE'
                                                ERROR STOP
                                    ELSE
                                          WRITE(KW(MSO),'(A,A8,A,I4,A)')' !!!!! ',siteName,&
                                                ' TILLAGE NO = ',JJ,' NOT IN TILL LIST FILE'
                                                ERROR STOP  ! Be careful here, it is different from the original code.
                                    END IF 
                              END IF
            END DO  searchFile
      
            IF(J==1)THEN
                  WDTE=WDT                    !  6  WDT  = WIDTH OF PASS(m)
                  SPDE=SPD                    !  7  SPD  = OPERATING SPEED(km/h)
                  EFM(NDT)=X1                 ! 14  EFM  = MACHINE harvestEffi = GRAZING RATE (KG/HD/D)currentOps=19
                  equipmentName(NDT)=TILX    !  equipmentName  = EQUIPMENT NAME
                  mixEfficiency(NDT)=XTP(1)        ! 16  mixEfficiency  = MIXING harvestEffi (0-1)     
                  tillRoughness(NDT)=XTP(2)
                  tillDepth(NDT)=XTP(3)*.001 ! 18  tillDepth  = TILLAGE DEPTH(MM)
                  ridgeHeight(NDT)=XTP(4)    ! 19  RHT  = RIDGE HEIGHT (MM)
                  ridgeInterval(NDT)=XTP(5)  ! 20  RIN  = RIDGE INTERVAL (m)
                  furrowHeight(NDT)=XTP(6)   ! 21  DKH  = HEIGHT OF FURROW DIKES (MM) (BLANK IF DIKES NOT USED)
                  furrowInterval(NDT)=XTP(7) ! 22  DKI  = distance BETWEEN FURROW DIKES (m)(BLANK IF DIKES NOT USED
                  currentOps(NDT)=INT(XTP(8))
                  harvestEffi(NDT)=XTP(9)      ! 24  harvestEffi  = HARVEST harvestEffi(0-1)= PESTICIDE APPLICATION harvestEffi
                  overRideHI(NDT)=XTP(10)    ! 25  overRideHI = OVER RIDES SIMULATED HI IF 0.<overRideHI<1.
                  soilCompactFrac(NDT)=XTP(11)! 26  soilCompactFrac = FRACTION OF SOIL COMPACTED(TIRE WIDTH/TILLAGE WIDTH)
                  PlantRedu_frac(NDT)=XTP(12) ! 27  PlantRedu_frac = FRACTION PLANT POPULATION REDUCED BY OPERATION
                  carbonEmiss(NDT)=XTP(13)   ! 28  carbonEmiss = CARBON EMISSION(kg/ha)
          
                  IF(currentOps(NDT)==operationCode(19))overRideHI(NDT)=.001*X1 ! X1: machine harvestEffi
                  JJ=JX(5)
                  IF(currentOps(NDT)==operationCode(17))THEN      ! operationCode(17): build furrow dikes
                        RidgeHeight2=ridgeHeight(NDT)
                        Ridge_IntervalX=ridgeInterval(NDT)
                  END IF
                  IF(currentOps(NDT)==operationCode(22)) HMO(NDT)=-tillDepth(NDT) ! operationCode(22): auto mow
                  IF(currentOps(NDT)==operationCode(23))THEN      ! operationCode(23): plastic cover
                        PALB=tillRoughness(NDT)
                        FCV=EFM(NDT)
                  END IF
            END IF
            
            IF(HRY<1.E-10)THEN
                  IF(PCD==PCAT(5))THEN
                        COTL(NDT)=PRIC-XLP
                        COOP(NDT)=PRIC-XLP
                        ICUS(NDT)=1
                  END IF
                  REWIND KR(3)
                  RETURN
            END IF
            ! economic analysis
            YR=MIN(30.,HRL/HRY)               ! The years of equiment can be used for.
            IF(numSimuYear>1.AND.DOY_realtime==0.AND.RTI>0.)THEN
                  X1=(1.+RTI)**YR
                  AMF=RTI*X1/(X1-1.+1.E-10)
            ELSE
            AMF=1./YR
            END IF
            SALV=RFV1*XLP*RFV2**YR
            CSTA=AMF*(PRIC-SALV)/HRY
            CSTI=RTI*SALV/HRY
            CSTO=CSTA+CSTI                 ! total cost?
            TAR=XLP*RC1*(.001*HRL)**RC2/(HRL+.0001)
            Fuel_Use(NDT)=FCM*PWR          ! FCM: fuel consumption multiplier, PWR: power of unit    
            CSFU=Fuel_Use(NDT)*fuelCost
            DO I=1,4                   ! Determine the power type 
                  IF(PCD==PCAT(I))EXIT
            END DO
            K=0
            SELECT CASE(I)       !"NON": the machine has no engine for power and it must be pulled by other machinery with engine power
                  CASE(1)
                        FLAB=1.1
                  CASE(2)
                        IF(NBT(NDT)==0)THEN
                              COTL(NDT)=CSTO+TAR
                              COOP(NDT)=TAR
                              K=1
                        ELSE
                              FLAB=.08
                        END IF
                  CASE(3)
                        FLAB=1.2
                  CASE(4)
                        COTL(NDT)=CSTO+TAR
                        COOP(NDT)=TAR
                        K=1 
            END SELECT
            IF(K==0)THEN
                  CSTW=laborCost*FLAB
                  TOCS=TAR+(1.+XLB)*CSFU+CSTW
                  CAPM=10./(SPDE*WDTE*EFM(NDT))
                  COOP(NDT)=TOCS*CAPM+COOP(NDT)
                  TOCS=TOCS+CSTO
                  COTL(NDT)=TOCS*CAPM+COTL(NDT)
                  Fuel_Use(NDT)=CAPM*Fuel_Use(NDT)
            END IF
            NN=JX(5)              ! Tractor ID number
            REWIND KR(3)          ! =============== the end of reading Tillcom file  
            IF(NN==0)RETURN
      END DO
   18 FORMAT(1X,I4,1X,A8,1X,A4,28F8.0)
   25 FORMAT(A8)
 
END SUBROUTINE Tillage_Table

! ----------------------- 5. Crop Table -----------------------------------
! ------------------ 1. Crop Table ------------------------------- 
SUBROUTINE Read_CropTable
      !     EPIC1102
      !     THIS SUBPROGRAM READS CROP TABLE TO DETERMINE CROP PARAMETERS      
      ! arguments list
 
      ! local variables  potential problems with local variables (SAVE attribute of XTP)
      CHARACTER(LEN = 12)::file_name
      REAL:: XTP(58), Y1, Y2, X1, X2, X3, X4, X5, G1, Z1, Z2, FU, &
             DFDU    
      INTEGER:: L, J2, IT
      ! ///////////////////////////////////////////////////////////////////
      
      IF(DOY_realtime>0)GO TO 174
      IF(LC==0)GO TO 174
      DO L=1,LC
            IF(cropID(L)==JX(6))GO TO 179
      END DO
  174 LC=LC+1
      varID(JX(6))=LC        ! How does  varID become 0?  varID is just be used here. 
      cropID(LC)=JX(6)                ! What is cropID?
      XMTU(LC)=JX(7)                      ! for trees, no PHU, but calcuated from XMTU
      ! READ CROP DATA
      READ(KR(4),1)file_name
      READ(KR(4),1)file_name
      J2=-1
      DO WHILE(J2/=JX(6))
            READ(KR(4),333,IOSTAT=NFL)J2,Crop_Name(LC),(XTP(L),L=1,58)  ! Crop_Name = NAMES OF CROPS IN CROP PARAMETER TABLE
            IF(NFL/=0)THEN
                  IF(runMode==0)THEN
                        WRITE(*,*)'CROP NO = ',JX(6),' NOT IN CROP LIST FILE'
                        ERROR STOP
                  ELSE
                        WRITE(KW(MSO),'(A,A8,A,I4,A)')' !!!!! ',siteName,&
	                        &'CROP NO = ',JX(6),' NOT IN CROP LIST FILE'
                        ERROR STOP  ! different from original code
                  END IF
              	              
            END IF       
      END DO
      RUE(LC)=XTP(1)                 !  1  WA   = BIOMASS/ENERGY (kg/ha/MJ)(FOR CO2 = 330 ppm
      HI(LC)=XTP(2)                  !  2  HI   = HARVEST INDEX (CROP YIELD/ABOVE GROUND BIOMASS)
      optimalTem(LC)=XTP(3)           !  3  optimalTem = OPTIMAL TEMP FOR PLANT GROWTH (C)
      plantMinT(LC)=XTP(4)          !  4  plantMinT = MIN TEMP FOR PLANT GROWTH (C)
      maxLAI(LC)=XTP(5)             !  5  maxLAI = MAX LEAF AREA INDEX
      fracGrowSeasonLAIDecline(LC)=XTP(6)       !  6  DLAI = FRACTION OF GROWING SEASON WHEN LEAF AREA STARTS DECLINING
      pointsLAIDevp(1,LC)=XTP(7) !  7 & 8  pointsLAIDevp(= LAI DEVELOPMENT PARMS--NUMBERS BEFORE DECIMAL = %
      pointsLAIDevp(2,LC)=XTP(8) !  1,2) OF GROWING SEASON. NUMBERS AFTER DDECIMAL = FRACTION OF maxLAI AT GIVEN %.
      factorLAIDecline(LC)=XTP(9)       !  9  factorLAIDecline = LAI DECLINE RATE FACTOR. 
      factorBioDecline(LC)=XTP(10)      ! 10  factorBioDecline = RUE  DECLINE RATE FACTOR.
      cropTolerantAl(LC)=XTP(11)   ! 11  ALT  = ALUMINUM TOLERANCE INDEX(1-5)1=SENSITIVE, 5=TOLERANT
      maxStomaCond(LC)=XTP(12)   ! 12  maxStomaCond  = MAX STOMATAL CONDUCTANCE (DROUGTH TOLERANT PLANTS HAVE LOW VALUES.)
      areationThreshold(LC)=XTP(13) ! 13  areationThreshold  = CRITICAL AERATION FACTOR (Tot_SoilWater/POR > areationThreshold CAUSES AIR STRESS)
      seedRate(LC)=XTP(14)          ! 14  seedRate  = SEEDING RATE (kg/ha)
      maxCropHeight(LC)=XTP(15)     ! 15  HMX  = MAX CROP HEIGHT (m)
      maxRootDep(LC)=XTP(16)        ! 16  maxRootDep = MAX ROOT DEPTH (m)
      CO2EffOnBio(2,LC)=XTP(17)          ! 17  NUMBER BEFORE DECIMAL = CO2 CONC IN FUTURE ATMOSPHERE (ppm).NUMBER AFTER DECIMAL = RESULTING RUE VALUE.
      fracYieldN(LC)=XTP(18)            ! 18  fracYieldN  = N FRACTION OF YIELD (KG/KG)
      fracYieldP(LC)=XTP(19)            ! 19  fracYieldP  = P FRACTION OF YIELD (KG/KG)
      fracYieldK(LC)=XTP(20)            ! 20  fracYieldK  = K FRACTION OF YIELD (KG/KG)
      lowerLimitHI(LC)=XTP(21)           ! 21  lowerLimitHI = LOWER LIMIT OF HARVEST INDEX
      pestDmagFactor(LC)=XTP(22)        ! 22  PST  = PEST(INSECTS,WEEDS,AND DISEASE)FACTOR (0-1)
      seedCost(LC)=XTP(23)          ! 23  seedCost = SEED COST ($/KG)
      seedYieldPrice(LC)=XTP(24)    ! 24  seedYieldPrice = PRICE FOR SEED YIELD ($/t)
      forageYieldPrice(LC)=XTP(25)  ! 25  forageYieldPrice = PRICE FOR FORAGE YIELD ($/t)
      waterFracYield(LC)=XTP(26)       ! 26  WCY  = FRACTION WATER IN YIELD.
      uptakeParaN(1,LC)=XTP(27)     ! 27-29 BN   = N FRACTION IN PLANT WHEN GROWTH IS 0.,.5,1.0
      uptakeParaN(2,LC)=XTP(28)
      uptakeParaN(3,LC)=XTP(29)
      uptakeParaP(1,LC)=XTP(30)     ! 30-32 BP   = P FRACTION IN PLANT WHEN GROWTH IS 0.,.5,1.
      uptakeParaP(2,LC)=XTP(31)   
      uptakeParaP(3,LC)=XTP(32)
      uptakeParaK(1,LC)=XTP(33)     ! 33-35 uptakeParaK   = K FRACTION IN PLANT WHEN GROWTH IS 0.,.5,1.
      uptakeParaK(2,LC)=XTP(34)
      uptakeParaK(3,LC)=XTP(35)
      windEroCrop(1,LC)=XTP(36) !36-38 BW  = WIND EROSION FACTORS FOR STANDING LIVE, STANDING DEAD, AND FLAT RESIDUE
      windEroCrop(2,LC)=XTP(37)
      windEroCrop(3,LC)=XTP(38)
      cropCode(LC)=INT(XTP(39))         ! 39  cropCode  = CROP ID NUMBERS. USED TO PLACE CROPS IN CATEGORIES
                                                         !   (1 FOR WARM SEASON ANNUAL LEGUME
                                                         !   2 FOR COLD SEASON ANNUAL LEGUME
                                                         !   3 FOR PERENNIAL LEGUME
                                                         !   4 FOR WARM SEASON ANNUAL
                                                         !   5 FOR COLD SEASON ANNUAL
                                                         !   6 FOR PERENNIAL
                                                         !   7 FOR EVERGREEN TREES
                                                         !   8 FOR DECIDUOUS TREES
                                                         !   9 FOR COTTON
                                                         !   10 FOR DECIDUOUS TREES (LEGUME)
      pointsFrostDamg(1,LC)=XTP(40)        ! 40 pointsFrostDamg(= FROST DAMAGE PARMS--NUMBERS BEFORE DECIMAL = MIN TEMP(deg C
      pointsFrostDamg(2,LC)=XTP(41)        ! 41 (1,2) NUMBERS AFTER DECIMAL = FRACTION YLD LOST WHEN GIVEN MIN TE IS EXPERIENCED.
      VPDPara(LC)=XTP(42)            ! 42  VPDPara = PARM RELATING VAPOR PRESSURE DEFFICIT TO WA
      VPDThreshold(LC)=XTP(43)       ! 43  VPDThreshold = THRESHOLD VPD (kpa)(F=1.
      VPDPara2(LC)=XTP(44)                ! 44  VPDPara2 = NUMBER BEFORE DECIMAL = VPD VALUE (kpa).  NUMBER AFTER DECIMAL = F2 < 1.
      rootPartition(1,LC)=XTP(45)  ! 45  RootWt_fracStage1= ROOT WEIGHT/BIOMASS PARTITIONING COEF
      rootPartition(2,LC)=XTP(46)  ! 46  RootWt_fracStage2= ROOT WEIGHT/BIOMASS PARTITIONING COEF
      germinateHeatUnit(LC)=XTP(47)       ! 47  germinateHeatUnit = HEAT UNITS REQUIRED FOR GERMINATION
      pointsPlantPopulation(1,LC)=XTP(48)    ! 48  PlantPu_points1= PLANT POP PARM--NUMBER BEFORE DECIMAL = # PLANTS
                                                            !  NO AFTER DEC = FRACTION OF MAX LAI
      pointsPlantPopulation(2,LC)=XTP(49)    ! 49  PlantPu_points2= 2ND POINT ON PLANT POP-LAI CURVE. PlantPu_points1<PlantPu_points2--PLANTS/M2
                                                            !   PlantPu_points1>PlantPu_points2--PLANTS/HA
      salinityThreshold(1,LC)=XTP(50)       ! 50  STX1 = YLD DECREASE/SALINITY INCREASE ((t/ha)/(MMHO/CM))
      salinityThreshold(2,LC)=XTP(51)       ! 51  STX2 = SALINITY THRESHOLD (MMHO/CM)
      ligninFrac(1,LC)=XTP(52)  ! 52  LigninFrac_Plant1 = LIGNIN FRACTION IN PLANT AT .5 MATURITY
      ligninFrac(2,LC)=XTP(53)  ! 53  LigninFrac_Plant2 = LIGNIN FRACTION IN PLANT AT MATURITY
      WUE(LC)=XTP(54)       ! 54  WUE  = WATER USE CONVERSION TO BIOMASS(T/MM)
      turnoutFracCottn(LC)=XTP(55)      ! 55  turnoutFracCottn  = FRACTION TURNOUT (COTTON LINT/STRIPPER YLD)
      lintFracCottn(LC)=XTP(56)      ! 56  lintFracCottn  = FRACTION LINT (COTTON LINT/PICKER YLD)
      carbonEmissionSeedWeight(LC)=XTP(57)     ! 57  carbonEmissionSeedWeight = CARBON EMISSION/SEED WEIGHT(KG/KG)
      leafWeightFrac(LC)=XTP(58)         ! 58  leafWeightFrac = FRACTION OF LEAF WEIGHT TO AboveGround_Biomass WEIGHT AT fracGrowSeasonLAIDecline
      REWIND KR(4) 
      !======================== End of reading CROPCOM.DAT file =====================================
  179 Crop_Num= varID(JX(6))
      IHU(Crop_Num)=IHU(Crop_Num)+1               ! What is IHU ?
      potentialHeatUnit(Crop_Num,IHU(Crop_Num))=OPV(1) ! Potential heat unit 
      IF(XMTU(Crop_Num)>0)potentialHeatUnit(Crop_Num,IHU(Crop_Num))=XMTU(Crop_Num)* Cal_Heat(1,365,plantMinT(Crop_Num),0) 
      Y1=pointsPlantPopulation(1,Crop_Num)
      Y2=pointsPlantPopulation(2,Crop_Num)
      IF(Y2>Y1)THEN    ! Crops 
          X4=Y1
          X5=Y2
      ELSE             ! Trees
          X4=Y2
          X5=Y1
      END IF   
      X1=split(X4)     ! first point on population curve X1: plants/ha for trees, plants/m2 for crops  
      X2=split(X5)     ! Second point on population curve X2: plants/ha for trees, plants/m2 for crops          
      CALL Coeff_Scurve(X4,X5,X1,X2)
      coeffPopuCurve(1,Crop_Num)=X4         ! Coefficient b1 in population curve
      coeffPopuCurve(2,Crop_Num)=X5         ! Coefficient b2 in population curve
      IF(OPV(5)>0.)THEN                     ! Plant population   if OPV5 >0 then population will be determined by .OPSC file 
            X3=OPV(5)
      ELSE               
            G1=X2
            DO IT=1,10
                  Z1=EXP(X4-X5*G1)
                  Z2=G1+Z1
                  FU=G1/Z2-.9
                  IF(ABS(FU)<1.E-5)EXIT
                  DFDU=Z1*(1.+X5*G1)/(Z2*Z2)
                  G1=G1-FU/DFDU
            END DO
            IF(IT>10)WRITE(KW(1),5)
            X3=G1
      END IF
      PPLA(Crop_Num,IHU(Crop_Num))=maxLAI(Crop_Num)*X3/(X3+EXP(X4-X5*X3))  ! ?
      POP(Crop_Num,IHU(Crop_Num))=X3
      IF(OPV(6)>0.) maxAnnualN=OPV(6)      ! if OPV6 >0 then population will be determined by .OPSC file
      FNMX(Crop_Num)= maxAnnualN            ! Maximum N Annualy  

    1 FORMAT(A12)
    5 FORMAT(5X,'!!!!! PLANT POP DID NOT CONVERGE')
  333 FORMAT(1X,I4,1X,A4,100F8.0)
END SUBROUTINE Read_CropTable



! ---------------------- 5. Operation inputs ------------------------------
SUBROUTINE Operation_Inputs(I3,II,JJ,JRT)
! EPIC1102
! THIS SUBPROGRAM ALLOWS INPUT TO OPERATION SCHEDULE FOR IRRIGATION,
! FERTILIZER, OR PESTICIDE APPLICATION
      IMPLICIT NONE
	! local variables
      INTEGER:: I1, I3, II, JJ, JRT, L
	
      I1=I3-6
      SELECT CASE(I1)
            CASE(1)
                  pestRate(II,JJ)=OPV(1)
                  pestKillEffi(II,JJ)=OPV(2)
                  CALL Pest_Table
                  LPC(II,JJ)=KDP1(JX(7))
                  numOfPest=numOfPest+1
            CASE(2)
                  irrVolumeArr(II,JJ)=OPV(1)
                  IF(OPV(4)>0.)IrrRunoffRatio=OPV(4)
                  numOfIrr=numOfIrr+1
            CASE(3)
                  WFA(II,JJ)=OPV(1)
                  CALL Fert_Table(L)    
                  LFT(II,JJ)=KDF1(JX(7))
                  numOfFert=numOfFert+1
            CASE(13)
                  IF(OPV(1)>0.)baseStockRate=OPV(1)
            CASE(21)
                  TLMA(II,JJ)=OPV(1)                          
            CASE DEFAULT
                  GO TO 169
      END SELECT
      irrTriggerArr(II,JJ)=irrTrigger
      CND(II,JJ)=CN2
      QIR(II,JJ)=IrrRunoffRatio
      baseStockRateArr(II,JJ)=baseStockRate
      JRT=1
      RETURN
  169 IF(OPV(2)<0.)THEN
            CN2=-OPV(2)
      ELSE
            IF(OPV(2)>0.)THEN
                  landUseCode=INT(OPV(2))          
                  CALL CN2_Table 
            END IF
      END IF
      CND(II,JJ)=CN2
      IF(ABS(OPV(3))>1.E-5)irrTrigger=OPV(3)
      irrTriggerArr(II,JJ)=irrTrigger
      IF(OPV(4)>0.)IrrRunoffRatio=OPV(4)
      QIR(II,JJ)=IrrRunoffRatio
      IF(OPV(8)>0.)minCFactor=OPV(8)
      fertCFactirArr(II,JJ)=minCFactor
      IF(OPV(9)>0.)HWC(II,JJ)=OPV(9)
      JRT=0
      RETURN
END SUBROUTINE Operation_Inputs
    
END MODULE ReadVarTable_Module