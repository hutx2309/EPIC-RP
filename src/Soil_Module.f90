MODULE Soil_Module
USE PARM
USE MisceFun_Module, ONLY: OPENV
IMPLICIT NONE
      CONTAINS
 
SUBROUTINE Read_Soil(FSOIL, soilID, yearTotPPT, soilName, soilHydGroup, fracFC, groundWaterResidenceTime0, &
                     splitedSoilLayer,soilCondCode, cultivaYrs, soilGroup, minThickSplit, fracSOC, fracHumus,&
                     XCC, NCC, sulfCon, maxLayer, fracVPip, fracHPip, IDSK)
      IMPLICIT NONE
      ! Argument list
      CHARACTER(LEN = 80), INTENT(IN):: FSOIL
      REAL, INTENT(IN):: yearTotPPT
      INTEGER, INTENT(IN):: soilID
      INTEGER, INTENT(OUT):: IDSK
      CHARACTER(LEN = 80), INTENT(OUT):: soilName
      REAL, DIMENSION(15), INTENT(OUT):: sulfCon, fracVPip, fracHPip
      REAL, INTENT(OUT):: soilHydGroup, fracFC, groundWaterResidenceTime0, splitedSoilLayer, &
                          soilCondCode, cultivaYrs,soilGroup, minThickSplit, fracSOC, fracHumus, &
                          XCC
      INTEGER, INTENT(OUT):: maxLayer, NCC
      ! Local variables declaration
      CHARACTER(LEN = 8):: FMT
      INTEGER:: I, J, II
      !=========================== Search Soil file ===================================
      CALL OPENV(KR(13),FSOIL,IDIR(3),KW(MSO) ) ! Catalog of soil data files
      II=-1 
      DO WHILE(II/=soilID) 
            READ(KR(13),*,IOSTAT=NFL)II,SOILFILE  
            IF(NFL/=0)THEN
                  IF(runMode==0)THEN
                        WRITE(*,'(T10,A,I8,A)')'SOIL NO = ',soilID,&
                              ' NOT IN SOIL LIST FILE'
                        ERROR STOP   
                  ELSE 
                        WRITE(KW(MSO),'(A,A8,A,I4,A)')' !!!!! ',siteName,&
                              'SOIL NO = ',soilID,' NOT IN SOIL LIST FILE'
                        ERROR STOP
                  END IF    
            END IF  
      END DO  
      REWIND KR(13)   

      !--------------------------  read soil data -----------------------------
      CALL OPENV(KR(14),SOILFILE,IDIR(3),KW(MSO))
      !   LINE 1    
      READ(KR(14),*)soilName, soilOrder   
      !   READ SOIL DATA  
      !   LINE 2/3   
      READ(KR(14),'(10F8.0)')soilAlbedo, soilHydGroup, fracFC, waterTableMinDep, &
            waterTableMaxDep, waterTableHigt, groundWaterStor, groundWaterMaxStor, &
            groundWaterResidenceTime0, returnFlowFrac, splitedSoilLayer, soilCondCode, &
            cultivaYrs,soilGroup,minThickMaxLayer,minThickProfile,minThickSplit,&
            fracSOC,fracHumus,XCC  

      NCC=INT(XCC)                  ! XCC: Code written automatically for *.sot (not user input)   
      
      IF(groundWaterStor<1.E-10) groundWaterStor=25.          ! mm  
      IF(groundWaterMaxStor<1.E-10) groundWaterMaxStor=50.    ! mm  
      IF(waterTableMaxDep<1.E-5)THEN  
            waterTableMinDep=50.                              ! m  
            waterTableMaxDep=100.                             ! m  
            waterTableHigt=75.                                ! m  
      END IF           
      
      IDSP=INT(soilCondCode+1.1)       ! P related control variable 

      IF(fracSOC<1.E-10) fracSOC=.04  
      ! Eq. 2.158 in EPIC-doc P38  cultivaYrs: period of cultivation before the simulation starts
      IF(fracHumus<1.E-10) fracHumus=.7-.3*EXP(-.0277*cultivaYrs) 

      IDSK=INT(MAX(soilGroup,1.))      ! K related control variable 

      IF(fracFC<1.E-10)fracFC=yearTotPPT/(yearTotPPT+EXP(9.043-.002135*yearTotPPT))  

      !------------------------- initialize the variables of each soil layers -------------------  
      DO I=1,MSL                 ! MSL =15 in allocate Para   
            NH3_Weight(I)=0.     ! NH3 weight     
            soilElement(I)=0.         ! Potential Soil Evaporation for a layer 
            wiltPoint(I)=0.      ! soil water content at wilting point (1500 Kpa) (m/m)  
            LORG(I)=I            ! ???  
            Layer_ID(I)=I                                                                 
            DO J=1,MPS           ! MPS =40                                                      
                  pestInSoil(J,I)=0.  ! ????                                                        
            END DO
      END DO   

      pestPlantIntercep=0.      ! P related   why is it here?
      maxLayer=INT(splitedSoilLayer)    
      IF(minThickMaxLayer<1.E-10) minThickMaxLayer=.1   
      IF(minThickProfile<1.E-10) minThickProfile=.1  
      IF(minThickSplit<1.E-3) minThickSplit=.15  
      groundWaterResidenceTime=groundWaterResidenceTime0           
      IF(groundWaterResidenceTime0<1.E-10) groundWaterResidenceTime=10.   
      ! ------------------------------- continue to read soil file -----------------------------
      ! THE SOIL IS DIVIDED VERTICALLY INTO LAYERS (MAX OF 10 LAYERS OF USER SPECIFIED THICKNESS)                                                   
      ! LINES 4/47 
      FMT='(15F8.2)'    
      READ(KR(14),FMT)(Z(I),I=1,MSL)                !  4  Z  = DEPTH TO BOTTOM OF LAYERS(m)  
      READ(KR(14),FMT)(bulkDensity(I),I=1,MSL)      !  5  bulkDensity = BULK DENSITY (T/M**3) 
      READ(KR(14),FMT)(wiltPoint(I),I=1,MSL)        !  6  wiltPoint  = SOIL WATER CONTENT AT WILTING POINT (1500 kpa)(m/m) ! (BIU) 
      READ(KR(14),FMT)(fieldCapacity(I),I=1,MSL)    !  7  fieldCapacity = WATER CONTENT AT FIELD CAPACITY (33kpa)(m/m) (BIU)  
      READ(KR(14),FMT)(sandFrac(I),I=1,MSL)         !  8  % sand
      READ(KR(14),FMT)(siltFrac(I),I=1,MSL)         !  9  % silt
      READ(KR(14),FMT)(SON(I),I=1,MSL)              !  10 INITIAL ORGANIC N CONC (g/t) (BIU)  
      READ(KR(14),FMT)(PH(I),I=1,MSL)               !  11 SOIL PH 
      READ(KR(14),FMT)(totBase(I),I=1,MSL)          !  12  SMB  = SUM OF BASES (cmol/kg) (BIU)   
      READ(KR(14),FMT)(SOC(I),I=1,MSL)              !  13  WOC  = ORGANIC CARBON CONC(%)
      READ(KR(14),FMT)(CaCO3(I),I=1,MSL)            !  14  CALCIUM CARBONATE (%) 
      READ(KR(14),FMT)(CEC(I),I=1,MSL)              !  15  CEC  = CATION EXCHANGE CAPACITY (cmol/kg) (BIU) 
      READ(KR(14),FMT)(rockFrac(I),I=1,MSL)         !  16  ROK  = COARSE FRAGMENTS (% VOL) (BIU)  
      READ(KR(14),FMT)(initialNO3(I),I=1,MSL)       !  17  INITIAL NO3 CONC (g/t)     (BIU) 
      READ(KR(14),FMT)(initialLabileP(I),I=1,MSL)   ! 18  INITIAL LABILE P CONC (g/t)(BIU) 
      READ(KR(14),FMT)(cropResidu(I),I=1,MSL)       ! CROP RESIDUE(t/ha)  
      READ(KR(14),FMT)(dryBulkDensity(I),I=1,MSL)   ! BULK DENSITY (OVEN DRY)(T/M**3) (BIU)
      READ(KR(14),FMT)(sorbRatioP(I),I=1,MSL)       ! = P SORPTION RATIO <1.          (BIU)  
                                                    ! = ACTIVE & STABLE MINERAL P (kg/ha) >1.     
      READ(KR(14),FMT)(satuateCondv(I),I=1,MSL)     ! SATURATED CONDUCTIVITY (mm/h)      (BIU)  
      READ(KR(14),FMT)(lateralFlowCondv(I),I=1,MSL) ! LATERAL HYDRAULIC CONDUCTIVITY(mm/h)(BIU) 
      READ(KR(14),FMT)(SOP(I),I=1,MSL)              ! INITIAL ORGANIC P CONC (g/t)       (BIU) 
      READ(KR(14),FMT)(exchangeK(I),I=1,MSL)        ! EXCHANGEABLE K CONC (g/t) 
      READ(KR(14),FMT)(elecCondv(I),I=1,MSL)        ! ELECTRICAL COND (mmHO/CM)
      READ(KR(14),FMT)(fracNO3Leach(I),I=1,MSL)     ! FRACTION OF STORAGE INTERACTING WITH NO3 LEACHING (BIU)
      READ(KR(14),FMT)(soilWater(I),I=1,MSL)        ! INITIAL SOIL WATER STORAGE (m/m)  
      READ(KR(14),FMT)(fracVPip(I),I=1,MSL)         ! FRACTION INFLOW PARTITIONED TO VERTICLE CRACK OR PIPE FLOW
      READ(KR(14),FMT)(fracHPip(I),I=1,MSL)         ! FRACTION INFLOW PARTITIONED TO HORIZONTAL CRACK OR PIPE FLOW 
      READ(KR(14),FMT)(structLitt(I),I=1,MSL)       ! STRUCTURAL LITTER(kg/ha)           (BIU)
      READ(KR(14),FMT)(metabLitt(I),I=1,MSL)        ! METABOLIC LITTER(kg/ha)            (BIU)
      READ(KR(14),FMT)(lgStructLitt(I),I=1,MSL)     ! LIGNIN CONTENT OF STRUCTURAL LITTER(kg/ha)(BIU)  
      READ(KR(14),FMT)(CStructLitt(I),I=1,MSL)      ! CARBON CONTENT OF STRUCTURAL LITTER(kg/ha)(BIU)
      READ(KR(14),FMT)(CMetabLitt(I),I=1,MSL)       ! C CONTENT OF METABOLIC LITTER(kg/ha)(BIU)
      READ(KR(14),FMT)(CLgStructLitt(I),I=1,MSL)    ! C CONTENT OF LIGNIN OF STRUCTURAL LITTER(kg/ha)(BIU)
      READ(KR(14),FMT)(NLgStructLitt(I),I=1,MSL)    ! N CONTENT OF LIGNIN OF STRUCTURAL LITTER(kg/ha)(BIU)
      READ(KR(14),FMT)(CBiomass(I),I=1,MSL)         ! C CONTENT OF BIOMASS(kg/ha)(BIU)  
      READ(KR(14),FMT)(CSlowHumus(I),I=1,MSL)       ! C CONTENT OF SLOW HUMUS(kg/ha)(BIU) 
      READ(KR(14),FMT)(CPassiveHumus(I),I=1,MSL)    ! C CONTENT OF PASSIVE HUMUS(kg/ha)(BIU) 
      READ(KR(14),FMT)(NStructLitt(I),I=1,MSL)      ! N CONTENT OF STRUCTURAL LITTER(kg/ha)(BIU)
      READ(KR(14),FMT)(NMetabLitt(I),I=1,MSL)       ! N CONTENT OF METABOLIC LITTER(kg/ha)(BIU) 
      READ(KR(14),FMT)(NBiomass(I),I=1,MSL)         ! N CONTENT OF BIOMASS(kg/ha)(BIU)
      READ(KR(14),FMT)(NSlowHumus(I),I=1,MSL)       ! N CONTENT OF SLOW HUMUS(kg/ha)(BIU)  
      READ(KR(14),FMT)(NPassiveHumus(I),I=1,MSL)    ! N CONTENT OF PASSIVE HUMUS(kg/ha)(BIU)
      READ(KR(14),FMT)(ironCon(I),I=1,MSL)          ! IRON CONTENT(%)  
      READ(KR(14),FMT)(sulfCon(I),I=1,MSL)          ! SULFUR CONTENT(%)
      ! LINE 48      
      READ(KR(14),'(15A8)')soilHorizon              ! SOIL HORIZON(A,B,C)
      ! LINES 49/51      
      IF(deNitri_Method>2)THEN
            READ(KR(14),FMT)(gasO2Con(J),J=1,MSL)       ! O2 CONC IN GAS PHASE  (g/m3 OF SOIL AIR) 
            READ(KR(14),FMT)(gasCO2Con(J),J=1,MSL)      ! CO2 CONC IN GAS PHASE (g/m3 OF SOIL AIR)
            READ(KR(14),FMT)(gasN2OCon(J),J=1,MSL)      ! N2O CONC IN GAS PHASE (g/m3 OF SOIL AIR)
      ELSE
            gasO2Con=0.
            gasCO2Con=0.
            gasN2OCon=0.
      END IF
      REWIND KR(14) !=============================== End of reading soil file ======================================== 

END SUBROUTINE Read_Soil

SUBROUTINE Soil_Process(NCC, IDSK, maxLayer, minThickSplit, fracFC, fracSOC, fracHumus, PZW, XX)
      ! initialzie soil variables layer by layer
      IMPLICIT NONE
      ! Argument lists
      REAL, INTENT(IN):: fracFC, minThickSplit
      INTEGER, INTENT(IN):: NCC, IDSK
      INTEGER, INTENT(INOUT):: maxLayer
      REAL, INTENT(INOUT):: fracSOC, fracHumus, PZW, XX
      ! local variables declaration
      INTEGER:: IT, J, K, KK, K1, L, L1, LZ, MXZ
      REAL:: DFDN, FU, DF, VG1, FU1, SUM, TPAW,  DG, DG1, F, RTO, VGA0, VGN0, &
             wiltPoint0, Porosity0, FC0, BSA, &
             WBM, WHP, WHS, & ! intermediate variables for cal C_con_Phumus and C_con_Shumus
             WT1, X1, X2, XY, XZ, ZD, ZMX, ZZ, XCB = 0.2
      ! ////////////////////////////////////////////////////////////////////
      L = 1
      SUM=0.  
      LZ=1
      K=1 
      TPAW=0. 
      XX=0.          ! XX: soil depth   
      
      Soil_Layers: DO J=1,MSL                               ! SOIL LOOP  
            IF(Z(J)<1.E-10)EXIT                             ! Z: from the surface to the bottom of the soil layer 1 in m
            clayFrac(J)=100.-sandFrac(J)-siltFrac(J)        ! Clay percentage  
            DG=1000.*(Z(J)-XX)                              ! m to mm?  
         
            CALL BulkDens_Effect(bulkDensity(J),modelPara(2),F,J,1 ) ! F will be calculated in the subroutine  
            CALL Distribute_NPK(cropResidu,DG,DG1,1.,.01,J, MSL)          ! DG1 is not initialized, but is output  
            CALL Distribute_NPK(initialLabileP,DG,DG,20.,.001,J, MSL)  
            CALL Distribute_NPK(initialNO3,DG,DG,10.,.001,J, MSL) 
            CALL Distribute_NPK(exchangeK,DG,DG,10.,.001,J, MSL) 
          
            ! (1) NO3 leach fraction 
            IF(fracNO3Leach(J)<1.E-10) fracNO3Leach(J)=storLeachFracNO3    ! SFT0: FRACTION OF STORAGE INTERACTING WITH NO3 LEACHING  
            ! (2) crop residual 
            totCropResi = totCropResi + cropResidu(J)                      ! RSD: CROP RESIDUE(t/ha)  
          
            ZD=.25*(XX+Z(J))                      ! XX = 0., for j =1 ; otherwise = Z(J-1)  
            F=ZD/(ZD+EXP(-.8669-2.0775*ZD))       ! Eq.269a in APEX-doc  F==FZ, ZD == X1   
            ! (3) Soil temperature 
            soilTem(J)=F*(yearAveTem-TX)+TX       ! Eq.269 in APEX-doc with LAG =0.0 
           
            ! (4) Initial Organic C concentration
            IF(SOC(J)<1.E-5) SOC(J)=XCB*EXP(-.001*DG) ! XCB =0.2                                                                    
          
            ! (5) Minearal Bulk Density by layer
            XZ=SOC(J)*.0172                                                                
            ZZ=1.-XZ                                                                       
            mineralBulkDensity(J)=ZZ/(1./bulkDensity(J)-XZ/.224) ! mineralBulkDensity: Minearal Bulk Density          
            IF(mineralBulkDensity(J)>2.65)THEN
                  mineralBulkDensity(J)=2.65
            ELSE                                                                 
                  IF(mineralBulkDensity(J)<1.)THEN     
                        mineralBulkDensity(J)=1.                                                                         
                        bulkDensity(J)=1./(ZZ+XZ/.224)  
                  END IF
            END IF           
            ! (6) Weight 
            WT(J)=bulkDensity(J)*DG*10.              ! Weight  
            DG1=DG                                                                         
            ! update Organic C content and organic N
            WT1=WT(J)/1000.                                                                
            X1=10.*SOC(J)*WT(J)                                                          
            SOC(J)=X1  
 
            IF(SON(J)>0.)THEN                                                              
                  SON(J)=WT1*SON(J)   
                  KK=0  
            ELSE  
                  SON(J)=.1*SOC(J)  
                  KK=1  
            END IF    
            
            ! -------------  XCC (NCC) control ---------- 
            ! NCC = 1: use .sot file created by EPIC model, no more divide;
            !     = 0: divide the soil into 10 layers           
            IF(NCC==0)THEN  
                  ! (7) C concentration in biomass (kg/ha)
                  WBM=fracSOC*X1             ! X1 here is organic C concentration                                                        
                  CBiomass(J)=WBM  
                  IF(KK==0)THEN  
                        RTO=SON(J)/SOC(J)  
                  ELSE    
                        RTO=.1   
                  END IF     
                  ! (8) N concentration in biomass (kg/ha)
                  NBiomass(J)=RTO*CBiomass(J)  
   
                  ! (9) CSlowHumus and NSlowHumus
                  WHP=fracHumus*(X1-WBM)                                                               
                  WHS=X1-WBM-WHP                                                                 
                  CSlowHumus(J)=WHS                                                                    
                  NSlowHumus(J)=RTO*CSlowHumus(J) 
   
                  ! (10) CPassiveHumus and NPassiveHumus 
                  CPassiveHumus(J)=WHP 
                  NPassiveHumus(J)=RTO*CPassiveHumus(J) 
   
                  ! (11) metabLitt and Structural Litt 
                  X1=cropResidu(J)   
                  IF(J==1)X1=X1+deadCropRes0  
                  metabLitt(J)=500.*X1  
                  structLitt(J)=metabLitt(J) 
                  ! (12) laging, C, N , C_laging, N_laging concentration
                  lgStructLitt(J) = .8*structLitt(J) 
                  CMetabLitt(J) =   .42*metabLitt(J) 
                  NMetabLitt(J) = .1*CMetabLitt(J)  
                  CStructLitt(J) = .42*structLitt(J) 
                  CLgStructLitt(J) = .8*CStructLitt(J) 
                  NLgStructLitt(J) = .2*CStructLitt(J) 
                  NStructLitt(J) = CStructLitt(J)/150. 

                  ! (13) total organic C and N concentration
                  SOC(J)=SOC(J)+CStructLitt(J)+CMetabLitt(J)  
                  SON(J)=SON(J)+NStructLitt(J)+NMetabLitt(J)
            END IF    
            ! (14.1) Fresh organic P
            FreshOrgP_Residu(J)=cropResidu(J)*1.1 
            soilElement(4)=soilElement(4)+FreshOrgP_Residu(J)  ! accumulate fresh P residual in each layer                                         
            ! (14.2) Active mineral P
            activeMineralP(J)=0.  
            ! (14.3) Initial organic P concentration 
            IF(SOP(J)>0.)THEN  
                  SOP(J)=WT1*SOP(J)  
            ELSE                                                                                
                  SOP(J)=.125*SON(J)                                                            
            END IF         
            ! (15) Porosity
            Porosity(J)=1.-bulkDensity(J)/2.65   ! Porosity: porosity of the soil (mm/mm)  SWAT2009-doc P109              
            ! (16) CEC
            ZZ=.5*(XX+Z(J))                                                                
            X2=.1*SOC(J)/WT(J)                                                             
            X1=MIN(.8*clayFrac(J),5.+2.428*X2+1.7*ZZ)                                           
            CEC(J)=MAX(CEC(J),X1)             
            ! (17) Cal soil water content and fieldCapacity using different models
            SELECT CASE(FC_Method+1) 
                  CASE(1,5)   
                        CALL Cal_SoilWaterW(clayFrac(J),sandFrac(J),X2,wiltPoint(J),fieldCapacity(J))  
                  CASE(2,6) 
                        CEM(J)=MAX(.1,CEC(J)-2.428*X2-1.7*ZZ)
                        CALL Cal_SoilWaterO(CEM(J),clayFrac(J),X2,sandFrac(J),wiltPoint(J),fieldCapacity(J),ZZ) 
                  CASE(7)
                        IF(wiltPoint(J)<1.E-10.OR.fieldCapacity(J)<1.E-10)THEN
                              CEM(J)=MAX(.1,CEC(J)-2.428*X2-1.7*ZZ)
                              CALL Cal_SoilWaterO(CEM(J),clayFrac(J),X2,sandFrac(J),wiltPoint(J),fieldCapacity(J),ZZ)
                        END IF  
                  CASE(8,9) 
                        CALL Cal_FC_WP(clayFrac(J),sandFrac(J),X2,wiltPoint(J),fieldCapacity(J)) 
                  CASE(10,11) 
                        CALL Cal_SoilWaterBNW(clayFrac(J),siltFrac(J),sandFrac(J),X2,bulkDensity(J),wiltPoint(J),fieldCapacity(J))
                  CASE DEFAULT
                        IF(FC_Method/=2.AND.FC_Method/=3)&
                        CALL Cal_SoilWaterW(clayFrac(J),sandFrac(J),X2,wiltPoint(J),fieldCapacity(J))
                        FC_Method=0  
            END SELECT

            IF(rockFrac(J)>=99.)rockFrac(J)=90. ! rockFrac: coarse fragments
            XY=1.-rockFrac(J)*.01  
            wiltPoint(J)=wiltPoint(J)*XY  
            XY=XY*DG  
            fieldCapacity(J)=fieldCapacity(J)*XY 
            ! CAP=CAP+fieldCapacity(J)    ! total field capacity   but not used by the program 

            ! (18) Wilt point cal by water content and DG
            Wilt_Point(J)=wiltPoint(J)*DG              !Wilting Point  
 
            ! update porosity
            Porosity(J)=Porosity(J)*XY
 
            ! saturated wilting point????
            CALL WiltPoint_Sat(J)                 ! recalculate wilt_point ?
            ! (19) Initial soil water 
            IF(soilWater(J)<1.E-10.AND.NCC==0)soilWater(J)=fracFC  
            soilWater(J)=soilWater(J)*(fieldCapacity(J)-Wilt_Point(J))+Wilt_Point(J) 
 
            ! (20) VGA?
            VGA(J)=EXP(-4.3003-.0097*clayFrac(J)+.0138*sandFrac(J)-.1706*X2)  !??? What is VGA?
            VGN0=1.                               ! Is this line necessary?
            VGA0=VGA(J)
            wiltPoint0=Wilt_Point(J)
            Porosity0=Porosity(J)
            FC0=fieldCapacity(J)

            DO IT=1,10           !?? calculate 10 times? what does this block mean? all variables are local
                  FU=wiltPoint0+(Porosity0-wiltPoint0)/(1.+(VGA0*336.27)**VGN0)**(1.-1./VGN0)-FC0
                  IF(ABS(FU)<1.E-5)EXIT
                  IF(IT==1)THEN
                        DF=.01
                  ELSE
                        DFDN=(FU-FU1)/(VGN0-VG1)
                        DF=FU/DFDN
                  END IF
                  VG1=VGN0
                  FU1=FU
                  VGN0=VGN0-DF
                  IF(ABS(DF)<1.E-5)EXIT
            END DO  

            ! (21) VGN ?
            VGN(J)=VGN0   
            soilElement(1)=soilElement(1)+XY             ! accumulate soil components except for coarse fragments  
            soilElement(3)=soilElement(3)+WT(J)          ! accumulate the weight of each layer 

            ! (22) Lateral hydraulic conductivity (mm/h)
            IF(lateralFlowCondv(J)<1.E-20) lateralFlowCondv(J)=uplandSteep*satuateCondv(J)  ! APEX-doc Eq. 82 
  
            ! (23) Al saturation 
            IF(CEC(J)>0.)THEN 
                  IF(CaCO3(J)>0.) totBase(J)=CEC(J)                                                     
                  IF(totBase(J)>CEC(J)) totBase(J)=CEC(J)       ! CEC: cation exchange capacity                                            
                  BSA=totBase(J)*100./(CEC(J)+1.E-20)! BSA: base saturation                                            
                  IF(PH(J)>5.6)THEN
                        ALS(J)=0.       ! ALS: AL saturation of soil layer 1 in % calculated as KCl-extractable Al divied by effective cation exchange capacity (ECEC)
                  ELSE
                        X1=.1*SOC(J)/WT(J)                                      ! SOC: organic carbon content in %        
                        ALS(J)=154.2-1.017*BSA-3.173*X1-14.23*PH(J)  ! Eq. 309 in APEX-doc   Lime  
                        IF(ALS(J)<0.)THEN
                              ALS(J)=0.
                        ELSE
                              IF(ALS(J)>95.) ALS(J)=95.
                        END IF
                  END IF
            ELSE
                  CEC(J)=PH(J) 
                  totBase(J)=CEC(J) 
                  ALS(J)=0.
            END IF  

            ! (24) P related varialbes: P absorbtion ratio, flow coefficient, stableMineralP
            SELECT CASE(IDSP)            ! soil condiations    P_related soil property 
                  CASE(1)                  ! in calcareous soils 
                        IF(CaCO3(J)>0.)THEN  
                              sorbRatioP(J)=.58-.0061*CaCO3(J) 
                              flowCoeffP(J)=.00076  
                        ELSE  
                              sorbRatioP(J)=.5 
                              flowCoeffP(J)=EXP(-1.77*sorbRatioP(J)-7.05)  
                        END IF 
                  CASE(2)                  ! in noncalcareous, slightly weathered soils  
                        sorbRatioP(J)=.02+.0104*initialLabileP(J) 
                        flowCoeffP(J)=EXP(-1.77*sorbRatioP(J)-7.05) 
                  CASE(3)                  ! in noncalcareous, moderately weathered soils 
                        sorbRatioP(J)=.0054*BSA+.116*PH(J)-.73 
                        flowCoeffP(J)=EXP(-1.77*sorbRatioP(J)-7.05)  
                  CASE(4)                  ! in noncalcareous 
                        sorbRatioP(J)=.46-.0916*LOG(clayFrac(J))    ! Eq.242 in APEX-doc    sorbRatioP: the P sorption coefficient (PSP) 
                        flowCoeffP(J)=EXP(-1.77*sorbRatioP(J)-7.05)! Eq. 244 in APEX-doc 
                  CASE(5)  
                        IF(sorbRatioP(J)<1.)THEN  
                              IF(CaCO3(J)>0.)THEN  
                                    sorbRatioP(J)=.58-.0061*CaCO3(J) 
                                    flowCoeffP(J)=.00076 
                              ELSE 
                                    sorbRatioP(J)=.5 
                                    flowCoeffP(J)=EXP(-1.77*sorbRatioP(J)-7.05) 
                              END IF 
                        ELSE  
                              activeMineralP(J)=.2*sorbRatioP(J)    
                              sorbRatioP(J)=initialLabileP(J)/(activeMineralP(J)+initialLabileP(J))  
                        END IF 
            END SELECT 
   
            IF(sorbRatioP(J)<.05) sorbRatioP(J)=.05      ! pg. 53 in APEX-doc                                                       
            IF(sorbRatioP(J)>.75) sorbRatioP(J)=.75      ! P absorption coefficient                                                  
            IF(activeMineralP(J)<1.E-5) activeMineralP(J)=initialLabileP(J)*(1.-sorbRatioP(J))/sorbRatioP(J)  
            stableMineralP(J)=4.*activeMineralP(J)       ! Eq. 243??  APEX-doc  
 
            ! (25) K related variables: soluble K, Fixed K   
            SELECT CASE(IDSK)    ! Soil group           K_related soil choice                                               
                  CASE(1)                                                                           
                        solubleK(J)=MAX(.05*exchangeK(J),.052*exchangeK(J)-.12)  
                        fixK(J)=374.+236.*clayFrac(J)  
                  CASE(2)  
                        solubleK(J)=.026*exchangeK(J)+.5                       !????  
                        fixK(J)=1781.+316.*clayFrac(J)  
                  CASE(3)  
                        solubleK(J)=.026*exchangeK(J)+.5  
                        fixK(J)=1781.+316.*clayFrac(J) 
            END SELECT  
              
            EQKS(J)=solubleK(J)/exchangeK(J) 
            EQKE(J)=exchangeK(J)/fixK(J) 

            activeMineralP(J)=activeMineralP(J)*WT1 
            stableMineralP(J)=stableMineralP(J)*WT1 
            labileP(J)=initialLabileP(J)*WT1 
            totLabileP=totLabileP+labileP(J)      ! labileP: Labile P concentration (1) Estimated current labile P concentration at start of simulation.                                           
            totActiveP=totActiveP+activeMineralP(J) 

            saltWeight(J)=6.4*elecCondv(J)*soilWater(J) ! Salt Weight in each layer  
            totSalt=totSalt+saltWeight(J)               ! totSalt 
 
            exchangeK(J)=exchangeK(J)*WT1 
            solubleK(J)=solubleK(J)*WT1  
            fixK(J)=fixK(J)*WT1 

            NO3_N_Soil(J)=initialNO3(J)*WT1  
            totNO3_N=totNO3_N+NO3_N_Soil(J)  
 
            IF(Z(J)<=RZ)THEN                            ! root zone  
                  soilWaterRZ=soilWaterRZ+soilWater(J)-Wilt_Point(J)  
                  potentialAvailWater=potentialAvailWater + fieldCapacity(J) - Wilt_Point(J) 
                  LZ=J                                  ! LZ: number of soil layers within root depth  
            END IF 
 
            IF(dryBulkDensity(J)<1.E-5) dryBulkDensity(J)=bulkDensity(J) 

            postTillBulkDensity(J) = bulkDensity(J)     ! postTillBulkDensity: post tillage bulk density -- Zhao  
            dryBulkDensity(J) = dryBulkDensity(J)/bulkDensity(J)  

            TOP=TOP+stableMineralP(J) 
            initialTotSOP=initialTotSOP+SOP(J) 

            totSolubleK=totSolubleK+solubleK(J) 
            totExchangeK=totExchangeK+exchangeK(J) 
            totFixK=totFixK+fixK(J) 
          
            totStructLitt=totStructLitt+structLitt(J) 
            totMetabLitt=totMetabLitt+metabLitt(J) 
            totLgStructLitt=totLgStructLitt+lgStructLitt(J) 
          
            totCStructLitt=totCStructLitt+CStructLitt(J) 
            totCMetabLitt=totCMetabLitt+CMetabLitt(J) 
            totCLgStructLitt=totCLgStructLitt+CLgStructLitt(J) 
            totNLgStructLitt=totNLgStructLitt+NLgStructLitt(J)  

            totCBiomass=totCBiomass+CBiomass(J) 
            totCSlowHumus=totCSlowHumus+CSlowHumus(J) 
            totCPassiveHumus=totCPassiveHumus+CPassiveHumus(J) 

            totNStructLitt=totNStructLitt+NStructLitt(J)  
            totNMetabLitt=totNMetabLitt+NMetabLitt(J)  
            totNBiomass=totNBiomass+NBiomass(J) 
            totNSlowHumus=totNSlowHumus+NSlowHumus(J)  
            totNPassiveHumus=totNPassiveHumus+NPassiveHumus(J)  
 
            IF(K<=3)THEN                                                               
                  TPAW=TPAW+fieldCapacity(J)-Wilt_Point(J)  ! total potential root available water?  
                  DO K1=K,3
                        IF(TPAW<WCS(K1))EXIT
                        ZCS(K1)=XX+(Z(J)-XX)*((WCS(K1)-PZW)/(TPAW-PZW)) 
                  END DO
                  K=K1
            END IF

            PZW=TPAW                                           ! potential root available water? 
            XX=Z(J)
      END DO   Soil_Layers        ! SOIL LOOP
      ! ----------------------------  soil loop end ------------------------------------
      IF(J>MSL)THEN
            Actual_SoilLayers=10        ! Number of allowed soil layers  
      ELSE
            L1=LZ+1                     ! round to the ceil layer for root depth 
            Actual_SoilLayers=J-1       ! actual soil layers  
            IF(L1/=J)THEN     
                  ZZ=RZ-Z(LZ)             ! the layers that has no roots  
                  RTO=ZZ/(Z(L1)-Z(LZ))    ! ratio that has no roots 
                  soilWaterRZ=soilWaterRZ+(soilWater(L1)-Wilt_Point(L1))*RTO  
                  potentialAvailWater=potentialAvailWater+RTO*(fieldCapacity(L1)-Wilt_Point(L1))
            END IF 
      END IF
    
      Layer_RD_Reach=Actual_SoilLayers  
      IF(maxLayer<Actual_SoilLayers) maxLayer=Actual_SoilLayers  
      LD1=1
      ! -------------XCC (NCC) control file ------------------------ 
      IF(NCC == 0)THEN   
            IF(Z(1)<.01)THEN
                  Z(1)=.01
            ELSE
                  IF(Z(1)>.01)THEN
                        Actual_SoilLayers=Actual_SoilLayers+1 
                        DO J=Actual_SoilLayers,2,-1  
                              Layer_ID(J)=Layer_ID(J-1)
                        END DO  
                        Layer_ID(1)=Actual_SoilLayers 
                        LD1=Actual_SoilLayers   
                        LORG(Actual_SoilLayers)=1  
                        RTO=.01/Z(1)  
                        CALL Soil_TopLayerMove(1,1,Actual_SoilLayers,0,RTO)     ! don't understand  
                        Z(Actual_SoilLayers)=.01 
                  END IF   
            END IF

            DO WHILE(Actual_SoilLayers<maxLayer)
                  L1=Layer_ID(1)  
                  ZMX=0.  
                  MXZ=2  
                  DO J=2,Actual_SoilLayers  
                        L=Layer_ID(J)  
                        ZZ=Z(L)-Z(L1)  
                        IF(ZZ>=minThickSplit)THEN
                              MXZ=J  
                              GO TO 130 
                        ELSE  
                              IF(ZZ>ZMX+.01)THEN
                                    ZMX=ZZ 
                                    MXZ=J
                              END IF 
                        END IF 
                        L1=L
                  END DO  
                  L=Layer_ID(MXZ) 
                  L1=Layer_ID(MXZ-1) 
              130 Actual_SoilLayers=Actual_SoilLayers+1  
                  CALL Soil_TopLayerMove(L,L1,Actual_SoilLayers,1,.5)   
                  DO J=Actual_SoilLayers,MXZ,-1  
                        Layer_ID(J)=Layer_ID(J-1) 
                  END DO  
                  Layer_ID(MXZ)=Actual_SoilLayers 
                  LORG(Actual_SoilLayers)=LORG(L) 
            END DO
      END IF   
END SUBROUTINE Soil_Process

!--------------------------- 1. Soil move------------------------------------------------
REAL FUNCTION Soil_Move(X,Y)
      ! EPIC1102
      ! THIS SUBPROGRAM IS CALLED BY ESLOS TO CALCULATE THE AMOUNT OF
      ! MATERIAL ADDED TO THE TOP LAYER AND REMOVED FROM THE SECOND LAYER.
      IMPLICIT NONE 
      REAL:: X, Y
     
      Soil_Move=X*Y
      X=X-Soil_Move
END FUNCTION Soil_Move
! 2 ------------------ Top layer movement --------------------------------------    
SUBROUTINE Soil_TopLayerMove(I,I1,K,L,RTO)
     
     IMPLICIT NONE
     ! local variables
     INTEGER:: I, I1, K, L
     REAL:: RTO
     
     IF(L>0)Z(K)=(Z(I)+Z(I1))*.5
      sorbRatioP(K)=sorbRatioP(I)
      bulkDensity(K)=bulkDensity(I)
      mineralBulkDensity(K)=mineralBulkDensity(I)
      clayFrac(K)=clayFrac(I)
      siltFrac(K)=siltFrac(I)
      sandFrac(K)=sandFrac(I)
      rockFrac(K)=rockFrac(I)
      satuateCondv(K)=satuateCondv(I)
      lateralFlowCondv(K)=lateralFlowCondv(I)
      fracNO3Leach(K)=fracNO3Leach(I)
      PH(K)=PH(I)
      soilTem(K)=soilTem(I)
      dryBulkDensity(K)=dryBulkDensity(I)
      flowCoeffP(K)=flowCoeffP(I)
      postTillBulkDensity(K)=postTillBulkDensity(I)
      totBase(K)=totBase(I)
      CEC(K)=CEC(I)
      gasO2Con(K)=gasO2Con(I)
      gasCO2Con(K)=gasCO2Con(I)
      gasN2OCon(K)=gasN2OCon(I)
      VGA(K)=VGA(I)
      VGN(K)=VGN(I)
      IF(FC_Method==1.OR.FC_Method==3.OR.FC_Method==5)CEM(K)=CEM(I)
      CaCO3(K)=CaCO3(I)
      ALS(K)=ALS(I)
      elecCondv(K)=elecCondv(I)
      WT(K)=Soil_Move(WT(I),RTO)
      NO3_N_Soil(K)=Soil_Move(NO3_N_Soil(I),RTO)
      SOP(K)=Soil_Move(SOP(I),RTO)
      NPassiveHumus(K)=Soil_Move(NPassiveHumus(I),RTO)
      NSlowHumus(K)=Soil_Move(NSlowHumus(I),RTO)
      NBiomass(K)=Soil_Move(NBiomass(I),RTO)
      NStructLitt(K)=Soil_Move(NStructLitt(I),RTO)
      NMetabLitt(K)=Soil_Move(NMetabLitt(I),RTO)
      CPassiveHumus(K)=Soil_Move(CPassiveHumus(I),RTO)
      CSlowHumus(K)=Soil_Move(CSlowHumus(I),RTO)
      CBiomass(K)=Soil_Move(CBiomass(I),RTO)
      structLitt(K)=Soil_Move(structLitt(I),RTO)
      metabLitt(K)=Soil_Move(metabLitt(I),RTO)
      lgStructLitt(K)=Soil_Move(lgStructLitt(I),RTO)
      CStructLitt(K)=.42*structLitt(K)
      CMetabLitt(K)=.42*metabLitt(K)
      CLgStructLitt(K)=.42*lgStructLitt(K)
      NLgStructLitt(K)=CStructLitt(K)-CLgStructLitt(K)
      cropResidu(K)=.001*(structLitt(K)+metabLitt(K))
      CStructLitt(I)=.42*structLitt(I)
      CMetabLitt(I)=.42*metabLitt(I)
      CLgStructLitt(I)=.42*lgStructLitt(I)
      NLgStructLitt(I)=CStructLitt(I)-CLgStructLitt(I)
      cropResidu(I)=.001*(structLitt(I)+metabLitt(I))
      SOC(I)=CBiomass(I)+CPassiveHumus(I)+CSlowHumus(I)+CMetabLitt(I)+CStructLitt(I)
      SOC(K)=CBiomass(K)+CPassiveHumus(K)+CSlowHumus(K)+CMetabLitt(K)+CStructLitt(K)
      SON(I)=NBiomass(I)+NPassiveHumus(I)+NSlowHumus(I)+NMetabLitt(I)+NStructLitt(I)
      SON(K)=NBiomass(K)+NPassiveHumus(K)+NSlowHumus(K)+NMetabLitt(K)+NStructLitt(K)
      labileP(K)=Soil_Move(labileP(I),RTO)
      activeMineralP(K)=Soil_Move(activeMineralP(I),RTO)
      FreshOrgP_Residu(K)=Soil_Move(FreshOrgP_Residu(I),RTO)
      stableMineralP(K)=Soil_Move(stableMineralP(I),RTO)
      solubleK(K)=Soil_Move(solubleK(I),RTO)
      exchangeK(K)=Soil_Move(exchangeK(I),RTO)
      fixK(K)=Soil_Move(fixK(I),RTO)
      EQKS(K)=Soil_Move(EQKS(I),RTO)
      EQKE(K)=Soil_Move(EQKE(I),RTO)
      saltWeight(K)=Soil_Move(saltWeight(I),RTO)
      Wilt_Point(K)=Soil_Move(Wilt_Point(I),RTO)
      fieldCapacity(K)=Soil_Move(fieldCapacity(I),RTO)
      Porosity(K)=Soil_Move(Porosity(I),RTO)
      soilWater(K)=Soil_Move(soilWater(I),RTO)
END SUBROUTINE Soil_TopLayerMove
 
!  3 --------  ESTIMATES DAILY AVEAGE TEMPERATURE AT THE CENTER OF EACH SOIL LAYER.     
SUBROUTINE Soil_Temperature   
      ! EPIC1102
      ! Simulate daily average soil temperature at the center of each layer, used for nutrient ccling and hydrology.
      ! Eqs reference: APEX doc - Eqs 269-273
      IMPLICIT NONE

      ! local variables
      REAL::  XLAG = 0.8, XLG1, F, DP, WW, B, WC, X4, XZ, X2, X3, TG, X1, X7, X5, ZZ,ZD, XX
      INTEGER:: J, ITMP = 1   ! questionable: if ITMP is initialized as 1 how could the section (=0) will be conducted?
 
      XLG1=1.-XLAG
      F=aveBulkDensity/(aveBulkDensity+686.*EXP(-5.63*aveBulkDensity)) ! Eq. 2.182 in EPIC-doc P43/ Eq. 270a APEX-doc P58  
      DP=1.+2.5*F                                       ! DP: damping depth
      WW=.356-.144*aveBulkDensity                      ! Eq. 270b in APEX-doc
      B=LOG(.5/DP)
      WC=.001*totSoilWater/(WW*Z(Layer_ID(Actual_SoilLayers))) ! totSoilWater: Water stored in the profile in mm (all layers)
      F=EXP(B*((1.-WC)/(1.+WC))**2)                     ! Eq. 2.184 in EPIC-doc/ Eq.270 in APEX-doc
      DD=F*DP
      IF(ITMP==0)THEN
            X4=.5*amplitudeT
            XZ=.5*(TMX-TMN)*soilRad/modelPara(94)
            X2=TX
            X3=(1.-coverLagingFactor)*X2+coverLagingFactor*soilTem(Layer_ID(2))
            TG=.5*(X2+X3)
            surfTem= yearAveTem +X4*COS((Date-200)/PIT)           ! surfTem: soil surface temperature
            X1=TX-surfTem
            DO J=1,Actual_SoilLayers
                  ISL=Layer_ID(J)
                  X7=modelPara(95)*Z(ISL)/DD
                  X5=EXP(-X7)
                  soilTem(ISL)=(yearAveTem +X4*COS((Date-200)/PIT-X7)*X5)+X1*X5
            END DO
      ELSE        
            XZ=.5*(TMX-TMN)*soilRad/modelPara(94)
             X2=TX+XZ  
            X3=(1.-coverLagingFactor)*X2+coverLagingFactor*soilTem(Layer_ID(2))  ! Eq.272
            surfTem=.5*(X2+X3)
            ZZ=2.*DD
            XX=0.
            X1=yearAveTem - surfTem
            DO J=1,Actual_SoilLayers
                  ISL=Layer_ID(J)
                  ZD=(XX+Z(ISL))/ZZ
                  F=ZD/(ZD+EXP(-.8669-2.0775*ZD))  ! F==FZ   Eq. 269a in APEX-doc 
                  soilTem(ISL)=XLAG*soilTem(ISL)+XLG1*(F*X1+surfTem)
                  soilElement(ISL)=0.
                  XX=Z(ISL)
            END DO
      END IF
END SUBROUTINE Soil_Temperature
                                
! 4 ------------------------ Percolation components----------------------
!  4.1 ------ cal percolation
SUBROUTINE cal_Percolation 
      ! EPIC1102
      ! THIS SUBPROGRAM COMPUTES PERCOLATION AND LATERAL SUBSURFACE FLOW FROM A SOIL LAYER WHEN FIELD CAPACITY IS EXCEEDED.
      ! Ref: APEX-doc  77-90
      IMPLICIT NONE
      !  local variables
      REAL:: AVW,X1 , X2, X3, ZZ, XZ
 
      Percolation=0.
      lateralFlow=0.
      AVW=soilWater(ISL)-fieldCapacity(ISL)
      
      IF(AVW<1.E-5)RETURN
      
      X1=24./(Porosity(ISL)-fieldCapacity(ISL)) ! porosity (mm)  , field capacity (mm)
      X2=lateralFlowCondv(ISL)*X1       
      ZZ=X1*satuateCondv(ISL)                  ! ZZ is TT in Eq. 2.30  SC: mm/h
      XZ=X2+ZZ
      IF(XZ>20.)THEN
            X3=AVW
      ELSE
            X3=AVW*(1.-EXP(-XZ))
      END IF
      Percolation=X3/(1.+X2/ZZ)                 ! Perlocation or vertical flow rate          mm/d
      lateralFlow=X3-Percolation               ! lateral subsurface or horizontal flow rate mm/d  
      IF(ISL/=IDR)lateralFlow=lateralFlow*groundWaterResidenceTime   ! IDR determine if it is a drainage system
END SUBROUTINE cal_Percolation
                                
! ---4.2  cal percolaltion  --------------                                
SUBROUTINE cal_Percolation2 
      ! EPIC1102
      ! THIS SUBPROGRAM COMPUTES PERCOLATION AND LATERAL SUBSURFACE FLOW FROM A SOIL LAYER BY &
      ! PASSING 4 mm SLUGS THROUGH AS A FUNCTION OF HYDRAULIC CONDUCTIVITY.
      IMPLICIT NONE
      ! local varialbes
      REAL:: ICW, AVW, POFC, X1, X5, X4, H, X2, ZZ, XZ, XX, X3, X6 

      Percolation=0.
      lateralFlow=0.
      ICW=0
      AVW=soilWater(ISL)-fieldCapacity(ISL)
      IF(AVW>0.)THEN
            POFC=Porosity(ISL)-fieldCapacity(ISL)
            X1=24./POFC
            DO WHILE(AVW>.01) 
                  X5=MIN(AVW/POFC,1.)
                  X4=MAX(1.E-5,X5**modelPara(93))
                  H=satuateCondv(ISL)*X4
                  X2=X1*lateralFlowCondv(ISL)*X4
                  ZZ=X1*H
                  XZ=X2+ZZ
                  XX=MIN(4.,AVW)
                  IF(XZ>20.)THEN
                        X3=XX
                  ELSE
                        X3=XX*(1.-EXP(-XZ))
                  END IF
                  X6=X3/(1.+X2/ZZ)
                  Percolation=Percolation+X6
                  lateralFlow=lateralFlow+X3-X6
                  AVW=AVW-4.
                  IF(ISL/=IDR)lateralFlow=lateralFlow*groundWaterResidenceTime
                  ICW=ICW+1
                  IF(ICW>50)EXIT
            END DO 
      END IF    
END SUBROUTINE cal_Percolation2
                                          
SUBROUTINE Soil_Percolation 
 
      ! EPIC1102
      ! THIS SUBPROGRAM IS THE MASTER PERCOLATION COMPONENT.  IT MANAGES THE ROUTING PROCESS
      IMPLICIT NONE
      ! local variables
      INTEGER:: KK, K, L1, J
      REAL, DIMENSION(5)::STPT = [.9,.75,.5,.25,.1]
      REAL:: POFC, ADD, XX, X1, X3, X2, RTO, T1, ZH, T2     
 
      ADD=0.
      Percolation=Rainfall-Runoff                   
      DO KK=1,Actual_SoilLayers
            ISL=Layer_ID(KK)
            soilWater(ISL)=soilWater(ISL)+Percolation      ! cause small problem when ISL = 4
            IF(waterTableHigt<=Z(ISL))THEN
                  subsurfLaterFlow(ISL)=0.
                  percolateFlow(ISL)=0.
                  Percolation=0.
            ELSE
                  IF(percolate_Method==0)THEN 
                        CALL cal_Percolation     
                  ELSE    
                        CALL cal_Percolation2 
                  END IF    
                  soilWater(ISL)=soilWater(ISL)-Percolation-lateralFlow   
                  subsurfLaterFlow(ISL)=lateralFlow
                  IF(ISL==IDR)THEN                                   ! IDR determines if it is a drainage system
                        SMM(18,MO)=SMM(18,MO)+lateralFlow             ! SMM(18,MO), VAR(18) : flow from drainage system
                        VAR(18)=lateralFlow
                  END IF
                  percolateFlow(ISL)=Percolation
                  ADD=ADD+lateralFlow
            END IF    
      END DO
      lateralFlow=ADD              ! totoal lateral flow
      L1=LD1
      DO K=Actual_SoilLayers,2,-1   ! back pass from bottom to surface to make sure porosity is not exceeded
            ISL=Layer_ID(K)           ! Ref: 90 a-d
            L1=Layer_ID(K-1)
            XX=soilWater(ISL)-fieldCapacity(ISL)
            IF(XX>0.)THEN
                  X1=VGN(L1)/(VGN(L1)-1.)
                  X2=soilWater(L1)-Wilt_Point(L1)
                  IF(X2>0.)THEN
                        RTO=(Porosity(L1)-Wilt_Point(L1))/X2
                        X3=RTO**X1-1.
                        IF(X3>0.)THEN
                              T1=(X3/VGA(L1)**VGN(L1))**(1./VGN(L1))
                              T1=T1/10.19
                        ELSE
                              T1=1.
                        END IF
                  ELSE
                        T1=1500.
                  END IF    
                  ZH=10.*(Z(ISL)-Z(L1))
                  X1=VGN(ISL)/(VGN(ISL)-1.)
                  X2=soilWater(ISL)-Wilt_Point(ISL)
                  IF(X2>0.)THEN
                        RTO=(Porosity(ISL)-Wilt_Point(ISL))/X2
                        X3=RTO**X1-1.
                        IF(X3>0.)THEN
                              T2=(X3/VGA(ISL)**VGN(ISL))**(1./VGN(ISL))
                              T2=T2/10.19
                        ELSE
                              T2=1.
                        END IF
                  ELSE
                        T2=1500.
                  END IF    
                  T2=T2+ZH
                  IF(T1<T2) CYCLE
                  X1=XX*MIN(.1,(T1-T2)/(T1+T2))
                  SLTP(ISL)=SLTP(ISL)+X1
                  soilWater(L1)=soilWater(L1)+X1
                  percolateFlow(L1)=percolateFlow(L1)-X1
                  soilWater(ISL)=soilWater(ISL)-X1
            END IF    
      END DO
      DO K=1,Actual_SoilLayers
            ISL=Layer_ID(K)
            IF(soilWater(ISL)>fieldCapacity(ISL))THEN
                  POFC=Porosity(ISL)-fieldCapacity(ISL)
                  RTO=(soilWater(ISL)-fieldCapacity(ISL))/POFC
                  DO J=1,5
                        IF(RTO>STPT(J))EXIT
                  END DO
                  NPT(J,ISL)=NPT(J,ISL)+1
            END IF    
      END DO
      IF(percolateFlow(L1)<0.)percolateFlow(L1)=0.
END SUBROUTINE Soil_Percolation

! ------------------------- 5. Add soil layer contents of water and nutrients ------------------------
SUBROUTINE BulkDensity_Change 
      ! EPIC1102
      ! THIS SUBPROGRAM SIMULATES THE CHANGE IN BULK DENSITY WITHIN THE
      ! PLOW LAYER CAUSED BY INFILTRATION SETTLING.
      IMPLICIT NONE 
      ! local variables
      INTEGER:: J
      REAL:: XX, F
      ! ////////////////////////////////////////////////////////////
      XX=Rainfall-Runoff
      DO J=1,Actual_SoilLayers
            ISL=Layer_ID(J)
            IF(XX>0..AND.postTillBulkDensity(ISL)<bulkDensity(ISL))THEN
                  ! Eq. 306b in APEX-doc
                  XX=XX*.2*(1.+2.*sandFrac(ISL)/(sandFrac(ISL)+EXP(8.597-.075*sandFrac(ISL))))/Z(ISL)**1.6
                  IF(XX<200.)THEN
                        F=XX/(XX+EXP(S_Curve(6,1)-S_Curve(6,2)*XX))
                  ELSE
                        F=1.
                  END IF
                  postTillBulkDensity(ISL)=postTillBulkDensity(ISL)+F*(bulkDensity(ISL)-postTillBulkDensity(ISL))
            END IF
            XX=percolateFlow(ISL)
            IF(Z(ISL)>BIG)EXIT
      END DO

END SUBROUTINE BulkDensity_Change
!------------------------------------------------- add soil water and nutrient in each layer--------------------------------------
SUBROUTINE Sum_WaterNutrient(KWX)
      ! THIS SUBPROGRAM ADDS SOIL LAYER CONTENTS OF WATER & NUTRIENTS
      IMPLICIT NONE
 
      ! local variables
      INTEGER:: J, K1, K2, L, LZ, L1, KWX, KK 
      REAL:: ATX(4),SWZ(4), XX, SUM, S1, S2, NTX, BL, RTO, SW10CM, ZZ  
      ! Initialize variables 
      ! Soil Water 
      totSoilWater=0.; soilWaterRZ=0.; potentialAvailWater=0.; Tot_Water_RootZone=0.; plawDepSoilWater=0.
      plawDepFieldCapa=0.
      ! N, P, K, C -- soil
      totNO3_N=0.; ZN2O=0.; ZN2OG=0.; ZN2OL=0.; ZNO2=0.; ZNH3=0.
      TOP=0.; totLabileP=0.; initialTotSOP=0.; Tot_FreshOrgP=0.; totActiveP=0.; plawDepLabileP=0. 
      totSolubleK=0.; totExchangeK=0.; totFixK=0.            
      TNOR=0.; totSalt=0.; Tot_Salt_RootZone=0.            
      ! Organic matterials 
      totCropResi=0.; totStructLitt=0.; totMetabLitt=0.; totLgStructLitt=0.
      ! Organic N and C
      totSOC=0.; totCStructLitt=0.; totCMetabLitt=0.; totCLgStructLitt=0.; totNLgStructLitt=0.
      totCBiomass=0.; totCSlowHumus=0.; totCPassiveHumus=0.; plawDepSOC=0.
      totSON=0.; totNBiomass=0.; totNStructLitt=0.; totNMetabLitt=0.; totNSlowHumus=0.; totNPassiveHumus=0.
 
      XX=0.
      IF(BIG>0.)CALL BulkDensity_Change 
      SUM=0.
      S1=0.
      S2=0.
      NTX=0 
      BL=ZCS(1)
      K1=1
      K2=LD1
      DO J=1,Actual_SoilLayers
            L=Layer_ID(J)
            ZN2O=ZN2O+WN2O(L)
            ZNO2=ZNO2+NO2Weight(L)
            ZN2OG=ZN2OG+WN2OG(L)
            ZN2OL=ZN2OL+Weight_N2OL(L)
            totNO3_N=totNO3_N+NO3_N_Soil(L)
            ZNH3=ZNH3+NH3_Weight(L)
            totExchangeK=totExchangeK+exchangeK(L)
            totFixK=totFixK+fixK(L)
            totSolubleK=totSolubleK+solubleK(L)
            totSalt=totSalt+saltWeight(L)
            SOC(L)=CBiomass(L)+CPassiveHumus(L)+CSlowHumus(L)+CMetabLitt(L)+CStructLitt(L)
            elecCondv(L)=.15625*saltWeight(L)/soilWater(L)
            IF(Z(L)<=plawDepth)THEN
                  plawDepSoilWater=plawDepSoilWater+soilWater(L)-Wilt_Point(L)
                  plawDepFieldCapa=plawDepFieldCapa+fieldCapacity(L)-Wilt_Point(L)
                  plawDepLabileP=plawDepLabileP+labileP(L)
                  plawDepSOC=plawDepSOC+SOC(L)
                  SUM=SUM+WT(L)
                  K1=J
                  K2=L
            END IF
            IF(Z(L)<=RZ)THEN
                  Tot_Water_RootZone=Tot_Water_RootZone+soilWater(L)
                  TNOR=TNOR+NO3_N_Soil(L)
                  Tot_Salt_RootZone=Tot_Salt_RootZone+saltWeight(L)
                  LZ=L
                  L1=J
                  soilWaterRZ=soilWaterRZ+soilWater(L)-Wilt_Point(L)
                  potentialAvailWater=potentialAvailWater+fieldCapacity(L)-Wilt_Point(L)
            END IF
            totLabileP=totLabileP+labileP(L)
            TOP=TOP+stableMineralP(L)
            totActiveP=totActiveP+activeMineralP(L)
            initialTotSOP=initialTotSOP+SOP(L)
            totCropResi=totCropResi+cropResidu(L)
            totSoilWater=totSoilWater+soilWater(L)
            Tot_FreshOrgP=Tot_FreshOrgP+FreshOrgP_Residu(L)
            totStructLitt=totStructLitt+structLitt(L)
            totMetabLitt=totMetabLitt+metabLitt(L)
            totLgStructLitt=totLgStructLitt+lgStructLitt(L)
            totCStructLitt=totCStructLitt+CStructLitt(L)
            totCMetabLitt=totCMetabLitt+CMetabLitt(L)
            totCLgStructLitt=totCLgStructLitt+CLgStructLitt(L)
            totNLgStructLitt=totNLgStructLitt+NLgStructLitt(L)
            totCBiomass=totCBiomass+CBiomass(L)
            totCSlowHumus=totCSlowHumus+CSlowHumus(L)
            totCPassiveHumus=totCPassiveHumus+CPassiveHumus(L)
            totNStructLitt=totNStructLitt+NStructLitt(L)
            totNMetabLitt=totNMetabLitt+NMetabLitt(L)
            totNBiomass=totNBiomass+NBiomass(L)
            totNSlowHumus=totNSlowHumus+NSlowHumus(L)
            totNPassiveHumus=totNPassiveHumus+NPassiveHumus(L)
            SON(L)=NBiomass(L)+NPassiveHumus(L)+NSlowHumus(L)+NMetabLitt(L)+NStructLitt(L)
            IF(KFL(18)==0.OR.KWX>0)CYCLE
            ! DETERMINE SOIL TEMP AT .1, .2, .5, 1. m
            RTO=.001*(soilWater(L)-Wilt_Point(L))/(Z(L)-XX)
            IF(XX<.1.AND.Z(L)>=.1)THEN
                  ATX(1)=soilTem(L)
                  SWZ(1)=RTO
                  SW10CM=.1*soilWater(L)/(Z(L)-XX)
            END IF
            IF(XX<.2.AND.Z(L)>=.2)THEN
                  ATX(2)=soilTem(L)
                  SWZ(2)=RTO
            END IF
            IF(XX<.5.AND.Z(L)>=.5)THEN
                  ATX(3)=soilTem(L)
                  SWZ(3)=RTO
            END IF
            IF(XX<1..AND.Z(L)>=1.)THEN
                  ATX(4)=soilTem(L)
                  SWZ(4)=RTO
            END IF 
            RTO=1000.*RTO
            ! COMPUTE CONTROL SECTION WATER CONTENT
            IF(Z(L)>ZCS(1).AND.BL<ZCS(3))THEN
                  IF(Z(L)<=ZCS(2))THEN
                        S1=S1+RTO*(Z(L)-BL)
                  ELSE
                        IF(BL<ZCS(2))THEN      
                        S1=S1+RTO*(ZCS(2)-BL)
                        BL=ZCS(2)
                        S2=S2+RTO*(Z(L)-BL)        
                        ELSE
                              IF(Z(L)<=ZCS(3))THEN
                                    S2=S2+RTO*(Z(L)-BL)
                              ELSE
                                    S2=S2+RTO*(ZCS(3)-BL)
                              END IF                                                    
                        END IF
                  END IF
                  BL=Z(L)
            END IF        
            XX=Z(L)
      END DO
      IF(KFL(18)>0.AND.KWX==0) WRITE(KW(18),554)IY,IYR,MO,DayOfMon,S1,S2,SWZ,ATX,SW10CM
      totSON=totNStructLitt+totNMetabLitt+totNBiomass+totNSlowHumus+totNPassiveHumus
      totSOC=totCStructLitt+totCMetabLitt+totCBiomass+totCSlowHumus+totCPassiveHumus
      VAR(85)=TWN0-totSON-sedimentLossN
      SMM(85,MO)=SMM(85,MO)+TWN0-totSON-sedimentLossN
      TWN0=totSON
      IF(LZ/=L)THEN
            ZZ=RZ-Z(LZ)
            L1=Layer_ID(L1+1)
            RTO=ZZ/(Z(L1)-Z(LZ))
            soilWaterRZ=soilWaterRZ+(soilWater(L1)-Wilt_Point(L1))*RTO
            potentialAvailWater=potentialAvailWater+RTO*(fieldCapacity(L1)-Wilt_Point(L1))
            Tot_Water_RootZone=Tot_Water_RootZone+soilWater(L1)*RTO
            TNOR=TNOR+NO3_N_Soil(L1)*RTO
            Tot_Salt_RootZone=Tot_Salt_RootZone+saltWeight(L1)*RTO
      END IF
      SSW=SSW+soilWaterRZ
      IF(K1/=Actual_SoilLayers)THEN
            KK=Layer_ID(K1+1)
            RTO=(plawDepth-Z(K2))/(Z(KK)-Z(K2))
            plawDepSoilWater=plawDepSoilWater+RTO*(soilWater(KK)-Wilt_Point(KK))
            plawDepFieldCapa=plawDepFieldCapa+RTO*(fieldCapacity(KK)-Wilt_Point(KK))
            plawDepLabileP=plawDepLabileP+RTO*labileP(KK)
            plawDepSOC=plawDepSOC+RTO*SOC(KK)
            SUM=SUM+RTO*WT(KK)
      END IF
      APBC=1000.*plawDepLabileP/SUM
      plawDepSOC=.001*plawDepSOC
 
  554 FORMAT(I4,1X,3I4,2F10.1,4F10.4,5F10.2)
END SUBROUTINE Sum_WaterNutrient
    
! ------------------------- 6. Wilt point saturate ---------------------------------
SUBROUTINE WiltPoint_Sat(I)
      !     EPIC1102
      IMPLICIT NONE
      ! local variables
      INTEGER:: I
      REAL:: X1, X2

      X1=.95*Porosity(I)
      IF(fieldCapacity(I)>X1)THEN
            X2=fieldCapacity(I)-Wilt_Point(I)
            fieldCapacity(I)=X1
            Wilt_Point(I)=fieldCapacity(I)-X2
            IF(Wilt_Point(I)<=0.)Wilt_Point(I)=.01*fieldCapacity(I)
      END IF    
END SUBROUTINE WiltPoint_Sat

! --------------------------- 7. Compute soil water contents ------------------------
! ----------------------------7.1 use Walter Rawls's method -------------------------
SUBROUTINE Cal_SoilWaterW(CL,SA,OC,OrgP_ini,FC_es)
      ! EPIC1102
      !THIS SUBPROGRAM USES WALTER RAWLS'S METHOD FOR ESTIMATING SOIL
      !WATER CONTENT AT 33 AND 1500 kpa.
      IMPLICIT NONE
      ! local variables
      REAL:: CL,SA,OC,OrgP_ini,FC_es
      OrgP_ini=.026+.005*CL+.0158*OC
      FC_es=.2576-.002*SA+.0036*CL+.0299*OC
END SUBROUTINE Cal_SoilWaterW
! ---------------------------- 7.2 use Otto Baumer's method --------------------------
SUBROUTINE Cal_SoilWaterO(CM,CL,OC,SA,WC15,WC3RD,ZZ)
      !EPIC1102
      !THIS SUBPROGRAM USES OTTO BAUMER'S METHOD FOR ESTIMATING SOIL
      !WATER CONTENT AT 33 AND 1500 kpa.
      IMPLICIT NONE

      ! local variables
      REAL:: SI, CL, SA, OC, OM, ZZ, BV, BW, BO, APD, CE, CM, CA, X1, X2, CAAF, VOMO, VODF, WC15G, &
       VFS, SICL, CF1, A1, BDX, SDF, CF3, SAF, W3RDG, OAIR, WSV, DBD, DBS, WSG, WC3RD, WC15 

      SI=100.-CL-SA
      OM=1.72*OC
      IF(ZZ<.25)THEN
          BV=1.085
          BW=1.035
          BO=1.9
      ELSE
            IF(ZZ<1.)THEN
                  BV=1.
                  BW=1.
                  BO=1.
            ELSE
                  BV=.915
                  BW=.96
                  BO=.1      
            END IF    
      END IF          
      APD=100./(37.74+.3366*OM)
      CE=CM+2.428*OC+1.7*ZZ
      CA=MIN(.8,CE/CL)
      X1=CL*CL
      X2=CA*CA
      IF(CL<10.)THEN
            CAAF=.099*X2*X1
      ELSE
            CAAF=9.9*X2
      END IF
      VOMO=BV*(42.84+1.37*OM+.00294*X1+CAAF+.0036*X1*X2-.0852*SA-.316*CL*CA)
      VODF=MAX(.01,BV*(.277+.16*OM-2.72*CA*X2+.268*CL*CA+.00546*X1*CA&
      -.00184*X1*X2))
      WC15G=.71+.45*OM+.336*CL+.117*CL*CA**1.5
      VFS=.1*SA
      SICL=CL+.3333*(SI+VFS)
      IF(SICL<15.)THEN
            CF1=1.
            GO TO 5
      END IF
      IF(SICL<30.)THEN
            CF1=2.-.0667*SICL
      ELSE
            CF1=0.
      END IF    
    5 A1=14.94+3.8*X2-.137*SA
      BDX=APD*(1.-.01*VOMO)
      SDF=SA-VFS
      CF3=37.74*SDF/((100.-SDF)/BDX+.3774*SDF)
      SAF=1.-.005*CF3*(CF1+1.)
      W3RDG=BW*(A1*SAF+WC15G+.746*OM)
      OAIR=BO*(3.8+.00018*X1-.03365*SA+.126*CL*CA+.000025*OM*SA*SA)
      WSV=VOMO*(1.-.01*OAIR)
      WSG=WSV/BDX
      DBD=APD*(1.-.01*(VOMO-VODF))
      DBS=(BDX-DBD)/WSG
      WC3RD=.01*W3RDG*(DBD+DBS*W3RDG)
      WC15=.01*WC15G*(DBD+DBS*WC15G)
END SUBROUTINE Cal_SoilWaterO
! ---------------------------- 7.3 -------------------------------------------------
SUBROUTINE Cal_SoilWaterBNW(CL,SI,SA,OC,bulkDen, WPs, FCs)
      !EPIC1102
      !THIS SUBPROGRAM USES THE BEHRMAN-NORFLEET-WILLIAMS (BNW) METHOD
      !FOR ESTIMATING SOIL WATER CONTENT AT 33 AND 1500 kpa.
      IMPLICIT NONE
 
      ! local variables
      REAL:: CL, SI, SA, OC, X1, AD1, BDO, bulkDen, PAO, RTO, PAWs, WPs, FCs

      WPs=.04285+.0001041*SI+.003958*CL+.00001555*CL*SI-.005606*LOG10(OC)
      X1=1.72*OC
      AD1=SA+SI+CL+X1
      BDO=(1.6*SA+1.3*SI+1.1*CL+.224*X1)/AD1
      PAO=(.05*SA+.26*SI+.08*CL+.9*X1)/AD1
      RTO=BDO/bulkDen
      PAWs=RTO*PAO
      FCs=WPs+PAWs
END SUBROUTINE Cal_SoilWaterBNW
! ---------------------------- 7.4 estimate soil water contents ---------------------
SUBROUTINE Cal_SoilWater 
      ! THIS SUBPROGRAM COMPUTES THE SOIL WATER CONTENT(m/m) AND NO3 
      !CONCENTRATION(G/M^3) OF A SOIL AT 0.15 AND 0.3 M DEPTHS. 
      IMPLICIT NONE
 
      ! local variables
      INTEGER:: I, L, L1, J
      REAL:: WT15, WT30, RTO, X1, Y1, Y2, W1

      SW15=0.
      SW30=0.
      SNN15=0.
      SNN30=0.
      SNA15=0.
      SNA30=0.
      WT15=0.
      WT30=0.
      DO I=1,Actual_SoilLayers
            L=Layer_ID(I)
            IF(Z(L)>.15)GO TO 2
            SW15=SW15+soilWater(L)
            SNN15=SNN15+NO3_N_Soil(L)
            SNA15=SNA15+NH3_Weight(L)
            WT15=WT15+WT(L)
            L1=L
      END DO
      GO TO 5
    2 RTO=(.15-Z(L1))/(Z(L)-Z(L1))   
      X1=RTO*soilWater(L)
      SW15=SW15+X1
      Y1=RTO*NO3_N_Soil(L)
      SNN15=SNN15+Y1
      Y2=RTO*NH3_Weight(L)
      SNA15=SNA15+Y2
      W1=RTO*WT(L)
      WT15=WT15+W1
      DO J=1,Actual_SoilLayers
            L=Layer_ID(J)
            IF(Z(L)>.3)GO TO 4
            SW30=SW30+soilWater(L)
            SNN30=SNN30+NO3_N_Soil(L)
            SNA30=SNA30+NH3_Weight(L)
            WT30=WT30+WT(L)
            L1=L
      END DO
      GO TO 6
    4 RTO=(.3-Z(L1))/(Z(L)-Z(L1))
      SW30=SW30+RTO*soilWater(L)
      SNN30=SNN30+RTO*NO3_N_Soil(L)
      SNA30=SNA30+RTO*NH3_Weight(L)
      WT30=WT30+RTO*WT(L)
    6 SW30=SW30/300.
      SNN30=1000.*SNN30/WT30
      SNA30=1000.*SNA30/WT30
    5 SW15=SW15/150.
      SNN15=1000.*SNN15/WT15
      SNA15=1000.*SNA15/WT15
 
END SUBROUTINE Cal_SoilWater

! --------------------------- 8. Estimate fieldCapacity and WP --------------------------------
! ------------------ 8.1 Sort numbers -----------------------------------
SUBROUTINE Shell_Sort(D, NXX, M1)
      IMPLICIT NONE
      ! local variables       potential problems with local variables (SAVE attribute)
      INTEGER:: M1, NXX(M1), M, K, I, J, N1
      REAL:: D(M1) 
      ! ///////////////////////////////////////////////////////////   
      M=1
      DO
            IF((2**M)>M1)EXIT
            M=M+1
      END DO
      M=M-1
      M=2**M
      
      DO
            K=M1-M
            DO I=1,K
                  DO J=I,1,-M
                        IF(D(NXX(J+M))>=D(NXX(J)))EXIT
                              N1=NXX(J)
                              NXX(J)=NXX(J+M)
                              NXX(J+M)=N1
                  END DO
            END DO
            M=M/2
            IF(M<=0)EXIT
      END DO
END SUBROUTINE Shell_Sort
! ---------------------------8.2 ------------------------------
SUBROUTINE Cal_FC_WP(CL,SA,OC,W1,F3)
      !     EPIC1102
      !     THIS SUBPROGRAM ESTIMATES fieldCapacity AND WP USING NEAREST NEIGHBOR 
      !     APPROACH.
      IMPLICIT NONE
 
      ! local variables potential problems with local variables (SAVE attribute XTP)
      INTEGER:: KX(NSX), I1, I, K, N1, J  
      REAL:: DXS(NSX), XTP(3), SUM, TOT, W1, F3, SA, CL, OC 

      XTP(1)=(SA-XAV(1))/XDV(1)
      XTP(2)=(CL-XAV(2))/XDV(2)
      XTP(3)=(OC-XAV(3))/XDV(3)
      XTP(1)=XTP(1)*BRNG/XRG(1)
      XTP(2)=XTP(2)*BRNG/XRG(2)
      XTP(3)=XTP(3)*BRNG/XRG(3)
      I1=1
      DO I=1,NSX
            DXS(I)=0.
            DO K=1,3
                  DXS(I)=DXS(I)+(XSP(I,K)-XTP(K))**2
            END DO
            DXS(I)=SQRT(DXS(I))
            KX(I1)=I
            I1=I1+1
      END DO
      I1=I1-1
      CALL Shell_Sort(DXS,KX,NSX)
      N1=MIN(I1,NSNN)
      SUM=0.
      DO I=1,N1
            J=KX(I)
            DO K=1,3
                  XTP(K)=XDV(K)*XSP(J,K)*XRG(K)/BRNG+XAV(K)
            END DO
            SUM=SUM+DXS(J)
      END DO
      TOT=0.
      DO I=1,N1
            J=KX(I)
            DXS(J)=(SUM/DXS(J))**EXNN
            TOT=TOT+DXS(J)
      END DO
      W1=0.
      F3=0.
      DO I=1,N1
            J=KX(I)
            DXS(J)=DXS(J)/TOT
            W1=W1+DXS(J)*XSP(J,4)
            F3=F3+DXS(J)*XSP(J,5)
      END DO
END SUBROUTINE Cal_FC_WP 
! ------------------------------9. distribute NPK ------------------------------
SUBROUTINE Distribute_NPK(X,DG,DG1,X1,X2,I, L)
      !     EPIC1102  potential problems with local variables (SAVE attribute X)
      IMPLICIT NONE
      INTEGER:: I, L
      REAL::X(L), DG, DG1, X1, X2, RTO, XX
      !//////////////////////////////////////////////////
      IF(X(I)>0.)RETURN
      IF(I==1)THEN
          X(1)=X1
          RETURN
      ELSE
            XX=X2*DG
            IF(XX>10.)THEN
                  RTO=.0001
            ELSE
                  RTO=DG*EXP(-XX)/DG1
            END IF
            X(I)=X(I-1)*RTO
      END IF 
END SUBROUTINE Distribute_NPK

! ----------------------- 10.1 Prepare Soil Variables for Output --------------------
SUBROUTINE Prep_SoilVar(YTP)
      !     EPIC1102
      !     THIS SUBPROGRAM PREPARES SOIL DATA TABLE FOR OUTPUT, AND
      !     CONVERTS WEIGHTS OF MATERIALS TO CONCENTRATION.
      IMPLICIT NONE
     
      ! local variabels  potential problems with local variables (SAVE attribute YTP)
      INTEGER:: I, J
      REAL::YTP(16), XX, WT1, DG, X1, X2 
      !/////////////////////////////////////////////////////////
      XX=0.
      YTP(1)=0.
      YTP(2)=0.
      YTP(3)=0.
      YTP(4)=0.
      YTP(5)=0.
      DO J=1,Actual_SoilLayers
            I=Layer_ID(J)
            WT1=WT(I)/1000.
            SOIL(1,I)=labileP(I)/WT1
            SOIL(2,I)=activeMineralP(I)/WT1
            SOIL(3,I)=stableMineralP(I)/WT1
            SOIL(4,I)=SOP(I)/WT1
            SOIL(5,I)=NO3_N_Soil(I)/WT1
            SOIL(14,I)=solubleK(I)/WT1
            SOIL(15,I)=exchangeK(I)/WT1
            SOIL(16,I)=fixK(I)/WT1
            SOIL(6,I)=SON(I)/WT1
            SOIL(7,I)=.1*SOC(I)/WT(I)
            DG=(Z(I)-XX)*1000.
            X1=soilWater(I)-Wilt_Point(I)
            X2=fieldCapacity(I)-Wilt_Point(I)
            SOIL(17,I)=X1/X2
            SOIL(18,I)=X1
            SOIL(19,I)=X2
            elecCondv(I)=.15625*saltWeight(I)/soilWater(I)
            SOIL(20,I)=Wilt_Point(I)/DG
            YTP(1)=YTP(1)+fieldCapacity(I)
            SOIL(9,I)=fieldCapacity(I)/DG
            SOIL(8,I)=Porosity(I)/DG
            SOIL(13,I)=dryBulkDensity(I)*bulkDensity(I)
            YTP(2)=YTP(2)+soilWater(I)
            SOIL(12,I)=soilWater(I)/DG
            YTP(3)=YTP(3)+Wilt_Point(I)
            YTP(4)=YTP(4)+Porosity(I)
            YTP(5)=YTP(5)+WT(I)
            XX=Z(I)
      END DO
      XX=Z(Layer_ID(Actual_SoilLayers))*1000.
      DO I=1,4
            YTP(I)=YTP(I)/XX
      END DO
 
END SUBROUTINE Prep_SoilVar
! -------------------------- 2.1 Bulk density effect on root growth ------------------------------
SUBROUTINE BulkDens_Effect(X1,X3,F,J,II )
      ! THIS SUBPROGRAM ESTIMATES THE EFFECT OF BULK DENSITY ON ROOT
      ! GROWTH AND SATURATED CONDUCTIVITY AS A FUNCTION OF TEXTURE      
      ! function: recalcuate Saturate_condv and return it

      ! local variables 
      INTEGER:: J, II
      REAL:: BD1, BD2, X2, X1, X3, B1, B2, F,XC

      BD1=X3+.00445*sandFrac(J)             ! Eq. 299b c d in APEX doc P66-67 / Eq. 2.240 EPIC-doc  P57
      BD2=X3+.35+.005*sandFrac(J)           ! X3: threshold buck density for root stress for a soil of zero sand constent
      X2=LOG(.01124*BD1)
      B2=(X2-LOG(8.*BD2))/(BD1-BD2)
      B1=X2-B2*BD1
      X2=B1+B2*X1
      IF(X2<-10.)THEN
            F=1.            
            GO TO 6
      END IF
      IF(X2>10.)THEN
            F=.0001
      ELSE
            F=X1/(X1+EXP(X2))
      END IF
    6 IF(II>2)RETURN
      IF(soilSaturate_Method==0.AND.satuateCondv(J)>0.)RETURN   
      XC=100.-clayFrac(J)
      satuateCondv(J)=12.7*XC/(XC+EXP(11.45-.097*XC))+1. !Eq. 92 in APEX-doc P27 estimate the saturated conductivity
      satuateCondv(J)=F*satuateCondv(J)                  ! Ref: Eq. 2.31 in EPIC-doc, F might be soil strength stress factor

END SUBROUTINE

END MODULE Soil_Module