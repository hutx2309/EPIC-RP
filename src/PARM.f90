MODULE PARM
      IMPLICIT NONE
      
      real, parameter:: PIT = 58.13
      real, parameter:: PI2 = 6.283185
      real, parameter:: CLT = 57.296   !CLT: 1/365/2pi for calcuate daylength See EPIC-doc P45
      integer, parameter:: MFT=20                                                                          
      integer, parameter:: MSL=15    ! number of maximum soil layers
      integer, parameter:: MSC=31                                                                              
      integer, parameter:: NSM=114   ! number of variables                                                                            
      integer, parameter:: MPS=40                                                                         
      integer, parameter:: MNC=30                                                                              
      integer, parameter:: MRO=150                                                                        
      integer, parameter:: MNT=100                                                                        
      integer, parameter:: NGF= 38   !origional model: MSO+5 control Id for gis file
 
      CHARACTER(4)::TITLE(60),SID(16)
      CHARACTER(LEN = 16)::HEDS(30),HEDC(20),HEDP(10)
      CHARACTER(9)::soilOrders(3)                     
      CHARACTER(80)::siteName,FWTH,OPSCFILE,SITEFILE,SOILFILE,soilOrder                                                                                   
      CHARACTER(4),DIMENSION(:),ALLOCATABLE::Crop_Name                             
      CHARACTER(8),DIMENSION(:),ALLOCATABLE::soilHorizon,fertName,equipmentName                                
      CHARACTER(LEN = 16),DIMENSION(:),ALLOCATABLE::pestName
      CHARACTER(LEN = 4),DIMENSION(:),ALLOCATABLE::varName
      Character(4), dimension(150)::simuYears
      
      ! ============================= Simulation Control =================================================
      integer, dimension(:),allocatable:: KFL   ! Control the output files
      integer:: KW(200), NFL, KR(30)            ! units of writing files
      integer:: IGIS                            ! GIS file ID
      integer:: NOP, IPYI, printInterval  ! print settings
      
      !  Date-related variables  
      integer:: NC(13)
      integer:: IYX, IY, IYER,& ! number of years with respect to 1880 (for dynamic CO2 cal)
                runMode, leapYr, IRO0, IRUN,  dateNum,  &
                IYR, IMON, IDAY, Hour, Mint, Sec, &  ! for timer: current year, month, day, hour, min, sec
                DayOfYear, MO, IBD, MSO, DayOfMon, cropRotationYrs, IGSD, IWIX(12),Date,NumOfDays,NWI,ISIX(6)  
      ! IBD = DOY
      !  loop control variables  
      integer:: JT1, & ! month with the minimum monthly average T observed (from stations)
                JT2, IAUF, IAUL,JD, ICCD,IRT, CIGZ,KTT, JCN
 
      ! ******************* 1. EPIC CONT.DAT: Simulation Control Variables *************************************************   
      integer:: numSimuYear,        & ! number of years of simulation 
                Year0,Month0,Day0,  & ! beginning year, month, and day of month of simulation 
                printChoice,        & ! = N1 FOR ANNUAL PRINTOUT                              | N YEAR INTERVA             
                                     !          = N2 FOR ANNUAL WITH SOIL TABLE              | N=0 SAME AS                
                                     !          = N3 FOR MONTHLY                             | N=1 EXCEPT                 
                                     !          = N4 FOR MONTHLY WITH SOIL TABLE             | N=0 PRINTS                 
                                     !          = N5 FOR MONTHLY WITH SOIL TABLE AT HARVEST  | OPERATIONS                 
                                     !          = N6 FOR N DAY INTERVAL                                                   
                                     !          = N7 FOR SOIL TABLE ONLY N DAY INTERVAL                                   
                                     !          = N8 FOR N DAY INTERVAL RAINFALL DAYS ONLY                                
                                     !          = N9 FOR N DAY INTERVAL DURING GROWING SEASON 
                weatherVarID,        & ! ID NUMBER OF WEATHER VARIABLES INPUT.  
                                     ! RAIN=1,  TEMP=2, RAD=3,  WIND SPEED=4, RH=5.  
                                     ! IF ANY VARIABLES ARE INP RAIN MUST BE INCLUDED. THUS, IT IS NOT NECESSARY TO SPECIF             
                                     ! ID=1 UNLESS RAIN IS THE ONLY INPUT VARIABLE.                            
                                     ! LEAVE BLANK IF ALL VARIABLES ARE GENERATED.  EXAMPLES: NGN=1 INPUTS RAIN.                                                      
                                                                                            ! NGN=23 INPUTS RAIN, TEMP, AND RAD.                                      
                                                                                            ! NGN=2345 INPUTS ALL 5 VARIABLES. 
                randomCycles,      & ! number times random number generation cycles before simulation starts 
                weatherInputMode,  & ! 0: normal operation of weather models; N Number years input weather before rewinding (used
                                     !                             for real time simulation)
                leapYrFlag,        & ! 0: leap year is considered; 1: leap year not considered
                ET_Method,         & ! Point Evaporation Method: 1, PENMAN-MONTEITH; 2, PENMAN; 3, RIESTLEY-TAYLOR; 4, HARGREAVES;
                                     !                           5, BAIER-ROBERTSON
                curveNumFlag,      & ! 0: Stochastic curve number estimator; >0 : Rigid curve number estimator
                peakRateMethod,    & ! 0: modified rational Eq. peak rate estimate; >0 SCS TR55 peak rate estimate(1: type 1 rainfall pattern;
                                     !                                   2: type 1A; 3: type 2; 4: type 3 ) 
                erosionMode,       & ! 0: normal erosion of soil profile; 1: static soil profile
                operationMode,     & ! 0: normal operation; 1: automatic heat unit shcedule (potentialHeatUnit must be input at planting)
                rootRespirFlag,    & ! 0: root respiration in NDNITCI; >0 : no root respirtation in NDNITCI
                CN_Method,         & ! 0: variable daily CN with depth soil water weigting; 
                                     ! 1: variable daily CN without depth weighting;
                                     ! 2: variable daily CN linear CN/totSoilWater, no depth weighting;
                                     ! 3: non-varying CN (CN2 used for all storms)
                                     ! 4: variable daily CN SMI (soil moisture index)  
                runoff_Method0,    & ! 0: Q estimated by CN; 1: for Green &AMPT estimate of Q, RF EXP DST, peak rainfall rate simualated;
                                     ! 2: G&A Q, rainfall Expotentional distribution, peak rainfall input; 
                                     ! 3: G&A Q, rainfall uniformally distribution, peak rainfall input  
                massPestFlag,      & ! <0: mass only, no pesticide in and out; =0: mass only with pesticide in and out; 
                                     ! >0: for pesticide and nutrients output in mass and concentration
                soilP_runoff_Method, & ! =0: for soluiable P runoff estimate using Gleams pesticide approach; >0: modified nonlinear appraoch
                DOY_realtime,      & ! real time day of year
                numSeedInitialized,  & ! number of times generator seeds are initialized for a site            
                Enrich_Method,  & ! 0: EPIC enrichment ratio method; 1: Gleams enrichemnt ratio method  
                biomsConver_Method,& ! 0: traditional EPIC radiation to biomass; >0: new experimental water use to biomass
                limeFlag,          & ! 0: apply lime; 1: don't apply lime
                erosionC_factor,   & ! 0: RUSLE C factor for all erosion Eqs.; >0: EPIC C factor for all erosion eqs
                FC_Method,         & != 0 FIELD CAP/WILTING PT EST RAWLS METHOD DYNAMIC.                        
                                     ! = 1 FIELD CAP/WILTING PT EST BAUMER METHOD DYNAMIC.                       
                                     ! = 2 FIELD CAP/WILTING PT INP RAWLS METHOD DYNAMIC.                        
                                     ! = 3 FIELD CAP/WILTING PT INP BAUMER METHOD DYNAMIC.                       
                                     ! = 4 FIELD CAP/WILTING PT EST RAWLS METHOD STATIC.                         
                                     ! = 5 FIELD CAP/WILTING PT EST BAUMER METHOD STATIC.                        
                                     ! = 6 FIELD CAP/WILTING PT INP STATIC.                                      
                                     ! = 7 FIELD CAP/WILTING PT NEAREST NEIGHBOR DYNAMIC                         
                                     ! = 8 FIELD CAP/WILTING PT NEAREST NEIGHBOR STATIC
                                     ! = 9 FIELD CAP/WILTING PT BEHRMAN-NORFLEET-WILLIAMS (BNW) DYNAMIC                         
                                     ! =10 FIELD CAP/WILTING PT BEHRMAN-NORFLEET-WILLIAMS (BNW) STATIC
                weatherVarForRuns, & ! 0: normal runs with daily weather input;                
                                     ! >0: continuous daily weather from run to run (no rewind)  
                 CO2_Method,       & ! 0: constant CO2; 1: dynamic CO2; 2: input CO2 
                 N_vol_Method,     & ! 0: original EPIC NITVOL Eqs; >0: revised NITVOL Eqs                 
                 correctionForRuns,& ! 0: normal runs; = DOY when weather correction to simulate input no means stops 
                 deNitri_Method,   & ! 1: orgional EPIC denitrification; 2: Armen Kemanian; 3: Cesar Izaurralde (origional DW); 4: Cesar Izaurralde (new DW)
                 NP_uptakeCurve,   & ! 0: Smith curve; >0: S curve
                 O2_Method,        & ! 0: original EPIC oxygen/depth function; 1: Armen Kemanian carbon/clay function; 2: New O2 ratio method
                 dirChoice,        & ! 0: reading data from working dir; 1: from weatdata dir; 
                                     ! 2: form working  dir plus 3 other 
                 soilSaturate_Method,& ! 0: reading saturated conductivity in soil file; >0: computing saturated conductivity with RAWLS method
                 LatChoice,        & ! 0: use input latitudes; >0: computing equivalent latitude based on azimuth orientation of land slope
                 P_apply_mode,     & ! 0: not auto p apply; >0: auto P apply
                 PAR_Method,       & ! 0: PAR driven by crop LAI development; >0: PAR driven by EVI from remote sensing
                 percolate_Method, & ! 0: for HPERC; 1: for HPERC1 (4mm slug flow) 
                 NC_Miner_Method         ! 0: PHOENIX; >0: CENTURY 
                 
      
      integer:: weatherVarID0, correctionForRuns0
      ! ------------------- add by TXH ----------------------------------
      real:: rainNcon0,            & ! Average concentration of N in rainfall (ppm)
             CO20,                 & ! CO2 concentration in the atmosphere (ppm)
             CNO30,                & ! Concentration of NO3 in irrigation water (ppm)
             saltConIrr,           & ! Concentration of Salt in irrigation water (ppm)
             pestDamgfactor,       & ! pest damage scaling factor(0- 10): --0. SHUTS OFF PEST                  
                                     ! DAMAGE FUNCTION. PEST DAMAGE FUNCTION CAN BE REGULATED FROM             
                                     ! VERY MILD(0.05-0.1) TO VERY SEVERE(1.-10.) 
             periodHalfHrRain,     & ! NO Y RECORD MAX .5H RAIN (BLANK IF WI IS NOT INPUT--LINE 19)
             wetDayCoeff,          & ! COEF (0-1)GOVERNING WET-DRY PROBABILITIES GIVEN DAYS OF RAIN (BIU OR IF W|D PROBS ARE                           
                                     ! INPUT--LINES 16 & 17)  
             rainPara,             & ! PARAMETER USED TO MODIFY EXPONENTIAL RAINFALL AMOUNT                    
                                     !  DISTRIBUTION (BIU OR IF ST DEV & SK CF ARE INPUT--LINES 14 & 15)  
             fieldLength,fieldWidth, & ! (km, BIU)   
             fieldAngle0,          & ! clockwise angle of field length from north (deg) (blank if unknown)      
             deadCropRes0,         & ! standing dead crop residue (t/ha) (blank if unknown)       
             windPowerPara,        & ! power parameter of modified expontential distribution of wind speed (bland if unknown)
             soilDiameter,         & ! soil particle diameter (um, blank if unknown )
             windErosionFactor,    & ! wind erosion control factor 0.0: no wind erosion; 1.0: normal simulation; >1.0: acceleates wind erosion        
             irrTrigger,           & ! 1. plant water stress factor(0-1); 2. soil water tension in top 200mm (>1 kpa); 3: plant available water deficit in root zone (-mm)     
             irrRunoffRatio,       & ! run off volume/ volume of IRR water applied (bland if IRR =0)
             irrAnnualMaxVol,      & ! max annual irrigation volume allowed(mm)        
             irrSingleMinVol,      & ! min single application volume allowed (mm) 
             irrSingleMaxVol,      & ! max single application volume allowed (mm) 
             fertTrigger0,         & ! auto fertilizer trigger 1. plant N stress factor (0-1); 
                                     ! 2. soil N concentration in root zone (g.t)
             fertRate,             & ! 1. application rate auto/fixed (kg/ha); 2. manure input to lagoon (kg/cow/d) IRR =4
             maxAnnualN,           & ! min annual N fertilizer appliaction (kg/ha) 
             dayReduceStress,      & ! time required for drainage sys to reduce plant stress (d) blank if drainage not used
             furrowDikeSafeFactor0,& ! furrow dike sagety factor(0-1)   
             conservFactor0,       & ! conservation practive factor (0.0 eliminates water erosion) 
             lagoonVol0,           & ! lagoon volume ratio (normal/maximum)
             lagoonInput,          & ! lagoon inputs from wash water (m3/cow/d)
             dayReduceLagoon,      & ! time to reduce lagoon storage from max to normal (d)
             liqManureRatio,       & ! ratio liqud/total manure applied  
             grazeLimit,           & ! above groud plant material grazing limit (t/ha)
             herdFrac,             & ! fraction of time herd is in feeding area 
             layerThickness,       & ! layer thickness for differential eqs soln to gas diffuse eqs (m)
             waterEroEq,           & ! water erosion driving eqs (0: MUST; 2: USLE; 3: MUSL; 6: RUSLE; 7: RUSL2)
             baseStockRate,        & ! base stock rate (ha/hd)
             storLeachFracNO3       ! fraction of storage interacting with NO3 leaching
      INTEGER, DIMENSION(4)::IDIR   ! data directory       
      ! ******************* 2. PARM1102.DAT: S-curve and other equation parameters  ***********************************************
      real:: S_Curve(30,2), modelPara(96)
      real:: irrCost,limeCost,fuelCost,laborCost  ! not declared in origional model
 
      ! --------------------- add by TXH --------------------------------
      real:: XKN50,XKN30,XKN10,CBVT0
      ! -----------------------------------------------------------------
      
      ! ******************* 3. TR55COM. DAT: TR55 Runoff estimation parameters ***************************************************
      real::coeffTR55(8,17,4) ! coefficients of 7th degree polynomial in TR55 peak runoff rate estimate

               
      ! ******************** 4. PRINT1102. DAT: print control *********************************************************************
      integer:: varID(800),       & ! an array to store all index of variables  
                stateVarID(100),  & ! an array to store index of variables print actually
                conVarID(4),           & ! an array to store index of variables related to concentration 
                monVarID(40),          & ! index of monthly variables
                annualVarID(40),         & ! index of annualy variables
                econVarID(40)            ! index of economic analysis variables: 0 and 1 
      integer:: numPrintVar, numConVar, numMonVar, numDayVar, numAnnuVar,numEconVar ! the numbers of variables
 
      ! ********************* 5. EPICRUN.DAT: Run settings **********************************************************************
      integer:: WP5_ID, dailyWeatherID ! daily weather station number
      real:: XKN1,XKN3,XKN5,CBVT
 
      
      ! ********************* 6. Site file: parameters of a particular site (e.g., umstead.sit) ********************************* 
      real:: YLAT, XLOG, ELEV, & ! site location (lat, lon) and elevation
             areaWshed,       & ! watershed area    
             uplandSlopLen,       & ! Upland slope length (m) 
             uplandSteep,     & ! Upland slope steepness (m/m)
             conservFactor,    & ! Conservation practice factor (0.0: elimates water erosion)
             timeIntervalGasDiff,  & ! Time interval for gas diffusion Eqs (h)
             CO2, CNO3I, rainNconX, & ! actual variables used for simulation, updated by readings
             ! Below are for day lengh cal according to latitude   
             SIN1, COS1, YLTS, YLTC, YTN,YTN1, HR0, Delta_Daylength, DayLength
 
      ! ------------------- add by TXH -----------------------------
      integer:: Num_GasDiff   ! related to timeIntervalGasDiff
      ! -------------------- OPS info from site file -------------------------------------
      integer:: irrCode, & ! = N0 FOR DRYLAND AREAS          | N = 0 APPLIES MINIMUM OF                
                            ! = N1 FROM SPRINKLER IRRIGATION  | VOLUME INPUT, ARMX, & fieldCapacity-totSoilWater             
                            ! = N2 FOR FURROW IRRIGATION      | N = 1 APPLIES INPUT VOLUME              
                            ! = N3 FOR IRR WITH FERT ADDED    | OR ARMX                                 
                            ! = N4 FOR IRR FROM LAGOON        |                                         
                            ! = N5 FOR DRIP IRR                   
                autoIrrInterval, & ! N DAY APPLICATION INTERVAL FOR AUTOMATIC IRRIGATION
                fertApplyMinInterval,    & ! min fertilizer application interval (blank if unknown)
                furrowFlag,    & ! 0: without furrow dikes; 1: with furrow dikes                
                fertTimes,     & ! Fertilizer # for auto N fertilizer and fertigation -- blank defaults to element N
                autoManureTrigger,     & ! >0: auto dry manure application without trigger
                mowMinInterval,   & ! min interval between auto mow
                numAutoP        ! fertilizer # for auto P fertilizer -- blank defaults t oelemental P
      ! ********************** 7. NCCLAYTO. WP1: read Weather file and weather variables*******************************************
      integer:: rainDistriFlag   ! 0: 
      
      real:: monAveTmax(6,12), monAveTmin(6,12), monTmaxStd(6,12), monTminStd(6,12), &
             halfHrRainInfo(6,12), monAveSoilRad(6,12), RH(6,12),PCF(6,12)
      
      real:: dailyRainInfo(3,6,12),probWetDay(2,6,12)
      
      real:: monAveT(12),                       & ! monthly average temperature
             UAVM(12),                      & ! rainy days in a month      |   monthly wind speed
             monMaxRad(12),                      & ! monthly radiation
             monDayLen(12),               & ! day length on the half of each month
             amplitudeT, & ! upper bound, lower bound, amplitude of average T during a year   
             rainfallVolum, rainfallStd, rainfallSkew,              & ! rainfall volume, std, and skew coefficient of distribution
             radCorrection,                & ! radiation correction 
             TX,                            & ! average temperature 
             V3,                            & ! random number
             BIG                              ! tempoary variable 
      ! ********************** 8. Soil Model and soil file *********************************************************************************
      real:: soilAlbedo,          & ! Soil Albedo
             waterTableMinDep,    & ! MIN DEPTH TO WATER TABLE(m)(BIU)  
             waterTableMaxDep,    & ! MAX DEPTH TO WATER TABLE(m)(BIU)
             waterTableHigt,    & ! INITIAL WATER TABLE HEIGHT(m) (BIU)
             groundWaterStor,     & ! GROUNDWATER STORAGE (mm)  
             groundWaterMaxStor,  & ! MAXIMUM GROUNDWATER STORAGE (mm)  
             returnFlowFrac,  & ! RETURN FLOW / (RETURN FLOW + DEEP PERCOLATION)   
             groundWaterResidenceTime,  & ! GROUNDWATER RESIDENCE TIME(d)(BIU) 
             totCropResi,totLabileP,totActiveP,&
             initialTotSOP, TOP,plawDepLabileP, APBC, APMU,totSalt,totNO3_N,soilWaterRZ,potentialAvailWater, &
             totSolubleK, totExchangeK,totFixK,totCPassiveHumus,totNPassiveHumus,totCSlowHumus,totNSlowHumus,totMetabLitt,&
             totCMetabLitt, totNMetabLitt,totStructLitt,totCStructLitt,totLgStructLitt,totCLgStructLitt,totNLgStructLitt,&
             totNStructLitt, totCBiomass,totNBiomass,totSON, TWN0, ZCS(3),plawDepSOC, plawDepSoilWater, &
             plawDepFieldCapa, ZMIX, aveBulkDensity, totSOC, DZ10,satuateCondv0, &
             BTN,BTK,BTP,AQV,Albedo,RCF,yearAveTem,DD, Tot_FreshOrgP,ZN2O,ZN2OG,ZN2OL,ZNO2, ZNH3,SSW, THK 
             ! BTCX,BTNX, !not used
      REAL::minThickMaxLayer,  & ! MINIMUM THICKNESS OF MAXIMUM LAYER (m)(SPLITTING STOPS WHEN ZQT IS REACHED) 
            minThickProfile      ! MINIMUM PROFILE THICKNESS(m)--STOPS SIMULATION.  
      integer:: IDSP,Actual_SoilLayers,Layer_RD_Reach, LD1, layersEqualThick, IUN
      ! CAP: total field capacity, but not used  -- TXH 
      ! *********************************  N P K in soil (by layers)******************************************
      real, dimension(:), allocatable:: soilElement, wiltPoint, satuateCondv, fieldCapacity, VFC, soilWater, &
                                        Z, bulkDensity, dryBulkDensity, postTillBulkDensity,PH, WT, &
                                        CEC, CEM, sandFrac, siltFrac, clayFrac, rockFrac, totBase, soilTem, mineralBulkDensity,& 
                                        SOC, CaCO3, NH3_Weight,SON, initialNO3,fracNO3Leach, NO3_N_Soil,&                                         
                                        initialLabileP, SOP,  sorbRatioP, flowCoeffP, stableMineralP, labileP,&
                                        exchangeK, solubleK,fixK, EQKE,EQKS, SoilWater_Fact, elecCondv, lateralFlowCondv, &
                                        structLitt, metabLitt, lgStructLitt, CBiomass, CMetabLitt,CStructLitt, &
                                        CLgStructLitt, CPassiveHumus, CSlowHumus, NBiomass, NMetabLitt,NStructLitt,&
                                        NLgStructLitt, NPassiveHumus,NSlowHumus, ironCon, &
                                        gasCO2Con,gasN2OCon,gasO2Con,FreshOrgP_Residu,activeMineralP,Porosity, &
                                        weightDenitriN, Wilt_Point, VGA,VGN ,ALS ,NO2Weight,WN2O,saltWeight,ZC, &
                                        XNetNMineralize,AWC, &
                                        UNA, &  ! N stress (use N stress function to set N application)
                                        ULYN    ! average annual N in crop yield
                                        
                  
      real,dimension(:,:),allocatable:: NPT, SOIL 
      integer, dimension(:), allocatable:: LORG, Layer_ID,numOfPestArr,JP
      
      ! ********************************** N P K in sediments ******************************************************
      real:: QNO3, RNO3, Sediment(8), totSoilErosion, rainEnergyFactor, mineralizedN, &
             sedimentLossN,sedimentLossP, runoffLossLabileP,SCN,dikeHeight, CYAV,CYMX,CYSD 
            
          
      real, dimension(:), allocatable:: VerticalFlow_CO2,verticalFlowN2O,VerticalFlow_NO2,verticalFlowNO3 ,VerticalFlow_O2, &
                                        SubSurfFlow_CO2L,SubSurfFlow_N2O,SubSurfFlow_O2L,SLTP,Tem_Fact, &
                                        HumusMine_RateCons,H2OF,WNO3F,CBNF 
      
      integer:: ISL
      ! ******************************** Nutrient (N P K cycling) **************************************************
      real:: totDenitriN,SGMN,SN2,totN2O,SMP,SIP, QCO2,QN2O,QNO2,QO2,SubSurfFlow_SolubleK,SubSurfFlow_Salt,SubSurfFlow_LaibleP,&
             TSFNO2,TSFNO3, TSFN2O,TSFK,TSFS,Vol_NH3,Nitri_NH3, Tot_SoilRes, SubSurfFlow_NO2,SubSurfFlow_NO3,VerticleFlow_SolubleK,&
             verticleFlowSalt, VerticleFlow_LaibleP,enrichRatio, LaibleP_Layer, SMNIM, denitriedN2O,layerThick,DN2, &
             FixedN_Final, PlantK_Demand, P_StressFactor,N_StressFactor, &
             Tot_Salt_RootZone, Tot_Water_RootZone,SSLT, SK, dailyFixedN
      ! ******************************** Soil Water Content **********************************************
      real:: SW15,SW30, SNA15,SNA30,SNN15,SNN30  
    
      ! ***************************** 9. OPS file: operation data ***********************************************************
      integer::landUseCode, autoIrrCode, hydroGroup,  runoff_Method, JX(7), &
               numOfIrr,numOfPest, numOfFert, operationCode(27)         
      ! 1 = KILL CROP                ! 2 = HARVEST WITHOUT KILL
      ! 3 = HARVEST ONCE DURING SIMULATION WITHOUT KILL     ! 4 =
      ! 5 = PLANT IN ROWS            ! 6 = PLANT WITH DRILL
      ! 7 = APPLY PESTICIDE          ! 8 = IRRIGATE
      ! 9 = FERTILIZE                ! 10 = BAGGING & TIES (COTTON)
      ! 11 = GINNING                 ! 12 = HAULING
      ! 13 = DRYING                  ! 14 = BURN
      ! 15 = PUDDLE                  ! 16 = DESTROY PUDDLE
      ! 17 = BUILDS FURROW DIKES     ! 18 = DESTROYS FURROW DIKES
      ! 19 = START GRAZING           ! 20 = STOP GRAZING
      ! 22 = AUTO MOW                ! 23 = PLASTIC COVER
      ! 24 = REMOVE PLASTIC COVER    ! 25 = STOP DRAINAGE SYSTEM FLOW
      ! 26 = RESUME DRAINAGE FLOW
      INTEGER, DIMENSION(:,:), ALLOCATABLE:: LPC,LFT         
      real:: CN2,CN0,totSoilWater, OPV(9), FCV
      real,dimension(:,:), allocatable:: pestKillEffi, pestRate, irrVolumeArr, WFA, TLMA,&
                                         CND, & ! CN2 
                                         QIR, & ! Runoff volume
                                         irrTriggerArr, & ! irrigation trigger
                                         baseStockRateArr, & ! Base Stock rate
                                         fertCFactirArr, & ! Mineral C Erosion Equation 
                                         HWC 
      ! --------------- irrigation -------------------------------------------
      real::  DALG, lagoonMaxVol, lagoonVol, fertCon,manureWeight,irrFromLagoon, &
              manureInputToLagoon, AFLG, aveWaterFromLagoon, furrowDikeSafeFactor, lateralCondvHyN, lateralCondvHyD,VIRT, &
              waterStress,WTN, NO3_IrrWater,AGP
      integer:: IDR,irrChoice,NII
      integer,dimension(:),allocatable:: KG, JPL
      real,dimension(:,:),allocatable:: VIL,VIR
      ! --------------  fertilization -------------------------------------------
      REAL::  TNOR,COST, CSFX, NStressTrigger
        ! auto fertilizer trigger (2 options); =1. plant N stress factor (0-1); 
        !                                       2. Soil N concentration in root zone g/T
      integer:: NDF, KDF1(8000),IDFT(2), NFA

      integer, dimension(:), allocatable::KDF,IAP,NCP,cropID
      real, dimension(:), allocatable:: FK,FN,FP, ANA, HUI, &
                                        FNH3,&         !AMMONIA N FRACTION(FNH3/FN)
                                        FNO,FOC,FPO,&  ! Organic N, P, C fraction
                                        FSLT, &        ! salt fraction 
                                        FCEM, &        ! C EMITTED/UNIT OF FERTILIZER(KG/KG)
                                        FCST           ! COST OF FERTILIZER($/KG) 
                                                     
      real,dimension(:,:),allocatable:: TotK_Frac,TotN_Frac,TotP_Frac
      ! ------------- pesticide table -------------------------------------------
      integer:: NDP, KDP1(8000), GrowingSeason_Days
      integer, dimension(:), allocatable:: KDP
      integer, dimension(:,:), allocatable:: NPSF 
      real:: pestGrowIndx, VQ(90), VY(90)     
      real, dimension(:), allocatable:: pestHfLifePlant, & !PESTICIDE HALF LIFE ON FOLIAGE (d)
                                        pestHfLifeSoil, & !PESTICIDE HALF LIFE IN SOIL (d)
                                        Coeff_OrgCarbon, & !PESTICIDE ORGANIC C ADSORPTION COEF
                                        pestSolubity, & !PESTICIDE SOLUBILITY (ppm)
                                        WashOff_Frac, & !PESTICIDE WASH OFF FRACTION
                                        Pest_C_emission, & !C EMMISSION/UNIT PESTICIDE(G/G)
                                        Pest_Cost, & !PESTICIDE COST ($/KG)
                                        pestPlantIntercep, PSTF, Pest_Leach,Pest_Drainged    
      
      real, dimension(:,:), allocatable::  totPestControl, pestInSoil, TPSF, SPQ, Pest_Sediment, PVQ, PVY, SPQC 
      real, dimension(:,:,:), allocatable::APQ,AQB,APY,AYB
      
      ! --------------------- add by TXH --------------------------------------
      integer:: IBGN, &
                JJK0   ! Crop ID
      ! ***************************** 10. Tillage file  ***********************************************************************
      integer:: NDT, MXT,cropCode(10),KHV, Mow_DayLap
      integer,dimension(:),allocatable:: NBE,NBT,currentOps, ICUS,IYH
      
      real:: RidgeHeight2, Ridge_IntervalX, PALB, limeRate, CN_Ratio, WBMX,standLiveDeadPlant, &
             GS_Rad ! Growing Season Solar Radiation (MJ/M^2)  
      real,dimension(:),allocatable::COOP,COTL, mixEfficiency, tillDepth,ridgeHeight,ridgeInterval, &
            furrowHeight,furrowInterval, harvestEffi, overRideHI,soilCompactFrac, PlantRedu_frac, &
            carbonEmiss,Fuel_Use, EFM, &        ! machine harvestEffi, grazing rate (Kg/HD/d)                                    
            tillRoughness, & ! random surface roughtness created by tillage operation (mm)
            HMO, abvGroundBiom, actualCropN,actualCropP,actualCropK, UNMX, UPMX, &
            soilSupplyRateN, soilSupplyRateP, soilSupplyRateK,TRA, WLV 
      real,dimension(:,:),allocatable::DMF, HIF,CSOF,CSTF,ET_GrowingSeason
      
      ! ****************************** 11. Crop Parameter *************************************
      ! -------------- read from crop database ---------------------
      integer:: LC
  
      real,dimension(:),allocatable::XMTU, RUE, HI, optimalTem, plantMinT, maxLAI, fracGrowSeasonLAIDecline, factorLAIDecline,&
            factorBioDecline, maxStomaCond, cropTolerantAl, areationThreshold, seedRate, maxCropHeight, &
            maxRootDep,fracYieldK,fracYieldN,fracYieldP,seedCost,pestDmagFactor, seedYieldPrice,forageYieldPrice,&
            waterFracYield, VPDPara, VPDThreshold, VPDPara2, germinateHeatUnit,WUE, &
            lintFracCottn,turnoutFracCottn,carbonEmissionSeedWeight, leafWeightFrac,FNMX,AboveGround_Biomass0
      
      
      real,dimension(:,:),allocatable::pointsLAIDevp, CO2EffOnBio, uptakeParaK,uptakeParaN,uptakeParaP,windEroCrop, &
            pointsFrostDamg, rootPartition,pointsPlantPopulation, salinityThreshold ,ligninFrac,coeffPopuCurve,&
            PPLA, POP
      
      integer,dimension(:,:),allocatable:: NBCX
      ! ------------------ derived crop-realted variables ------------  
      ! -------------------- crop residule -------------------
      real:: standDeadResOrganism,standDeadResiP, StandDeadRes_OrgP, standDeadResiOrgN, standDeadLignin,&
             standDeadResiK, StandDeadRes_OrgK, groundCover,coverFactor,coverLagingFactor
      real, dimension(:), allocatable:: standDeadResiN, standCropResi, cropResidu 
      integer:: plantCategoryCode(10)
      ! ------------------- crop -----------------------------
      real,dimension(:,:),allocatable::YLD1, YLD2,YLCF, YLKF,YLNF,YLPF,SFMO,growStressFactor,RWTX  
      real, dimension(:), allocatable:: sumActuralEP,sumPotentialEP,ACET,XDLA0,XDLAI, totRootWeight,CPHT, &
                                        XLAI, DMLX, LfWidth, PPL0, YLD, AJHI, lowerLimitHI, AJWA, Min_Stress_Factor, WCHT,&
                                        BLYN,potentialBioIncrese,UK2,Optimal_N,Optimal_P, RWX 
      
      real:: GSEP, GSVP, GS_PlantTranspiration, &  ! Growing Season Plant Transpiration (mm)  
             GS_VPD, SRA, WaterFrac_Yield, ROSP, SM99, Yield_N, YLP, YLC, YLK, CLG, Var_GS(8), TYC, TYK, &
             TYN, TYP, SHRL, SAT, RootZone_Depth, plantEP, Poten_WaterUse_Rate, potentialDemandN, &
             potentialDemandP, SUN, SUK, SUP, HIX
          
      real,dimension(:,:,:),allocatable::Growing_Stress 
                                    !     Growing_Stress(1,     = WATER STRESS (d)                                                   
                                    !     Growing_Stress(2,     = N STRESS (d)                                                       
                                    !     Growing_Stress(3,     = P STRESS (d)                                                       
                                    !     Growing_Stress(4,     = K STRESS (d)                                                       
                                    !     Growing_Stress(5,     = TEMPERATURE STRESS (d)                                             
                                    !     Growing_Stress(6,     = AERATION STRESS (d)                                                
                                    !     Growing_Stress(7,     = SALT STRESS (d)
      
      ! ********************************** Root development **************************************
      real::RZ,& ! root zone depth/ minimum of soil profile depth/soil profile depth?
            UB1,UOB
      real, dimension(:), allocatable:: RD, RDF           
      
      ! *************************************** Erosion Model ************************************  
      real:: XAV(3), XDV(3), XRG(3), BRNG, EXNN, SX, WSX, plawDepth, minCFactor, slowHumusTranRate, &
             passivHumusTranRate, dampDepAdj, runoffAdj, USL
      REAL, ALLOCATABLE:: XSP(:,:)  
      integer:: NSX, waterEroModel,NDVSS, NSNN
      ! ******************************************************************************************* 
      ! *************************************** Water erosion *************************************
      real:: RLF, RSF, SL, REK, EK, RSK, RSLK, abvGroundResidue, GroundCover_Frac, &
             GroundCover_StandLiveBio, coverFactorPlantPopu, CropManage_Factor, BETA, slopeLength, &
             RUSM(3), CX(12), totHfHourRain(12), maxRH, DR
      integer:: IRO
      
      real,dimension(:,:),allocatable:: RWT, RWF
 
      ! *************************************** Wind erosion ***************************************
      real:: RNCF(12), U10MX(12), TMNF(12), windDirection(12,16), &
             soilRad, surfTem,Wind_Erod, windErodSoil, fieldAngle, RandRough, Roughtness_Factor, &
             vegCoverFactor, TLMF, enrichRatioFinal  ! add by TXH
      integer:: NWDA ! why could not declare enrichRatio ?   
 
      ! ***************************************** EP model ****************************************** 
      real:: barometricPressure, psychrometerConst, EO, HGX, soilEvp, VPD
      integer:: NEV
      real, dimension(:), allocatable:: Current_LAI, EP
 
      ! *************************************** Heat accumulation ***********************************
      real:: yearTotHeat ,HSM, XHSM
      real,dimension(:),allocatable:: HU,  FLF0, totCropBio,totCropBioX
      real,dimension(:,:),allocatable:: potentialHeatUnit, fracHeatUnit
      
      integer:: Crop_Num,IGO, JDHU,IPL,IHV,KC,KT  
      integer,dimension(:),allocatable:: numOfTill, NHU,NBC,KOMP,IHU, JE
      integer,dimension(:,:),allocatable::LY,tillDOY,IPLD, IHVD, LYR, LT, JH
 
      ! ********************************* Gas diffusion model ************************************
      real, dimension(:), allocatable:: soilResp, CLO2, AO2C, AFP, HKPO, VWC, VWP, TPOR, CLCO2,&
                                        ACO2C, AN2OC, HKPC, HKPN, CLN2O, DRWX, Weight_CO2L, WCO2G, &
                                        Weight_N2OL, WN2OG, WO2G, Weight_O2L, SOT, &
                                        DPRC,DPRN,DPRO,DCO2GEN,DN2G,DN2OG,DO2CONS,SMEA,SMES, RWTZ,EAR
          
      real:: DFCO2T, DFO2T, DFN2OT, DFCO2S, DFO2S, DFN2OS 
      
      ! ******************************** Daily Weather ********************************
      real:: TMXM,TMNM,TXSD,TNSD,SRAM,TMX,TMN,Rainfall,RHD, EVI, TXXM, &
             WX(3), XIM(3)
      integer::KGN(5) = 0,LW
      ! --------------------------- Solar Radiation ------------------------
      real:: maxRad,RM,SRAD
      ! --------------------------- Wind Speed (10m height)-----------------
      real:: U10
      ! --------------------------- RH -------------------------------------
      real:: RHM
      ! --------------------------- add by TXH ----------------------------- 
      integer:: ICV
      
      ! ********************************** Hydrology (Rain, Snow, and Runoff) **************************
      real::SRD(12),Inflow(12), ARF, snowPackAge,totRunoff, snowMelting, snowWater, RainfallIntesity_5hr, &
            Rainfall_Duration,DKHL,dikeInterval, Dike_Volum, IDRL, Flow_Rate, SCI, SMX,&
            accumRainfall30, PMORF,Percolation, lateralFlow,peakRunoffRate, Runoff, alpha05L, &
            crackFlow  
      
      real, dimension(:), allocatable:: GWPS, RSPS, subsurfLaterFlow, &
                                        percolateFlow ! Added by TXH, percolateFlow is separated from 
                                                         ! initialLabileP
                                
      real, dimension(:,:), allocatable:: availWaterCrop,CQV,CRF
         
      !------------------------- runoff ----------------------------------
      real:: Tot_Rainfall, PRAV, PRSD, PRB,peakRainRate, QPQB, QPS, TCAV,TC, TCS, TCC, TCMX, TCMN
      integer:: NQP
 
      ! ************************** random numbers ********************************
      INTEGER, DIMENSION(21):: IX, IX0
      integer:: randID(21), randID0(21)  ! V1, V3 ,V4 are also random numbers, but not declaired

      ! ************************** Variables *************************************
      real, dimension(:), allocatable:: VAR, &    !  daily weather variables
                                        SMY, &    !  store monthly variables temporarily
                                        SM        !  yearly output
      
      ! SMY(3) = SOLAR RADIATION (MJ/M^2)           SMY(4) = PRECIPITATION (mm)  
      ! SMY(5) = SNOWFALL (WATER EQUIVALENT mm)     SMY(6) = SNOWMELT (WATER EQUIVALENT mm)  
      ! SMY(10)   = POTENTIAL ET (mm)  
      ! SMY(11)= ACTUAL ET (mm)                     SMY(12)= POTENTIAL PLANT TRANSPIRATION (mm) 
      ! SMY(13)= ACTUAL PLANT TRANSPIRATION (mm)    SMY(14)= RUNOFF (mm) 
      ! SMY(15)= SCS CURVE NUMBER                   SMY(16)= LATERAL SUBSURFACE FLOW (mm)  
      ! SMY(17)= PERCOLATION (mm)                   SMY(18)= DRAINAGE SYSTEM FLOW (mm)   
      ! SMY(19)= IRRIGATION VOLUME (mm)             SMY(20)= EXTERNAL INFLOW TO MAINTAIN WATER TABLE (mm) 
      ! SMY(30,31,33,34,35,36)= results for different water erosion models (t/ha)                                                                                               
      ! SMY(42)= WIND EROSION (t/ha)                SMY(43)= ORGANIC N LOSS WITH SEDIMENT (kg/ha)                                             
      ! SMY(46)= N LEACHED (kg/ha)                  SMY(47)= NET MINERALIZATION (kg/ha)                               
      ! SMY(49)= DENITRIFICATION (kg/ha)            SMY(50)= N FIXATION (kg/ha)                                                                                                                                       
      ! SMY(51)= NITRIFICATION (kg/ha)              SMY(52)= VOLATILIZATION (kg/ha)                                             
      ! SMY(53)= N LOSS IN DRAINAGE SYSTEM (kg/ha)  SMY(54)= P LOSS WITH SEDIMENT (kg/ha)                                 
      ! SMY(56)= P MINERALIZATION (kg/ha)           SMY(57)= P LEACHED (kg/ha)                                           
      ! SMY(58)= ENRICHMENT RATIO                   SMY(59)= N FERT ORGANIC FORM (kg/ha)                                         
      ! SMY(60)= N FERT NO3 FORM (kg/ha)            SMY(61)= N FERT NH3 FORM(kg/ha) 
      ! SMY(62)= P FERT ORGANIC FORM (kg/ha)        SMY(63)= P FERT MINERAL FORM (kg/ha) 
      ! SMY(68)= SOIL WATER IN TOP 10mm OF SOIL(mm) SMY(77)= ORGANIC C LOSS WITH SEDIMENT(kg/ha) 
      ! SMY(85)= N MINERALIZATION FROM HUMUS (kg/ha)  
 
      real, dimension(:,:), allocatable:: SMM, &      ! Monthly output variables
                                          SMS, &      ! Soil factors for microbiology process for all layers
                                          cropVar, &  ! Crop related variables
                                          pestVar, &  ! Pesticide related variables
                                          SMAP,SMYP,SOL
      
      real, dimension(:,:,:), allocatable:: SMMC, & 
                                            SMMP      ! Pestiside related daily variables 
      
      INTEGER::IPY,IRL,ISTP,ISX,JCN0,JCN1            
      
      INTEGER::K21,M21,MO1,MYR,NQP0,NQP1,NWD0              
      
      INTEGER::KDT2(100),dailyVarID(40),IHRL(12),IYS(8)                            
      
      INTEGER,DIMENSION(:),ALLOCATABLE::NCR,NHV,NYLN,IHT, NX                                                                   
      
      INTEGER,DIMENSION(:,:),ALLOCATABLE::IGMD
      
      REAL::BARF,CN,CN1,CN3,CSO1,CST1,DARF,solubleNConc 
      
      REAL::PBWM,RGSM,SARF,SNOF,S3                                                                
       
      REAL::VALF1, VALF2, V1, winterDormancy, YLDX                      
                                               
      REAL::VARS(30),ASW(12),RCM(12),SET(12),RSY(12),TAMX(12),TEI(12),TET(12),TMXF(12),TQ(12),TR(12),&             
            TRHT(12),TSN(12),TSR(12),TSTL(12),TSY(12),TXMN(12),TXMX(12),TYW(12),W(12),WCS(3)                                                  
      REAL::STV(30,12),wetDayFreq(6,12)
      
      REAL,DIMENSION(:), ALLOCATABLE::HUF,PSTM,TFTN,TFTP,TCAW,TCQV,TCRF,TCSO,TCST,TDM                                             
      REAL,DIMENSION(:),ALLOCATABLE::TETG,TFTK,THU,TIRL,TRD,TVAL,TVIR,TYLC,TYLK,TYLN,TYLP,TYL1,TYL2
      REAL,DIMENSION(:,:),ALLOCATABLE::STDA, TSFC                                                
      REAL,DIMENSION(:,:,:),ALLOCATABLE::APQC  
      
END MODULE PARM                                                                           
