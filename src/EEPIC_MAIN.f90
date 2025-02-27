! **************************************************************************************
!  This is a modified version of the EPIC model (v.1102) by TXH (March 2020). 
!  Extended functionbilities including radiative transfer model and photosynthesis process.
!  note: BIU: Blank If Unknown   
! ************************************************************************************** 
    
PROGRAM EEPIC                                                                        
      USE PARM   
      USE ReadVarTable_Module
      USE MisceFun_Module
      USE Initialize
      USE Weather_Generator
      USE Soil_Module
      USE GasDiffuse_Module
      USE Crop_Module
      USE Hydrology
      USE Nutrient_Module
      USE Print_Module
      USE ReadWeather_Module 
      IMPLICIT NONE

      ! /////////////////////Declaration of local variables /////////////////////////////////// 
      
      CHARACTER(LEN = 1), DIMENSION(4)::ASG = ['A','B','C','D']                                                           
      CHARACTER(LEN = 4)::ANM                                                          
      CHARACTER(LEN = 80)::fileName,CM(9),FCMOD,FCROP,FFERT,FMLRN,FOPSC,&
                        FPARM,FPEST,FPRNT,FSITE,FSOIL,FTILL,FTR55,FWIDX,FWIND,FWLST,FWPM1,&
                        FWPM5,soilName,TITW5(2),WINDFILE,WPM1FILE                                                                                                             
      REAL(KIND = 8), DIMENSION(13):: PPX                                                                    
                                                
      INTEGER:: NXX(3,5),NYY(3,5), JZ(4), NY(3), &           
                IIX(7)             ! an array for store site info for multi runs  
      
      REAL, DIMENSION(30,2):: SCRX                      ! copy of s curve coeff
      REAL, DIMENSION(13,16):: XZP
      REAL, DIMENSION(200):: XTP,XYP
      REAL, DIMENSION(16):: YTP
      REAL, DIMENSION(15):: sulfCon, fracVPip, fracHPip
      REAL, DIMENSION(2)::otherCost 

      INTEGER:: I, I1, II, IY1, IRTC, J, J1, JRT, L, L1, maxLayer, N1, NN, N2, &
                KRX, KK, K, K1, K2, XCN, XLC, XYR, XQP, IWRT, LRG, &
                siteID, soilID, operationID, windID, WP1_ID, IDSK

      REAL::X1, X2, X3, X4, X6, XX, FTN, FTC, FTP, FTK, SUM, TOT, maxRainfall, meanDayLength, &
            RXM, SWW, SRF2, SM1, YLAZ, PZW, lagoonVolX, BTC, AD1, AD2, AD3, AD4, ADD 

      REAL::vernalizeDay,         & ! added by TXH    
            soilHydGroup,         & ! HYDROLOGIC SOIL GROUP--1.=A; 2.=B; 3.=C; 4.=D  
            fracFC,               & ! FRACTION OF FIELD CAP FOR INITAL WATER STORAGE (BIU) 
            groundWaterResidenceTime0,  & ! Initial GROUNDWATER RESIDENCE TIME(d)(BIU) 
            splitedSoilLayer,           & ! NUMBER OF SOIL LAYERS AFTER SPLITTING (3-15); 0 NO SPLITTING OCCURS INITIALLY
            soilCondCode,               & ! 0. FOR CALCAREOUS SOILS AND NON CALCAREOUS WITHOUT WEATHERING INFORMATION;
                                          ! 1. FOR NON CACO3 SLIGHTLY WEATHERED; 2. NON CACO3 MODERATELY WEATHERED  
                                          ! 3. NON CACO3 HIGHLY WEATHERED; 4. INPUT PSP OR ACTIVE + STABLE MINERAL P (kg/ha) 
            cultivaYrs,            & ! NUMBER YEARS OF CULTIVATION AT START OF SIMULATION.
            soilGroup,             & ! 1 FOR KAOLINITIC SOIL GROUP; 2 FOR MIXED SOIL GROUP; 3 FOR SMECTITIC SOIL GROUP         
            fracSOC,               & ! FRACTION OF ORG C IN BIOMASS POOL(.03-.05) 
            fracHumus,             & ! FRACTION OF HUMUS IN PASSIVE POOL(.3-.7)
            XCC,                   & ! CODE WRITTEN AUTOMATICALLY BY .SOT (NOT USER INPUT)
            adjustCN2,             & ! slop adjusted CN2 (added by TXH)
            manureWeightX
      !======================= site information ==================     
      REAL::peakRateFactor,        & ! peak rate - EI adjustment factor (blank if unkown)
            snowCont0,             & ! water con of snow on ground at start of simulation (mm) 
            siteAzimuth,           & ! azimuth orientation of land slope (degrees clockwise from north)
            channelLen,            & ! Mainstream channel length (km, blank if unknown) 
            channelDep,            & ! Channel depth (m)    
            channelSlp,            & ! Mainstream channel slope (m/m, blank if unknown) 
            surfNForChannel,       & ! Surface N for channel (blank if unknown)
            nManningCoeff            ! Manning's roughness coefficient (APEX doc Eq. 37 and Eq. 76 )
      INTEGER:: NCC, drainageDepth   ! 0: no drainage  ;  1: depth of drainage system (mm)    
 
            
      REAL, DIMENSION(6)::coeffChannelGeometry  ! coefficients for channel geometry: 
                                                ! X = Channel_geometry(N)*areaWshed** channel geometry(N+1) 
                                   ! X = depth for Channel_geometry(3) and (4); X = length for Channel_geometry(5) and (6) 
      !=================weather data================================
      REAL, DIMENSION(6,12)::monAvePPT 
      REAL::windEroFactor,  &
            yearTotPPT           ! total annual precipitation
      !=================Soil data from .sol file====================                         
      REAL::minThickSplit,     & ! MINIMUM LAYER THICKNESS FOR BEGINNING SIMULATION LAYER SPLITTING--
                                 ! MODEL SPLITS FIRST LAYER WITH THICKNESS GREATER THAN ZTK(m); 
                                 ! IF NONE EXIST THE THICKEST LAYER IS SPLIT. 
            totSalt0
         
      ! ////////////////////////////////////////////////////////////////////////////////////////////////
      
      ! Begin
      CALL Timer(0)       ! record the date and time
      CALL Setting_Var    ! Setting values for some parameters: e.g., control variables                                                                               
      CALL Setting_Random
      
      IRUN=0
      MSO=33              ! output file ID
      
      KW = [(i + 50, i = 1, 200)] ! KW: 51-250, units for writting files
    
      OPEN(KW(MSO+6),FILE='G:\EEPIC_PilotStudy\LM\WORKSPACE.DAT') 
      ! runMode = 0 RUNS FROM EPICRUN.DAT
      !         > 0 RUNS BATCH MODE FROM RUNALL.BAT                                   	  
      READ(KW(MSO+6),'(I4)') runMode
      
      ! open a file to save error info 
      fileName='EPICERR.DAT'            
      CALL OPENV(KW(MSO),fileName,0,KW(MSO))   
      ! if run in batch mode
      IF(runMode>0)THEN
            I=1
            DO
                  READ(KW(MSO),'(A80)',IOSTAT=NFL)TITW5   
                  IF(NFL/=0)EXIT
                  I=I+1
            END DO
            IF(I==1)THEN
                  REWIND KW(MSO)
                  WRITE(KW(MSO),621)IYER,IMON,IDAY,Hour,Mint,Sec                               
                  WRITE(*,621)IYER,IMON,IDAY,Hour,Mint,Sec                               
            END IF
      END IF
       
      ! ================ 1. read input file lists and model settings============ 

      fileName='EPICFILE.DAT'           ! input file lists 
      CALL OPENV(KR(12),fileName,0,KW(MSO))  ! KR = 11...40 

      READ(KR(12),509)FSITE,FWPM1,FWPM5,FWIND,FWIDX,FCROP,FTILL,FPEST,&              
                  FFERT,FSOIL,FOPSC,FTR55,FPARM,FMLRN,FPRNT,FCMOD,FWLST  
      CLOSE(KR(12))

      CALL Model_Setting(FPARM, FPRNT, FTR55, WINDFILE, otherCost, SCRX, coeffChannelGeometry)

      weatherVarID0=weatherVarID
      correctionForRuns0 = correctionForRuns  
      ! If operationCode and plantCategoryCode were not provided values, they will be initialized here  
      ! operationCode has already been initialized with 1..27    
      IF(operationCode(1)==0) operationCode = [(I, I = 1, 27)] 
	  
      ! plantCategoryCode has already been initialized in AHEAD file with 1..10   
      IF(plantCategoryCode(1)==0) plantCategoryCode = [(I, I = 1, 10)]
      
      CALL Initialize_Var(sulfCon)   ! initialzies variables    
      
      randID = [(i, i = 1, 21)]

      CALL split_int(printChoice,printInterval)  ! printInterval : the flag of print out interval: monthly daily annually
      
      NOP=1                                                                          
      IF(printChoice<=5.AND.printInterval>0)NOP=0 ! NOP controls to write IYR, MO, DayOfMon to KW(1)  
      IF(printInterval==0)printInterval=1         ! =1: 1 day interval  
      IF(printChoice<=5)IPYI=printInterval        ! What is IPYI? 
      
      ! These parameters can be find in the user manual
      plawDepth=modelPara(16)            ! plaw layer depth (m) used to track soluble phosphorus concentration or weight  
      minCFactor=modelPara(32)           ! Minimum C factor values in soil erosion equation 0.0001-0.8
      slowHumusTranRate=modelPara(47)    ! century slow humus transformation rate 
      passivHumusTranRate=modelPara(48)  ! century passive humus transfromation rate
      dampDepAdj=modelPara(82)           ! damping depth adjustment for soil temperature
      runoffAdj=modelPara(83)            ! runoff volume adjustment for direct link

      ! ========================= 2. read EPICRUN.DAT ============================== 
      IF(runMode>0)THEN
            CALL Get_CommendLine(CM)
            siteName=CM(2)
      ELSE
            fileName='EPICRUN.DAT'
            CALL OPENV(KR(11),fileName,0,KW(MSO) ) ! Open EPICRUN file                                            
      END IF 
     
      ! ------------------------ Site info from EPICRUN.DAT file ------------------
      ! EPICRUN.DAT include site number, soil number, weather station, wind station...  
      IF(runMode==0)THEN  
            ! read line 1                                                                        
            READ(KR(11),*,IOSTAT=NFL)siteName,siteID,WP1_ID,WP5_ID,windID,soilID,operationID,dailyWeatherID, &   
            XKN5,XKN3,XKN1,CBVT                  ! these are zero values
            IF(NFL/=0)GO TO 219 
            IF(siteID==0)GO TO 219
      ELSE
            WRITE(KW(MSO+6),'(A80,A8,1X,7A8)')CM  ! for batch run mode
            siteName=CM(2)
            REWIND KW(MSO+6)
            READ(KW(MSO+6),'()')
            READ(KW(MSO+6),'(88X,7I8)')IIX
            siteID=IIX(1)
            WP1_ID=IIX(2)
            WP5_ID=IIX(3)
            windID=IIX(4)
            soilID=IIX(5)
            operationID=IIX(6)
            dailyWeatherID=IIX(7)  
      END IF
      
      siteName=TRIM(ADJUSTR(siteName))   ! right aligned
      IF(XKN5<1.E-10)XKN5=XKN50     ! parameters in PARM file
      IF(XKN3<1.E-10)XKN3=XKN30
      IF(XKN1<1.E-10)XKN1=XKN10
      IF(CBVT<1.E-10)CBVT=CBVT0       
      
      WRITE(*,679) siteName,siteID,WP1_ID,soilID,operationID,dailyWeatherID   ! print to the screen   


      CALL Open_OutputFiles         ! open output files 

      ! ==================== 3. write to the summary file (.out) ========================= 
      WRITE(KW(1),621)IYER,IMON,IDAY,Hour,Mint,Sec        ! KW(1) write date and time to the .out file
      WRITE(KW(1),508)'FSITE',FSITE,'FWPM1',FWPM1,'FWPM5',FWPM5,'FWIND',&     
            FWIND,'FWIDX',FWIDX,'FCROP',FCROP,'FTILL',FTILL,'FPEST',FPEST,&    
            'FFERT',FFERT,'FSOIL',FSOIL,'FOPSC',FOPSC,'FTR55',FTR55,'FPARM',&  
            FPARM,'FMLRN',FMLRN,'FPRNT',FPRNT,'FCMOD',FCMOD,'FWLST',FWLST  
      ! print out variable names
      WRITE(KW(1),'(/1X,A)')'-----VARIABLE NAMES' 
      WRITE(KW(1),229)  
      WRITE(KW(1),287)varName
      ! Print out S curve, coeffecients of the S curve, and other parameters     
      WRITE(KW(1),245)                                                                       
      WRITE(KW(1),242) (J,modelPara(J),(SCRX(J,I),I=1,2),(S_Curve(J,I),I=1,2), J = 1, 30) 
      WRITE(KW(1),244) (J,modelPara(J), J = 31, SIZE(modelPara))  

      WPM1FILE=' '  
      WINDFILE=' ' 
      
      leapYr=1       ! default leap marker  
      IWRT=0         ! control to print head info to each output file (only need to write once for each file)    
      IRO0=1         ! years of crop rotation   
      IRUN=IRUN+1 
      ! =========================== 4. Site info from .SIT file =======================
      CALL Site_Info(FSITE, siteID, coeffChannelGeometry, peakRateFactor, snowCont0, siteAzimuth,  &
                        channelLen, channelSlp, channelDep, surfNForChannel, nManningCoeff, &
                        drainageDepth, meanDayLength, RXM)  

      
      CALL Print_Page(0)   
      
      IF(weatherVarForRuns == 0 .OR. IRUN==1) IYR=Year0  ! readed from EPICCONT.DAT 
      dateNum=Day0+100*Month0+10000*IYR  
      IF(WP5_ID>0)THEN
            READ(KR(20),'()')  
            READ(KR(20),300)II  
            II=IYR-II  
            DO I=1,II  
                  READ(KR(20),'()')  
            END DO
      END IF 
 
      !==================== 5. Read and print weather info =======================
      CALL OPENV(KR(27),FWLST,IDIR(1),KW(MSO) ) ! Catalog of weather stations with daily weather data
      IF(weatherVarID>0)THEN                            ! ID numbers of weather variables
            IF(weatherVarForRuns == 0 .OR.IRUN==1)THEN  ! weatherVarForRuns ==0: normal runs with daily weather input
                  CALL SearchWeatherList(IDIR(1))          !IDIR : directory  
                  ! read daily weather or locate the nearest station  
                  IYR=Year0                                                                       
            END IF
      END IF      
      ! leapYrFlag: 0 if leap year is considered, 1 if leap year is ignored
      ! leapYr    : 0 for leap year, 1 for common year
      CALL Leap_Yr_Check(IYR,leapYr,leapYrFlag)
      CALL Print_Page(0)   
      WRITE(KW(1),'(//1X,A,A/)')'____________________WEATHER DATA________________________'

      ! -------------- cal and print CO2 concentration -----------------------
      IYX=Year0-1880  
      IF(CO2_Method==0)THEN  
            WRITE(KW(1),'(T10,A)')'STATIC ATMOSPHERIC CO2'
      ELSE   
            IF(CO2_Method==1)THEN 
                  IF(IYX<25)THEN  
                        CO2=280.33  
                  ELSE    
                        X1=IYX  
                        CO2=280.33-X1*(.1879-X1*.0077)  
                  END IF 
                  WRITE(KW(1),'(T10,A)')'DYNAMIC ATMOSPHERIC CO2' 
            ELSE                                                                           
                  WRITE(KW(1),'(T10,A)')'INPUT ATMOSPHERIC CO2' 
            END IF                                                                              
      END IF
      WRITE(KW(1),'(T10,A,F6.0,A)')'CO2 CONC ATMOSPHERE = ',CO2,' ppm'       
      WRITE(KW(1),'(T10,A,F4.0,A)')'PERIOD OF RECORD P5MX =',periodHalfHrRain,' Y'   
      ! -------------------- print climate variables ----------------------------
      ! climate variables control
      KGN = 0  ! Precip = 1; Temp = 2; SolarRad = 3; WindSpd = 4; RelHum = 5
      IF(weatherVarID>0)THEN                                                                  
            CALL Print_WeatherVar   ! write the names of input weather variables                                                         
      ELSE                                                                           
            WRITE(KW(1),'(/T10,A,A)')'**********RAIN, TEMP, RAD, WIND SPEED,',&    
                  ' & REL HUM ARE GENERATED**********'                                         
      END IF    
     
      CALL MonthlyWeatherInfo(FWPM5, FWPM1, FWIND, WP1_ID, windID, siteAzimuth, peakRateFactor, WPM1FILE, WINDFILE, monAvePPT, &
                              meanDayLength, maxRainfall, vernalizeDay, yearTotPPT, windEroFactor, YLAZ)


      ! ============= 6. Print out details of models/methods used in the simulation ============ 
      CALL Print_Page(0)    
      CALL Print_GeneralInfo(channelLen, channelSlp, nManningCoeff, channelDep, siteAzimuth,&
                              surfNForChannel, RXM, peakRateFactor, snowCont0, otherCost, KRX)     
      
      ! ============================== End of printing general info =====================
      
      ! =============================== Initialize data again ===========================
  531 CALL Initialize_Var(sulfCon)                                           
      ! initialize random number                                                       
      randID = [(i, i = 1, 21)]                                                                     
      CALL Setting_Random     ! initialize IX                                                        
   
      ! ============================== 7. Random numbers================================ 
      CALL Print_Page(0)       
      CALL Print_RandomNums                                                           

      V3=Generate_Random(randID(3))  
      V1=Generate_Random(randID(2)) 
      BIG=MAX(.2,modelPara(24))     ! maximum depth for biological mixing (m) Range: 0.1-0.3 
      IBD= Cal_DOY(Month0, Day0, leapYr)     ! IBD: User decided date when a simulation begin from, by setting day0  
      DayOfYear=IBD              
      MO=1                         
      CALL Cal_Mon(IBD,MO)          ! decide the month a simulation begin             
      MO1=MO  
      ! ============================ 8. read soil data =================================================== 
      CALL Read_Soil(FSOIL, soilID, yearTotPPT, soilName, soilHydGroup, fracFC, groundWaterResidenceTime0, &
                  splitedSoilLayer,soilCondCode, cultivaYrs, soilGroup, minThickSplit, fracSOC, fracHumus,&
                  XCC, NCC, sulfCon, maxLayer, fracVPip, fracHPip, IDSK) 

      ! =========================== 9. Soil process by layers (constants and initialization) ============== 
      PZW=0. 
      RZ=2.         ! RZ: root zone depth/minimum of soil profile depth/soil profile depth?  
      ! -------------------------- soil variables related to soils layer by layer ----------------------- 
      CALL Soil_Process(NCC, IDSK, maxLayer, minThickSplit, fracFC, fracSOC, fracHumus, PZW, XX)

      ! ============================= 10. C, N, P, K  in actual soil layers ===============================     
      ! only C and N concentration for each layer    13 variables                                        
      ! The first 15 colums is for each soil layer, the last is summary                                                
      XZP = 0.0
      ! Store variables of each layer into arraies: XZP and SOL
      CALL Nutrient_Process(NCC, XX, XZP)

      ! ============================= 11. Field operation schedule ==========================
      ! ---------------------  Check file existence ------------------------------ 
      CALL OPENV(KR(3),FTILL,IDIR(1),KW(MSO) )  ! Tillage database
      CALL OPENV(KR(4),FCROP,IDIR(1),KW(MSO) )  ! crop parameters database      
      CALL OPENV(KR(8),FPEST,IDIR(1),KW(MSO) )  ! pesticide database 
      CALL OPENV(KR(9),FFERT,IDIR(1),KW(MSO) )  ! fertilizer database   
      CALL OPENV(KR(15),FOPSC,IDIR(4),KW(MSO) ) ! Catalog of available operation schedules                                                                          
      !  CALL OPENV(KR(20),FWIDX,IDIR(1),KW(MSO) ) ! Southern oscillation coefficients file 
      IY1=1 
      CALL Management_Info(IY1, KRX,  operationID, JRT, soilHydGroup, adjustCN2)
      K2=1
      ! ===================================== 12. Print soil information by layers ====================

      CALL Prep_SoilVar(YTP)   ! PREPARES SOIL DATA TABLE FOR OUTPUT, AND CONVERTS WEIGHTS OF MATERIALS TO CONCENTRATION.  
      CALL Print_SoilVar1(soilName, ASG, adjustCN2, splitedSoilLayer, minThickSplit, &
                              fracSOC, fracHumus, cultivaYrs, soilCondCode, soilGroup)
      CALL Print_SoilVar2(YTP,1)            ! print out the soil table  YTP: the total of parameters for all layer    
                                                                 
      YTP(2)=YTP(2)*XX*1000.               ! XX is the depth  
      totSoilWater=YTP(2)                  ! totSoilWater: Soil Water  
      SWW=YTP(2)+snowWater                 ! SWW: soil water + snow water  
      totSalt0=totSalt                     ! SLT0: Initial Soil Salt  
      SCI=MAX(3.,SMX*(1.-soilWaterRZ/potentialAvailWater))                                                       
      
      CALL Print_Page(1)
      CALL Print_ManagementInfo(maxRainfall, yearTotPPT, drainageDepth, lagoonVolx, manureWeightX)

      XNetNMineralize=0.       ! Initialize
      !=================================== 13. Write head info to each output file ======================================== 
      
      IF(IWRT==0)THEN          ! ??? What is IWRT? 
            CALL Write_HeadInfo1(soilName, WPM1FILE, WINDFILE, XCC, soilHydGroup, fracFC, groundWaterResidenceTime0, &
            splitedSoilLayer, soilCondCode, soilGroup, minThickSplit, fracSOC, fracHumus, YLAZ) 
            IWRT=1   
      END IF  

      ! What does this block mean?
      DO I=1,cropRotationYrs                        ! cropRotationYrs is the number of years in .OPC file                                             
            IF(NBC(I)==0)NBC(I)=1   
            IF(LY(I,1)>0)CYCLE   
            I1=I-1      
            IF(I1==0)I1=cropRotationYrs   
            LY(I,1)=LY(I1,NBC(I1))  
      END DO     
      ! ======================================15 Crop growth =========================================== 
      IF(IGO>0)THEN
            NBCX(1,LY(cropRotationYrs,NBC(cropRotationYrs)))=NBCX(1,LY(cropRotationYrs,NBC(cropRotationYrs)))+1                                 
            NN=NBC(1)                                                                      
            DO J=1,MNC                                                                 
                  IF(NHU(J)==0)CYCLE
                  DO I=1,NN                                                                      
                        IF(LY(1,I)==J)EXIT
                  END DO
                  IF(I<=NN)CYCLE                                                                         
                  NBC(1)=NBC(1)+1                                                                
                  DO L=NBC(1),2,-1                                                                       
                        L1=L-1                                                                         
                        LY(1,L)=LY(1,L1)                                                               
                  END DO                                                               
                  LY(1,1)=NHU(J)
            END DO
      END IF      
      
     ! ========================= 16. update the price of crop using CMOD1102.DAT =============== 
      CALL OPENV(KR(26),FCMOD,0,KW(MSO) )       ! database of crop prices                                                    
      DO                                                                      
            READ(KR(26),630,IOSTAT=NFL)ANM,  &  ! ANM = CROP NAME  
                                       X1,   &  ! X1  = GRAIN PRICE ($/t)  
                                       X2       ! X2  = FORAGE PRICE ($/t)                                               
            IF(NFL/=0)EXIT
            IF(LEN_TRIM(ANM)==0)EXIT   ! could not find the crop in the crop list file, so the CMOD1102 is useless for the OAK example
            DO J=1,LC                                    ! LC was given in CPTBL                                           
                  IF(ANM==Crop_Name(J))GO TO 684                                                    
            END DO                                                                              
            CYCLE
        684 seedYieldPrice(J)=X1                                                                     
            forageYieldPrice(J)=X2                                                                     
      END DO
      
      REWIND KR(26)     
      ! ============================== 17. print crop parameters ===============================
      ! These parameters are read from the Crop table file, CROPCOM.DAT
      ! PRINTOUT CROP PARAMETERS                                                       
      CALL Print_Page(1)                                                                  
      CALL Print_CropParms1
      
      IPL=0        
      LRG=0         
      !------------------- update some crop-related variables -----------------------
      
      CALL Update_CropParms(LRG)
      CALL Print_CropParms2(LRG)
      
      ! ==========================18. write pesticide data ==========================
      IF(NDP>0)THEN              ! Number of pesticide applictions
            WRITE(KW(1),'(//1X,A/)')'____________________PESTICIDE DATA____________________'  
            WRITE(KW(1),382)    
            ! PRINTOUT PESTICIDE DATA                                                        
            DO I=1,NDP                                                                     
                  WRITE(KW(1),380)pestName(I),pestSolubity(I),pestHfLifeSoil(I),pestHfLifePlant(I),&
                        WashOff_Frac(I),Coeff_OrgCarbon(I),Pest_Cost(I)  
                  pestSolubity(I)=pestSolubity(I)*10.   
                  pestHfLifeSoil(I)=1.-EXP(-.693/pestHfLifeSoil(I))   
                  pestHfLifePlant(I)=1.-EXP(-.693/pestHfLifePlant(I))  
            END DO 
      END IF     
      
      JD=LY(cropRotationYrs,NBC(cropRotationYrs))               ! What is NBC?                                                
      ICCD=0                                                                         
      IRTC=1                                                                         
      I=1                                                                            
      IF(NBC(cropRotationYrs)>1)GO TO 634                                                        
      IF(IHVD(1,JD)==0)GO TO 677                                                     
      IF(IPLD(1,JD)<IHVD(1,JD))GO TO 677                                             
  634 N1=DOY_realtime+365*(cropRotationYrs-1)  
      
      DO I=1,NBC(cropRotationYrs)                                                                
            IF(N1<IHVD(1,LY(cropRotationYrs,I)))GO TO 677                                            
      END DO     
      
      GO TO 678                                                                      
  677 ICCD=1                                                                         
      IRTC=I 
  678 XLC=LC
      X1=standDeadResOrganism/XLC
      SUM=0.      
      
      ! ---------------- N, P, K in the dead crop residual ---------
      DO J=1,LC
            standDeadResiN(J)=8.29*deadCropRes0
            SUM=SUM+standDeadResiN(J)                                                                  
            standDeadResiP=1.04*deadCropRes0                                                                  
            standDeadResiK=8.29*deadCropRes0                                                                  
            standDeadLignin=.1*deadCropRes0                                                                    
            standCropResi(J)=X1
      END DO
      
      IF(RZ>XX)RZ=XX                ! RZ: soil profile depth in m                                                          
      IF(BIG>XX)BIG=XX              ! BIG ?                                                           
      BTN=totSON+totNO3_N+SUM       ! Total N                                                          
      ! BTNX=BTN                                                                            
      BTC=totSOC                    ! Total C
      ! BTCX=BTC        ! not used                                                                   
      BTP=totLabileP+totActiveP+TOP+initialTotSOP+ soilElement(4)+ standDeadResiP     ! Phrausaus                                            
      BTK=totSolubleK+totExchangeK+totFixK+standDeadResiK                                        ! K                                                    
      KK=numOfTill(1)                                                          
      KC=1   
      
      DO KT=1,KK                                                                     
            IF(tillDOY(1,KT)>=IBD)GO TO 200                               
            II=currentOps(LT(1,KT))                                                             
            IF(II==operationCode(1))KC=KC+1                                                        
      END DO      
      
      KT=numOfTill(1)                                                                      
  200 JT2=LT(1,KT)                                                                   
      ! IF(Climate_ID>0)CALL WREAD     
      ! Fert_trig(BFT0): auto fertilizer trigger (2 options); =1. plant N stress factor (0-1); 
      !                                                       =2. Soil N concentration in root zone g/T                                                   
      NStressTrigger=fertTrigger0
      ! aveBulkDensity: Average soil bulk density of the profile in t/m^3  
      IF(fertTrigger0>1.)NStressTrigger=10.*fertTrigger0*aveBulkDensity*RZ  
      ! APEX-doc Eq 283a    ! RZ: root zone depth/ minimum of soil profile depth/soil profile depth?   
      UB1=RZ*modelPara(54)  ! para(54): exponential coeff in potential water use root growth distribution equation
      UOB=1.-EXP(-UB1)                                                               
      AWC=0.                ! available water content ?                                                                   
      AQV=0.                ! runoff volume ?                                                        
      ARF=0.                ! rainfall      ?                                                  
      Albedo=soilAlbedo                                                                            
      IF(snowWater>5.)Albedo=.6   ! if snow cover exists with 5 mm or greater water content, albdo is set to 0.6 ----APEX-doc P30                                                           
      wiltPoint=0.                                                                                
      subsurfLaterFlow=0.
      soilElement=0.
      IGO=0  
      Crop_Num = varID(JJK0) 
      KC=0  
      IPLD=0  
      MO=MO1  
      DayOfYear=IBD-1  
      IF(DayOfYear<=0) DayOfYear=365  
      call Cal_DL_MaxRd 
      HR0 = DayLength

      ! BEGIN ANNUAL SIMULATION LOOP                                                   
      IRO=IRO0-1                          !IRO: times of planting ?      
      
      !************************************19. Daily simulation ********************************** 
  533 CALL Simulation_Daily(JRT, IRTC, PZW)         ! Major simulations          
      
      !********************************************************************************************
      IF(ISTP==1)GO TO 219                                                           
      IY=MAX(1,IY-1)                                                                 
      XYR=IY                                                                         
      CALL Print_Page(1)      
      ! print soil table after simulation
      
      IF(KFL(9)>0)THEN
            CALL Print_SoilTable(sulfCon, fracVPip, fracHPip)
            
      END IF
      
      WRITE(KW(1),'(//1X,A/)')'____________________FINAL SOIL DATA_____________________'        
          
      CALL Print_SoilVar2(YTP,1)
      
      CALL Sum_WaterNutrient(1)
      
      WRITE(KW(1),'(/T10,A,F7.1,A)')'ERODED SOIL THICKNESS = ',THK,' mm'             
      WRITE(KW(1),'(/T10,A,F7.2,A)')'FINAL WATER CONTENT OF SNOW = ',snowWater,' mm'               
      
      !  PRINTOUT WATER BALANCE                                                         
      CALL Print_WaterVar(SM(4),SM(14),SM(11),SM(16),SM(17),SM(19),SWW,SM(20),SM(84))
      
      IF(irrCode==4)THEN  
            lagoonVol=.1*lagoonVol/(DALG+1.E-10)                                                     
            CALL Check_LagoonWater(SM(23),SM(21),SM(24),SM(78),lagoonVol,lagoonVolX,SM(22))               
            CALL Check_LagoonMaure(SM(79),SM(80),manureWeightX,manureWeight)                                    
      END IF   
      
      RNO3=SM(4)*rainNconX   
      SUM=0.  
      TOT=0.   
      ADD=0.
      AD1=0.
      AD2=0.
      AD3=0.
      AD4=0.   
      
      DO J=1,LC                                                                     
            TOT=TOT+actualCropP(J)                                                               
            ADD=ADD+actualCropK(J)                                                               
            SUM=SUM+actualCropN(J)                                                               
            AD1=AD1+standDeadResiN(J)
            AD2=AD2+1000.*TYLC(J)
            AD3=AD3+420.*totCropBio(J)
            AD4=AD4+420.*standCropResi(J)
      END DO          
      
      FTN=ZN2OG+ZN2OL+ZNO2+totNO3_N+totSON+AD1+standDeadResiOrgN+SUM+ZNH3
      
      X1=SM(49)
      
      IF(deNitri_Method>2) X1=SM(89)
      
      CALL Cal_NP_Balance(BTN,RNO3,SM(43),SM(44),SM(45),SM(46),SM(106),SM(102),SM(105),SM(107),SM(103),&
                          SM(104),X1,SM(59),TYN,SM(52),SM(60),SM(61),SM(50),SM(98),SM(101),FTN,1)                                      
      WRITE(KW(1),636)ZN2OG,ZN2OL,ZNO2,totNO3_N,ZNH3,totSON,AD1,standDeadResiOrgN,SUM
      
      FTC=totSOC+AD3+AD4+standDeadResOrganism*420.
      
      CALL Cal_C_Balance(BTC,SM(77),SM(75),SM(76),SM(74),SM(65),SM(99),AD2,SM(97),FTC)
      WRITE(KW(1),637)totCStructLitt,totCMetabLitt,totCBiomass,totCSlowHumus,totCPassiveHumus,totSOC,AD3,AD4                                     
      
      FTP=totLabileP+TOP+totActiveP+initialTotSOP+Tot_FreshOrgP+standDeadResiP+StandDeadRes_OrgP+TOT    
      
      CALL Cal_NP_Balance(BTP,0.,SM(54),SM(55),0.,SM(57),0.,0.,0.,0.,0.,0.,0.,SM(62),&
           TYP,0.,SM(63),0.,0.,0.,0.,FTP,2)    
      
      WRITE(KW(1),638)totLabileP,TOP,totActiveP,initialTotSOP,Tot_FreshOrgP,standDeadResiP,StandDeadRes_OrgP,TOT                             
      
      FTK=totSolubleK+totExchangeK+totFixK+standDeadResiK+StandDeadRes_OrgK+ADD            
      
      CALL Cal_NP_Balance(BTK,0.,SM(78),SM(79),SM(80),SM(81),0.,0.,0.,0.,0.,0.,0.,&
                         0.,TYK,0.,SM(64),0.,0.,0.,0.,FTK,3)                                                            
      WRITE(KW(1),639)totSolubleK,totExchangeK,totFixK,standDeadResiK,StandDeadRes_OrgK,ADD    
      
      CALL Cal_Salt_Balance(SM(69),SM(72),SM(82),SM(71),SM(70),totSalt0,totSalt)                 
      XQP=NQP                                                                        
      XCN=JCN                                                                        
      XX=XQP+.01                                                                     
      PRSD=PRSD-PRAV*PRAV/XX                                                         
      PRAV=PRAV/XX       
      
      IF(PRSD>0.)PRSD=SQRT(PRSD/XX)                                                  
      CYSD=CYSD-CYAV*CYAV/XX                                                         
      CYAV=CYAV/XX                                                                   
      IF(CYSD>0.)CYSD=SQRT(CYSD/XX)                                                  
      QPS=QPS/XX                                                                     
      TCAV=TCAV/XX                                                                   
      SM(58)=SM(58)/XX                                                               
      !============================ annual values===========================
      !     DETERMINE AVE ANNUAL VALUES    
      SMY(13)=0.  
      SUM=0.                                                                         
      TOT=0.  
      ADD=0.  
      AD2=0.
      
      DO I=1,12  
            X1=IHRL(I)
            IF(X1>0.)THEN
                  monDayLen(I)=monDayLen(I)/X1  
                  monMaxRad(I)=monMaxRad(I)/X1
                  TXMX(I)=TXMX(I)/X1   
                  TXMN(I)=TXMN(I)/X1  
                  TSR(I)=TSR(I)/X1
            END IF  
            TR(I)=TR(I)/XYR  
            TSN(I)=TSN(I)/XYR  
            TSY(I)=TSY(I)/XYR   
            RSY(I)=RSY(I)/XYR  
            TYW(I)=TYW(I)/XYR   
            TQ(I)=TQ(I)/XYR   
            W(I)=W(I)/(TEI(I)+1.E-20)  
            RCM(I)=RCM(I)/(TEI(I)+1.E-20)  
            totHfHourRain(I)=totHfHourRain(I)/CX(I)  
            CX(I)=CX(I)/XYR  
            TEI(I)=TEI(I)/XYR   
            SMY(I)=SRD(I)/XYR   
            SET(I)=SET(I)/XYR    
            TET(I)=TET(I)/XYR   
            ASW(I)=ASW(I)/XYR   
            Inflow(I)=Inflow(I)/XYR   
            TSTL(I)=TSTL(I)/XYR  
            TRHT(I)=TRHT(I)/XYR  
            TAMX(I)=TAMX(I)/XYR   
            TOT=TOT+Inflow(I)    
            ADD=ADD+TAMX(I)   
            SUM=SUM+ASW(I)  
            SMY(13)=SMY(13)+SMY(I)   
            AD2=AD2+CX(I)
      END DO  
      
      SM(15)=SM(15)/(XCN+1.E-10)                                                     
      X1=SM(28)+1.E-10  
      SM(29)=SM(29)/X1   
      SM(37)=SM(37)/X1   
      SM(90)=SM(90)/X1  
      SM(91)=SM(91)/X1  
      SRF2=SM(4)*SM(4)   
      SUM=SUM/12.  
      
      DO K=1,14    
            SM(K)=SM(K)/XYR  
      END DO      
      
      DO K=16,28   
            SM(K)=SM(K)/XYR  
      END DO     
      
      DO K=30,36                                                                     
            SM(K)=SM(K)/XYR  
      END DO     
      
      DO K=38,57   
            SM(K)=SM(K)/XYR   
      END DO      
      
      DO K=59,NSM   
            SM(K)=SM(K)/XYR    
      END DO    
      
      IF(NDP==0.OR. massPestFlag<0.OR.numSimuYear==1)GO TO 212   
      
      CALL Print_Page(1)  
      
      WRITE(KW(1),'(//1X,A/)')'______________PESTICIDE SUMMARY TABLE________________'  
      WRITE(KW(1),460)  
      NY(1)=1 
      NY(2)=INT(.1*XYR+1.5)  ! ??? why
      NY(3)=INT(.5*XYR+1.5)    
      
      ! NDP: times of applying pesticide?
      DO K=1,NDP                                                                 
            WRITE(KW(1),462)pestName(K)                                                        
            DO I=1,5                                                                   
                  DO J=1,numSimuYear                                                                    
                        NX(J)=J                                                                      
                        XYP(J)=APY(I,K,J)                                                            
                        IF(XYP(J)<=1.E-4)THEN                                                        
                        APY(I,K,J)=0.                                                            
                        AYB(I,K,J)=0.                                                            
                        END IF                                                                            
                        XTP(J)=APQ(I,K,J)                                                            
                        IF(XTP(J)>1.E-4)CYCLE                                                        
                        APQ(I,K,J)=0.                                                                
                        AQB(I,K,J)=0.                                                                
                  END DO                                                                              
                  CALL Sort_RealNum(XTP,NX,numSimuYear)                                                       
                  NXX(1,I)=NX(NY(1))                                                             
                  NXX(2,I)=NX(NY(2))                                                             
                  NXX(3,I)=NX(NY(3))                                                             
                  CALL Sort_RealNum(XYP,NX,numSimuYear)                                                       
                  NYY(1,I)=NX(NY(1))                                                             
                  NYY(2,I)=NX(NY(2))                                                             
                  NYY(3,I)=NX(NY(3))
            END DO             
          
            !     PRINTOUT PESTICIDE FREQ SUMMARY                                                
            DO N2=1,3                                                                      
                  WRITE(KW(1),463)(APQ(I,K,NXX(N2,I)),I=1,5)                                   
                  WRITE(KW(1),464)(AQB(I,K,NXX(N2,I)),I=1,5)                                   
                  WRITE(KW(1),465)(APY(I,K,NYY(N2,I)),I=1,5)                                   
                  WRITE(KW(1),466)(AYB(I,K,NYY(N2,I)),I=1,5)                                   
                  IF(N2==1)WRITE(KW(1),474)                                                         
                  IF(N2==2)WRITE(KW(1),473)                                                         
            END DO    
      END DO   

      WRITE(KW(1),'(/1X,A)')'-----AVE ANNUAL VALUES (g/ha)'                          
      DO K=1,NDP                                                                     
            DO L=1,7                                                                     
                  SMAP(L,K)=SMAP(L,K)/XYR                                                      
            END DO                                                                       
            SMAP(10,K)=SMAP(10,K)/XYR                                                    
      END DO   
      
      I1=0                                                                           
      K1=0                                                                           
  418 I1=I1+10                                                                       
      N1=MIN(I1,NDP)                                                                 
      K2=K1+1                                                                        
      N2=MIN(10,NDP-K1)       
      
      !     PRINTOUT PESTICIDE SUMMARY                                                     
      WRITE(KW(1),383)(pestName(K),K=K2,N1)                                              
      IF(massPestFlag==0)GO TO 467                                                           
      WRITE(KW(1),384)HEDP(1),(SMAP(1,K),K=K2,N1)  
      
      I=1     
      
      DO K=K2,N1                                                                     
            PPX(I)=SMAP(2,K)                                                             
            CALL ACOUT(PPX(I),SM(14),1.)                                                 
            I=I+1                                                                        
      END DO    
      
      WRITE(KW(1),387)HEDP(2),(PPX(I),I=1,N2)    
      
      I=1                                                                            
      DO K=K2,N1                                                                     
            PPX(I)=SMAP(3,K)                                                             
            CALL ACOUT(PPX(I),SM(17),1.)                                                 
            I=I+1                                                                        
      END DO       
      
      WRITE(KW(1),387)HEDP(3),(PPX(I),I=1,N2)    
      
      I=1                                                                            
      DO K=K2,N1                                                                     
            PPX(I)=SMAP(4,K)                                                             
            CALL ACOUT(PPX(I),SM(16),1.)                                                 
            I=I+1                                                                        
      END DO   
      
      WRITE(KW(1),387)HEDP(4),(PPX(I),I=1,N2)    
      
      DO L=5,7                                                                       
            WRITE(KW(1),384)HEDP(L),(SMAP(L,K),K=K2,N1)                                  
      END DO          
      
      I=1                                                                            
      DO K=K2,N1                                                                     
            PPX(I)=SMAP(10,K)                                                            
            CALL ACOUT(PPX(I),SM(18),1.)                                                 
            I=I+1                                                                        
      END DO                                                                         
      WRITE(KW(1),387)HEDP(10),(PPX(I),I=1,N2)                                       
      GO TO 469      
      
  467 DO L=1,7                                                                       
            WRITE(KW(1),470)HEDP(L),(SMAP(L,K),K=K2,N1)                                  
      END DO                                                                         
      WRITE(KW(1),470)HEDP(10),(SMAP(10,K),K=K2,N1)                                  
  469 IF(N1>=NDP)GO TO 212                                                           
      K1=I1                                                                          
      GO TO 418                                                                      
  212 PPX(1)=SM(44)                                                                  
      PPX(2)=SM(45)                                                                  
      PPX(3)=SM(46)                                                                  
      PPX(4)=SM(55)    
      
      IF(massPestFlag>0)THEN                                                                 
            CALL ACOUT(PPX(1),SM(14),1000.)                                              
            CALL ACOUT(PPX(2),SM(16),1000.)                                              
            CALL ACOUT(PPX(3),SM(17),1000.)                                              
            CALL ACOUT(PPX(4),SM(14),1000.)                                              
      END IF          
      
      TYN=TYN/XYR                                                                    
      TYP=TYP/XYR                                                                    
      TYC=TYC/XYR                                                                    
      X6=otherCost(1)+otherCost(2)                                                             
      X1=CST1/XYR+X6                                                                 
      X4=CSO1/XYR                                                                    
      X2=VALF1/XYR                                                                   
      X3=X2-X1                                                                       
      SM1=0.      
      
      DO K=1,6                                                                       
            XTP(K)=ISIX(K)                                                               
            SM1=SM1+XTP(K)                                                               
      END DO       
      
      DO K=1,6                                                                       
            XTP(K)=XTP(K)/SM1                                                            
      END DO     
      
      CALL Print_Page(1)                                                                  
      WRITE(KW(1),350)                                                               
      WRITE(KW(1),661)(XTP(K),K=1,6)   
      
      IF(DARF>0..AND.XYR>1.)THEN                                                     
            DARF=DARF-SRF2/XYR                                                                
            IF(DARF>0.) DARF=SQRT(DARF/(XYR-1.))                                          
      ELSE                                                                                
          DARF=0.                                                                      
      END IF     
      
      WRITE(KW(1),662)BARF,SARF,DARF                                                 
      WRITE(KW(1),323)PRB,PRAV,PRSD,QPQB,QPS,NQP                                     
      IF(TCMN>1.E10)TCMN=0.                                                          
      WRITE(KW(1),417)TCAV,TCMN,TCMX                                                 
      WRITE(KW(1),448)CYAV,CYMX,CYSD   
      
      AD1=0.                                                                         
      DO K=1,10                                                                      
            AD1=AD1+initialNO3(K)                                                              
      END DO    
      
      DO K=1,10                                                                      
            initialNO3(K)=initialNO3(K)/AD1                                                          
      END DO  
      
      WRITE(KW(1),234)(initialNO3(K),K=1,10)
      WRITE(KW(1),'(/1X,A)')'-----SOIL LAYER (soilWater-fieldCapacity)/(Porosity-fieldCapacity)'
      WRITE(KW(1),'(T17,A)')'<.1      .1-.25   .25-.50   .50-.75   .75-.90     >.9'
      XX=XYR*365.25
      DO K=1,Actual_SoilLayers
            ISL=Layer_ID(K)
            J1=1
            DO J=6,1,-1
                  YTP(J1)=REAL(NPT(J,ISL))/XX
                  J1=J1+1
            END DO    
            WRITE(KW(1),'(1X,I3,F7.3,6F10.4)')K,Z(ISL),(YTP(I),I=1,6)
      END DO
      
      WRITE(KW(1),'(/1X,A)')'UPWARD FLOW BY LAYER'
      WRITE(KW(1),731)(Z(Layer_ID(K)),SLTP(Layer_ID(K)),K=1,Actual_SoilLayers)
  731 FORMAT(T10,30F10.2)      
      WRITE(KW(1),'(10F8.2)')RUSM     

      !below write .SCN file: summary soil organic C and N table 
      IF(KFL(16)>0)THEN   
            CALL Print_SCN(XYR, XZP)               
      END IF
      
      CALL Print_Page(1)    
      
      WRITE(KW(1),350)                                                               
      WRITE(KW(1),295)                                                               
      WRITE(KW(1),319)IY    
      
      !     PRINTOUT SUMMARY MONTHLY                                                       
      WRITE(KW(1),321)varName(1),TXMX,SM(1),varName(1)                                       
      WRITE(KW(1),321)varName(2),TXMN,SM(2),varName(2)                                       
      WRITE(KW(1),243)varName(4),TR,SM(4),varName(4)                                         
      WRITE(KW(1),321)'DAYP',(SMY(I),I=1,13),'DAYP'                                  
      WRITE(KW(1),243)varName(17),TSN,SM(17),varName(17)                                     
      WRITE(KW(1),243)varName(14),TQ,SM(14),varName(14)                                      
      WRITE(KW(1),321)'soilWaterRZ',ASW,SUM,'soilWaterRZ'                                          
      WRITE(KW(1),311)varName(28),TEI,SM(28),varName(28)                                     
      WRITE(KW(1),224)'ALPH',totHfHourRain,'ALPH'                                              
      WRITE(KW(1),321)varName(29),W,SM(29),varName(29)                                       
      WRITE(KW(1),325)varName(37),RCM,SM(37),varName(37)                                     
      WRITE(KW(1),321)varName(NDVSS),TSY,SM(NDVSS),varName(NDVSS)                            
      WRITE(KW(1),321)varName(36),RSY,SM(36),varName(36)                                     
      WRITE(KW(1),321)varName(7),TET,SM(7),varName(7)                                        
      WRITE(KW(1),224)varName(7),U10MX,varName(7)                                            
      WRITE(KW(1),243)'DAYW',TAMX,ADD,'DAYW'                                         
      WRITE(KW(1),243)varName(39),TRHT,SM(39),varName(39)                                    
      WRITE(KW(1),243)'DAYQ',CX,AD2,'DAYQ'                                           
      WRITE(KW(1),224)HEDC(7),TSTL,HEDC(7)                                           
      WRITE(KW(1),311)varName(42),TYW,SM(42),varName(42)                                     
      WRITE(KW(1),321)varName(20),Inflow,TOT,varName(20)                                        
      WRITE(KW(1),311)varName(10),SET,SM(10),varName(10)                                     
      WRITE(KW(1),311)varName(3),TSR,SM(3),varName(3)                                        
      WRITE(KW(1),224)'DayLength',monDayLen,'DayLength'                                             
      WRITE(KW(1),'(/1X,A)')'-----AVE ANNUAL VALUES'                                 
      WRITE(KW(1),293)IY,(varName(stateVarID(K)),SM(stateVarID(K)),K=5,numPrintVar),'COST',X1,'RTRN',X2
      WRITE(KW(1),369)(varName(conVarID(K)),PPX(K),K=1,numConVar)                                    
      SMY(3)=0.                                                                      
      SMY(4)=0.                                                                      
      SMY(1)=SM(59)+SM(60)+SM(61)                                                    
      SMY(2)=SM(62)+SM(63)                                                           
      SMY(5)=1000.*labileP(LD1)/WT(LD1)                                                   
      X1=LC                                                                          
      X6=X6/X1    
      ! print crop related values
      DO K=1,LC                                                                  
            IF(NCR(K)==0)CYCLE
            XX=MIN(NCR(K),IY)
            X3=TETG(K)+.01                                                                   
            TETG(K)=1000.*(TYL1(K)+TYL2(K))/X3                                       
            TYL1(K)=TYL1(K)/XX                                                                  
            TYL2(K)=TYL2(K)/XX                                                             
            TYLN(K)=TYLN(K)/XX                                                             
            TYLP(K)=TYLP(K)/XX
            TYLC(K)=TYLC(K)/XX                    
            TYLK(K)=TYLK(K)/XX                                                             
            SMY(3)=SMY(3)+TYLN(K)                                                          
            SMY(4)=SMY(4)+TYLP(K)                                                          
            TDM(K)=TDM(K)/XX                                                               
            THU(K)=THU(K)/XX                                                                    
            IF(cropCode(K)==plantCategoryCode(3).OR.cropCode(K)==plantCategoryCode(6).OR.cropCode(K)==&
            plantCategoryCode(7).OR.cropCode(K)==plantCategoryCode(8).OR.cropCode(K)==plantCategoryCode(10))XX=IY                                           
            PSTM(K)=PSTM(K)/XX                                                                  
            TCSO(K)=TCSO(K)/XX                                                             
            TCST(K)=TCST(K)/XX+X6                                                          
            TVAL(K)=TVAL(K)/XX                                                                  
            TFTK(K)=TFTK(K)/XX                                                                  
            TFTN(K)=TFTN(K)/XX                              
            TFTP(K)=TFTP(K)/XX                                                             
            TRA(K)=TRA(K)/XX                                                               
            TRD(K)=TRD(K)/XX                                                               
            TVIR(K)=TVIR(K)/XX                                                             
            TIRL(K)=TIRL(K)/XX                                                             
            TCAW(K)=TCAW(K)/XX                                                             
            TCRF(K)=TCRF(K)/XX                                                             
            TCQV(K)=TCQV(K)/XX 
            X3=X3/XX                                                            
            DO L=1,4                                                                       
                  TSFC(L,K)=TSFC(L,K)/XX                                                       
                  STDA(L,K)=STDA(L,K)/XX                                                       
            END DO                                                                         
            DO L=5,7                                                                       
                  TSFC(L,K)=TSFC(L,K)/XX                                                       
            END DO                                                                         
            X1=TVAL(K)-TCST(K)                                                             
            X2=TVAL(K)-TCSO(K)    
          
          !     PRINTOUT CROP SUMMARY                                                          
            WRITE(KW(1),326)Crop_Name(K),TYL1(K),TYL2(K),TDM(K),TYLN(K),TYLP(K),&               
                        TYLK(K),TYLC(K),TFTN(K),TFTP(K),TFTK(K),TVIR(K),TIRL(K),TCAW(K),&
                        X3,TETG(K),TRA(K),THU(K),PSTM(K),TCST(K),TCSO(K),TVAL(K),X1,X2                           
            WRITE(KW(1),577)(TSFC(L,K),L=1,7),(STDA(L,K),L=1,3)  
          
      END DO
      
      !     PRINTOUT SUMMARY                                                               
      IF(WP5_ID>0)REWIND KR(20)                                                        
      IF(KFL(3)>0)THEN                                                               
            WRITE(KW(3),578)(TITLE(I),I=21,35),IYER,IMON,IDAY,IY                              
            WRITE(KW(3),560) varName(4),varName(10),varName(11),varName(14),varName(16),varName(17),&             
                           varName(29),varName(NDVSS),varName(42),varName(48),varName(47),varName(50), &
                           varName(51),varName(52),varName(49),varName(43),varName(44),varName(45),&
                           varName(46),varName(56),varName(54),varName(55),varName(57),varName(66),varName(77)
          
            WRITE(KW(3),498) SM(4), SM(10), SM(11), SM(14), SM(16), SM(17), SM(29), SM(NDVSS), SM(42), SM(48), &
                           SM(47),SM(50),SM(51),SM(52),SM(49),SM(43),SM(44),SM(45),SM(46),SM(56),SM(54),SM(55), &
                           SM(57),SM(66),SM(77),&
                           (pestName(K),SMAP(1,K),K=1,10),&
                           (Crop_Name(K),TYL1(K),TYL2(K),TDM(K),TYLN(K),TYLP(K),TYLC(K),TFTN(K),TFTP(K),TVIR(K),TIRL(K),TETG(K),&
                           TCAW(K),TCRF(K),TCQV(K),THU(K),potentialHeatUnit(K,IHU(K)),TCST(K),TCSO(K),TVAL(K),PSTM(K),&
                           (TSFC(L,K),L=1,7),(STDA(L,K),L=1,3),K=1,LC)
      END IF  
      
      IF(KFL(MSO+1)>0)WRITE(KW(MSO+1),694)siteName,Crop_Name(1),TYL1(1),(SM(annualVarID(J)),J=1,numAnnuVar)      
      
      IF(KFL(NGF)>0)THEN
            X2=SM(49)+SM(52)                                                                    
            X3=SM(43)+SM(44)+SM(45)+SM(46)                                                      
            DO K=1,8                                                                            
                  Var_GS(K)=Var_GS(K)/XYR                                                          
            END DO                                                                         
            Var_GS(2)=Var_GS(2)/365.25                                                         
            WRITE(KW(NGF-1),512)XLOG,YLAT,Var_GS(1),TETG(1),Var_GS(3),Var_GS(4),SM(19),Var_GS(5), &
                              SM(11),Var_GS(6),Var_GS(7),TRA(1),TYLN(1),Var_GS(2),X2,X3,SM(42), &
                              SM(NDVSS),SM(3),TDM(1),(TSFC(L,1),L=1,7),SM(46),SM(43),SM(47),SM(49),&
                              SM(50),SM(85),SM(51),SM(52),SM(53),SM(54),SM(56),SM(57),SM(58),SM(59),&
                              SM(60),SM(61),SM(62),SM(63),TYLP(1),TYLK(1),SM(77),SM(4),SM(17),SM(14),&
                              SM(5),SM(6),SM(10),SM(12),SM(13),SM(16),SM(18),SM(15),SM(20),Var_GS(8),&
                              SM(68)
      END IF     
      
      IF(Actual_SoilLayers<3.OR.minThickProfile<1.E-10)GO TO 507                                               
      !  1  JZ(1) = NUMBER OF Y FOR SECOND THRU LAST SIMULATION.                           
      !  2  JZ(2) = 0 FOR NORMAL EROSION OF SOIL PROFILE                                   
      !           = 1 FOR STATIC SOIL PROFILE                                              
      !  3  JZ(3) = ID NUMBER OF WEATHER VARIABLES INPUT.  RAIN=1,  TEMP=2,                
      !            RAD=3,  WIND SPEED=4,  REL HUM=5.  IF ANY VARIABLES ARE printInterval             
      !            RAIN MUST BE INCLUDED.  THUS, IT IS NOT NECESSARY TO SPECIF             
      !            ID=1 UNLESS RAIN IS THE ONLY INPUT VARIABLE.                            
      !            LEAVE BLANK IF ALL VARIABLES ARE GENERATED.  EXAMPLES                   
      !            Climate_ID=1 INPUTS RAIN.                                                      
      !            Climate_ID=23 INPUTS RAIN, TEMP, AND RAD.                                      
      !            Climate_ID=2345 INPUTS ALL 5 VARIABLES.                                        
      !  4  JZ(4) = DAILY WEATHER STA # FROM KR(27) WTHCOM.DAT   
      CALL OPENV(KR(6),FMLRN,0,KW(MSO) )        ! multi run application                          
      READ(KR(6),300)JZ                                                              
      IF(JZ(1)==0)GO TO 507                                                          
      numSimuYear=JZ(1)                                                                         
      erosionMode=JZ(2)                                                                     
      weatherVarID=JZ(3)                                                                      
      dailyWeatherID=JZ(4)                                                                     
      IPY=1
      
      IF(numSimuYear>MYR)THEN
            DEALLOCATE(APQ,APQC,APY,AQB,AYB)
            MYR=numSimuYear
            ALLOCATE(APQ(5,MPS,MYR),APQC(5,MPS,MYR),APY(5,MPS,MYR),AQB(5,MPS,MYR),AYB(5,MPS,MYR))             
      END IF       
      
      IF(correctionForRuns ==0)GO TO 203                                                           
      DO K=1,12                                                                      
            IF(TR(K)>0.)RNCF(K)=monAvePPT(1,K)/TR(K)                                                
            TMNF(K)=(monAveTmax(1,K)-monAveTmin(1,K))/(TXMX(K)-TXMN(K))                                   
            TMXF(K)=monAveTmax(1,K)-TXMX(K)                                                    
            IF(weatherVarID>0)CYCLE                                                               
            IF(correctionForRuns <=NC(K+1)-leapYr)EXIT                                                    
      END DO       
      
      correctionForRuns =0                                                                         
 203  CALL Reset   
      
      IF(weatherVarID==0)THEN
            WRITE(KW(1),'(/T10,A)')'**********RAIN, TEMP, RAD, WIND SPEED, REL HUM ARE GENERATED**********'                                                
            GO TO 533
      END IF     
      
      IF(weatherVarForRuns==0)THEN
            REWIND KR(7)                                                                   
            CALL SearchWeatherList(IDIR(1))                                                       
            !     IYR=Year0
      END IF        
      
      WRITE(KW(1),293)                                                               
      CALL Print_WeatherVar                                                                  
      GO TO 533                                                                      
  507 randomCycles=randomCycles+100                                                                    
      REWIND KR(6)                                                                   
      IF(randomCycles<numSeedInitialized*100)GO TO 538                                                      
      IRO0=IRO0+1                                                                    
      randomCycles=0                                                                          
      IF(IRO0>cropRotationYrs.OR.numSeedInitialized==0)GO TO 532                                               
  538 REWIND KR(1)                                                                   
      IF(weatherVarID0>0)REWIND KR(7)                                                         
      GO TO 531     
! ====================== close all files =============================
  532 CLOSE(KR(1))                                                                   
      IF(weatherVarForRuns==0)CLOSE(KR(7))                                                         
      DO I=2,32                                                                           
          CLOSE(KW(I))                                                                 
      END DO                                                                              
      CALL Timer(1)                                                                      
      CLOSE(KW(1))
      IF(runMode>0)GO TO 219
 
  219 DO I=1,SIZE(KR)                                                                      
            CLOSE(KR(I))                                                                 
      END DO
      
      DO I=2,SIZE(KW)                                                                         
            CLOSE(KW(I))                                                                 
      END DO   
      
      STOP    ! //////////////////////////   end ////////////////////////////////                                                                                                       
  224 FORMAT(1X,A4,12F9.2,11X,A4)                                                                                                   
  229 FORMAT(3X,'1',4X,'2',4X,'3',4X,'4',4X,'5',4X,'6',4X,'7',4X,'8',&               
      4X,'9',3X,'10',3X,'11',3X,'12',3X,'13',3X,'14',3X,'15',3X,'16',3X,&            
      '17',3X,'18',3X,'19',3X,'20')                           
  234 FORMAT(/1X,'-----CURVE NUMBER DISTRIBUTION'/T10,'>95=',F6.2,3X,'>90=',F6.2,3X,'>85=',F6.2,3X,&
             '>80=',F6.2,3X,'>75=',F6.2,3X,'>70=',F6.2,3X,'>65=',F6.2,3X,'>60=',F6.2,3X,'>55=',F6.2,3X,'<55=',F6.2)                                                                                                                
  242 FORMAT(I5,E13.5,2F10.3,2E13.5)                                                 
  243 FORMAT(1X,A4,13F9.1,2X,A4)
  244 FORMAT(I5, E13.5)                                                     
  245 FORMAT(////1X,'-----MISCELLANEOUS PARAMETERS'//T10,'modelPara',9X,'S_Curve1I',4X,'S_Curve2I',4X,&
            'S_Curve1C',7X,'S_Curve2C')                                                           
  287 FORMAT(20(1X,A4))                   
  293 FORMAT(//I5,9(2X,A4,F10.2)/(5X,9(2X,A4,F10.2)))                                  
  295 FORMAT(/1X,'-----AVE MO VALUES')                                                                             
  300 FORMAT(20I4)                                                             
  311 FORMAT(1X,A4,13F9.0,2X,A4)                                                                                                 
  319 FORMAT(45X,'YR=',1X,I4,1X,I4/T12,'JAN',6X,'FEB',6X,'MAR',6X,'APR',6X,'MAY',6X,'JUN',6X,'JUL',6X,'AUG',6X,&
            'SEP',6X,'OCT',6X,'NOV',6X,'DEC',6X,' YR')                                                       
  321 FORMAT(1X,A4,13F9.2,2X,A4)                                                     
  323 FORMAT(/1X,'-----PEAK FLOW RATE STATS(mm/h)',T70,'UNIT PEAKS(1/h)'/T10,'MAX = ',F7.2,5X,'MEAN = ',F6.2,5X,&
            'soilWater Dike_Volum = ',F6.2,5X,'MAX = ',F7.4,5X,'MEAN = ',F7.4,5X,'NO PKS = ',I6)                             
  325 FORMAT(1X,A4,13F9.4,2X,A4)                                                     
  326 FORMAT(/2X,A4,1X,'YLD=',F5.1,'/',F5.1,2X,'BIOM=',F5.1,'t/ha',2X,&              
      'Yield_N=',F5.0,2X,'YLP=',F5.0,2X,'YLK=',F5.0,2X,'YLC=',F5.0,2X,'FN=',&
      F5.0,2X,'FP=',F5.0,2X,'FK=',F5.0,'kg/ha'/T7,'IRGA=',F5.0,2X,'IRDL=',&
      F5.0,2X,'Crop_AvailWater=',F6.0,'mm',2X,'GSET=',F6.0,'mm',2X,'WUEF=',F5.2,&
      'kg/mm',2X,'RAD=',F7.0,'Mj/m2',2X,'HU=',F7.0,2X,'PSTF=',F6.2/T7,&
      'COST=',F7.2,2X,'COOP=',F7.2,2X,'RTRN=',F7.0,2X,'NTRN=',F7.0,2X,&
      'NTRO=',F7.0,'$/ha')                                                                                                         
  350 FORMAT(/1X,'____________________SUMMARY TABLE____________________')                                                          
  369 FORMAT(6X,9(1X,A4,F9.4))                                         
  380 FORMAT(10X,A16,3F10.0,F10.2,F10.0,F10.2)                                       
  382 FORMAT(T43,'HALF LIFE(DAYS)  WASH OFF',T83,'COST'/T13,'NAME',11X,&             
      'SOLUBILITY',5X,'SOIL',4X,'FOLIAGE',3X,'FRACTION',3X,'KOC',7X,'($/KG)')                                                                          
  383 FORMAT(/13X,10(A8,4X))                                                         
  384 FORMAT(5X,A4,10F10.0)                                                          
  387 FORMAT(5X,A4,10F10.4)                                                         
  417 FORMAT(/1X,'-----TIME OF CONCENTRATION STATS(h)'/T10,'MEAN = ',&               
             F6.2,5X,'MIN = ',F6.2,5X,'MAX = ',F6.2)                                                 
  448 FORMAT(/1X,'-----SEDIMENT CONCENTRATION STATS(g/m3)'/T10,'MEAN = ',F10.0,5X,'MAX = ',F10.0,5X,'STDV = ',F10.0)                                
  460 FORMAT(/1X,'-----FREQUENCY & DURATION OF PESTICIDE LOSS(g/ha)'/T35, &            
             'DURATION(d)'/20X,'1',12X,'4',11X,'21',11X,'60',11X,'90')                     
  462 FORMAT(5X,A16,T35,'MAXIMUM')                                                   
  463 FORMAT(8X,'SOL  ',5E13.5)                                                      
  464 FORMAT(8X,'Q+subsurfLaterFlow',5E13.5)                                                      
  465 FORMAT(8X,'ADSRB',5E13.5)                                                      
  466 FORMAT(8X,'SED Y',5E13.5)                                                      
  470 FORMAT(5X,A4,10E12.4)                                                          
  473 FORMAT(T35,'50 % EXCEED')                                                      
  474 FORMAT(T35,'10 % EXCEED')
  498 FORMAT(3F8.0,3F8.1,F8.3,18F8.2,10(1X,A16,F8.0),10(4X,A4,19F8.2,11F8.1)) 
  508 FORMAT(1X,A5,4X,A80)                                                           
  509 FORMAT(10X,A80)                  
  512 FORMAT(1X,4F10.2,5F10.1,2F10.3,F10.0,50F10.2)                                    
  560 FORMAT(25(4X,A4),10(9X,'pestName',8X,'APRT'),10(4X,'Crop_Name',4X,'YLDG',4X&            
      ,'YLDF',4X,'BIOM',4X,' Yield_N',4X,' YLP',4X,' YLC',4X,' FTN',4X,' FTP',&
      4X,'IRGA',4X,'IRDL',4X,'WUEF',4X,' Crop_AvailWater',4X,' CRF',4X,' CQV',4X,' THU'&            
      ,4X,' potentialHeatUnit',4X,'COST',4X,'COOP',4X,'RETN',4X,'PSTF',4X,'  Water_Stress',4X,&              
      '  NS',4X,'  PS',4X,'  KS',4X,'  TS',4X,'  AS',4X,'  SS',4X,'  bulkDensity'&            
      ,4X,'ALSA',4X,' SRT'))                                                        
  577 FORMAT(T7,'STRESS (BIOM) WATER=',F5.1,2X,'N=',F5.1,2X,'P=',F5.1,2X&            
      ,'K=',F5.1,2X,'TEMP=',F5.1,2X,'AIR=',F5.1,2X,'SALT=',F5.1,5X,&                 
      '(ROOT) bulkDensity=',F5.1,2X,'ALSAT=',F5.1,2X,'TEMP=',F5.1,'D')                        
  578 FORMAT(1X,15A4,1X,3I4,1X,I3)                                                                                                                                                                                
  621 FORMAT(/T5,'EEPIC1102_Ph',2X,3I4,2X,2(I2,':'),I2/) 
  630 FORMAT(4X,A4,2F8.0)             
  636 FORMAT(5X,'FINAL CONTENTS:'/5X,'N2OG=',E13.6,2X,'N2OL=',E13.6,2X,&
      'NO2 =',E13.6,2X,'NO3 =',E13.6,2X,'NH3 =',E13.6,2X,'ORGN=',E13.6/5X,&
      'standDeadResiN=',E13.6,2X,'SDON=',E13.6,2X,'actualCropN =',E13.6)              
  637 FORMAT(5X,'FINAL CONTENTS:'/5X,'LSC =',E13.6,2X,'LMC =',E13.6,2X,&             
      'BMC =',E13.6,2X,'HSC =',E13.6,2X,'HPC =',E13.6,2X,'totSOC =',E13.6/&
      5X,'BIOC=',E13.6,2X,'standCropResi =',E13.6)                               
  638 FORMAT(5X,'FINAL CONTENTS:'/5X,'PLAB=',E13.6,2X,'PMS =',E13.6,2X,&             
      'PMA =',E13.6,2X,'PHUM=',E13.6,2X,'FreshOrgP_Residu =',E13.6,2X,'standDeadResiP=',E13.6/&             
      5X,'STOP=',E13.6,2X,'actualCropP =',E13.6)                                             
  639 FORMAT(5X,'FINAL CONTENTS:'/5X,'totSolubleK =',E13.6,2X,'totExchangeK =',E13.6,2X,&             
      'totFixK =',E13.6,2X,'standDeadResiK=',E13.6,2X,'STOK=',E13.6,2X,'actualCropK =',E13.6)              
  661 FORMAT(/1X,'-----SOI STAGE DISTRIBUTION'/T10,'1 = ',F6.3,2X,'2 = '&            
      ,F6.3,2X,'3 = ',F6.3,2X,'4 = ',F6.3,2X,'5 = ',F6.3,2X,'6 = ',F6.3)             
  662 FORMAT(/1X,'-----ANNUAL RAINFALL DISTRIBUTION'/T10,'MAX = ',F6.0,&             
      ' mm',2X,'MIN = ',F6.0,' mm',2X,'STDV = ',F6.0,'mm')                                
  679 FORMAT(/5X,'RUN= ',A8,2X,'SIT#= ',I8,2X,'WP1#= ',I8,2X,'SOL#= ',&             
      I8,2X,'OPS#= ',I8,2X,'WTH#= ',I8)                                                                                                                     
  694 FORMAT(1X,A8,4X,A4,F8.2,F8.0,F8.1,6F8.2,F8.0,F8.4,2F8.3,10(4X,A4,&             
      10F8.2,5F8.0,4F8.2,10F8.1))                
END Program  EEPIC                                                                          