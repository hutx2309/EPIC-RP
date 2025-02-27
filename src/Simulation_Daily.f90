SUBROUTINE Simulation_Daily(JRT, IRTC, PZW)                                                                
      ! EPIC1102                                                                       
      ! THIS SUBPROGRAM DRIVES THE DAILY SIMULATION FOR ANY NUMBER OF YEARS   
      ! ********************************************************************************
      ! A module to calculate biomass accumulation involving photosynthesis process
      ! was added by TXH   04-26-2020
      ! ********************************************************************************
      USE PARM
      USE ReadVarTable_Module
      USE MisceFun_Module
      USE Soil_Module
      USE Hydrology
      USE Erosion_Module
      USE Management
      USE Nutrient_Module
      USE GasDiffuse_Module
      USE Crop_Module
      USE Bio_Acc
      USE Print_Module
      USE ReadWeather_Module 
      
      IMPLICIT NONE
      ! ============================ local variables ====================================                                                                                                                     
      CHARACTER(LEN = 4):: ANMX !, strYear                                                             
      REAL(KIND = 8), DIMENSION(13):: PPX    
      ! --------------- Add by TXH -----------------
      integer:: JBG = 0, &
                Growth_Model = 0 ! 0: traditional RUE method; 1 : photosynthesis method
      
      real:: Mix_Efficiency= 0.0, &  !====== Mix_Efficiency is XXXX, add by TXH =======   
             YGIS=0.,   &    ! CROP YLD                   (t/ha) 
             BGIS=0.,   &    ! Biomass Production         (t/ha)  
             WGIS=0.,   &    ! CROP YLD/Growing Season ET (kg/mm)  
             FGIS= 0.0, &    ! N fertilizer applied       (kg/ha) 
             GS_ET    , &    ! Growing Season ET          (mm)  
             HGIS=0.,   &    ! Crop Harvest Index 
             N_YLD=0.,  &    ! N content of crop yield    (kg/ha)
             P_YLD=0. , &    ! P content of crop yield    (kg/ha)  
             K_YLD=0.,  &    ! K content of crop yield    (kg/ha)
             GS_SoilWater,&  ! root zone soil water (mm) 
             WSGS=0., SNGS=0.,SPGS=0., STGS=0., SAGS=0., SSGS=0.                                               
 
      integer:: JRT, KDT(12,12),KTP(7), JGO(MNC), KTF = 0, FMO1, FMO2,  &
                I, II, I1, I2, I3, IN1, IN2, IPC, IP15, IRLX, IRGX, IRY, IGZ, ISM, IVR, &
                IRTC, IYR1, J1, J2, J, JJ, JJX, K, K1, KK, KTMX, L, L1, actualLayer, LD2, MOP, &
                N1, N2, NN, NN1, NT1, NB1, XDA, XDA1
      real:: XTP(4,12,12),YTP(16), IrrWater, accumIrrWater, TSMY, TSMQ, actualET, CY, groundwaterNConc, CNV, CNH, &
             CSTX, COX, WindDir, &
             BDX, groundWaterLeachRate, X1, X2, X3, X4, X5, XX, XY, &
             percolateRate, PMTE, PZW, QPQ, RFRR, RTO, returnFlowRate, returnFlowNRate, SUM, DG, &
             waterFromLagoon, fertFromLagoon, N_Soil, & ! lagoon related variables
             groundwaterTotN, TPAW, DSNO, WB1, WB2, XPR, XZ, Y1,Z5, ZR, ZZ
      real:: totWaterFromLagoon

      integer, dimension(7):: MNST = [(i, i = 1,7)]
      DATA IrrWater,accumIrrWater,TSMQ,TSMY/4*0./                                             
       
      ! function declare
      REAL, EXTERNAL:: PMAV     
      ! PMAV(FMO1,FMO2,X1,X2)=(FMO1*X1+FMO2*X2)/30.
 
      ! *******************************  Testing NEE *****************************
      ! WRITE(filename,"(I4)") IYR

      ! **************************************************************************  
      ! prepare daily/montly output files
      CALL Write_HeadInfo2                                  
      J2=0                                                                           
      J1=1                                                                           
      IGIS=NGF                             ! NGF= MSO+5  GIS file control 
      JP(Crop_Num)=0                       ! P related control variables 
      waterFromLagoon=0.  
      totWaterFromLagoon=aveWaterFromLagoon        

      ! *********************************************************************************************
      ! Simulations Begin here .....   outer loop: yearly
      !                                inter loop: daily
      ! *********************************************************************************************

yearLoop: DO IY=1,numSimuYear           
         
      IF(runMode==0)WRITE(*, "(T5,'YEAR ',I4,' OF ',I4)") IY, numSimuYear      
      
      ! *********************** for testing ***********************************************
      !WRITE(strYear,"(I4)") IYR
      !OPEN(KW(51), file = 'G:\\EEPIC_TestData\\output\\'//TRIM(strYear)//'RP_LMNEE.txt')
      ! OPEN(KW(51), file = 'F:\Diag_EEPIC\Simu_corn_650_rue.txt') 
      !WRITE(KW(51), "(1X, 'Year', 1X, 'Mon', 1X,'Date', 3X,'  RSPC', 3X,'  NPPC', 3X, ' NEE')")    
      ! ***********************************************************************************

      DayOfYear=IBD                  ! IBD from Main
      KDT=0                          ! initialize KDT                                                                       
      GS_SoilWater=0.                                                                             
      HSM=0.                         ! simulated heat unit
      IAP=0                          ! ???
      FLF0=0.                        ! ???
      AboveGround_Biomass0=0.        
      
      !-------------------- ! CO2 concentration calculation method -------------------
      IF(CO2_Method==1)THEN                                                                     
            IF(IYX<25)THEN       ! IYX: number of years relative to 1880                                                          
                  CO2=280.33                                                               
            ELSE                                                                         
                  X1=IYX                                                                   
                  CO2=280.33-X1*(.1879-X1*.0077)                                           
            END IF                                                                            
      END IF
      CALL Print_Page(1)                                                                  
      WRITE(KW(1),"(T10,'ATMOS CO2 = ',F5.0)") CO2   
      !-------------------------------------------------------------------------------
      
      IRO=IRO+1                      ! Crop Rotation                                             
      IF(IRO>cropRotationYrs)IRO=1   ! cropRotationYrs (yr) : 1 ~ 10                                     
      IF(IRO<cropRotationYrs)THEN                                                                
            IRLX=0                                                                       
      ELSE                                                                           
            IRL=IRL+1                                                                    
            IRLX=IRL                                                                     
      END IF                                                                         
      
      NN=NBC(IRO)                     ! What is NBC?                                                              
      NN1=numOfTill(IRO)              ! number of tillage in a year     
      
      DO I=1,LC                       ! LC: crops?                                                           
      RUE(I)=100.*CO2/(CO2+EXP(CO2EffOnBio(1,I)-CO2*CO2EffOnBio(2,I))) ! Eq. 275c potential growth in APEX-doc                  
      END DO                                                                         
      
      IF(KFL(7)>0)WRITE(KW(7),'(4X,4I4)')IRUN,IRO0,randomCycles,IY   
      
      DO I=1,MNC                                               
            JE(I)=MNC+1                ! ? JE = = 31 ?                                                           
            IF(KG(I)==0)CYCLE          ! KG seems to save crop ID                                             
            DO J=1,NN                                                                    
                  IF(I==LY(IRO,J))EXIT                                                     
            END DO                                                                       
            JE(J)=I                                                                      
      END DO          
      
      totWaterFromLagoon = waterFromLagoon + totWaterFromLagoon                                                                 
      aveWaterFromLagoon = totWaterFromLagoon/IY                                                                   
      waterFromLagoon=0.                                                                        
      fertFromLagoon=0.                                                                        
      KOMP=0                                                                         
      IGZ=0                                                                          
      KTT=0                                                                          
      KTMX=1                                                                         
      JJ=IBD-1                         ! Yesterday                                                   
      JT1=LT(IRO,KT)                                                                 
                                                 
      WRITE(KW(1),'(//T5,A)')'Y M D  OPERATION'                                                                                                                
      
      NumOfDays=366-leapYr                         ! Total number of days in a year                                                 
      JP(1)=0                                                                        
      JP(2)=0                                                                        
      IPC=1                                        ! Counting dates for printing       
      
      IF(P_apply_mode>0)THEN
            IF(APBC<20.)THEN
                  APMU=2.25*(30.-APBC)
                  Crop_Num=LY(IRO,1)
                  IF(APMU>45.)CALL Apply_Fert(5,JT1)
            END IF
      END IF    
      
      ! =================== Soil Erodibility Factor: EK = X1*X2*X3*X4 for topsoil layer ============
      ! ------------------- at the start of each year of simulation -----------------------
      
      X1=.2+.3*EXP(-.0256*sandFrac(LD1)*(1.-.01*siltFrac(LD1))) ! Eq. 115a in APEX-doc Pg.35                        
      X2=(siltFrac(LD1)/(clayFrac(LD1)+siltFrac(LD1)))**.3                              
      X5=.1*SOC(LD1)/WT(LD1)                                                          
      IF(X5>5.)THEN                                                                  
            X3=.75                                                                       
      ELSE                                                                           
            X3=1.-.25*X5/(X5+EXP(3.718-2.947*X5))                   ! Eq.115c                                     
      END IF                                                                         
      XX=1.-.01*sandFrac(LD1)                                                             
      X4=1.-.7*XX/(XX+EXP(-5.509+22.899*XX))                                     
      EK=X1*X2*X3*X4                           
      
      ! SL: slope length and steepness factor
      ! conservFactor: conservation practice factor (=0 eliminates water erosion) 
      USL=EK*SL*conservFactor                                               
      
      ! -------------------- ????? cann't find referene for following five lines  ????? -----------------
      SUM=(sandFrac(LD1)*.0247-siltFrac(LD1)*3.65-clayFrac(LD1)*6.908)/100.                        
      DG=EXP(SUM)                                                                    
      REK=9.811*(.0034+.0405*EXP(-.5*((LOG10(DG)+1.659)/.7101)**2))                    
      RSK=REK*conservFactor*RSF                                                                     
      RSLK=RSK*RLF     
      !  Soil Erodibility by Wind  
      IF(erosionMode==0.OR.IY==1) Wind_Erod = Wind_Erodibility(sandFrac(LD1),siltFrac(LD1), clayFrac(LD1), CaCO3(LD1))        
      ! =======================================================================================
      NT1=0                                  ! ????                                                                        
      IF(WP5_ID==0)THEN                      ! Weather station # from WPM5US
            IWIX=1                           ! IWIX: fist day of each month                               
      ELSE
            READ(KR(20),533,IOSTAT=NFL)IWIX
            IF(NFL/=0)THEN
                  IWIX=1
            ELSE                                       
                  DO I=1,12                                                                      
                        IF(IWIX(I)/=0)THEN                                                           
                              IWIX(I)=IWIX(I)+1                                                        
                              CYCLE                                                                         
                        END IF                                                                            
                        IWIX(I)=1                                                                    
                  END DO
            END IF                                                                                       
      END IF
      
      ! *************************************************************************************************      
      ! Daily Simulation.............
      ! *************************************************************************************************    
   dayLoop: DO Date=IBD,NumOfDays   ! Date: Day of month of fertilizer application (1-31) ? explained by EPIC-doc P61                                             
      
      DayOfMon = Cal_DOM(DayOfYear, MO, leapYr)          ! Cal day of the month, DayOfMon          
      XDA=DayOfMon                                                                 
      Mow_dayLap=Mow_dayLap+1               ! Mow_DayLap : day laps from last mow                                       
      NWI=IWIX(MO)                          ! NWI: looks like calculation beginning of a month                                  
      ISIX(NWI)=ISIX(NWI)+1                 ! ISIX(6):      ????
      
      ! ------------------- Prepare for weather simulation: Weather of 15 days before --------------- 
      IP15=Date+15                                                                  
      CALL Cal_Mon(IP15, MOP)               ! MOP is not initialized, which is very risky                                          
      IF(IP15>NumOfDays)IP15=IP15-NumOfDays ! next year                                                       
      FMO2=MIN(30,IP15-NC(MOP))             ! FMO1: half month 1                                        
      FMO1=30-FMO2                          ! FMO2: half month 2                                         
      I1=MOP-1                              ! I1: last month                                         
      IF(I1<1)I1=12                                                                  
      I2=IWIX(I1)                           ! First day of last month                                           
      I3=IWIX(MOP)                          ! First day of this month                                                   
      rainfallVolum=PMAV(FMO1,FMO2,dailyRainInfo(1,I2,I1),dailyRainInfo(1,I3,MOP))     ! rainfall                                  
      rainfallStd=PMAV(FMO1,FMO2,dailyRainInfo(2,I2,I1),dailyRainInfo(2,I3,MOP))       ! std of rainfall                          
      rainfallSkew=PMAV(FMO1,FMO2,dailyRainInfo(3,I2,I1),dailyRainInfo(3,I3,MOP))      ! Skew of rainfall                             
      PBWM=PMAV(FMO1,FMO2,probWetDay(LW,I2,I1),probWetDay(LW,I3,MOP)) ! Wet day probability , what is LW?    LW was initialized as 1                            
      TMXM=PMAV(FMO1,FMO2,monAveTmax(I2,I1),monAveTmax(I3,MOP))       ! Tmax                             
      TMNM=PMAV(FMO1,FMO2,monAveTmin(I2,I1),monAveTmin(I3,MOP))       ! Tmin                             
      TXSD=PMAV(FMO1,FMO2,monTmaxStd(I2,I1),monTmaxStd(I3,MOP))       ! std of Tmax                          
      TNSD=PMAV(FMO1,FMO2,monTminStd(I2,I1),monTminStd(I3,MOP))       ! std of Tmin                              
      SRAM=PMAV(FMO1,FMO2,monAveSoilRad(I2,I1),monAveSoilRad(I3,MOP)) ! Soil radiation 
     ! ------------------------------------------------------------------------     
      IPC=MAX(IPC,Date)                        ! What is IPC?                                           
      XDA1=INT(31.-XDA)                             ! XDA: the day of the month

      ! initial value of variables 
      actualLayer=Layer_ID(Actual_SoilLayers)   
      windErodSoil=0.                          ! (t/ha)                                                                
      peakRainRate=0.                                                                                      
      enrichRatio=1.                                                                          
      VAR=0.                                   ! Store output variables  
      growStressFactor=0.                  ! Why does it have a dimension of (7,30)?    7 stress for 30 crops ?                                                                   
      coverLagingFactor=1.                                                                         
      RCF=.9997*RCF   
      dailyFixedN=0.    
      
      ! APEX-doc P58  Eq 273 a and b 
      ! coverLagingFactor: lagging factor for simulating residue and snow cover effects on surface temperature
      ! groundCover: sum of above ground biomass and crop residue in t/ha  
      
      IF(groundCover<10.)coverLagingFactor=groundCover/(groundCover+EXP(S_Curve(5,1)-S_Curve(5,2)*groundCover))  !  S_Curve(5)  Estimate soil cover factor                  
      SNOF=0.                                                                        
      IF(snowWater>0.)THEN                                             !  Snow                                                            
            SNOF=snowWater/(snowWater+EXP(S_Curve(17,1)-S_Curve(17,2)*snowWater))    !  S_Curve(17) Snow cover factor                            
            coverLagingFactor=MAX(SNOF,coverLagingFactor)                                                            
      END IF                                                                         
      coverLagingFactor=coverLagingFactor*.85      
      
      ! *****************************************************************************************
      !  1. Read or generate daily climate data: Tmax, Tmin, Radiation, Rainfall, RH, CO2, U10
     
      IF(weatherVarID>0)THEN
            CALL ReadDailyWeather1
      ELSE    
            ! GENERATE WEATHER VARIABLES totally                                                                
            U10=0.                                                                         
            RHD=0.                                                                         
            CALL Generate_Rainfall(0 )                                                                   
            Rainfall=Rainfall*RNCF(MO)
          
            CALL Para_SimuWeather  ! return RHM, RM, WX, XIM                                                                         
            CALL Generate_AirTem                                                                      
            X1=TMX-TMN                                                                     
            TMX=TMX+TMXF(MO)                                                               
            TMN=TMX-TMNF(MO)*X1                                                            
            CALL Generate_SolarRad                                                                    
            CALL Generate_RH      ! return RHD                                                                      
            CALL Generate_WindSpeed 
      END IF   
       
      ! ----------------------- Store to daily weather variables ---------------------
      SMM(1,MO)=SMM(1,MO)+TMX         ! SMM(1,MO), VAR(1) : Daily Tmax                                                                                                  
      VAR(1)=TMX                                                                     
      SMM(2,MO)=SMM(2,MO)+TMN         ! SMM(2,MO), VAR(2) : Daily Tmean                                                
      VAR(2)=TMN                                                                     
      SMM(3,MO)=SMM(3,MO)+SRAD        ! SMM(3,MO), VAR(3) : Daily Solar Radiation                                                   
      VAR(3)=SRAD     
                                       ! what do SMM(5,MO) and SMM(6,MO) store?
      SMM(7,MO)=SMM(7,MO)+U10          ! SMM(7,MO), VAR(7) : Daily Wind Speed                                               
      VAR(7)=U10                                                                     
      SMM(8,MO)=SMM(8,MO)+RHD          ! SMM(8,MO), VAR(8) : Daily RH                                             
      VAR(8)=RHD      
                                       ! SMM(9,MO), VAR(9) : Daily VPD
      ! *************************************************************************************** 
      ! ********************************* 2. Wind Erosion *************************************   
      IF(windErosionFactor>0.AND.snowWater<10.)THEN
             WindDir = Wind_Dir()       ! simulate wind direction
            CALL SoiLoss_WindEro(JRT, WindDir)! return: windErodSoil (kg/m) and TLMF                                                                
          
            IF(JRT==0)THEN                 ! JRT: 0, wind erosion calculated ; 1: no wind erosion
                  windErodSoil=windErodSoil*windErosionFactor                                                                      
                  SMM(41,MO)=SMM(41,MO)+Roughtness_Factor ! SMM(41, MO), VAR(41): Roughtness Factor                                      
                  VAR(41)=Roughtness_Factor
            END IF
      END IF
      ! ****************************************************************************************
      ! ****************** 3. Theoretical Daylength and Maximum Radiation **********************  
 
      call Cal_DL_MaxRd                            ! return variables                                                                  
      IHRL(MO)=IHRL(MO)+1                          ! seems to label a month that daylength was calcuated                                           
      monDayLen(MO)=monDayLen(MO)+DayLength    ! accumulate day length                                                  
      monMaxRad(MO)=monMaxRad(MO)+maxRad                       ! accumulate maximum solar radiation                                                                                                                     
      IF(LatChoice >0)THEN                                                       
            SRAM=maxRad*MAX(.8,.21*SQRT(TMXM-TMNM))
            CALL Para_SimuWeather 
            CALL Generate_SolarRad 
      END IF
      ! ****************************************************************************************** 
      ! ******************** 4. Estimate soil temperature based on air temperature ***************      
      TX=(TMN+TMX)/2.
      IF(TX>0.)HSM=HSM+TX                        ! Accumulate the average T, for HUI of trees
       
      CALL Soil_Temperature                      ! returned variables                                                                  
      
      LD2=Layer_ID(2)                                                                   
      SMM(67,MO)=SMM(67,MO)+soilTem(LD2)       ! SMM(67,MO), VAR(67): Soil Temperature in second layer                                 
      VAR(67)=soilTem(LD2)             
      ! *******************************************************************************************
      ! ****************** 5. Hydrology Processes and Water Erosion (sediment) ********************  
      IF(Rainfall>0.)THEN                                                                 
            SRD(MO)=SRD(MO)+1.           ! Radiation : how many wet days in a month for radiation                                             
            SMM(4,MO)=SMM(4,MO)+Rainfall ! SMM(4,MO), VAR(4) : Daily Rainfall                                                      
            VAR(4)=Rainfall                                                                   
            ARF=ARF+Rainfall             ! ARF: Annual Rainfall?                                                            
            snowPackAge=0.                                                                      
      END IF  
      verticalFlowNO3(LD1)=rainNconX*Rainfall! NO3 in the first layer; Ave_N_rainX: Average N concentration in rainfall                                                     
      RNO3=verticalFlowNO3(LD1)                ! NO3 from rain to surface layer of soil     
      
      ! ------------------------------  5.1  Snow fall and melting --------------------------------    
      totRunoff=0.                             ! rainfall + snow melting    
      snowMelting=0.                                    ! Snow melting
      snowPackAge=snowPackAge+1.              ! Snow Pack Age In Days  
      IF(TX<=0.)THEN
            DSNO=Rainfall                                                                    
            snowWater=snowWater+DSNO                                                                   
            SMM(5,MO)=SMM(5,MO)+Rainfall          ! SMM(5,MO), VAR(5): Rainfall (Snowfall) below 0c                                              
            VAR(5)=Rainfall                       ! This is snowfall                                                 
            Rainfall=0.                                                                         
      ELSE
            IF(snowWater>0..AND.SRAD>10.)CALL Snow_Melt 
            !PREDICTS DAILY SNOW MELT WHEN THE AVERAGE AIR TEMPERATURE EXCEEDS 0 DEGREES C                                              
            totRunoff=Rainfall+snowMelting                                                                    
            SMM(6,MO)=SMM(6,MO)+snowMelting               ! SMM(6,MO), VAR(6): snowMelting--snow melting rate  (mm/d)                                           
            VAR(6)=snowMelting
      END IF       
      ! --------------------------------- 5.2 Runoff -----------------------------------------
      ! Hydrology
      Runoff=0.    
      peakRunoffRate=0.                         ! mm/h    APEX-doc P34    
      rainEnergyFactor=0.  
      alpha05L=.02083                           ! Alpha_0.5L: SWAT-doc P70                                               
      RFRR=0.                                    ! What is RFRR?                                  
      crackFlow=0.
      subsurfLaterFlow=0.   
      totRunoff=totRunoff+accumIrrWater   
      
      ! Nutrients  
      QNO3=0.                                   ! QNO3: amount of NO3-N in runoff 
      mineralizedN=0.                                                                             
      verticleFlowSalt=0.                      ! Salt Volume                                                                                
      sedimentLossP=0.                                                                          
      sedimentLossN=0.                                                                          
      runoffLossLabileP=0.  
 
      ! Erosion
      Sediment=0.                                                                         
      totSoilErosion=0.                                                                        
      coverFactor=0.                                                                         
 
      IVR=1                                      ! IVR determine how peak rainfall rate was calcualted       
      IF(totRunoff>0.)THEN
            IF(Rainfall>accumIrrWater)IVR=0         
            Rainfall=totRunoff                                                                        
            CALL RainEnergy_Factor(IVR)  
            accumIrrWater=0.
            VAR(28)=rainEnergyFactor             ! VAR(28): Rainfall Energer factor                                                             
            X1=clayFrac(LD1)                                                                    
            X5=.1*SOC(LD1)/WT(LD1)                                                         
            X2=MAX(50.,63.+62.7*LOG(X5)+15.7*X1-.25*X1*X1)       ! ?????                         
            RFRR=(Rainfall/X2)**.6                               ! ?????
 
            CALL Predict_Runoff(IVR)      
          
            Tot_Rainfall=Tot_Rainfall+Rainfall               ! Rainfall accumulation                                     
         
            JCN=JCN+1                                                                      
            SMM(15,MO)=SMM(15,MO)+CN                         ! SMM(15, MO), VAR(15):  Curve number                                    
            VAR(15)=CN                                                                     
          
            IF(Runoff>1.)THEN
                  NQP=NQP+1                                    ! counting numbers of runoff > 1.                                                    
                  PRAV=PRAV+peakRunoffRate                    ! accumulating peak runoff rate                                                        
                  PRSD=PRSD+peakRunoffRate*peakRunoffRate                                                                
                  IF(peakRunoffRate>PRB)PRB=peakRunoffRate   ! maximum peak runoff rate                                                            
                  QPQ=peakRunoffRate/Runoff                   ! ratio = peak/runoff                                                   
                  IF(QPQ>QPQB)QPQB=QPQ                         ! maximum ratio                                  
                  QPS=QPS+QPQ                                  ! accumulate ratios                                 
                  TCAV=TCAV+TC                                 ! cal average concention of time                                    
                  IF(TC>TCMX)THEN
                        TCMX=TC                                  ! maximum concention of time
                  ELSE           
                        IF(TC<TCMN)TCMN=TC                       ! minimum concention of time                                      
                  END IF
            END IF
            ! ----------------------  5.3 Water Erosion (SEDIMENT)  ------------------------------------   
            ! 35: pavement and urban area
            IF(erosionMode==0.AND.landUseCode/=35)  CALL Water_Erosion    
          
      END IF    ! end runoff > 0
      
      XX=EXP(-RFRR)                                                                  
      RandRough=MAX(1.,RandRough*XX)                ! ????                                           
      SMM(40,MO)=SMM(40,MO)+RandRough               ! SMM(40,MO), VAR(40):  RandRough: random roughness of soil                                      
      VAR(40)=RandRough   
      
      XX=EXP(-.1*Sediment(3))                  ! erosio                                    
      RidgeHeight2=MAX(.001,RidgeHeight2*XX)         ! ridge height and dike are altered by erosion                                             
      dikeHeight=MAX(.001,dikeHeight*XX)                                                             
      IF(dikeHeight/(DKHL+1.E-5)<.7)dikeHeight=DKHL                                      
     
      SMM(39,MO)=SMM(39,MO)+RidgeHeight2             ! SMM(39,MO), VAR(39): Ridge height                                         
      VAR(39)=RidgeHeight2                                                                   
	  
      totSoilErosion=Sediment(waterEroModel)+windErodSoil                                                                   
      X1=.9*WT(LD1)                                  ! WT: weight                                                       
      IF(totSoilErosion>X1)THEN                      ! erosided soil is not allow to exceed 90% of top soil layer?                                               
            RTO=X1/totSoilErosion                                                                           
            Sediment(waterEroModel)=Sediment(waterEroModel)*RTO                                                      
            windErodSoil=windErodSoil*RTO                                                                    
            totSoilErosion=X1                                                                      
      END IF
      CY=1.E5*Sediment(waterEroModel)/(Runoff+1.E-5)    ! weight of erosided soil/runoff                                                
      CYAV=CYAV+CY                                                                   
      CYSD=CYSD+CY*CY                                   ! What are CYAV, CYSD?                                 
      IF(CY>CYMX)CYMX=CY       
      
      soilWaterRZ=soilWaterRZ+accumIrrWater      
      
      IrrWater=0.         
      ! ----------------------------------- 5.4 Soil Percolation --------------------------------------
      CALL Soil_Percolation 
          
      ! Question 1: what is return flow? how is it connected with groundwater storage and soil layers?
                    ! What are VGN and VGA ??  
      percolateFlow(actualLayer)=percolateFlow(actualLayer)+crackFlow     
      XX=percolateFlow(actualLayer)
      SMM(17,MO)=SMM(17,MO)+XX                    ! SMM(17,MO), VAR(17):  Percolation/vetical flow                  
      VAR(17)=XX                                                  
      SMM(16,MO)=SMM(16,MO)+lateralFlow          ! SMM(16,MO), VAR(16):  Lateral/horizontal flow               
      VAR(16)=lateralFlow
      ! --------------------------------- 5.5 Groundwater Component -----------------------------------
      ! P92 in APEX-doc 
      groundWaterStor=groundWaterStor+XX        ! XX is percolation from the bottom layer to ground water storage  
      X1=groundWaterResidenceTime*groundWaterStor     ! Eq. 368  X1 is percolation rate + return flow rate
      returnFlowRate=X1*returnFlowFrac                                    
      SMM(110,MO)=SMM(110,MO)+returnFlowRate     ! SMM(110,MO), VAR(110): Return flow rate  (mm/d)
      VAR(110)=returnFlowRate                     
      
      percolateRate=X1-returnFlowRate           ! SMM(111,MO), VAR(111): Percolation rate to ground water  (mm/d) 
      SMM(111,MO)=SMM(111,MO)+percolateRate                         
      VAR(111)=percolateRate                   
 
      !  Soluble N content in groundwater 
      groundwaterNConc=solubleNConc/(groundWaterStor+1.E-10)   ! groundwaterNConc: N concentration in ground waterkg/mm
      groundwaterTotN=groundwaterNConc*X1
      CNV=groundwaterTotN/(returnFlowRate*modelPara(84)+percolateRate)
      CNH=CNV*modelPara(84)                    ! ???? Para(84)
      returnFlowNRate=CNH*returnFlowRate                  ! CNH==CGWN: Concentration of N in the groundwater Kg/mm
      SMM(112,MO)=SMM(112,MO)+returnFlowNRate             ! SMM(112,MO), VAR(112): return flow N rate        Kg/ha/d
      VAR(112)=returnFlowNRate                               
      
      groundWaterLeachRate=CNV*percolateRate
      SMM(113,MO)=SMM(113,MO)+groundWaterLeachRate        ! SMM(113,MO), VAR(113): Leaching N rates           Kg/ha/d 
      VAR(113)=groundWaterLeachRate
      
      solubleNConc=MAX(1.E-10,solubleNConc-returnFlowNRate-groundWaterLeachRate) ! solubleNConc
      GWPS = 0.0            ! added by TXH (whether it is initialized daily or yearly needs to be confirmed)
      DO K=1,NDP                                  ! NDP : times of pesticide applied
            X3=GWPS(K)/groundWaterStor            ! ??? questionable: how GWPS was inilized?
            RSPS(K)=MIN(GWPS(K),returnFlowRate*X3)
            GWPS(K)=GWPS(K)-RSPS(K)
            pestVar(12,K)=MIN(GWPS(K),percolateRate*X3)
            GWPS(K)=GWPS(K)-pestVar(12,K)
            pestVar(11,K)=RSPS(K)
            SMMP(11,K,MO)=SMMP(11,K,MO)+RSPS(K)
            SMMP(12,K,MO)=SMMP(12,K,MO)+pestVar(12,K)
      END DO
      groundWaterStor=MAX(1.E-10,groundWaterStor-percolateRate-returnFlowRate)
     
      ! ---------------------------- 5.6 Evaportranspiration ----------------------------------
      CALL Cal_EVP 
                                          
      SMM(10,MO)=SMM(10,MO)+EO          ! SMM(10,MO), VAR(10): EO: potential evaporation in mm/d                                               
      VAR(10)=EO                                                                     
                                        ! SMM(12,MO), VAR(12): Potential plant evaporation   
      actualET=soilEvp                      ! soilEvp: potential soil water evaporation for all layers                       
 
      ! accumRainfall30: accumulated rainfall in mm for 30 days preceding the daily estimate day i in mm
      accumRainfall30=(XDA*(SMM(4,MO)-SMM(14,MO))+XDA1*PMORF)/31.   
      
      !------------------------------ 5.7 Watertable Dynamic ---------------------------------------
      IF(waterTableMinDep<Z(actualLayer))  CALL WaterTable_Dyn 
      
      ! ************************************************************************************************
      ! ************************************ 6. lagoon *************************************************
      IF(irrCode==4)THEN                      ! irrigation from lagoon   no ref equestions found
            manureWeight=manureWeight+manureInputToLagoon                                                                 
            SMM(26,MO)=SMM(26,MO)+manureInputToLagoon   ! SMM(26,MO), VAR(26): manure input to lagoon        kg/ha                                          
            VAR(26)=manureInputToLagoon                                                                   
            CALL Lagoon(JRT)  ! return or modified variables
          
            IF(soilWaterRZ-potentialAvailWater<irrTrigger .AND. JRT==0)THEN  ! JRT =0 : indicate irr from lagoon
                  fertRate=MIN(manureWeight,fertCon*irrFromLagoon)                                                        
                  SMM(27,MO)=SMM(27,MO)+fertRate   ! SMM(27,MO), VAR(27): manure output from lagoon   kg/ha                                                   
                  VAR(27)=fertRate                                                                    
                  manureWeight=manureWeight-fertRate                                                                  
                  fertRate=fertRate/areaWshed                                                                    
                  fertFromLagoon=fertFromLagoon+fertRate                                                                  
                  waterFromLagoon=waterFromLagoon+irrFromLagoon                                                                 
                  fertCon=10.*manureWeight*areaWshed/lagoonVol  

                  CALL Apply_Fert(2,JT1)                  ! no ref equestions found                                                       
            
                  IrrWater=irrFromLagoon          
             
                  CALL Cal_IrrWater(IrrWater,1.,0.,JRT,JT1,1)                                           
                  accumIrrWater=accumIrrWater+IrrWater
            END IF
      END IF    
      ! *************************************************************************************************
      ! **************************************** 7. Management ******************************************
      
      ! ---------------------------------------- 7.1 Scheduled Fertilizer -------------------------------
      NFA=NFA+1              ! NFA: days since last fertilizer appliction                                                     
      NII=NII+1              ! NII: days since last irrigation                                                         
      IF(autoManureTrigger>0.AND.NFA>=fertApplyMinInterval) CALL Apply_Fert(2,JT1 )                                       
      IRGX=0                                                                         
      IF(leapYr==0.AND.Date==60)THEN      ! Two months?
            NT1=1                                                                          
            GO TO 30
      END IF      
      ! ----------------------------------------- Management continued ----------------------------------
      NB1=KT                                                                         
      KHV=0                                                                          
      DO KT=NB1,NN1                                ! NN1 is number of operations                                    
            IF(KOMP(KT)>0)CYCLE                      ! What is KOMP?                                         
            XHSM=HSM/yearTotHeat                ! HSM: current accumulated heat                                                        
            IF(KTF==0)KTF=KT                                                               
            DO K=1,LC                                                                      
                  IF(JH(IRO,KT)==cropID(K))EXIT    ! KDC: Crop identification number (1-100)
            END DO                                                                         
            IF(K>LC)K=1                                                                            
            Crop_Num=K                                                                          
            IF(operationMode >0 .AND. IY>1 )GO TO 443                                                   
            IF(Date<tillDOY(IRO,KT)+NT1)GO TO 25       
          
        443 IF(KG(Crop_Num)==0.AND.JPL(Crop_Num)==0)THEN                                             
                  IF(IY==1.AND.currentOps(KT)==operationCode(2))GO TO 590                                                           
            ELSE                                                                           
                  XHSM=HU(Crop_Num)/potentialHeatUnit(Crop_Num,IHU(Crop_Num))                                               
            END IF                                                                                                                              
            IF(XHSM>=fracHeatUnit(IRO,KT))GO TO 590   ! fracHeatUnit: fraction of potential heat units at which operation takes place
            IF(totCropBio(Crop_Num)/(totCropBioX(Crop_Num)+1.E-5)<.99.AND.&
                               rootPartition(1,Crop_Num)<rootPartition(2,Crop_Num))GO TO 590                                                     
            IF(cropCode(Crop_Num)==plantCategoryCode(7).OR.cropCode(Crop_Num)==plantCategoryCode(8).OR.&
                               cropCode(Crop_Num)==plantCategoryCode(10))GO TO 25                                                        
            IF(YLAT>0.)THEN
                  IF(MO==12)GO TO 121
            ELSE
                  IF((MO==7.AND.plantMinT(Crop_Num)>7.).OR.(MO==12.AND.DayOfMon==30))GO TO 121                                                                               
            END IF    
            GO TO 25
        590 IF(plawDepSoilWater/plawDepFieldCapa<modelPara(59))GO TO 121   ! EPIC0810 manual P57  para(59): Soil water to delay tillage                                             
	                  !WRITE(KW(1),589)IY,MO,DayOfMon,plawDepSoilWater,plawDepFieldCapa                                                 
            IF(MO<12)GO TO 25                                                            
        121 JT1=LT(IRO,KT)                                                                 
            KOMP(KT)=1                                                                     
            IF(KT>KTMX)KTMX=KT                                                             
            CSTX=COTL(JT1)                                                                 
            COX=COOP(JT1)    
            ! ---------------------------------- 7.2 Scheduled Irrigation -----------------------------------
            IF(currentOps(JT1)==operationCode(8))THEN
                  IF(IRGX>0)THEN
                        KOMP(KT)=0                                                                     
                        CYCLE
                  END IF                                                             
                  IrrWater=irrVolumeArr(IRO,KT)                                                               
                  IRY=0                                                                          
                  IF(IrrWater>0.)IRY=1                                                                
                  CALL Cal_IrrWater(IrrWater,EFM(JT1),tillDepth(JT1),JRT,JT1,IRY)                                 
                  IF(JRT/=0)THEN                                                                 
                        KOMP(KT)=0                                                                        
                        CYCLE                                                                          
                  ELSE                                                                                
                        IRGX=1                                                                       
                        accumIrrWater=accumIrrWater+IrrWater                                                                
                        GO TO 469                                                                    
                  END IF 
            END IF                                                                        
            ! ------------------------------------ 7.3 Scheduled Fertilization ------------------------------
            IF(currentOps(JT1)==operationCode(9))THEN
                  CALL Apply_Fert(1, JT1)
                  IF(Date/=J1.OR.NBT(JT1)>0.OR.NBE(JT1)/=J2)THEN
                        J1=Date                                                                         
                        J2=NBE(JT1)                                                                    
                        GO TO 469
                  ELSE                                                                                                   
                        CSTX=0.                                                                        
                        COX=0.                                                                         
                  END IF
            END IF     
            ! ------------------------------------ 7.4 Scheduled Pesticide -----------------------------------
            IF(currentOps(JT1)==operationCode(7))THEN
                  IF(Runoff>modelPara(58))THEN
                        WRITE(KW(1),588)IY,MO,DayOfMon,Runoff
                        KOMP(KT)=0                                                                     
                        CYCLE                    
                  ELSE                                                       
                        CALL Apply_Pesticide                                                                        
                        IF(Date==J1.AND.NBT(JT1)==0.AND.NBE(JT1)==J2)THEN
                              CSTX=0.                                                                        
                              COX=0.
                        END IF
                  END IF
            END IF
            ! --------------------------------------- 7.5 Scheduled Lime----------------------------------
            limeRate=0.
            IF(currentOps(JT1)==operationCode(27))THEN
                  limeRate=TLMA(IRO,KT)
                  CALL Apply_Lime(limeRate ) 
                  SMM(66,MO)=SMM(66,MO)+limeRate   ! Limestone applied (CaCO3 equivalent)    T/ha                                                  
                  VAR(66)=limeRate                                                                    
                  X3=limeRate*limeCost                                                                 
                  COST=COST+X3
                  X1=COTL(JT1)                                                              
                  X2=X1-COOP(JT1)                                                           
                  COST=COST+X1                                                               
                  CSFX=CSFX+X2                                                               
                  SMM(92,MO1)=SMM(92,MO1)+Fuel_Use(JT1)                                              
                  SMY(92)=SMY(92)+Fuel_Use(JT1)
                  IF(NOP>0)THEN 
                        WRITE(KW(1),63)IYR,MO1,DayOfMon,equipmentName(JT1),limeRate,Crop_Name(Crop_Num),&
                        X1,X2
                  END IF                                                         
                  IF(KFL(20)>0)THEN                                                    
                        WRITE(KW(20),567)IYR,MO1,DayOfMon,cropID(Crop_Num),X3,X3,limeRate                            
                        WRITE(KW(20),666)IYR,MO1,DayOfMon,equipmentName(JT1),cropID(Crop_Num),&
                                   currentOps(JT1),NBE(JT1),NBT(JT1),X1,X2,Fuel_Use(JT1)
                  END IF                                          
            END IF                                                                                                                                                                
            IF(currentOps(JT1)==operationCode(20))THEN
                  KOMP(NB1)=1                                                                    
                  KTF=KTF+1                                                                           
                  CYCLE
            END IF                                                                     
            J1=Date                                                                         
            J2=NBE(JT1)      
            ! ----------------------------- 7.6 Tillage & Operations ----------------------------------
            ! 35: pavement and urban area          
       469  IF(landUseCode/=35) CALL Till_OPS(CSTX,COX,JRT)
          
            IF(JRT>0)GO TO 25
            SMM(96,MO)=SMM(96,MO)+carbonEmiss(JT1)                                                              
            IF(furrowFlag/=0.AND.DKHL>0.)THEN                                                     
                  dikeHeight=DKHL                                                                     
                  IF(NOP>0)WRITE(KW(1),92)IYR,MO,DayOfMon,dikeHeight,dikeInterval,XHSM                             
            END IF                                                                         
            COST=COST+CSTX                                                                 
            CSFX=CSFX+COX                                                                  
            JT2=JT1
      END DO                                                                                  
      KTF=NB1                                                                        
   25 IF(KTT>0)IGZ=IGZ+1                                                             
      KT=KTF                                                                         
      irrTrigger=irrTriggerArr(IRO,KTMX)                                                              
      irrRunoffRatio=QIR(IRO,KTMX)                                                              
      minCFactor=fertCFactirArr(IRO,KTMX)
      JT1=LT(IRO,KT)                                                                 
      KTF=0                                                                          
      IF(ABS(irrTrigger)<1.E-5)GO TO 30                                                     
      IF(irrTrigger>0.)GO TO 122                                                            
      IF(soilWaterRZ-potentialAvailWater>irrTrigger)GO TO 30                                                       
      GO TO 29                                                                       
  122 IF(irrTrigger<1.)THEN
            IF(waterStress>irrTrigger)GO TO 30                                                             
      ELSE                                                             
            WTN = Cal_WTN()                                                                    
            IF(WTN<irrTrigger)GO TO 30                                                            
      END IF
   29 IF(irrCode>=3.AND.irrCode/=5)THEN
            IF(irrCode==4.OR.NFA<fertApplyMinInterval)GO TO 30                                                  
            CALL Apply_Fert(2, JT1)                                                             
            IF(JRT>0)GO TO 30
      END IF                                                               
      CALL Cal_IrrWater(IrrWater,EFM(autoIrrCode),tillDepth(autoIrrCode),JRT,autoIrrCode,0) 
     
      accumIrrWater=accumIrrWater+IrrWater              
      
      ! *******************************************************************************************************
      ! ***************************************** 8. N, P, K cycling ******************************************
      ! 35 : pavement and urban area
   30 IF(landUseCode/=35)CALL Nutrient_Cyc 
      
      ! ???? Here Mix_Efficiency is not assigned a value ?????
      CALL Tillage_Mix(Mix_Efficiency, ZMIX, 1, 1)
      
      VAR(104)=verticalFlowN2O(actualLayer)                   ! SMM(104, MO), VAR(104): N2O in vertical flow 
      SMM(104,MO)=SMM(104,MO)+verticalFlowN2O(actualLayer)
      SMNIM=0.
      
      IF(deNitri_Method>2)  CALL Solve_GasDiff

      XX = 0.0
      soilResp=0.      ! soil respiration every day at each layer  
      
      DO J=1,Actual_SoilLayers
            ISL=Layer_ID(J)
            IF(soilTem(ISL)>0.)THEN
                  Z5=.5*(Z(ISL)+XX)*1000    
                  IF(NC_Miner_Method==0)THEN            
                        CALL Phoenix(Z5, EAR(ISL))
                  ELSE    
                        CALL CENTURY(Z5)
                  END IF
                  IF(deNitri_Method<3)THEN    ! ------------------- Denitrification ------------------
                        layerThick=1000.*(Z(ISL)-XX)     ! Thickness of a layer
                        IF(deNitri_Method==2)THEN
                              CALL DenitriN1  
                        ELSE
                              weightDenitriN(ISL)=0.
                              DN2=0.
                              denitriedN2O=0.
                              CALL DenitriN2 
                              SN2=SN2+DN2
                        END IF  
                        totDenitriN=totDenitriN+weightDenitriN(ISL)
                        totN2O=totN2O+denitriedN2O
                  END IF
            END IF
            XX=Z(ISL)
      END DO
      
      ! .DBG file --- Daily Denitrification 
      IF(KFL(32)>0)CALL Print_NO3Loss    
      
      SMM(87,MO)=SMM(87,MO)+SMNIM           ! ????
      VAR(87)=SMNIM
      
      IF(irrCode==1)totRunoff=totRunoff+accumIrrWater
      ! ********************************************************************************************************   
      ! ************************ 9. Soil and nutrient loss by water and wind erosion *************************** 
    
      SMM(30,MO)=SMM(30,MO)+Sediment(3)          ! SMM(30,MO), VAR(30): Water Erosion by USLE                                         
      VAR(30)=Sediment(3) 
      SMM(31,MO)=SMM(31,MO)+Sediment(5)          ! SMM(31,MO), VAR(31): Water Erosion by MUSLE                                          
      VAR(31)=Sediment(5)      
      SMM(32,MO)=SMM(32,MO)+Sediment(2)          ! SMM(32,MO), VAR(32): Water Erosion by Onstad-Foster                                      
      VAR(32)=Sediment(2)                                                                 
      SMM(33,MO)=SMM(33,MO)+Sediment(4)          ! SMM(33,MO), VAR(33): Water Erosion by MUSS                                             
      VAR(33)=Sediment(4)                                                                 
      SMM(35,MO)=SMM(35,MO)+Sediment(8)          ! SMM(35,MO), VAR(35): Water Erosion by RUSLE2                                               
      VAR(35)=Sediment(8)                                                                 
      SMM(36,MO)=SMM(36,MO)+Sediment(7)          ! SMM(36,MO), VAR(36): Water Erosion by RUSLE                                         
      VAR(36)=Sediment(7)                                                                 
      SMM(42,MO)=SMM(42,MO)+windErodSoil         ! SMM(42,MO), VAR(42): Wind Erosion                                                     
      VAR(42)=windErodSoil                            
 
      SMM(83,MO)=SMM(83,MO)+Sediment(6)          ! SMM(83,MO), VAR(83): Water Erosion by Modified MUSLE                                           
      VAR(83)=Sediment(6)                                                                 
                                                                
      VAR(29)=coverFactor                                                                           
      VAR(58)=enrichRatio                                                                                                                         
      
      enrichRatioFinal=MIN(enrichRatio*(Sediment(waterEroModel)+windErodSoil)/WT(LD1),.9) 
      
      ! -------------------------- Organic N Loss ----------------------------
      CALL OrgN_Loss  
      
      SMM(43,MO)=SMM(43,MO)+sedimentLossN             ! SMM(43,MO), VAR(43): N loss due to soil erosion                                            
      VAR(43)=sedimentLossN                                                                     
      
      CALL C_Loss                                                 
      
      IF(NDP>0)CALL Pest_Trans 
      
      standLiveDeadPlant=0.            ! standing live and dead plant material   t/ha
      vegCoverFactor=0.                 ! wind ersoin vegative cover factor = C1 * BIOM + C2 * STD + C3 * RSD
                                         ! C1, C2, C3 are coeffcients, 
                                         ! BIOM: above ground bio; STD: standing dead plant residue
                                         ! RSD: flat residue    
      coverFactorPlantPopu=0.                                                                         
      groundCover=cropResidu(LD1)+standDeadResOrganism                                                        
      abvGroundResidue=cropResidu(LD1)+standDeadResOrganism 
      waterStress=1.                                                                          
      N1=1
      VARS(9)=standDeadResOrganism 
      STV(9,MO1)=standDeadResOrganism
      plantEP=0.
      
      ! *************************************************************************************************
      ! ***********************************10. crop growth simulation ************************************** 
      IF(IGO>0)THEN             ! IGO: number of crops growed   JRT: whether a management practice is conducted
            JGO=0                 ! ?????  JGO, JBG, KG, JPL
            JBG=JBG+1        
            IF(JBG>IGO)JBG=1
            I=JBG
            DO J=1,LC
                  IF(KG(J)>0)THEN
                        JGO(I)=KG(J)
                        I=I+1
                        IF(I>IGO)I=1
                  END IF
            END DO
            IN1=0
            DO IN2=1,IGO                                                                    
                  IN1=IN1+1
                  Crop_Num=JGO(IN1)
                  N1=MAX(1,NCP(J))                                                               
                  IF(JPL(Crop_Num)>0)THEN
                        HU(Crop_Num)=HU(Crop_Num)+MAX(0.,surfTem-plantMinT(Crop_Num))                                         
                        IF(plawDepSoilWater/plawDepFieldCapa<modelPara(11).OR.&
                        (HU(Crop_Num)<germinateHeatUnit(Crop_Num).AND.MO<12))CYCLE
                        JPL(Crop_Num)=0                                                                     
                        IF(NOP>0)WRITE(KW(1),89)IYR,MO,DayOfMon,plawDepSoilWater,HU(Crop_Num),XHSM,&
                                          Crop_Name(Crop_Num)                           
                        HU(Crop_Num)=0.                                                                     
                        IGMD(N1,Crop_Num)=IYR*10000+MO*100+DayOfMon
                  END IF                                              
                  standLiveDeadPlant  = standLiveDeadPlant+standCropResi(Crop_Num)     ! standing live and dead plant  t/ha   ref: APEX-doc                                                  
                  abvGroundResidue = abvGroundResidue+standCropResi(Crop_Num)                                                                       
                  vegCoverFactor=vegCoverFactor+windEroCrop(2,Crop_Num)*standCropResi(Crop_Num)                                                         
                  groundCover=groundCover+totCropBio(Crop_Num)-totRootWeight(Crop_Num)+standCropResi(Crop_Num)   !CV: ground cover  including both live bio and residue                                                                       
                  XX=PPL0(Crop_Num)              ! population of plants    plants/m2                                                     
                  coverFactorPlantPopu=MAX(coverFactorPlantPopu,XX/(XX+EXP(S_Curve(15,1)-S_Curve(15,2)*XX)))    ! para(15): express plant population effect on EPIC water erosion cover factor                         
                  vegCoverFactor=vegCoverFactor+windEroCrop(1,Crop_Num)*abvGroundBiom(Crop_Num)                                                    
                  AWC(Crop_Num)=AWC(Crop_Num)+Rainfall-Runoff                                                                 
                  AQV=AQV+Runoff                                                                     
             
                  wiltPoint=0.                     ! U: change to wiltPoint ?                                                    
                  soilSupplyRateN=0.               ! UN: the rate of N supplied by the soil in kg/ha/d                                                       
                  soilSupplyRateK=0.               ! UK: the rate of K supplied by the soil in kg/ha/d                                                      
                  soilSupplyRateP=0.               ! UP: the rate of P supplied by the soil in kg/ha/d                                                         
                  potentialDemandN=0.              ! UNO3                                                           
                  potentialDemandP=0.
                  plantEP=0.                       ! plant evaporation
              
                  potentialBioIncrese(Crop_Num)=0.  ! DDM: the potential increase in biomass in t/ha/d                                                                  
             
                  CALL LAI_HU_RootDep(JRT )
                       
                  ! print daily LAI and LeafWidth
                  ! WRITE(KW(51),"(1X, I4, 4F8.3 )") Date, Current_LAI(Crop_Num), CPHT(Crop_Num), LfWidth(Crop_Num), totCropBio(Crop_Num)   
                  IF(JRT==0)THEN
                        SUN=0.                                                                         
                        SUP=0.                                                                         
                        SUK=0.                                                                         
                        SAT=0.                                                                         
                        CALL Water_Nutrient(JRT )

                        ! Potential Growth
                        IF (Growth_Model == 0 .OR. CPHT(Crop_Num) < 1.E-5 ) THEN
                              ! daily biomass accumulation use traditional methods: RUE or Water Stress
                              CALL Poten_Growth0 
                        ELSE
                              ! Added by TXH    
                              ! daily biomass accumulation use hourly photosynthesis model    
                              CALL Poten_Growth
                  
                        END IF                     
 
                        CALL N_Demand                                                            
                  
                        CALL P_Demand                                                                    
                  
                        CALL K_Demand                                                                  
                  
                        CALL NPK_Supply 
                  
                        CALL Overall_Stress(JRT) 
                  
                        VAR(50)=FixedN_Final                                                                    
                        SMM(13,MO)=SMM(13,MO)+plantEP                                                              
                        actualET=plantEP + actualET                                                                     
                        IF(HUI(Crop_Num)>modelPara(3))THEN                                                       
                              sumActuralEP(Crop_Num)=sumActuralEP(Crop_Num) + plantEP                                                       
                              sumPotentialEP(Crop_Num)=sumPotentialEP(Crop_Num)+EP(Crop_Num)                                                  
                        END IF                                                                         
                        VAR(13)= plantEP   !  plant evaporation                                                               
                        GSEP=GSEP+ plantEP                                                                 
                        ACET(Crop_Num)=ACET(Crop_Num)+ plantEP + soilEvp
                  END IF                                                           
              
                  CALL Actual_Growth 
              
                  cropVar(1,Crop_Num)=HUI(Crop_Num)                                                           
                  cropVar(2,Crop_Num)=Current_LAI(Crop_Num)                                                          
                  cropVar(3,Crop_Num)=RD(Crop_Num)                                                            
                  cropVar(4,Crop_Num)=totRootWeight(Crop_Num)                                                            
                  cropVar(5,Crop_Num)=totCropBio(Crop_Num)                                                            
                  cropVar(6,Crop_Num)=.42*totCropBio(Crop_Num)                                                            
                  cropVar(7,Crop_Num)=abvGroundBiom(Crop_Num)                                                           
                  cropVar(8,Crop_Num)=CPHT(Crop_Num)                                                           
                  cropVar(9,Crop_Num)=standCropResi(Crop_Num)                                                           
                  cropVar(10,Crop_Num)=actualCropN(Crop_Num)                                                           
                  cropVar(11,Crop_Num)=actualCropP(Crop_Num)                                                           
                  cropVar(12,Crop_Num)=actualCropK(Crop_Num)                                                          
              
                  IF(cropCode(Crop_Num)/=plantCategoryCode(7).AND.cropCode(Crop_Num)/=plantCategoryCode(8).AND.&
                        cropCode(Crop_Num)/=plantCategoryCode(10))&
                        standLiveDeadPlant=standLiveDeadPlant+abvGroundBiom(Crop_Num)                                                        
              
                  IF(DOY_realtime==Date.AND.IY==IGSD) CALL Update_CropVar 
                  VAR(100)=VAR(74)-VAR(99)
                  IF(KFL(12)>0)THEN
                        II=IY                                                                          
                        IF(DOY_realtime>0)II=IRLX 
                        X2=1000.*YLD1(N1,Crop_Num)                                                                                                                                          
                        WRITE(KW(12),513)IYR,MO,DayOfMon,II,Crop_Name(Crop_Num),HUI(Crop_Num),AJHI(Crop_Num),&
                                 Current_LAI(Crop_Num),totCropBio(Crop_Num),totRootWeight(Crop_Num),&
                                 abvGroundBiom(Crop_Num),HIX,YLDX,X2,actualCropN(Crop_Num),VAR(99),VAR(100),&
                                (Growing_Stress(J,N1,Crop_Num),J=1,7),(RWT(Layer_ID(K),Crop_Num),K=1,Actual_SoilLayers)                                              
                  END IF            
            END DO
      END IF
      ! *******************************End Crop Growth Simulation*******************************
      ! ****************************************************************************************
      CALL Sum_WaterNutrient(0)      
      
      GS_SoilWater=GS_SoilWater+totSoilWater      ! totSoilWater : soil water   
      ! =================================== N, P, K in EP ======================================
      
      X1=MAX(actualET,EO*EXP(-modelPara(42)*SCI/SMX))                                          
      !     X1=EO*EXP(-modelPara(42)*SCI/SMX)                                          
      ! 	  SCI=SCI+X1-Rainfall+Runoff+LabileP_Con_Ini(actualLayer)+lateralFlow                                                    
      SCI=SCI+X1-Rainfall+Runoff
      
      SCI=MIN(SCI,modelPara(73)*SMX)           ! para(73)  
      
      X1=(accumRainfall30-modelPara(9))/100. ! para(9): pest damage moisture threshold (mm), previous 30-day rainfall minus runoff (25-150)                                                     
      X2=groundCover-modelPara(10)            ! para(10): pest damage cover threshold (t/ha), crop residue + above ground biomass.                                                   
      IF(IGO>0.AND.NFA>=fertApplyMinInterval)THEN                    ! This is the amount of cover reqiured for pests begin to grow 
            IF(NStressTrigger>=1..AND.TNOR<NStressTrigger) CALL Apply_Fert(4, JT1)
      END IF       
      GrowingSeason_Days=GrowingSeason_Days+1                                                                    
      IF(TMN>0.)THEN                                                                 
            IF(X1>0..AND.X2>0.)pestGrowIndx=pestGrowIndx+RHD*TX                                          
      ELSE                                                                           
            pestGrowIndx=pestGrowIndx+TMN                                                                
      END IF           
       
      IF(DOY_realtime/=Date.OR.IY/=IGSD)GO TO 200                                             
      IF(ICCD>0)GO TO 559                                                            
      IGSD=IGSD+1                                                                    
      IYR=IYR-1                                                                      
      ICCD=1                                                                         
      GO TO 200                                                                      
      559 CALL Update_SoilVar                                                                    
      IF(ISTP==1)GO TO 88                                                            
  200 CALL NO3_EP                            ! Upward N by soil evaporation                                                    
      CALL PUpward_EP                        ! Upward P                                                    
      CALL SaltMove_EP                       ! salt  movement                                                  
      VAR(47)=mineralizedN                                                                   
      VAR(48)=SGMN                                                                   
      VAR(49)=totDenitriN                                                                    
      VAR(89)=SN2
      VAR(93)=totN2O
      VAR(94)=DFO2S
      VAR(95)=DFCO2S
      VAR(101)=DFN2OT
      VAR(56)=SMP                                                                    
      VAR(52)=Vol_NH3                                                                   
      VAR(51)=Nitri_NH3                                                                   
      VAR(11)=actualET                                                                    
      VAR(14)=Runoff                                                                     
      VAR(34)=Sediment(1)                                                                 
      VAR(44)=QNO3
      VAR(106)=QNO2
      VAR(107)=QN2O
      VAR(105)=VerticalFlow_NO2(actualLayer)                                                                   
      VAR(46)=verticalFlowNO3(actualLayer)
      VAR(55)=runoffLossLabileP                                                                    
      VAR(45)=TSFNO3                                                                   
      VAR(71)=TSFS                                                                   
      VAR(80)=TSFK
      SMM(100,MO)=SMM(100,MO)+VAR(100)
      VAR(102)=TSFNO2 
      VAR(103)=TSFN2O                                                                  
      VARS(1)=ZNH3                                                                   
      VARS(2)=totNO3_N                                                                   
      VARS(3)=totLabileP                                                                    
      VARS(4)=totSolubleK                                                                    
      VARS(5)=snowWater                                                                    
      VARS(6)=soilWaterRZ                                                                   
      VARS(7)=waterTableHigt                                                                   
      VARS(8)=groundWaterStor                                                                   
      VARS(11)=plawDepSOC                                                                  
      VARS(12)=totSOC                                                                   
      VARS(13)=totStructLitt                                                                   
      VARS(14)=totMetabLitt                                                                   
      VARS(15)=totLgStructLitt                                                                  
      VARS(16)=totCStructLitt                                                                  
      VARS(17)=totCMetabLitt                                                                  
      VARS(18)=totCLgStructLitt                                                                 
      VARS(19)=totNLgStructLitt                                                                
      VARS(20)=totCBiomass                                                                  
      VARS(21)=totCSlowHumus                                                                  
      VARS(22)=totCPassiveHumus                                                                  
      VARS(23)=totNStructLitt                                                                  
      VARS(24)=totNMetabLitt                                                                  
      VARS(25)=totNBiomass                                                                  
      VARS(26)=totNSlowHumus                                                                  
      VARS(27)=totNPassiveHumus                                                                  
      VARS(28)=totSON                                                                   
      VARS(29)=totSalt
      VARS(30)=ZNO2                                                                  
      CALL Prep_SoilVar(YTP)                                                                
      DO I=1,Actual_SoilLayers                                                                    
          J=Layer_ID(I)                                                                     
          SMS(2,J)=SMS(2,J)+soilTem(J)                                                    
      END DO        
      
      ! ========================================== write daily simulation results ==============================
      ! Write to .DCN file: daily soil organic C and N table 
      IF(KFL(15)>0)CALL Print_OrgNC(DayOfMon)        
      
      ! .DPS file : daily pesticide 
      IF(KFL(5)>0.AND.NDP>0)THEN
          II=IY                                                                          
          IF(DOY_realtime/=0)THEN
              II=IRLX                                                                        
              IF(IY/=IGSD)GO TO 561
          END IF                                                          
          DO L=1,NDP                                                                     
              X1=100.*(pestVar(2,L)+pestVar(4,L))/(Runoff+lateralFlow+1.E-5)                               
              WRITE(KW(5),468)IYR,MO,DayOfMon,II,pestName(L),(pestVar(K,L),K=1,10),Runoff,&              
              lateralFlow,percolateFlow(actualLayer),X1                                                           
          END DO                                                                         
      561 pestVar=0.
      END IF            
      ! .DTP file: daily soil temperature
      IF(KFL(10)>0)WRITE(KW(10),107)IYR,MO,DayOfMon,DD,surfTem,(soilTem(Layer_ID(K)),&               
      K=1,Actual_SoilLayers)                                                                      
      X1=plawDepSoilWater/plawDepFieldCapa
      IF(cropCode(Crop_Num)/=plantCategoryCode(7).AND.cropCode(Crop_Num)/=plantCategoryCode(8)&
         .AND.cropCode(Crop_Num)/=plantCategoryCode(10))THEN
          X2=1000.*YLD(Crop_Num)                                                                   
      ELSE
          X2=YLD(Crop_Num)
      END IF    
      ! .DGN file: daily general table
      IF(KFL(17)>0)WRITE(KW(17),107)IYR,MO,DayOfMon,X1,(VAR(dailyVarID(K)),K=1,numDayVar),&
                         ZNH3,totNO3_N,NO3_N_Soil(LD1),percolateFlow(LD1),verticalFlowNO3(LD1),&
                         Albedo,HUI(Crop_Num),AJHI(Crop_Num),Current_LAI(Crop_Num),totCropBio(Crop_Num),&
                         totRootWeight(Crop_Num),abvGroundBiom(Crop_Num), HIX,YLDX,X2,actualCropN(Crop_Num),&
                         VAR(99),VAR(100),(VARS(K),K=23,28)
                  
      ! .DWC file: daily water cycle
      IF(KFL(27)>0)WRITE(KW(27),91)IYR,MO,DayOfMon,(VAR(dailyVarID(K)),K=1,numDayVar),&
      (VARS(monVarID(K)),K=1,numMonVar)      
      
      CALL Cal_SoilWater  
      ! .DHS file
      IF(KFL(28)>0)WRITE(KW(28),682)IYR,MO,DayOfMon,SW15,SW30,SNN15,SNN30,SNA15,SNA30,VAR(4),VAR(10),VAR(11),&
                        (VAR(K),K=13,20),(VARS(K),K=6,8),(Z(Layer_ID(K)),K=1,Actual_SoilLayers),&
                        (soilWater(Layer_ID(K)),K=1,Actual_SoilLayers),&
                        (wiltPoint(Layer_ID(K)),K=1,Actual_SoilLayers),&
                        (soilElement(Layer_ID(K)),K=1,Actual_SoilLayers),&
                        (percolateFlow(Layer_ID(K)),K=1,Actual_SoilLayers),&
                        (subsurfLaterFlow(Layer_ID(K)),K=1,Actual_SoilLayers),&
                        (NO3_N_Soil(Layer_ID(K)),K=1,Actual_SoilLayers),&
                        (soilSupplyRateN(Layer_ID(K)),K=1,Actual_SoilLayers),&
                        (verticalFlowNO3(Layer_ID(K)),K=1,Actual_SoilLayers),&
                        (soilTem(Layer_ID(K)),K=1,Actual_SoilLayers)
      
      N_Soil=ZNO2+totNO3_N+totSON+ZNH3      ! N CONTENT OF SOIL (kg/ha)                                                   
      Var_GS(2)=Var_GS(2)+N_Soil
      ! .DSL file : daily soil table
      IF(KFL(21)>0)THEN                                                              
            WRITE(KW(21),107)IYR,MO,DayOfMon                                                
            CALL Prep_SoilVar(YTP)                                                                 
            CALL Print_SoilVar2(YTP,21)                                                         
      END IF     
      ! --------------------------- Print by different time intervals ----------------------
      IF(printChoice<6)GO TO 44      
      
      IF(printChoice==9)THEN
            IF(IGO==0)GO TO 44                                                             
      END IF
      IF(Date/=IPC)GO TO 44                                                           
      IF(printChoice==8)THEN
            IF(Rainfall<1.)GO TO 44
      END IF          
      IF(printChoice/=7)THEN                                                             
            !    PRINTOUT DAILY                                                                 
            WRITE(KW(1),532)IYR,MO,DayOfMon,(varName(dailyVarID(K)),VAR(dailyVarID(K)),K=1,numDayVar)                     
            DO I=1,NN                                                                      
                  K1=LY(IRO,I)                                                               
                  IF(KG(K1)==0)CYCLE                                                         
                  WRITE(KW(1),105)IYR,MO,DayOfMon,Crop_Name(K1),(HEDC(K),cropVar(K,K1),K=1,19)                                                                        
            END DO                                                                         
            WRITE(KW(1),532)IYR,MO,DayOfMon,(HEDS(K),VARS(K),K=1,13)        
 
      ELSE  
            CALL Prep_SoilVar(YTP)                                                                
            WRITE(KW(1),107)IYR,MO,DayOfMon                                                     
            WRITE(KW(1),101)                                                               
            CALL Print_SoilVar2(YTP,1)                                                              
      END IF      
      IPC=IPC+printInterval             !printInterval was set to 1 at the MAIN  
      !-------------------------------------------------------------------------------------------
   44 SMM(11,MO)=SMM(11,MO)+actualET                                                      
      SMM(47,MO)=SMM(47,MO)+mineralizedN                                                     
      SMM(48,MO)=SMM(48,MO)+SGMN                                                     
      SMM(49,MO)=SMM(49,MO)+totDenitriN                                                      
      SMM(89,MO)=SMM(89,MO)+SN2                                                           
      SMM(56,MO)=SMM(56,MO)+SMP                                                      
      SMM(52,MO)=SMM(52,MO)+Vol_NH3                                                     
      SMM(51,MO)=SMM(51,MO)+Nitri_NH3                                                     
      SMM(14,MO)=SMM(14,MO)+Runoff                                                       
      SMM(34,MO)=SMM(34,MO)+Sediment(1)                                                   
      SMM(44,MO)=SMM(44,MO)+QNO3
      SMM(106,MO)=SMM(106,MO)+QNO2
      SMM(107,MO)=SMM(107,MO)+QN2O
      SMM(105,MO)=SMM(105,MO)+VerticalFlow_NO2(actualLayer)                                                     
      SMM(46,MO)=SMM(46,MO)+verticalFlowNO3(actualLayer)                                                
      SMM(55,MO)=SMM(55,MO)+runoffLossLabileP                                                      
      SMM(45,MO)=SMM(45,MO)+TSFNO3                                                     
      SMM(80,MO)=SMM(80,MO)+TSFK                                                     
      SMM(71,MO)=SMM(71,MO)+TSFS
      SMM(93,MO)=SMM(93,MO)+totN2O
      SMM(94,MO)=SMM(94,MO)+DFO2S
      SMM(95,MO)=SMM(95,MO)+DFCO2S 
      SMM(101,MO)=SMM(101,MO)+DFN2OT                                                          
      SMM(102,MO)=SMM(102,MO)+TSFNO2
      SMM(103,MO)=SMM(103,MO)+TSFN2O
      !CALL NCONT
      IF(Date==DOY_realtime+NT1.AND. weatherVarID==0)THEN
            randomCycles=randomCycles+100                                                                    
            DO KK=1,randomCycles                                                                    
                  DO J=1,20                                                                    
                        XX=Generate_Random(21)                                                             
                        IX(J)=IX(21)                                                             
                  END DO                                                                       
            END DO                                                                         
      END IF    
      
      JRT=0                                                                          
      IF(totSoilErosion>1.E-5.AND. erosionMode==0.AND.Z(actualLayer)>minThickProfile.AND.Actual_SoilLayers>=3)&
            CALL Thickness_Ero(JRT)                                    
      
      IF(JRT>0)GO TO 88           

      ! ////////////////// test NEE /////////////////////////////
      ! WRITE(KW(51), "(3(1X, I4), 10(3X,F12.6))") IYR, MO, Date, (CStructLitt(K), K = 1, 10)    
      ! WRITE(KW(51), "(3(1X, I4),3(3X,F10.4))") IYR, MO1, Date, VAR(74), VAR(99), VAR(100)  
      ! /////////////////////////////////////////////////////////

      DayOfYear=Date+1         ! day increase by 1                                  
      
      CALL Cal_Mon(DayOfYear,MO)                                                              
      
      IF(MO==MO1) CYCLE   
      ! **********************************************************************************************
      ! Continue to the next day, if it is the end of a month, begin to summarize monthly variables
      ! ***********************************************************************************************                                                                   
      PMORF=SMM(4,MO1)-SMM(14,MO1)                                                   
      XX=Date-JJ      ! JJ: store the last day of a month     XX: the number of days in a month, used to averge data over a month.                                                 
      STV(1,MO1)=ZNH3                                                                
      STV(2,MO1)=totNO3_N                                                                
      STV(4,MO1)=totSolubleK                                                                 
      STV(5,MO1)=snowWater                                                                 
      STV(6,MO1)=soilWaterRZ                                                                
      STV(7,MO1)=waterTableHigt                                                                
      STV(8,MO1)=groundWaterStor                                                                
      STV(11,MO1)=plawDepSOC                                                               
      STV(12,MO1)=.001*totSOC                                                           
      STV(13,MO1)=.001*totStructLitt                                                           
      STV(14,MO1)=.001*totMetabLitt                                                           
      STV(15,MO1)=.001*totLgStructLitt                                                          
      STV(16,MO1)=.001*totCStructLitt                                                          
      STV(17,MO1)=.001*totCMetabLitt                                                          
      STV(18,MO1)=.001*totCLgStructLitt                                                         
      STV(19,MO1)=.001*totNLgStructLitt                                                        
      STV(20,MO1)=.001*totCBiomass                                                          
      STV(21,MO1)=.001*totCSlowHumus                                                          
      STV(22,MO1)=.001*totCPassiveHumus                                                          
      STV(23,MO1)=.001*totNStructLitt                                                          
      STV(24,MO1)=.001*totNMetabLitt                                                          
      STV(25,MO1)=.001*totNBiomass                                                          
      STV(26,MO1)=.001*totNSlowHumus                                                          
      STV(27,MO1)=.001*totNPassiveHumus                                                          
      STV(28,MO1)=.001*totSON                                                           
      STV(29,MO1)=totSalt
      STV(30,MO1)=ZNO2                                                               
      IF(KFL(6)==0)GO TO 472                                                              
      I1=IY                              ! I1 is the year calcuator                                            
      IF(DOY_realtime/=0)THEN
            I1=IRLX                                                                        
            IF(IY/=IGSD)GO TO 472
      END IF                    
      !======================== print monthly variables ==============================   
      WRITE(KW(6),475)IYR,MO1,I1,(SMM(econVarID(J),MO1),J=1,numEconVar),(STV(K,MO1),&          
      K=4,6)                                                                         
  472 DO 118 K=1,NN                                                                  
      K1=LY(IRO,K)                                                                   
      IF(KG(K1)==0)GO TO 50                                                          
      SMMC(1,K1,MO1)=HUI(K1)                                                         
      SMMC(2,K1,MO1)=Current_LAI(K1)                                                        
      SMMC(3,K1,MO1)=RD(K1)                                                          
      SMMC(4,K1,MO1)=totRootWeight(K1)                                                          
      SMMC(5,K1,MO1)=totCropBio(K1)                                                          
      SMMC(6,K1,MO1)=.42*totCropBio(K1)                                                          
      SMMC(7,K1,MO1)=abvGroundBiom(K1)                                                         
      SMMC(8,K1,MO1)=CPHT(K1)                                                         
      SMMC(10,K1,MO1)=actualCropN(K1)                                                         
      SMMC(11,K1,MO1)=actualCropP(K1)                                                         
      SMMC(12,K1,MO1)=actualCropK(K1)                                                        
      TSTL(MO1)=TSTL(MO1)+abvGroundBiom(K1)                                                    
      IF(KFL(11)==0)GO TO 50                                                         
      II=IY                                                                          
      IF(DOY_realtime/=0)THEN
            II=IRLX                                                                        
            IF(IY/=IGSD)GO TO 50
      END IF                                                           
      WRITE(KW(11),523)IYR,MO1,II,Crop_Name(K1),(SFMO(J,K1),J=1,7),soilWaterRZ,SMM(4,MO1),SMM(11,MO1),&
            SMM(14,MO1),SMM(17,MO1),SMM(16,MO1)   
      
   50 IF(KFL(22)>0)THEN                                                              
            RNO3=SMM(4,MO1)*rainNconX                                                       
            SUM=ZNO2+totNO3_N+ZNH3+totSON+standDeadResiN(K1)+standDeadResiOrgN+actualCropN(K1)                                       
            WRITE(KW(22),565)IYR,MO1,SMM(4,MO1),SMM(10,MO1),SMM(11,MO1), SMM(13,MO1),SMM(14,MO1), &
                (SMM(J,MO1),J=16,20),(STV(J,MO1), J=6,8),RNO3,(SMM(J,MO1),J=43,46),SMM(49,MO1), &
                 SMM(89,MO1),SMM(52,MO1),SMM(85,MO1),SMM(50,MO1),(SMM(J,MO1),J=59,61),actualCropN(K1),&            
                 YLNF(1,K1),Crop_Name(K1),YLD1(1,K1),SUM                                         
      END IF                                                                         
      
      SMMC(9,K1,MO1)=standCropResi(K1)                                                         
      ISM=0                                                                          
      
      DO J=1,7                                 ! What does this block mean?                                      
            KTP(J)=INT(MIN(31.,SFMO(J,K1)+.5))                                                   
            ISM=ISM+KTP(J)                                                             
            SFMO(J,K1)=0.                                                              
      END DO  
      
      IF(ISM==0)GO TO 118                                                            
      
      CALL Sort_IntNum(KTP,MNST,7)                                                        
      KDT(MO1,K1)=KTP(MNST(5))+100*MNST(5)+1000*KTP(MNST(6))+100000*&                 
      MNST(6)+1000000*KTP(MNST(7))+100000000*MNST(7)                                 
118   CONTINUE                                                                       
      
      IF(KFL(25)>0)THEN                                                              
            IF(IY==1.AND.Date==1)WRITE(KW(25),680)areaWshed                                   
            YTP(1)=SMM(14,MO1)+SMM(16,MO1)+SMM(18,MO1)                                 
            YTP(2)=SMM(NDVSS,MO1)                                                      
            YTP(3)=SMM(43,MO1)                                                         
            YTP(4)=SMM(54,MO1)                                                         
            YTP(5)=SMM(44,MO1)+SMM(45,MO1)+SMM(53,MO1)                                 
            YTP(6)=SMM(55,MO1)                                                         
            WRITE(KW(25),679)IY,MO1,YTP                                                
      END IF
      
      STV(10,MO1)=cropResidu(LD1)                                                                        
      VARS(10)=cropResidu(LD1)                                                                        
      STV(3,MO1)=1000.*labileP(LD1)/WT(LD1)                                               
      TXMX(MO1)=TXMX(MO1)+SMM(1,MO1)                                                 
      TXMN(MO1)=TXMN(MO1)+SMM(2,MO1)                                                 
      TSR(MO1)=TSR(MO1)+SMM(3,MO1)                                                   
      SMM(2,MO1)=SMM(2,MO1)/XX                                                       
      SMM(1,MO1)=SMM(1,MO1)/XX                                                       
      SMM(3,MO1)=SMM(3,MO1)/XX                                                       
      SMM(67,MO1)=SMM(67,MO1)/XX                                                     
      SMM(68,MO1)=SMM(68,MO1)/XX                                                          
      SMM(7,MO1)=SMM(7,MO1)/XX                                                       
      TET(MO1)=TET(MO1)+SMM(7,MO1)                                                        
      SMM(8,MO1)=SMM(8,MO1)/XX                                                       
      SMM(9,MO1)=SMM(9,MO1)/XX                                                       
      SMM(39,MO1)=SMM(39,MO1)/XX                                                     
      SMM(40,MO1)=SMM(40,MO1)/XX                                                     
      X1=NWDA-NWD0+1.E-5                                                             
      TAMX(MO1)=TAMX(MO1)+X1                                                         
      SMM(38,MO1)=SMM(38,MO1)/X1                                                     
      SMM(41,MO1)=SMM(41,MO1)/X1                                                     
      NWD0=NWDA                                                                      
      SSW=SSW/XX                                                                     
      ASW(MO1)=ASW(MO1)+SSW                                                          
      SSW=0.                                                                         
      TR(MO1)=TR(MO1)+SMM(4,MO1)                                                     
      TSN(MO1)=TSN(MO1)+SMM(17,MO1)                                                  
      TSY(MO1)=TSY(MO1)+SMM(NDVSS,MO1)                                               
      RSY(MO1)=RSY(MO1)+SMM(36,MO1)                                                  
      TYW(MO1)=TYW(MO1)+SMM(42,MO1)                                                  
      TQ(MO1)=TQ(MO1)+SMM(14,MO1)                                                    
      SET(MO1)=SET(MO1)+SMM(10,MO1)                                                  
      TRHT(MO1)=TRHT(MO1)+SMM(39,MO1)                                                
      JJ=Date                                                                         
      DO K=1,NSM                                                                     
            SMY(K)=SMY(K)+SMM(K,MO1)                                                     
      END DO                                                                         
      IF(NDP>0)THEN                                                                  
            DO K=1,NDP                                                                   
                  SMMP(8,K,MO1)=pestPlantIntercep(K)                                                    
                  DO K1=1,7                                                                
                        SMYP(K1,K)=SMYP(K1,K)+SMMP(K1,K,MO1)                                 
                  END DO                                                                   
                  SMYP(10,K)=SMYP(10,K)+SMMP(10,K,MO1)                                     
            END DO                                                                       
      END IF                                                                         
      W(MO1)=W(MO1)+SMM(29,MO1)                                                      
      RCM(MO1)=RCM(MO1)+SMM(37,MO1)                                                  
      TEI(MO1)=TEI(MO1)+SMM(28,MO1)                                                  
      X1=SMM(28,MO1)+1.E-10                                                               
      SMM(29,MO1)=SMM(29,MO1)/X1                                                     
      SMM(37,MO1)=SMM(37,MO1)/X1                                                     
      SMM(83,MO1)=SMM(83,MO1)/X1                                                     
      SMM(90,MO1)=SMM(90,MO1)/X1                                                          
      SMM(91,MO1)=SMM(91,MO1)/X1                                                          
      X1=JCN-JCN0                                                                                      
      IF(X1>0.)SMM(15,MO1)=SMM(15,MO1)/X1                                                 
      X1=NQP-NQP0                                                                    
      IF(X1>0.)SMM(58,MO1)=SMM(58,MO1)/X1                                                 
      NQP0=NQP                                                                       
      JCN0=JCN                                                                       
      IF(massPestFlag>0)THEN
            PPX(MO1)=SMM(44,MO1)                                                           
            CALL ACOUT(PPX(MO1),SMM(14,MO1),1000.)                                         
            SMM(44,MO1)=PPX(MO1)                                                           
            PPX(MO1)=SMM(45,MO1)                                                           
            CALL ACOUT(PPX(MO1),SMM(16,MO1),1000.)                                         
            SMM(45,MO1)=PPX(MO1)                                                           
            PPX(MO1)=SMM(46,MO1)                                                           
            CALL ACOUT(PPX(MO1),SMM(17,MO1),1000.)                                         
            SMM(46,MO1)=PPX(MO1)                                                           
            PPX(MO1)=SMM(55,MO1)                                                           
            CALL ACOUT(PPX(MO1),SMM(14,MO1),1000.)                                         
            SMM(55,MO1)=PPX(MO1)                 
      END IF                                                    
      ! ------------------------ WRITE MO VALUES AND SUM YEARLY VALUES -------------------------------                                         
      IF(MO>MO1)GO TO 83          ! go to next month                                                       
      IF(limeFlag==0.AND.limeRate<1.E-5)THEN
            CALL Apply_Lime(0.0 )                                                                   
            IF(limeRate>0.)THEN                                                                 
                  X3=limeRate*limeCost                                                                 
                  COST=COST+X3                                                               
                  X1=COTL(IAUL)                                                              
                  X2=X1-COOP(IAUL)                                                           
                  COST=COST+X1                                                               
                  CSFX=CSFX+X2                                                               
                  SMM(92,MO1)=SMM(92,MO1)+Fuel_Use(IAUL)                                              
                  SMY(92)=SMY(92)+Fuel_Use(IAUL)
                        IF(NOP>0)THEN 
                        WRITE(KW(1),63)IYR,MO1,DayOfMon,equipmentName(IAUL),limeRate,Crop_Name(Crop_Num),&
                        X1,X2
                  END IF                                                         
                  IF(KFL(20)>0)THEN                                                    
                        WRITE(KW(20),567)IYR,MO1,DayOfMon,cropID(Crop_Num),X3,X3,limeRate                            
                        WRITE(KW(20),666)IYR,MO1,DayOfMon,equipmentName(IAUL),cropID(Crop_Num),&
                        currentOps(IAUL),NBE(IAUL),NBT(IAUL),X1,X2,Fuel_Use(IAUL)
                  END IF                                          
            END IF
      END IF                                                                         
      SMM(66,MO1)=SMM(66,MO1)+limeRate
      SMY(66)=SMY(66)+limeRate                                                    
      VAR(66)=limeRate                                                                    
      SMY(1)=SMY(1)/12.                                                              
      SMY(2)=SMY(2)/12.                                                              
      SMY(3)=SMY(3)/12.                                                              
      SMY(7)=SMY(7)/12.                                                              
      SMY(8)=SMY(8)/12.                                                              
      SMY(9)=SMY(9)/12.                                                              
      SMY(67)=SMY(67)/12.                                                            
      SMY(68)=SMY(68)/12.                                                                 
      SMY(39)=SMY(39)/12.                                                            
      SMY(40)=SMY(40)/12.                                                            
      SMY(41)=SMY(41)/12.                                                            
      DO K=1,NSM                                                                     
            SM(K)=SM(K)+SMY(K)                                                           
      END DO                                                                         
      
      IF(NDP>0)THEN                                                                  
            DO K=1,NDP                                                                   
                  DO L=1,7                                                                      
                        SMAP(L,K)=SMAP(L,K)+SMYP(L,K)                                        
                  END DO                                                                   
                  SMAP(10,K)=SMAP(10,K)+SMYP(10,K)                                         
            END DO                                                                       
      END IF                                                                         
      
      X1=SMY(28)+1.E-10                                                                   
      SMY(29)=SMY(29)/X1                                                             
      SMY(37)=SMY(37)/X1                                                             
      SMY(90)=SMY(90)/X1                                                                  
      SMY(91)=SMY(91)/X1                                                                  
      
      IF(massPestFlag>0)THEN
            PPX(1)=SMY(44)                                                                 
            PPX(2)=SMY(45)                                                                 
            PPX(3)=SMY(46)                                                                 
            PPX(4)=SMY(49)                                                                 
            CALL ACOUT(PPX(1),SMY(14),1000.)                                               
            CALL ACOUT(PPX(2),SMY(16),1000.)                                               
            CALL ACOUT(PPX(3),SMY(17),1000.)                                               
            CALL ACOUT(PPX(4),SMY(14),1000.)                                               
            SMY(44)=PPX(1)     ! Why kind = 8 for PPX                                                             
            SMY(45)=PPX(2)                                                                 
            SMY(46)=PPX(3)                                                                 
            SMY(49)=PPX(4)                                                                 
      END IF          
      
      X1=JCN-JCN1                                                                    
      SMY(15)=SMY(15)/(X1+1.E-20)                                                    
      JCN1=JCN                                                                       
      X1=NQP-NQP1                                                                    
      SMY(58)=SMY(58)/(X1+1.E-20)                                                    
      NQP1=NQP                                                                       
      
      IF(KFL(7)>0)THEN                                                               
            DO K=1,NDP                                                                   
                  WRITE(KW(7),462)pestName(K)                                                  
                  DO L=1,5                                                                 
                        WRITE(KW(7),143)HEDP(L),(SMMP(L,K,J),J=1,12),SMYP(L,K),HEDP(L)                                                              
                  END DO                                                                   
            END DO                                                                       
      END IF     
      
      IF(KFL(26)>0.AND.NDP>0)THEN
            DO J=1,NN
                  JJX=LY(IRO,J)                                                                    
                  IF(totPestControl(JJX,1)>0.)EXIT                                                     
            END DO
            IF(J>NN)THEN
                  XPR=0.
                  ANMX=' '
            ELSE
                  XPR=totPestControl(JJX,1)
                  ANMX=Crop_Name(JJX)
            END IF
            WRITE(KW(26),681)IYR,IY,SMY(14),SMY(16),SMY(17),SMY(18),SMY(NDVSS),SMY(77),pestName(1),ANMX,XPR,&
                           (SMYP(J,1),J=2,7),SMYP(10,1),APQC(2,1,IY)
            DO K=2,NDP
                  DO J=1,NN
                        JJX=LY(IRO,J)                                                                    
                        IF(totPestControl(JJX,K)>0.)EXIT                                                     
                  END DO
                  IF(J>NN)THEN
                        XPR=0.
                        ANMX=' '
                  ELSE
                        XPR=totPestControl(JJX,K)
                        ANMX=Crop_Name(JJX)
                  END IF
                  WRITE(KW(26),683)pestName(K),ANMX,XPR,(SMYP(J,K),J=2,7),&
                        SMYP(10,K),APQC(2,K,IY)
            END DO
      END IF                                                              
      
      II=0  
      
      IF(IY==IPY)THEN
            II=IPYI                                                                        
            IF(printChoice<=2)THEN
                  WRITE(KW(1),96)IYR,(varName(annualVarID(K)),SMY(annualVarID(K)),K=1,numAnnuVar)                          
            ELSE
                  IF(NDP/=0.AND. massPestFlag>=0)THEN
                        DO K=1,NDP                                                                  
                              IF(K==6.OR.K==1)THEN                                                           
                                    CALL Print_Page(1)                                                              
                                    WRITE(KW(1),112)                                                           
                                    WRITE(KW(1),99)IYR,IY                                                      
                              END IF                                                                         
                              WRITE(KW(1),111)pestName(K)                                                        
                              IF(massPestFlag>0)THEN
                                    WRITE(KW(1),113)HEDP(1),(SMMP(1,K,J),J=1,12),&
                                    SMYP(1,K),HEDP(1)                 
                                    DO L=1,12                                                                      
                                          PPX(L)=SMMP(2,K,L)                                                         
                                          CALL ACOUT(PPX(L),SMM(14,L),1.)                                            
                                    END DO                                                                         
                                    PPX(13)=SMYP(2,K)                                                              
                                    CALL ACOUT(PPX(13),SMY(14),1.)                                                 
                                    WRITE(KW(1),109)HEDP(2),(PPX(J),J=1,13),HEDP(2)                                
                                    DO L=1,12                                                                      
                                          PPX(L)=SMMP(3,K,L)                                                         
                                          CALL ACOUT(PPX(L),SMM(17,L),1.)                                            
                                    END DO                                                                         
                                    PPX(13)=SMYP(3,K)                                                              
                                    CALL ACOUT(PPX(13),SMY(17),1.)                                                 
                                    WRITE(KW(1),109)HEDP(3),(PPX(J),J=1,13),HEDP(3)                                
                                    DO L=1,12                                                                      
                                          PPX(L)=SMMP(4,K,L)                                                         
                                          CALL ACOUT(PPX(L),SMM(16,L),1.)                                            
                                    END DO                                                                         
                                    PPX(13)=SMYP(4,K)                                                              
                                    CALL ACOUT(PPX(13),SMY(16),1.)                                                 
                                    WRITE(KW(1),109)HEDP(4),(PPX(J),J=1,13),HEDP(4)                                
                                    DO L=5,7                                                                       
                                          WRITE(KW(1),113)HEDP(L),(SMMP(L,K,J),J=1,12),&
                                          SMYP(L,K),HEDP(L)               
                                    END DO                                                                         
                                    DO L=8,9                                                                       
                                          WRITE(KW(1),114)HEDP(L),(SMMP(L,K,J),J=1,12),&
                                          HEDP(L)                         
                                    END DO                                                                         
                                    DO L=1,12                                                                      
                                          PPX(L)=SMMP(10,K,L)                                                        
                                          CALL ACOUT(PPX(L),SMM(18,L),1.)                                            
                                    END DO                                                                         
                                    PPX(13)=SMYP(10,K)                                                             
                                    CALL ACOUT(PPX(13),SMY(18),1.)                                                 
                                    WRITE(KW(1),109)HEDP(10),(PPX(J),J=1,13),HEDP(10)                              
                                    CYCLE
                              END IF                                                                           
                              DO L=1,7                                                                       
                                    ! PRINTOUT PESTICIDE MONTHLY                                                 
                                    WRITE(KW(1),143)HEDP(L),(SMMP(L,K,J),J=1,12),SMYP(L,K),HEDP(L)             
                              END DO                                                                         
                              DO L=8,9                                                                       
                                    WRITE(KW(1),146)HEDP(L),(SMMP(L,K,J),J=1,12),HEDP(L)                         
                              END DO                                                                         
                              WRITE(KW(1),143)HEDP(10),(SMMP(10,K,J),J=1,12),SMYP(10,K),HEDP(10)
                        END DO 

                        DO K=1,NDP                                                                     
                              WRITE(KW(1),462)pestName(K)                                                    
                              WRITE(KW(1),463)(APQ(I,K,IY),I=1,5)                                        
                              WRITE(KW(1),464)(AQB(I,K,IY),I=1,5)                                        
                              WRITE(KW(1),465)(APY(I,K,IY),I=1,5)                                        
                              WRITE(KW(1),466)(AYB(I,K,IY),I=1,5)                                        
                        END DO 
                  END IF                                                                            
                  ! ************** Print Yearly Simulation Results ***************************************
                  CALL Print_Page(1)                                                                  
                  WRITE(KW(1),102)IYR,IY                                                         

                  IF(numPrintVar>0)THEN                                                                  
                        DO J=1,numPrintVar                                                                   
                              K=stateVarID(J)                                                                  
                              ! ---------------------  PRINTOUT MONTHLY --------------------------------                                                         
                              WRITE(KW(1),104)varName(K),(SMM(K,L),L=1,12),SMY(K),varName(K)                   
                        END DO                                                                       
                        WRITE(KW(1),108)'WIDX',(IWIX(L),L=1,12),'WIDX'                               
                  END IF                                                                         
                  IF(numConVar>0)THEN                                                                  
                        DO J=1,numConVar                                                                   
                              K=conVarID(J)                                                                  
                              WRITE(KW(1),109)varName(K),(SMM(K,L),L=1,12),SMY(K),varName(K)                   
                        END DO                                                                       
                  END IF                                                                         
                  IF(numMonVar>0)THEN
                        DO J=1,numMonVar                                                                     
                              K=monVarID(J)                                                                      
                              WRITE(KW(1),100)HEDS(K),(STV(K,L),L=1,12),HEDS(K)                             
                        END DO 
                  END IF                                                                            
            END IF
      END IF
      
      K=1                                                                            
      N2=NN     
      
      DO WHILE(K<=N2)
            IF(K==1.OR.K==6)THEN                                                           
                  CALL Print_Page(1)                                                                
                  WRITE(KW(1),102)IYR,IY                                                       
            END IF                                                                         
            J=LY(IRO,K)
            IYH(J)=IYH(J)+1                                                                
            N1=MAX(1,NCP(J))                                                               
            
            DO L=1,20                                                                      
                  !  ------------------- PRINTOUT CROP MONTHLY  -----------------------------                                                     
                  WRITE(KW(1),95)HEDC(L),(SMMC(L,J,K1),K1=1,12),HEDC(L)                        
                  DO K1=1,12                                                                   
                        SMMC(L,J,K1)=0.                                                          
                  END DO                                                                       
            END DO                                                                         
            
            DO K1=1,N1                                                                     
                  DO L=1,7                                                                     
                        TSFC(L,J)=TSFC(L,J)+Growing_Stress(L,K1,J)                                           
                  END DO                                                                       
            END DO                                                                         
            
            WRITE(KW(1),108)'STRS',(KDT(L,J),L=1,12),'STRS'                                
            DO J1=N1,1,-1                                                                       
                  IF(availWaterCrop(J1,J)<1.E-10)THEN                                                     
                        availWaterCrop(J1,J)=AWC(J)                                                            
                        AWC(J)=0.                                                                   
                        CRF(J1,J)=ARF                                                            
                        ARF=0.                                                                   
                        CQV(J1,J)=AQV                                                            
                        AQV=0.                                                                   
                        JP(J)=0                                                                  
                  END IF                                                                       
                  TCAW(J)=TCAW(J)+availWaterCrop(J1,J)                                                    
                  TCQV(J)=TCQV(J)+CQV(J1,J)                                                    
                  TCRF(J)=TCRF(J)+CRF(J1,J)                                                    
            END DO
            PMTE=0.    
          
            DO J1=1,N1                                                                 
                  IF(ET_GrowingSeason(J1,J)>0.)THEN                                                           
                        GS_ET=ET_GrowingSeason(J1,J)                                                             
                        Var_GS(3)=Var_GS(3)+GS_ET                                                       
                  ELSE                                                                           
                        ET_GrowingSeason(J1,J)=ACET(J)                                                          
                  END IF
                  PMTE=PMTE+totCropBio(J1)+standCropResi(J1)                                                                                                                                                     
                  XTP(1,J1,J)=YLD1(J1,J)*seedYieldPrice(J)                                                 
                  XTP(2,J1,J)=YLD2(J1,J)*forageYieldPrice(J)
                  XTP(4,J1,J)=WaterFrac_Yield
                  VALF1=XTP(1,J1,J)+VALF1+XTP(2,J1,J)                                            
                  IF(ET_GrowingSeason(J1,J)>.1)THEN                                                                
                       XTP(3,J1,J)=MIN(20.,1000.*(YLD1(J1,J)+YLD2(J1,J))/ET_GrowingSeason(J1,J))                               
                  ELSE                                                                                
                        XTP(3,J1,J)=0.                                                                    
                  END IF                                                                              
                  XX=NPSF(J1,J)                                                                  
                  TPSF(J1,J)=TPSF(J1,J)/(XX+1.E-10)                                              
                  YGIS=MAX(YGIS,YLD1(J1,J))                                                           
                  Var_GS(1)=Var_GS(1)+YGIS                                                                
                  BGIS=MAX(BGIS,DMF(J1,J))                                                            
                  WGIS=MAX(WGIS,XTP(3,J1,J))                                                          
                  FGIS=MAX(FGIS,TotN_Frac(J1,J))                                                           
                  Var_GS(5)=Var_GS(5)+FGIS                                                                
                  N_YLD=MAX(N_YLD,YLNF(J1,J))                                                           
                  P_YLD=MAX(P_YLD,YLPF(J1,J))                                                           
                  K_YLD=MAX(K_YLD,YLKF(J1,J))                                                           
                  HGIS=MAX(HGIS,HIF(J1,J))                                                            
                  Var_GS(7)=Var_GS(7)+HGIS                                                                
                  WSGS=MAX(WSGS,Growing_Stress(1,J1,J))                                                           
                  SNGS=MAX(SNGS,Growing_Stress(2,J1,J))                                                           
                  SPGS=MAX(SPGS,Growing_Stress(3,J1,J))                                                           
                  STGS=MAX(STGS,Growing_Stress(5,J1,J))                                                           
                  SAGS=MAX(SAGS,Growing_Stress(6,J1,J))                                                           
                  SSGS=MAX(SSGS,Growing_Stress(7,J1,J))    
              
                  IF(IY==IPY)THEN                                                                
                        IF(CSTF(J1,J)<=0.)THEN                                                       
                              CSTF(J1,J)=COST                                                          
                              CSOF(J1,J)=COST-CSFX                                                     
                              COST=0.                                                                  
                              CSFX=0.                                                                  
                        END IF                                                                       
                        ! ---------------------------  PRINTOUT CROP ANNUAL -----------------------------------
                        IF(cropCode(J)==plantCategoryCode(7).OR.cropCode(J)==plantCategoryCode(8).OR.&
                              cropCode(J)==plantCategoryCode(10))THEN
                              X2=.0001*PPL0(J)
                              X1=YLCF(J1,J)
                              XTP(3,J1,J)=0.                                                                    
                        ELSE
                              X2=PPL0(J)
                              X1=1000.*YLCF(J1,J)
                        END IF
                        WRITE(KW(1),106)Crop_Name(J),YLD1(J1,J),YLD2(J1,J),totCropBioX(J),&
                              YLNF(J1,J),YLPF(J1,J),YLKF(J1,J),X1,TotN_Frac(J1,J),TotP_Frac(J1,J),&
                              TotK_Frac(J1,J),VIR(J1,J),VIL(J1,J),availWaterCrop(J1,J),ET_GrowingSeason(J1,J), &
                              XTP(3,J1,J),X2,TPSF(J1,J),CSTF(J1,J),CSOF(J1,J),XTP(1,J1,J),&
                              XTP(2,J1,J),EK,REK,Wind_Erod
                        WRITE(KW(1),97)(Growing_Stress(L,J1,J),L=1,7)                                            
                        TSMQ=0.                                                                      
                        TSMY=0.                                                                      
                  END IF      
              
                  TDM(J)=TDM(J)+DMF(J1,J)                                                        
                  TYL1(J)=TYL1(J)+YLD1(J1,J)                                                     
                  TYL2(J)=TYL2(J)+YLD2(J1,J)                                                     
                  TYLN(J)=TYLN(J)+YLNF(J1,J)                                                     
                  TYLP(J)=TYLP(J)+YLPF(J1,J)                                                     
                  TYLK(J)=TYLK(J)+YLKF(J1,J)
                  TYLC(J)=TYLC(J)+YLCF(J1,J)            
                  IF(YLNF(J1,J)>0.)THEN                                                          
                        NYLN(J)=NYLN(J)+1                                                          
                        XX=NYLN(J)                                                                 
                        X1=XX+1.                                                                   
                        X2=MAX(1.,BLYN(J),1.1*Growing_Stress(2,J1,J))                                          
                        UNA(J)=MIN(1000.,UNA(J)*(modelPara(28)/(1.-Growing_Stress(2,J1,J)/X2))**2)                  
                        X3=UNA(J)                                                                       
                        ULYN(J)=(ULYN(J)*XX+YLNF(J1,J))/X1                                              
                        UNA(J)=modelPara(46)*UNA(J)+(1.-modelPara(46))*ULYN(J)                                    
                  END IF                                                                              
                  TRD(J)=TRD(J)+RDF(J)                                                           
                  THU(J)=THU(J)+HUF(J)                                                           
                  TETG(J)=TETG(J)+ET_GrowingSeason(J1,J)                                                      
                  TVIR(J)=TVIR(J)+VIR(J1,J)                                                      
                  TIRL(J)=TIRL(J)+VIL(J1,J)                                                      
                  CST1=CST1+CSTF(J1,J)                                                           
                  CSO1=CSO1+CSOF(J1,J)                                                           
                  TFTN(J)=TFTN(J)+TotN_Frac(J1,J)                                                     
                  TFTP(J)=TFTP(J)+TotP_Frac(J1,J)                                                     
                  TFTK(J)=TFTK(J)+TotK_Frac(J1,J)                                                     
                  TCST(J)=TCST(J)+CSTF(J1,J)                                                     
                  TCSO(J)=TCSO(J)+CSOF(J1,J)                                                     
                  TVAL(J)=TVAL(J)+XTP(1,J1,J)+XTP(2,J1,J)                                        
                  PSTM(J)=PSTM(J)+TPSF(J1,J)
            END DO   
          
            IF(N1>1)N2=N2-1                                                                
            K=K+1
          
      END DO                                                                          
      VIRT=0.                                                                        
      DARF=DARF+SMY(4)*SMY(4)                                                        
      IF(SMY(4)>BARF)BARF=SMY(4)                                                                  
      IF(SMY(4)<SARF)SARF=SMY(4)                                                   
      !------------------------------------ PRINTOUT ANNUAL FILE ---------------------------------------                                                          
      X1=.001*totSOC                                                                    
      IF(DOY_realtime==0)GO TO 556                                                           
      IF(IY/=IGSD)GO TO 555                                                               
      K=IRTC                                                                         
      L1=LY(IRO,K)                                                                   
      J=MAX(1,NHV(L1))                                                               
      WRITE(KW(2),498)IYR,IRLX,SMY(4),SMY(10),SMY(11),SMY(14),SMY(16),SMY(17),SMY(29),SMY(33),SMY(42),&
            SMY(48),SMY(47),SMY(50),SMY(51),SMY(52),SMY(49),SMY(43),SMY(44),SMY(45),SMY(46),SMY(56),SMY(54),&              
            SMY(55),SMY(57),limeRate,plawDepSOC,X1,APBC,totLabileP,totNO3_N                                      
      
      IF(KFL(19)==0)GO TO 555      
      
      !------------------------- Print out Crop Yield --------------------------------------------
      
      WRITE(KW(19),558)IYR,IRLX,Crop_Name(L1),YLD1(J,L1),YLD2(J,L1),XTP(4,J,L1),HIF(J,L1),DMF(J,L1),&
            RWF(J,L1),YLNF(J,L1),YLPF(J,L1),YLCF(J,L1),TotN_Frac(J,L1),TotP_Frac(J,L1),TotK_Frac(J,L1),&
            VIR(J,L1),VIL(J,L1),XTP(3,J,L1),ET_GrowingSeason(J,L1),availWaterCrop(J,L1),CRF(J,L1),CQV(J,L1),&
            CSTF(J,L1),CSOF(J,L1),XTP(1,J,L1),XTP(2,J,L1),TPSF(J,L1),(Growing_Stress(L,J,L1),L=1,7),PPL0(L1),IPLD(J,L1),&
            IGMD(J,L1),IHVD(J,L1)                                                                     
      ! *************** annual yield file for testing ***************************
      ! added by TXH

      ! WRITE(KW(50), "(1X,2I4,1X,A4,3F8.3)")IYR, IRLX, Crop_Name(L1), YLD1(J,L1), HIF(J,L1), DMF(J,L1)
      ! *************************************************************************
      IF(J==1)GO TO 555                                                              
      J1=J-1                                                                         
      IPLD(J1,L1)=IPLD(J,L1)                                                         
      IGMD(J1,L1)=IGMD(J,L1)                                                         
      GO TO 555                                                                      
 556 IF(KFL(2)>0)WRITE(KW(2),498)IYR,IY,SMY(4),SMY(10),SMY(11),SMY(14),SMY(16),SMY(17),SMY(29),SMY(33),&
                            SMY(42),SMY(48),SMY(47),SMY(50),SMY(51),SMY(52),SMY(49),SMY(43),SMY(44),SMY(45),&
                            SMY(46),SMY(56),SMY(54),SMY(55),SMY(57),limeRate,plawDepSOC,X1,APBC,totLabileP,&
                            totNO3_N                         
      IF(KFL(19)==0)GO TO 555                                                        
      K=1                                                                                 
      N2=NN                                                                               
  686 J=LY(IRO,K)                                                                    
      N1=MAX(1,NCP(J))                                                               
      
      DO J1=1,N1                                                                     
            WRITE(KW(19),558)IYR,IY,Crop_Name(J),YLD1(J1,J),YLD2(J1,J),XTP(4,J1,J),HIF(J1,J),DMF(J1,J),&
               RWF(J1,J),YLNF(J1,J),YLPF(J1,J),YLCF(J1,J),TotN_Frac(J1,J),TotP_Frac(J1,J),TotK_Frac(J1,J),&
               VIR(J1,J),VIL(J1,J),XTP(3,J1,J),ET_GrowingSeason(J1,J),availWaterCrop(J1,J),CRF(J1,J),&
               CQV(J1,J),CSTF(J1,J),CSOF(J1,J),XTP(1,J1,J),XTP(2,J1,J),TPSF(J1,J),(Growing_Stress(L,J1,J),L=1,7),&
               PPL0(J),IPLD(J1,J),IGMD(J1,J),IHVD(J1,J)
            ! *************** annual yield file for testing ***************************
            ! added by TXH

            ! WRITE(KW(50), "(1X,2I4,1X,A4,3F8.3)")IYR, IY, Crop_Name(J), YLD1(J1,J), HIF(J1,J), DMF(J1,J)
            ! *************************************************************************
            IF(J1==1)CYCLE                                                             
            L1=J1-1                                                                    
            IPLD(L1,J)=IPLD(J1,J)                                                      
            IGMD(L1,J)=IGMD(J1,J)                                                      
      END DO 
      
      IF(N1>1)N2=N2-1                                                                
      K=K+1                                                                               
      IF(K<=N2)GO TO 686                                                                  
  555 IF(KFL(23)>0)THEN                                                              
            DO K=1,NN                                                                    
                  L1=LY(IRO,K)                                                             
                  DO J=1,NCP(L1)                                                           
                        WRITE(KW(23),668)IYR,IY,M21,K21,Crop_Name(L1),totCropBioX(L1),(RWTX&              
                        (Layer_ID(I),L1),I=1,Actual_SoilLayers),RWX(L1)                                        
                  END DO                                                                   
            END DO                                                                       
      END IF  

      IF(KFL(22)>0)THEN                                                              
            RNO3=SMY(4)*rainNconX                                                       
            WRITE(KW(22),566)IYR,SMY(4),SMY(10),SMY(11),SMY(13),SMY(14),&
                  (SMY(J),J=16,20),RNO3,(SMY(J),J=43,46),SMY(49),SMY(89),SMY(52),&
                  SMY(85),SMY(50),(SMY(J),J=59,61)
      END IF                             
      APBC=.5*(SMY(1)+SMY(2))                                                        
      IF(KFL(8)>0)WRITE(KW(8),293)IRO0,IYR,APBC,PMTE,(SMY(annualVarID(K)),K=1,numAnnuVar)                                                                          
      RWTX=0.
      DO K=1,LC                                                                       
            J=LY(IRO,K)
            IF(J==0)CYCLE
            IF(KFL(24)>0)THEN                                                              
                  IF(cropCode(J)==plantCategoryCode(7).OR.cropCode(J)==plantCategoryCode(8))WRITE(KW(24),669)IYR,IY,&                   
                  Crop_Name(J),YLD2(K,J),totCropBio(J),totRootWeight(J),Current_LAI(J),standCropResi(J)                                  
            END IF                                                                         
            NCP(J)=1
            IF(KG(J)==0)NCP(J)=0                                                           
            NHV(J)=0                                                                       
            RWX(J)=0.                                                                           
            DO I=N1,1,-1                                                                        
                  NPSF(I,J)=0                                                                
                  TotN_Frac(I,J)=0.                                                               
                  TotP_Frac(I,J)=0.                                                               
                  TotK_Frac(I,J)=0.                                                               
                  YLD1(I,J)=0.                                                               
                  YLD2(I,J)=0.                                                               
                  YLNF(I,J)=0.                                                               
                  YLPF(I,J)=0.                                                               
                  YLKF(I,J)=0.
                  YLCF(I,J)=0.                                                               
                  DMF(I,J)=0.
                  RWF(I,J)=0.                                                                     
                  VIL(I,J)=0.                                                                
                  VIR(I,J)=0.                                                                
                  availWaterCrop(I,J)=0.                                                                
                  CQV(I,J)=0.                                                                
                  CRF(I,J)=0.                                                                
                  TPSF(I,J)=0.                                                                    
                  CSOF(I,J)=0.                                                               
                  CSTF(I,J)=0.                                                               
                  ET_GrowingSeason(I,J)=0.                                                                
            END DO                                                                         
            HUF(J)=0.                                                                           
            RDF(J)=0.                                                                      
            BLYN(J)=0.                                                                     
            IF(cropCode(J)==plantCategoryCode(3).OR.cropCode(J)==plantCategoryCode(6).OR.cropCode(J)==&
                  plantCategoryCode(7).OR.cropCode(J)==plantCategoryCode(8).OR.cropCode(J)==plantCategoryCode(10))THEN                                            
                  ANA(J)=0.                                                                  
                  sumActuralEP(J)=0.                                                                       
                  sumPotentialEP(J)=0.
                  ACET(J)=0.
                  AWC(J)=0.                                                                   
                  ARF=0.                                                                   
                  AQV=0.                                                                                                                                          
            END IF
            IF(YLAT>0.)THEN                                                                              
                  IF(cropCode(J)/=plantCategoryCode(2).AND.cropCode(J)/=plantCategoryCode(5))totCropBioX(J)=0.
            ELSE
                  IF(cropCode(J)/=plantCategoryCode(1).AND.cropCode(J)/=plantCategoryCode(4))totCropBioX(J)=0.
            END IF
      END DO

      IF(KFL(NGF)==0)GO TO 499                                                       
            X2=SMY(49)+SMY(52)                 ! N LOST TO ATMOSPHERE (kg/ha)                                                           
            X3=SMY(43)+SMY(44)+SMY(45)+SMY(46) ! N LOST WITH RUNOFF, SUBSURFACE FLOW, AND EROSION                                                      
            GS_SoilWater=GS_SoilWater/REAL(NumOfDays)                                                                                                                                                                                                                   
 
            WRITE(KW(IGIS),685)XLOG,YLAT,YGIS,WGIS,GS_ET,GS_PlantTranspiration,SMY(19),FGIS,SMY(11),GS_VPD,HGIS,GS_Rad,&
                N_YLD,N_Soil,X2,X3,SMY(42),SMY(NDVSS),SMY(3),BGIS,(Growing_Stress(L,1,1),L=1,7),SMY(46),SMY(43),&
                SMY(47),SMY(49),SMY(50),SMY(85),SMY(51),SMY(52),SMY(53),SMY(54),SMY(56),SMY(57),&
                SMY(58),SMY(59),SMY(60),SMY(61),SMY(62),SMY(63),P_YLD,K_YLD,SMY(77),SMY(4),SMY(17),&             
                SMY(14),SMY(5),SMY(6),SMY(10),SMY(12),SMY(13),SMY(16),SMY(18),SMY(15),SMY(20),&
                GS_SoilWater,SMY(68)                                                     
      
            Var_GS(8)=Var_GS(8)+GS_SoilWater  
        499 WB1=SMY(34)                                                                    
            WB2=SMY(42)
            Growing_Stress=0.                                                                               
            SMM=0.                                                                         
            IF(NDP==0)GO TO 83                                                             
            SMYP=0.                                                                        
            SMMP=0.                                                                        
            totPestControl=0.     
         83 MO1=MO                                                                         
      END DO dayLoop  ! ----------- Day Loop --------------------------
      ! ********************************************************************************************
      ! End Daily Simulation     go to next day
      ! ********************************************************************************************
      IBD=1                                                                          
      MO=1                                                                           
      KT=1                                                                           
      KC=0                                                                           
      IYR1=IYR                                                                            
      IYR=IYR+1                                                                      
      IYX=IYX+1                                                                      
      leapYr=1                                                                          
      IPY=IPY+II                                                                     
      CALL Sum_WaterNutrient(1)
      
      IF(correctionForRuns0 >0)THEN                                                                     
            DO J=1,21                                                                         
                  IX(J)=IX0(J)                                                             
                  randID(J)=randID0(J)                                                                 
            END DO
      V3=Generate_Random(randID(3))                                                                                                                                               
      END IF
      
      IF(FC_Method<=3.OR.FC_Method==7.OR.FC_Method==9)THEN
            XX=0.                                                                          
            K=1                                                                            
            TPAW=0.                                                                        
            DO I=1,Actual_SoilLayers                                                                
                  J=Layer_ID(I)                                                                       
                   Y1=.1*SOC(J)/WT(J)                                                                  
                  XZ=.0172*Y1                                                                         
                  ZZ=1.-XZ
                  bulkDensity(J)=1./(XZ/.224+ZZ/mineralBulkDensity(J))
                  BDX=modelPara(2)+.35+.005*sandFrac(J)
                  bulkDensity(J)=MIN(bulkDensity(J),BDX)
                  ZR=1.E-4*WT(J)/bulkDensity(J)
                  IF(I/=Actual_SoilLayers)Z(J)=XX+ZR                                                                                                                                            
                  X1=1000.*(Z(J)-XX)
                  SELECT CASE(FC_Method+1)
                        CASE(1,3)
                              CALL Cal_SoilWaterW(clayFrac(J),sandFrac(J),Y1,X2,X3)                                            
                        CASE(2,4)
                              ZZ=.5*(XX+Z(J))                                                              
                              CALL Cal_SoilWaterO(CEM(J),clayFrac(J),Y1,sandFrac(J),X2,X3,ZZ)
                        CASE(8)
                              CALL Cal_FC_WP(clayFrac(J),sandFrac(J),Y1,X2,X3) 
                        CASE(10)         
                              CALL Cal_SoilWaterBNW(clayFrac(J),siltFrac(J),sandFrac(J),Y1,bulkDensity(J),X2,X3)
                  END SELECT                                                                     
                  XY=1.-rockFrac(J)*.01                                                               
                  XZ=XY*X1                                                                       
                  Wilt_Point(J)=X2*XZ                                                                   
                  fieldCapacity(J)=X3*XZ                                                                    
                  CALL WiltPoint_Sat(J)
                  IF(K<=3)THEN
                        TPAW=TPAW+fieldCapacity(J)-Wilt_Point(J)
                        DO                                                          
                              IF(TPAW<WCS(K))EXIT
                              ZCS(K)=XX+(Z(J)-XX)*((WCS(K)-PZW)/(TPAW-PZW))    ! where would the PZW come?                              
                              K=K+1
                              IF(K>3)EXIT
                        END DO
                  END IF    
                  PZW=TPAW                                                                       
                  XX=Z(J)                                                                  
            END DO
            layersEqualThick= INT(Z(Layer_ID(Actual_SoilLayers))/layerThickness+.999) !ZC(layersEqualThick)=layersEqualThick*LayerThickness
      END IF                                                                                            
      IF(erosionMode>0)THEN                                                                 
            DO J=1,Actual_SoilLayers                                                                  
                  L=Layer_ID(J)                                                                 
                  CSlowHumus(L)=SOL(1,L)                                                              
                  CPassiveHumus(L)=SOL(2,L)                                                         
                  CStructLitt(L)=SOL(3,L)                                                         
                  CMetabLitt(L)=SOL(4,L)                                                         
                  CBiomass(L)=SOL(5,L)                                                         
                  SOC(L)=SOL(6,L)                                                          
                  NSlowHumus(L)=SOL(7,L)                                                         
                  NPassiveHumus(L)=SOL(8,L)                                                         
                  NStructLitt(L)=SOL(9,L)                                                         
                  NMetabLitt(L)=SOL(10,L)                                                        
                  NBiomass(L)=SOL(11,L)                                                        
                  SOC(L)=SOL(12,L)                                                              
                  activeMineralP(L)=SOL(13,L)                                                              
                  SOP(L)=SOL(14,L)                                                          
                  stableMineralP(L)=SOL(15,L)                                                          
                  exchangeK(L)=SOL(16,L)                                                        
                  fixK(L)=SOL(17,L)                                                        
                  ! soilWater(L)=SOL(18,L)                                                          
                  structLitt(L)=SOL(19,L)                                                              
                  metabLitt(L)=SOL(20,L)                                                              
                  lgStructLitt(L)=SOL(21,L)                                                             
                  CLgStructLitt(L)=SOL(22,L)                                                       
                  NLgStructLitt(L)=SOL(23,L)                                                      
            END DO                                                                       
            ! CALL SCONT(1)                                                                   
      END IF                                                                         
      
      CALL Prep_SoilVar(YTP)                                                                
	  
      IF(KFL(14)>0)CALL Print_OrgNCvar(IYR1,12,DayOfMon)

      IF(KFL(31)>0)CALL Print_SoilVar2(YTP,31)
	  
      SMY=0.                                                                         
      IF(printChoice==2 .OR. printChoice==4)THEN                                                            
            WRITE(KW(1),101)                                                           
            CALL Print_SoilVar2(YTP,1)                                                          
      END IF                                                                         
      IF(IGSD/=0.AND.IY==IGSD)THEN                                                   
            REWIND KR(7)                                                               
            CALL ReadDailyWeather0                                                             
            weatherVarID=weatherVarID0                                                                   
            IYR=Year0                                                                   
            IGSD=IGSD+cropRotationYrs                                                              
      END IF                                                                         
      CALL Leap_Yr_Check(IYR,leapYr,leapYrFlag)                                                       
      IGIS=IGIS+1    
      !CLOSE(KW(51)) ! for testing                                                                    
END DO  yearLoop               ! year loop   

! *************************************************************************
!  End Yearly Simultion ..... go to next year
! *************************************************************************
     ! IY=numSimuYear+1                                                                      
   88 RETURN      
   63 FORMAT(1X,3I4,2X,A8,2X,'LIME RATE=',F8.2,' t/ha',1X,'CROP=',A4,&
            1X,'COTL=',F7.2,' $',1X,'COOP=',1X,F7.2,' $')          
   89 FORMAT(1X,3I4,2X,'GERMINATION--0.2 m totSoilWater = ',F7.1,'mm',2X,'HU = ',&
            F4.0,'c',2X,'fracHeatUnit = ',F6.2,2X,A4)       
      !90 FORMAT(1X,3I4,2X,'LIME',2X,'RATE=',F6.0,'t/ha',1X,'fracHeatUnit=',F6.2,2X,&
      !   'COST=',F7.0,'$/ha')      
   91 FORMAT(1X,3I4,F12.1,400F12.3)         
   92 FORMAT(1X,3I4,2X,'RB FR DK',20X,'furrowHeight=',F5.0,'mm',2X,'furrowInterval=',F6.2&            
            ,'m',2X,'fracHeatUnit=',F6.2)         
   95 FORMAT(A16,12F10.2,5X,A16)         
   96 FORMAT(//I5,9(2X,A4,F8.2)/(5X,9(2X,A4,F8.2)))        
   97 FORMAT(T10,'STRESS DAYS (BIOM)     WATER=',F5.1,2X,'N=',F5.1,2X,&              
            'P=',F5.1,2X,'K=',F5.1,2X,'TEMP=',F5.1,2X,'IrrWater=',F5.1,2X,'SALT=',F5.1)                      
   99 FORMAT(45X,'YR=',I4,2X,'YR#=',I4/T22,'JAN',9X,'FEB',9X,'MAR',9X,&              
            'APR',9X,'MAY',9X,'JUN',9X,'JUL',9X,'AUG',9X,'SEP',9X,'OCT',9X,&               
            'NOV',9X,'DEC',9X,' YR')         
   100 FORMAT(A16,12F9.2,5X,A16)     
   101 FORMAT(T5,'SOIL DATA')            
   102 FORMAT(60X,'YR=',I4,2X,'YR#=',I4/T24,'JAN',7X,'FEB',7X,'MAR',7X,&               
            'APR',7X,'MAY',7X,'JUN',7X,'JUL',7X,'AUG',7X,'SEP',7X,'OCT',7X,&              
            'NOV',7X,'DEC',10X,' YR')                
   104 FORMAT(1X,A15,13F9.2,2X,A15)          
   105 FORMAT(1X,3I4,1X,A4,5(2X,A4,F8.3)/(10X,5(2X,A4,F8.3)))          
   106 FORMAT(2X,A4,1X,'YLD=',F6.1,'/',F6.1,2X,'BIOM=',F6.1,'t/ha',2X,&               
            'Yield_N=',F5.0,2X,'YLP=',F5.0,2X,'YLK=',F5.0,2X,'YLC=',F6.0,2X,'FN=',&
            F5.0,2X,'FP=',F5.0,2X,'FK=',F5.0,'kg/ha'/T7,'IRGA=',F5.0,2X,'IRDL=',&
            F5.0,2X,'availWaterCrop=',F6.0,'mm',2X,'GSET=',F6.0,'mm',2X,'WUEF=',F10.2,&
            'kg/mm',2X,'POP=',F9.4,'p/m2',2X,'PSTF=',F6.2/T7,'COST=',F7.0,2X,&
            'COOP=',F7.0,2X,'RTRN=',F7.0,'/',F7.0,'$/ha',2X,'EK=',F6.3,2X,&
            'REK=',F6.3,2X,'Wind_Erod=',F6.3)
   107 FORMAT(1X,I4,2X,I2,2X,I2,100F10.3)         
   108 FORMAT(1X,A15,12I10,5X,A15)           
   109 FORMAT(1X,A15,13F9.4,2X,A15)           
   111 FORMAT(34X,'-------------------------',A8,'-------------------------')          
   112 FORMAT(47X,'PESTICIDE SIMULATION (G/HA)')           
   113 FORMAT(1X,A4,13F9.0,2X,A4)             
   114 FORMAT(1X,A4,12F9.0,11X,A4)            
   !119 FORMAT(10I8)                             
   143 FORMAT(1X,A4,13F12.5,2X,A4)              
   146 FORMAT(1X,A4,12F12.5,14X,A4)                 
   293 FORMAT(1X,2I4,42F8.2)          
   465 FORMAT(8X,'ADSRB',5E13.5)           
   462 FORMAT(5X,A16)                        
   464 FORMAT(8X,'Q+SubSurf_LateralFlow',5E13.5)         
   463 FORMAT(8X,'SOL  ',5E13.5)           
   466 FORMAT(8X,'SED Y',5E13.5)             
   468 FORMAT(1X,3I4,1X,I4,1X,A16,14F10.4)         
   475 FORMAT(1X,I4,I2,I4,2X,40F8.1)        
   498 FORMAT(1X,2I4,6F8.1,F8.3,17F8.1,F8.2,F8.1,F8.2,2F8.1,10(A16,F8.0))          
   !505 FORMAT(A80)                            
   513 FORMAT(1X,3I4,1X,I4,4X,A4,100F10.3)                            
   523 FORMAT(1X,I4,I2,I4,2X,A4,13F8.1)                                                   
   532 FORMAT(1X,3I4,1X,5(2X,A4,F10.3)/(10X,5(2X,A4,F10.3)))                         
   533 FORMAT(5X,12I4)                                                                
   558 FORMAT(1X,2I4,1X,A4,7F8.3,6F8.1,2F8.3,4F8.1,4F8.2,9F8.3,3I9,&                  
            10(1X,A16,F8.0))                                                               
 ! 562 FORMAT(1X,I4,I2,2X,60F8.1)                                                     
   565 FORMAT(1X,2I4,28F8.2,4X,A4,F8.2,F8.0)                                          
   566 FORMAT(1X,I4,4X,10F8.2,24X,15F8.2,4X,A4,F8.2,F8.0)                                          
   567 FORMAT(1X,3I4,2X,'LIME',14X,I4,6X,'   9',8X,F10.2,10X,2F10.2)                                                
   588 FORMAT(5X,'!!!!!',3I4,' Q= ',F6.1)                                          
 ! 589 FORMAT(5X,'^^^^^',3I4,' plawDepSoilWater = ',F6.1,1X,'plawDepFieldCapa = ',F6.1)                    
   666 FORMAT(1X,3I4,2X,A8,8X,I6,6X,3I4,2F10.2,20X,F10.2)                             
   668 FORMAT(1X,4I4,1X,A4,20F8.2)                                                
   669 FORMAT(1X,2I4,1X,A4,5F8.2)                                                     
   679 FORMAT(1X,I4,1X,I4,1X,13(1X,F10.3))                                            
   680 FORMAT(10X,'WATERSHED AREA = ',F10.2,' HA', /18X,'Q',10X,'Y',9X,'sedimentLossN', &            
             9X,'sedimentLossP',9X,'QN',9X,'peakRunoffRate', /2X,'YR    MO',7X,'(mm)',6X,'(t/ha)',3X,&              
            '|-----------------(kg/ha)---------------|')                                   
   681 FORMAT(1X,I4,1X,I4,6F8.2,1X,A16,1X,A4,8F8.2,E12.5)                               
   682 FORMAT(1X,3I4,400E12.5)
   683 FORMAT(T60,A16,1X,A4,8F8.2,E12.5)                                                                             
   685 FORMAT(1X,4F10.2,5F10.1,2F10.3,F10.0,50F10.2)
            
END SUBROUTINE Simulation_Daily                                                                           

! ==================================================================
REAL FUNCTION PMAV(FMO1, FMO2, X1, X2)
      IMPLICIT NONE
      REAL, INTENT(IN)::  X1, X2
      INTEGER, INTENT(IN):: FMO1, FMO2
      PMAV = (FMO1*X1+FMO2*X2)/30.
END FUNCTION PMAV