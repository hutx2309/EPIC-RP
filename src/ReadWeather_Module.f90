MODULE ReadWeather_Module
      USE PARM
      USE MisceFun_Module
      USE Weather_Generator 
      USE Hydrology, ONLY: Sat_Vapor
      IMPLICIT NONE
      CONTAINS

! ------------------------ 1. Read montly weather data in the Main program -------------------------------
SUBROUTINE MonthlyWeatherInfo(FWPM5, FWPM1, FWIND, WP1_ID, windID, siteAzimuth, peakRateFactor, WPM1FILE, WINDFILE,& 
                              monAvePPT, meanDayLength, maxRainfall, vernalizeDay, yearTotPPT, windEroFactor, YLAZ)
      ! added by TXH: this is to read weather data measured at stations, including montly average data and daily data


      ! Argument lists
      CHARACTER(LEN = 80), INTENT(IN):: FWPM5, FWPM1, FWIND
      INTEGER, INTENT(IN):: WP1_ID, windID
      REAL, INTENT(IN):: siteAzimuth, peakRateFactor, meanDayLength
      REAL, DIMENSION(6,12), INTENT(OUT)::monAvePPT
      REAL, DIMENSION(12):: XYP, XTP  ! the XYP and XTP in main program have a size of 200, but here only 6 or 12 is needed
      REAL, INTENT(INOUT):: maxRainfall, vernalizeDay, yearTotPPT, windEroFactor, YLAZ
      CHARACTER(LEN = 80), INTENT(OUT):: WPM1FILE, WINDFILE
      
      ! Local variables declaration
      CHARACTER(LEN = 80), DIMENSION(2):: TITW5
      INTEGER:: I, I1, II, IW, XM, J
      REAL:: W0, W1, Y, X, ELEX, RY, XX, D, E, F, upperT, lowerT, & 
             ! upper and lower bounds of montly average temperature of weather station 
             SUM, REXP, ADD, X1, X2, X3, XL, XY2, RAMM,R1, R6, RTO, V4, extremeRain, &
             HGXW, HGX0, &                           ! for Potential ET model
             AWV, windEnergy, RN2, WV, EV, CLF       ! for reading wind data
      REAL, DIMENSION(12):: monWindSpeed             ! monthly average wind speed  
      ! //////////////////////////////////////////////////////////////////////////////////////////////////
      ! ========================== 18. Weather station choice ====================================== 
      ! WPM5US.DAT    is monthly weather data of various stations
      ! WPM1USEL.DAT  is the list of stations.
      CALL OPENV(KR(18),FWPM5,IDIR(1),KW(MSO) ) ! Alternate weather station catalog (used with fWIDX) 
      IF(WP5_ID>0)THEN                          ! WP5_ID: Weather station #  KR(18) WPM5US.DAT
            I1=WP5_ID-1    
            II=67*I1     
            IF(II>0)THEN
                  DO I=1,II
                        READ(KR(18),'()',IOSTAT=NFL)
                        IF(NFL/=0)EXIT  
                  END DO
                   IF(runMode==0)THEN  
                        WRITE(*,*)'WPM5 NO = ',WP5_ID,' NOT IN WPM5 FILE'
                        ERROR STOP  
                  ELSE
                        WRITE(KW(MSO),'(A,A8,A,I4,A)')' !!!!! ',siteName,&
                              'WPM5 NO = ',WP5_ID,' NOT IN WPM5 FILE'
                        ERROR STOP
                  END IF
            END IF 
            READ(KR(18),'(A80)')TITW5      ! title of WPM5US.DAT   !     LINE 1/2    
      END IF     

      ! WP1_ID: WEATHER STA # FROM KR(17) WPM1USEL.DAT   
      ! search for the nearest weather station from WPM1USEL.DAT if WPM1FILE is not provided
      CALL OPENV(KR(17),FWPM1,IDIR(1),KW(MSO) ) ! Catalog of weather stations with monthly weather data 
      IF(WP1_ID==0)THEN   
            W0=1.E20   
            DO     
                  READ(KR(17),*,IOSTAT=NFL)II,OPSCFILE,Y,X,ELEX   
                  IF(NFL/=0)EXIT   
                  RY=Y/CLT
                  XX=SIN1*SIN(RY)+COS1*COS(RY)*COS((X-XLOG)/CLT)  ! solar zenith angle 
                  D=6378.8*ACOS(XX)                               ! cal the distance  rad(\theta) = r/l for a sphere
                  E=ABS(ELEV-ELEX)  
                  ! modelPara(79) = 1 weighting factor for locating appropriate weather stations (==1 gives strictly distance)
                  W1=modelPara(79)*D+(1.-modelPara(79))*E   
                  IF(W1>=W0)CYCLE  
                  W0=W1  
                  WPM1FILE=OPSCFILE  
            END DO   
      ELSE
            II=-1  
            DO WHILE(II/=WP1_ID)  
                  READ(KR(17),*,IOSTAT=NFL)II,WPM1FILE             ! WPM1FILE: weather list file  
                  IF(NFL/=0)THEN
                        IF(runMode==0)THEN                                                                    
                              WRITE(*,*)'WPM1 NO = ',WP1_ID,' NOT IN MO WEATHER LIST FILE'
                              ERROR STOP                                                      
                        ELSE
                              WRITE(KW(MSO),'(A,A8,A,I4,A)')' !!!!! ',siteName,&
	                              ' WPM1 NO = ',WP1_ID,' NOT IN MO WEATHER LIST FILE'
                              ERROR STOP
                        END IF 
                  END IF                                                                            
            END DO
      END IF
      
      REWIND KR(17)
      ! Finally we find the a WP1 file that is the nearest one
      ! ====================== 20. read data from the nearest stations NCCLAYTO.WP1 ==========================    
      
      CALL OPENV(KR(24), WPM1FILE, IDIR(1), KW(MSO))   
      !     LINE 1/2                                                                       
      READ(KR(24),'()')
      READ(KR(24),'()')
      !----------------continue to read NCCLAYTO.WP1-------------------  
      readData: DO IW=1,6          ! Why it is 6?   what is IW?        
      !     LINES 3/15 
            readMonthlyData: IF(IW==1)THEN
                  READ(KR(24),1)(monAveTmax(IW,I),I=1,12)  
                  READ(KR(24),1)(monAveTmin(IW,I),I=1,12)  
                  READ(KR(24),1)(monTmaxStd(IW,I),I=1,12) 
                  READ(KR(24),1)(monTminStd(IW,I),I=1,12) 
                  READ(KR(24),1)(monAvePPT(IW,I),I=1,12) 
                  READ(KR(24),1)(dailyRainInfo(2,IW,I),I=1,12) 
                  READ(KR(24),1)(dailyRainInfo(3,IW,I),I=1,12) 
                  READ(KR(24),1)(probWetDay(1,IW,I),I=1,12) 
                  READ(KR(24),1)(probWetDay(2,IW,I),I=1,12) 
                  READ(KR(24),1)(UAVM(I),I=1,12)                ! Average days of precipitation in each month  
                  READ(KR(24),1)(halfHrRainInfo(IW,I),I=1,12)  
                  READ(KR(24),1)(monAveSoilRad(IW,I),I=1,12)   
                  READ(KR(24),1)(RH(IW,I),I=1,12)               ! This is actual temperature for calculating for cal RH   
                  READ(KR(24),1)(monWindSpeed(I),I=1,12)        ! UAV0: Here is wind speed                                  
                  REWIND KR(24)     
            ELSE         !IW = 2..6   read from WPM5US data    
                  READ(KR(18),1)(monAveTmax(IW,I),I=1,12)    
                  READ(KR(18),1)(monAveTmin(IW,I),I=1,12)   
                  READ(KR(18),1)(monTmaxStd(IW,I),I=1,12)  
                  READ(KR(18),1)(monTminStd(IW,I),I=1,12) 
                  READ(KR(18),1)(monAvePPT(IW,I),I=1,12)  
                  READ(KR(18),1)(dailyRainInfo(2,IW,I),I=1,12) ! dailyRainInfo(1,), (2,), (3,) are the montly storm mean, std, and skew coefficient  
                  READ(KR(18),1)(dailyRainInfo(3,IW,I),I=1,12) 
                  READ(KR(18),1)(probWetDay(1,IW,I),I=1,12)    ! probWetDay(1,): Probability of wet day after wet day                                      
                  READ(KR(18),1)(probWetDay(2,IW,I),I=1,12)    ! probWetDay(2,): Probability of wet day after dry day                                      
                  READ(KR(18),1)(UAVM(I),I=1,12)               ! Average number of precipitation in each month                             
                  READ(KR(18),1)(halfHrRainInfo(IW,I),I=1,12)                                              
                  READ(KR(18),1)(monAveSoilRad(IW,I),I=1,12)                                            
                  READ(KR(18),1)(RH(IW,I),I=1,12)                                              
                  READ(KR(18),1)(monWindSpeed(I),I=1,12)       ! UAV0: Average wind speed of each month                                             
            END IF readMonthlyData
                      
            WRITE(KW(1),2)                                   ! write the average month values  

            ! ------------------- get the montly T, minimum monthly T -------------
            monAveT(1)=.25*(monAveTmax(IW,12)+monAveTmin(IW,12)+monAveTmax(IW,1)+monAveTmin(IW,1))                      
            upperT=monAveT(1)
            lowerT=upperT   
          
            JT1=1                                                                          
            IF(monAveTmin(IW,1)>monAveTmin(IW,12)) JT1=12                                               
            TMN=monAveTmin(IW,JT1)               ! TMN: Global variable:  the mininum temperature                                    
          
            DO I=2,12                                                                      
                  I1=I-1                                                                       
                  monAveT(I)=.25*(monAveTmax(IW,I)+monAveTmin(IW,I)+monAveTmax(IW,I1)+monAveTmin(IW,I1))    
              
                  IF(monAveT(I)>upperT)THEN      ! get the lower (lowerT) and upper (upperT) bounds of average T
                        upperT=monAveT(I)                   
                  ELSE
                        IF(monAveT(I)<lowerT)lowerT=monAveT(I)
                  END IF
              
                  IF(monAveTmin(IW,I)>TMN)CYCLE                                                      
                  JT1=I                                                                        
                  TMN=monAveTmin(IW,I)       ! get the minimum Tmin of the year                                                             
            END DO      
      
            amplitudeT=upperT-lowerT             ! Amplitude of average T for a year
     
            ! --------------- monthly rainfall ----------------- 
            DO I=1,12         !check if the std and skew coefficient are available                                                  
                  IF(dailyRainInfo(2,IW,I)<1.E-5.OR.dailyRainInfo(3,IW,I)<1.E-5)EXIT
            END DO
      
            ! Calcualte rainfall from two different methods 
            IF(I>12)THEN                                                                         
                  rainDistriFlag=0    
            ELSE                                                                        
                  rainDistriFlag=1    
                  SUM=0.                                                                         
                  DO I=1,10000                                                                   
                        XX=Generate_Random(randID(3))     !? Generate_Random                                                        
                        SUM=SUM+(-LOG(XX))**rainPara                                                    
                  END DO                                                                         
                  REXP=10100./SUM                                                                
            END IF
      
            BIG=0.
            ADD=0.                                                                     
            V3=Generate_Random(randID(3))  
      
            !  Using probability to cal rainfall of each month 
            monthLoop1: DO I=1,12                                                                   
                  I1=I+1                                                                         
                  XM=NC(I1)-NC(I)                ! NC:accumulated days by month                                                     
                  DayOfYear=INT((NC(I1)+NC(I))*.5) ! use the half month day to cal Daylength and radiation
                  CALL Cal_DL_MaxRd              ! compute day length and maximum radiation                                                       
                  monMaxRad(I)=maxRad                                                                   
                  monDayLen(I)=DayLength         ! DayLength: day length in hrs                                                            
          
                  IF(DayLength>BIG) BIG=DayLength! the largest DL in a year                                                     
          
                  XYP(I)=0.                      ! XYP(200) : 1:12                                                               
                  XX=monTmaxStd(IW,I)-monTminStd(IW,I)  ! XX: std of Tmin and Tmax of each month                                      
          
                  IF(XX>10.)THEN                 ! no more than 10 degree?                                                             
                        monTmaxStd(IW,I)=(monTmaxStd(IW,I)-monAveTmax(IW,I))*.25  
                        monTminStd(IW,I)=(monAveTmin(IW,I)-monTminStd(IW,I))*.25  
                  ENDIF      
          
                  ! point rainfall:  Eq. 1-5 in APEX doc on page 10 and 11 
                  IF(probWetDay(1,IW,I)>0.)THEN  ! probWetDay(1,): Probability of wet day after wet day   
                        UAVM(I)=XM*probWetDay(1,IW,I)/(1.-probWetDay(2,IW,I)+probWetDay(1,IW,I))           
                  ELSE ! if the wet-dry probabilities are not aviailable, the average monthly number of rainy days may be subsituted                                                                   
                        probWetDay(1,IW,I)=wetDayCoeff *(UAVM(I)+.0001)/XM        ! XM: the number of days in a month                                   
                        probWetDay(2,IW,I)=1.- wetDayCoeff + probWetDay(1,IW,I)   ! Eq.5 in APEX-doc on page 11                                           
                  ENDIF         
          
                  dailyRainInfo(1,IW,I)=monAvePPT(IW,I)/(UAVM(I)+.01)     ! Daily rainfall stored in dailyRainInfo(1,:,:,)
          
                  X2=SQRT(monAveTmax(IW,I)-monAveTmin(IW,I))  
          
                  IF(monAveSoilRad(IW,I)<=0.)THEN                      ! If Ave month soil Radiation is not provided                                                    
                        X1=MAX(.8,.21*X2)                                        
                        monAveSoilRad(IW,I)=X1*maxRad  
                  END IF 
          
                  !----------------  point evaporation -----------------
                  ! seems the mix use of Penman and Hargreaves methods
                  TX=.5*(monAveTmax(IW,I)+monAveTmin(IW,I))      ! monthly average T
                  XL=2.501-2.2E-3*TX                             ! Eq. 93c in APEX doc  Penman on page 29-30
                  RAMM=maxRad/XL                                 ! Clear day radiation MJ/m^2/d 
                  EO=modelPara(38)*RAMM*(TX+17.8)*X2             ! Eq. 97 Hargreaves method: point evaporation  
                  ADD=ADD+EO                                     ! Total evaporation mm/day
          
                  !--------------------------- point rainfall ?------------------------------
                  ! APEX-doc Eq.1 on page 10
                  IF(rainDistriFlag==0)THEN                      ! rainDistriFlag = 0: rainfall distribution is skewed-normal
                        SUM=0. 
                        rainfallVolum=dailyRainInfo(1,IW,I)      ! montly rainfall volume  
                        rainfallStd=dailyRainInfo(2,IW,I)        ! monthly sd of rainfall 
                        rainfallSkew=dailyRainInfo(3,IW,I)       ! monthly Skew coefficient of rainfall distribution  
                        R6=rainfallSkew/6.                       ! Eq. 1c  
                        DO J=1,1000   
                              V4=Generate_Random(randID(3))   
                              XX=SQRT(-2.*LOG(V3))*COS(6.283185*V4) ! COMPUTES A STANDARD NORMAL DEVIATE GIVEN TWO RANDOM NUMBERS  
                              V3=V4  
                              R1=Sim_Rain(R6, XX )                ! calculate daily rainfall from skewed distribution  
                              SUM=SUM+R1                          ! 1000 iteraction 
                        END DO                                                                         
                        PCF(IW,I)=1010.*dailyRainInfo(1,IW,I)/SUM  
                  ELSE                                            ! uniform distribution  
                        dailyRainInfo(1,IW,I)=dailyRainInfo(1,IW,I)*REXP 
                        PCF(IW,I)=1.  
                  ENDIF
          
            END DO monthLoop1 ! end month loop
            ! ================================================================================================
      
            XYP(1)=monAveTmax(IW,1) 
            BIG=monAveSoilRad(IW,1) 
            maxRH=RH(IW,1)                       ! relative humudity  
            maxRainfall=monAvePPT(IW,1)          ! Average monthly rainfall         
            extremeRain=halfHrRainInfo(IW,1)     ! extreme rain within a month
      
            DO I=2,12   
                  IF(monAveSoilRad(IW,I)>BIG)BIG=monAveSoilRad(IW,I)                  ! Store the maximum monthly soil radiation  
                  IF(monAvePPT(IW,I)>maxRainfall)maxRainfall=monAvePPT(IW,I)          ! Store the maximum monthly average T 
                  IF(RH(IW,I)>maxRH)maxRH=RH(IW,I)                                    ! Store the maximum RH   
                  IF(halfHrRainInfo(IW,I)>extremeRain)extremeRain=halfHrRainInfo(IW,I)! Store the maximum rain at 0.5h  
                  XYP(1)=XYP(1)+monAveTmax(IW,I)                                      ! Store accumulated montly Tmax  
            END DO  
      
            XYP(1)=XYP(1)/12.                    ! Average monthly Tmax   
            radCorrection=1.
            IF(BIG>100.)radCorrection=.04184     ! BIG is the maximum soil radiation  
            X3=.3725/(XYP(1)+20.)                !???   
	    
            !----------------------------- Rainfall and humuidity --------------------------------
            ! APEX doc p15-16   p35-36    
            vernalizeDay=0.    
            monthLoop2: DO I=1,12   
                  XM=NC(I+1)-NC(I)                  ! days in one month    
                  wetDayFreq(IW,I)=UAVM(I)/XM       ! wet day frequency?   
                  XYP(2)=XYP(2)+monAveTmin(IW,I)    ! accumulation of Tmin
                  XYP(3)=XYP(3)+monAvePPT(IW,I)     ! accumulation of PPT  
                  XYP(4)=XYP(4)+UAVM(I)             ! accumulation of precipitation days  
                  monAveSoilRad(IW,I)=radCorrection*monAveSoilRad(IW,I)  
                  XYP(5)=XYP(5)+monAveSoilRad(IW,I)  
                  X1=MAX(monAvePPT(IW,I),12.7)      ! at least  12.7 for simulations?   
                  TX=.5*(monAveTmax(IW,I)+monAveTmin(IW,I)) ! monthly average T 
                  IF(maxRH>1.)THEN                  ! >1.0 : the input is actual temperature for cal RH
                        RH(IW,I)=Sat_Vapor(RH(IW,I)+273.)/Sat_Vapor(TX+273.) ! Saturated Vapor calcuated given T   
                  ELSE   
                        IF(RH(IW,I)<1.E-10)THEN
                              XX=monAveTmax(IW,I)-monAveTmin(IW,I)  
                              RH(IW,I)=.9-.8*XX/(XX+EXP(5.122-.1269*XX)) ! Eq. 26 in APEX-doc  
                        END IF   
                  END IF
                  X2=MAX(TX,-1.7)                                  ! at least -1.7 ?   
                  XYP(6)=XYP(6)+((X1/25.4)/(1.8*X2+22.))**1.111    ! ????  
                  X1=monAvePPT(IW,I)/(UAVM(I)+1.E-10)            ! average monthly rainfall/average number of days of rainfall  
                  IF(extremeRain<1.)THEN
                        IF(extremeRain<1.E-10)  halfHrRainInfo(IW,I)=MAX(.05,peakRateFactor*X3*(monAveTmax(IW,I)+20.))  
                        XTP(I)=5.3*X1*halfHrRainInfo(IW,I)  
                        CYCLE
                  END IF
                  ! Eq. 125-126 in APEX-doc P36  
                  XY2=.5/periodHalfHrRain  
                  F=XY2/(UAVM(I)+.01)     
                  XTP(I)=halfHrRainInfo(IW,I)  
                  halfHrRainInfo(IW,I)=-XTP(I)/LOG(F)  
                  halfHrRainInfo(IW,I)=peakRateFactor*halfHrRainInfo(IW,I)/(X1+1.)  !?? questional why halfHrRainInfor was assgined twice? 
                  IF(halfHrRainInfo(IW,I)<.1)  halfHrRainInfo(IW,I)=.1   
                  IF(halfHrRainInfo(IW,I)>.95) halfHrRainInfo(IW,I)=.95  
                  X1=1.4-.0778*TX  
                  X2=.5+.37*TX  
                  X1=MIN(1.,X1,X2)  
                  IF(X1<=0.)CYCLE
                  vernalizeDay=vernalizeDay+X1*XM    ! Vernalization time (Days) 
            END DO monthLoop2   
      
            XYP(2)=XYP(2)/12.            ! Average Tmin (yearly)  
            XYP(5)=XYP(5)/12.            ! Average soil radiation (yearly)     
            yearTotPPT=XYP(3)            ! Total ppt of a year
      
            IF(baseStockRate<1.E-10)THEN   
                  ADD=30.44*ADD          ! ADD: total Evaporation
                  RTO=ADD/yearTotPPT     ! tot EP/ tot rainfall for a year
                  F=RTO/(RTO+EXP(S_Curve(25,1)-S_Curve(25,2)*RTO)) ! S(25): regulate denitrification as a function of soil water content
                  baseStockRate=10.*F+1.
            END IF  
      
            IF(IW>1)THEN   
                  II=IW-1  
                  WRITE(KW(1),3)II,TITW5(1) 
            ELSE  
                  WRITE(KW(1),4)WPM1FILE  ! write the title of each month 
            END IF  
      
            ! ------------------------ write to KW(1) sitename.out ------------------------
      
            WRITE(KW(1),5)varName(1),(monAveTmax(IW,I),I=1,12),XYP(1),varName(1)  
            WRITE(KW(1),5)varName(2),(monAveTmin(IW,I),I=1,12),XYP(2),varName(2) 
            WRITE(KW(1),6)'SDMX',(monTmaxStd(IW,I),I=1,12),'SDMX' 
            WRITE(KW(1),6)'SDMN',(monTminStd(IW,I),I=1,12),'SDMN' 
            WRITE(KW(1),9)varName(4),(monAvePPT(IW,I),I=1,12),yearTotPPT,varName(4) 
            WRITE(KW(1),8)'SDRF',(dailyRainInfo(2,IW,I),I=1,12),'SDRF'  
            WRITE(KW(1),6)'SKRF',(dailyRainInfo(3,IW,I),I=1,12),'SKRF' 
            WRITE(KW(1),7)'PW/D',(probWetDay(1,IW,I),I=1,12),'PW/D'  
            WRITE(KW(1),7)'PW/W',(probWetDay(2,IW,I),I=1,12),'PW/W' 
            WRITE(KW(1),5)'DAYP',(UAVM(I),I=1,12),XYP(4),'DAYP'  
            WRITE(KW(1),8)'P5MX',(XTP(I),I=1,12),'P5MX' 
            WRITE(KW(1),9)varName(3),(monAveSoilRad(IW,I),I=1,12),XYP(5),varName(3) 
            WRITE(KW(1),8)'maxRad',monMaxRad,'maxRad'  
            WRITE(KW(1),6)'DayLength',monDayLen,'DayLength'  
            WRITE(KW(1),6)'RHUM',(RH(IW,I),I=1,12),'RHUM'  
            WRITE(KW(1),6)'ALPH',(halfHrRainInfo(IW,I),I=1,12),'ALPH'  
            WRITE(KW(1),6)' PCF',(PCF(IW,I),I=1,12),' PCF' 
            IF(WP5_ID==0)EXIT      ! no more reading IW = 2...6  

      END DO readdata  
    
      REWIND KR(18) ! ========================== The end of read from montlhy weather (.WP1) data ===========================
      XYP(6)=115.*XYP(6)                           !???  
      yearAveTem=(XYP(1)+XYP(2))*.5                ! average T over a year 

      ! LatChoice: 0 for using input latitudes for subareas; >0 computing equalent latitude based on azimuth
      IF(LatChoice >0)THEN
            X1=ASIN(uplandSteep)
            YLAZ=YLAT/CLT
            X2=siteAzimuth/CLT
            YLAZ=CLT*ASIN(uplandSteep*COS(X2)*COS(YLAZ)+COS(X1)*SIN(YLAZ))
            WRITE(KW(1),'(T10,A,F8.3)')'EQUIVALENT LATITUDE = ',YLAZ
      ELSE
            YLAZ=YLAT          
      END IF
      
      XX=YLAZ/CLT
      YLTS=SIN(XX)   
      YLTC=COS(XX)
      YTN=TAN(XX)
      JDHU=400                  ! What is JDHU?      JDHU is DOY and it can only be 1-366 ?                                                    
      winterDormancy=meanDayLength  ! Winter dormancy(h) causes dormancy in winter growth crops. no growth when day length < winter_dormancy                                                                    
      IF(meanDayLength<11.)THEN   
            JDHU = Cal_DOY(JT1,15, leapYr)                 ! JT1 is the month with minimum Tmin  
            winterDormancy=modelPara(6)+meanDayLength  ! Annual minimum day length + PARM(6)  
      END IF   

      ! ===============================  read wind data =======================================
      ! ------------------ search for the wind file name-------------------
      ! If...ELSE...END IF: choose a wind file from WINDUSEL.DAT if the wind file is not provided
      CALL OPENV(KR(19),FWIND,IDIR(1),KW(MSO) ) ! Catalog of weather stations and monthly wind data  
      IF(windID==0)THEN                           ! Wind_ID: wind station #
            W0=1.E20     
            DO    
                  READ(KR(19),*,IOSTAT=NFL)II,OPSCFILE,Y,X,ELEX   
                  IF(NFL/=0)EXIT  
                  RY=Y/CLT
                  XX=SIN1*SIN(RY)+COS1*COS(RY)*COS((X-XLOG)/CLT)
                  D=6378.8*ACOS(XX)
                  E=ABS(ELEV-ELEX)
                  W1=modelPara(79)*D+(1.-modelPara(79))*E                                                                  
                  IF(W1>=W0)CYCLE                                                                
                  W0=W1           
                  WINDFILE=OPSCFILE                     ! NCCLAYTO.WND is the nearest station                                                        
            END DO  
      ELSE   
            II=-1  
            DO WHILE(II/=windID)  
                  READ(KR(19),*,IOSTAT=NFL)II,WINDFILE   
                  IF(NFL/=0)THEN
                        IF(runMode==0)THEN
                              WRITE(*,'(T10,A,I8,A)')'WIND NO = ',windID,' NOT IN'&
	                              ' MO WIND LIST FILE'
                              ERROR STOP                                                                  
                        ELSE
                              WRITE(KW(MSO),'(A,A8,A,I4,A)')' !!!!! ',siteName,&
	                              ' WIND NO = ',windID,' NOT IN MO WIND LIST FILE'
                              ERROR STOP
                        END IF 
                  END IF  
            END DO  
      END IF   
      REWIND KR(19)  
      ! ------------------------  read windfile  .WND ------------------------
      CALL OPENV(KR(25),WINDFILE,IDIR(1),KW(MSO))   
      !   LINE 1/2                                                                       
      READ(KR(25),'()')
      READ(KR(25),'()')
                                            
      !   UAVM = AV MO WIND SPEED (M/S)(REQUIRED TO SIMULATE WIND                        
      !   EROSION--Wind_ero_fact>0 LINE 24  AND POTENTIAL ET IF PENMAN OR                   
      !            PENMAN-MONTEITH EQS ARE USED--LINE 4)                                   
      !   LINE 3                                                                         
      READ(KR(25),1) UAVM    
      ! ------------------------------- Wind erosion  ------------------------------- 
      AWV=0.    
      windEnergy=0.   
      
      monthLoop3: DO I=1,12   
            RNCF(I)=1.   
            TMNF(I)=1.  
            IF(UAVM(I)<1.E-5) UAVM(I)=monWindSpeed(I)   
            SMY(I)=0.                                                                         
            DO J=1,100                                        ! 100 simulations ?  
                  RN2=Generate_Random(randID(5))                                                        
                  WV=UAVM(I)*(-LOG(RN2))**windPowerPara       ! WV: daily average wind velocity (m/s)  
                  IF(WV<modelPara(67))CYCLE                   ! modelPara: wind erosion threshold wind speed, normal value :6   
                  EV=193.*EXP(1.103*(WV-30.)/(WV+1.))         ! Eq. 2.131 in EPIC-doc P32   EV: daily wind energy   
                  SMY(I)=SMY(I)+EV    
            END DO                                                                       
            windEnergy=windEnergy+SMY(I)                      ! total wind energy of a year  
            AWV=AWV+UAVM(I)  
      END DO monthLoop3
                     
      AWV=AWV/12.                                             ! Average wind speed of a year  
      windEroFactor = (3.86*AWV**3/XYP(6)**2)**.3484          ! ???  
      WRITE(KW(1),'(T10,A,F7.3)')'WIND EROS CLIMATIC FACTOR = ',windEroFactor 
      X1=MAX(1.,yearAveTem)                                        ! Yearly average T   
      CLF=SQRT(yearTotPPT/X1)                                      ! Total ppt / yearly average T ?  
      WRITE(KW(1),'(T10,A,F7.3)')'CLIMATIC FACTOR = ',CLF  
      
      ! DIR  = AV MO FRACTION OF WIND FROM 16 DIRECTIONS (BLANK UNLESS                 
      !        WIND EROSION IS SIMULATED--Wind_ero_fact>0 LINE 24).        
      windDir: DO J=1,16    
            ! LINES 4/19  
            READ(KR(25),1)(windDirection(I,J),I=1,12)  
            IF(windDirection(1,J)>0.)CYCLE  
            DO I=1,12  
                  windDirection(I,J)=1.  
            END DO   
      END DO  windDir  
      REWIND KR(25)  
      WRITE(KW(1),'(/T20,A,A12)')'WIND = ',WINDFILE   
      WRITE(KW(1),5)varName(7),(UAVM(I),I=1,12),AWV,varName(7)        ! AWV: average wind speed
     ! ------------------------- normalize wind direction fraction ----------------------    
      monthLoop4: DO I=1,12                                                                   
            IF(UAVM(1)>0.)THEN                ! UAVM here is average monthly wind speed  
            CALL Integrate_ExpotentialEqs(windPowerPara, SUM) ! SUM is still vernalization time (days) here
            UAVM(I)=UAVM(I)/SUM                               ! ???  What does it mean?  
            END IF   

            DO J=2,16   
                  windDirection(I,J)=windDirection(I,J)+windDirection(I,J-1)        ! accumulate wind fractions 
            END DO 
 
            DO J=1,16  
                  windDirection(I,J)=windDirection(I,J)/windDirection(I,16)         ! normalized to the 16th direction  
            END DO
      END DO  monthLoop4

      ! --------------------------------- Rainfall distribution ---------------------------
      IF(rainDistriFlag==0)THEN  
            WRITE(KW(1),'(/T10,A)')'RAINFALL DIST IS SKEWED NORMAL'  
      ELSE   
            WRITE(KW(1),'(/T10,A,F5.2)')'RAINFALL DIST IS EXP--PARM = ', rainPara  
            WRITE(KW(1),'(/T10,A,F6.3)')'WET-DRY PROB COEF = ', wetDayCoeff  
      END IF  
 
      ! ---------------------------------  Theroretical Heat accumulation --------------------- 
      yearTotHeat = Cal_Heat(1,365,0.,1)         ! theoretical heat unit accumulation for use in CPTHU
                                                 ! MO was modified to 12 during calcuating heat
      TX=(monAveTmax(1,MO)+monAveTmin(1,MO))/2.  ! MO is 12 here.  Does it make sense?
      soilRad=monAveSoilRad(1,MO) 
      surfTem=TX 

      WRITE(KW(1),10)vernalizeDay,yearTotHeat ! SUM: vernalization time 
      WRITE(KW(1),'(//1X,A)')'-----WIND EROSION DATA' 
      WRITE(KW(1),11)fieldLength,fieldWidth,fieldAngle0,windPowerPara, soilDiameter,windErosionFactor  
      
      ! ---------------------------------- print methods to Estimate ET --------------------------
      ! EP_Meth = 4 for Hargreaves mothod
      IF(ET_Method==4.AND.modelPara(13)>0.)THEN  
            HGX=modelPara(13)              ! Hargreaves PET equation exponent: range 0.5-0.6
      ELSE
            IF(windEroFactor<.5)THEN
                  HGXW=.5
            ELSE
                  HGXW=MIN(.6,windEroFactor)
            END IF    
            IF(CLF>7.)THEN
                  HGX0=.5
            ELSE
                  HGX0=MIN(.6,1.2-.1*CLF)
            END IF 
            HGX=.5*(HGXW+HGX0)  
      END IF 
      
      SELECT CASE(ET_Method) 
            CASE(1)   
                  WRITE(KW(1),'(/T10,A,A)')'**********PENMAN-MONTEITH  EQ USED TO EST POT ET**********'   
            CASE(2) 
                  WRITE(KW(1),'(/T10,A,A)')'**********PENMAN EQ USED TO EST POT ET**********'  
            CASE(3) 
                  WRITE(KW(1),'(/T10,A,A)')'**********PRIESTLEY-TAYLOR EQ USED TO EST POT ET**********'  
            CASE(4)
                  WRITE(KW(1),'(/T10,A,F6.3,A)')'**********HARGREAVES EQ USED TO EST POT ET HGX = ', &
                            HGX,'**********'  
            CASE(5) 
                  WRITE(KW(1),'(/T10,A,A)')'********** BAIER-ROBERTSON EQ USED TO EST POT ET **********' 
            CASE DEFAULT
                  WRITE(KW(1),'(/T10,A,F6.3,A)')'**********HARGREAVES EQ USED TO EST POT ET HGX = ', &
                            HGX,'**********'   
      END SELECT  

      ! ======================= End of reading from wind speed file ================================
      1 FORMAT(12F6.2,8X)  
      2 FORMAT(/1X,'-----AVE MO VALUES')
      3 FORMAT(/5X,'STAGE =',I2,5X,'WPM5 = ',A80/T11,'JAN',6X,'FEB',6X,'MAR',6X,'APR',6X,'MAY',6X,'JUN',6X,'JUL',&
            6X,'AUG',6X,'SEP',6X,'OCT',6X,'NOV',6X,'DEC',6X,' YR')
      4 FORMAT(19X,'WPM1 = ',A12/T11,'JAN',6X,'FEB',6X,'MAR',6X,'APR',6X,'MAY',6X,'JUN',6X,'JUL',6X,'AUG',6X,'SEP',&
            6X,'OCT',6X,'NOV',6X,'DEC',6X,' YR')
      5 FORMAT(1X,A4,13F9.2,2X,A4)        
      6 FORMAT(1X,A4,12F9.2,11X,A4)
      7 FORMAT(1X,A4,12F9.3,11X,A4)           
      8 FORMAT(1X,A4,12F9.1,11X,A4)    
      9 FORMAT(1X,A4,13F9.1,2X,A4)
     10 FORMAT (T10,'VERNALIZATION TIME = ',F5.0,' D'/T10,'ANNUAL HEAT UNITS = ',F8.0,' C')       
     11 FORMAT(T10,'FIELD LENGTH = ',F6.2,' km'/T10,'FIELD WIDTH = ',F6.2,' km'/T10,'FIELD ANGLE = ',F4.0,' deg'/T10,&
            'WIND SPEED MOD EXP POWER modelPara = ',F6.2/T10,&
            'SOIL PARTICLE SoilPartl_dia = ',F5.0,' UM'/T10,'ACCELERATE WIND EROS FACTOR = ',F7.3)    
END SUBROUTINE MonthlyWeatherInfo

! --------------------------2. Search Daily Weather File List -------------------------------------------------------
SUBROUTINE SearchWeatherList(ID) 
      !     EPIC1102
      !     THIS SUBPROGRAM READS THE DAILY WEATHER LIST AND LOCATES THE 
      !     SPECIFIED STATION (WP_Daily) OR THE NEAREST STATION IF WP_Daily=0.
      IMPLICIT NONE
 
      ! local variables
      INTEGER:: II, ID
      REAL:: W0, W1, Y, X, ELEX, RY, D, E, XX  
 
 
      IF(dailyWeatherID == 0)THEN
          W0=1.E20
          DO 
              READ(KR(27),*,IOSTAT=NFL)II,OPSCFILE,Y,X,ELEX
                IF(NFL/=0)EXIT
              RY=Y/CLT
              XX=SIN1*SIN(RY)+COS1*COS(RY)*COS((X-XLOG)/CLT)
              D=6378.8*ACOS(XX)
              E=ABS(ELEV-ELEX)
              W1=modelPara(79) *D+(1.- modelPara(79))*E                                                                  
              IF(W1>=W0)CYCLE                                                                
              W0=W1        
              FWTH=OPSCFILE
            END DO
      ELSE
          II=-1    
          DO WHILE(II/=dailyWeatherID)
              READ(KR(27),*,IOSTAT=NFL)II,FWTH
                IF(NFL/=0)THEN
                    IF(runMode==0)THEN
                        WRITE(*,*)'FWTH NO = ',dailyWeatherID,' NOT IN DAILY WEATHER LIST FILE'
                        STOP
                    ELSE
                      WRITE(KW(MSO),'(A,A8,A,I4,A)')'!!!!! ',siteName,' FWTH NO = ',&
                      &dailyWeatherID,' NOT IN DAILY WEATHER LIST FILE'
                  END IF
                  STOP
              END IF                  
            END DO
        END IF 
      REWIND KR(27)
      CALL OPENV(KR(7),FWTH,ID,KW(MSO))
      CALL ReadDailyWeather0 
 
END SUBROUTINE SearchWeatherList
               
! ---------------------------- 3. read daily weather file to the day before simulation -------------
SUBROUTINE ReadDailyWeather0 
      !     EPIC1102
      !     THIS SUBPROGRAM READS THE DAILY WEATHER FILE TO THE DAY BEFORE THE
      !     SIMULATION BEGINS.
      IMPLICIT NONE
 
      ! local variables potential problems with local variables (SAVE attribute XTP)
      INTEGER, DIMENSION(12):: MOFD = [31,29,31,30,31,30,31,31,30,31,30,31]
      INTEGER:: I1, I2, I3, J1, J2, J3, II, K1, yr, N1, L 
      REAL:: XTP(6)
      
      READ(KR(7),2)I3,I2,I1,(XTP(L),L=1,6)
      J3=10000*I3
      J1=100*I2+J3
      II=I1+J1
      K1=I2
      IF(II>= dateNum)THEN           ! if the beginning year is larger than the start year of simulation return
            dateNum=II
            yr=I3
            REWIND KR(7)
            RETURN
      ELSE
            DO
                  CALL Leap_Yr_Check(I3, leapYr, leapYrFlag) 
                  DO I2=K1,12
                        N1=MOFD(I2)
                        IF(I2==2)N1=N1-leapYr
                        J2=100*I2
                        J1=J2+J3
                        DO WHILE(I1<N1)
                              I1=I1+1
                              II=J1+I1
                              IF(II==dateNum)RETURN
                              READ(KR(7),1,IOSTAT=NFL)(XTP(L),L=1,6)! Weather of every day SAVE to XTP. but how XTP will be used
                              IF(NFL/=0)THEN
                                    WRITE(KW(MSO),'(T10,A)')'START DATE EXCEEDS&
                                          & WEATHER FILE--WEATHER GENERATED'
                                    weatherVarID=0
                                    RETURN
                              END IF
                        END DO
                        I1=0
                  END DO
                  K1=1
                  I3=I3+1
                  J3=10000*I3
            END DO
      END IF

    1 FORMAT(14X,6F7.2)
    2 FORMAT(2X,3I4,6F7.2)
    !1 FORMAT(14X, 7F6.0)
    !2 FORMAT(2X, 3I4, 7F6.0)
END SUBROUTINE ReadDailyWeather0

! --------------------------- Read or simulate daily weather in the Simulation Program---------------------------------
SUBROUTINE ReadDailyWeather1 
      IMPLICIT NONE
      ! arguments list
 
      ! local variables
      INTEGER:: III, IRH=0      
      REAL:: X1      
      ! READ WEATHER VARIABLES INCLUDED IN Climate_ID
      IF(IGSD==0.OR.IY/=IGSD)THEN                              ! what is IGSD??
            DO                                                 ! DO...Exit...END DO
                  !     READ DAILY WEATHER IF NOT GENERATED                                            
                  !  1  SRAD   = SOLAR RADIAION(MJ/M2 OR LY)(BLANK TO GENERATE                           
                  !  2  TMX  = MAX TEMP(C)                                                             
                  !  3  monTmin  = MIN TEMP(C)                                                             
                  !  4  Rainfall  = RAINFALL(mm)(999. TO GENERATE OCCURRENCE & AMOUNT;                      
                  !                               -1. TO GENERATE AMOUNT GIVEN OCCURRENCE)                                
                  !  5  RHD  = RELATIVE HUMIDITY (FRACTION)(BLANK TO GENERATE                          
                  !  6  U10  = WIND VELOCITY (M/S)(BLANK TO GENERATE                                   
                  !  7  X1   = ATMOSPHERIC CO2 CONC (ppm)
                  !  8  EVI  = VEGETATION INDEX                                              
                  READ(KR(7),103,IOSTAT=NFL)SRAD,TMX,TMN,Rainfall,RHD,U10,X1,EVI,peakRainRate   
                  !RHD = RHD/100.0   ! Check RHD when reading it from your data files
                  
                  IF(NFL==0.AND.(Date<=DOY_realtime.OR.DOY_realtime==0))THEN
                        IF(X1>0..AND.CO2_Method==2)CO2=X1     ! CO2_mod_flag: = 0 for constant atmospheric CO2                                        
                        SRAD=SRAD*radCorrection                              ! = 1 for dynamic atmospheric CO2  
                        III=0                                                 ! = 2 for inputing CO2 
                  ELSE 
                        weatherVarID=0                                                                       
                        WRITE(KW(1),'(/T10,3I4,A,A)')IYR,MO,DayOfMon,' RAIN, TEMP, RAD, WIND'&            
                              ,' SPEED, & REL HUM ARE GENERATED'                                             
                        EXIT
                  END IF    
                  
                  ! ------------------- if rainfall is missing, generate it---------------------------
                  IF(Rainfall<0.)THEN                                                   
                        CALL Generate_Rainfall(1 )     
                  ELSE
                        IF(Rainfall>900.)CALL Generate_Rainfall(0)
                  END IF        
                  Rainfall=Rainfall*RNCF(MO)  ! RNCF: rainfall frequency of each month
                  
                  ! ------------------- if temperature is missing, generate it-------------------------
                  IF(KGN(2)==0.OR.TMX<=TMN)THEN      ! KGN: control temperature
                        CALL Para_SimuWeather          ! return TXXM, RHM, RM, WX, XIM                                                                         
                        CALL Generate_AirTem           ! return TMX, TMN
                        X1=TMX-TMN                                                                     
                        TMX=TMX+TMXF(MO)                                                               
                        TMN=TMX-TMNF(MO)*X1 
                  ELSE                                                                                                                                    
                        IF(TMX>100..OR.TMN>100.)THEN
                              CALL Para_SimuWeather      ! return TXXM, RHM, RM, WX, XIM                                                                         
                              CALL Generate_Tem          ! return: TMX, TMN                                                                    
                              X1=TMX-TMN                                                                     
                              TMX=TMX+TMXF(MO)                                                               
                              TMN=TMX-TMNF(MO)*X1
                        END IF
                  END IF 
                  !-------------------- if solar radiation is missing, generate it ------------------------
                  IF(SRAD<1.E-5.OR.SRAD>900..OR.KGN(3)==0)THEN               ! solar radiation                                             
                        IF(III==0)CALL Para_SimuWeather        ! return RHM, RM, WX, XIM                                                             
                        CALL Generate_SolarRad                                                               
                  END IF
                  !--------------------- if relative humudity is missing, generate it-----------------------
                  IF(RHD<0..OR.RHD>1..OR.IRH>0.AND.RHD<900.)THEN    ! IRH is not initialized here?
                        ! RHD=DEW POINT TEMP--CONVERT TO RELATIVE HUM
                        TX=.5*(TMX+TMN) 
                        RHD=MIN(.99,Sat_Vapor(RHD+273.)/Sat_Vapor(TX+273.))
                        IRH=1                                                                    
                  ELSE
                        IF(RHD<1.E-5.OR.RHD>900..OR.KGN(5)==0)THEN 
                              ! RHD IS MISSING DATA (0. OR 999.) OR NOT INCLUDED IN Climate_ID                                                
                              IF(III==0)CALL Para_SimuWeather     ! return RHM, RM, WX, XIM                                                             
                              CALL Generate_RH                    ! return RHD                                                                  
                        END IF                                                                         
                  END IF    
                  !------------------------if wind speed is missing, generate it-----------------------------
                  IF(U10<1.E-5.OR.U10>900..OR.KGN(4)==0)CALL Generate_WindSpeed   ! return U10
                  EXIT
            END DO                                                     
      END IF
  103 FORMAT(14X,9F7.2)
  !103 FORMAT(14X,9F6.0) 
END SUBROUTINE ReadDailyWeather1
 
END MODULE ReadWeather_Module