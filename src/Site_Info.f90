SUBROUTINE Site_Info(FSITE, siteID, coeffChannelGeometry, peakRateFactor, snowCont0, siteAzimuth,  &
                     channelLen, channelSlp, channelDep, surfNForChannel, nManningCoeff, &
                     drainageDepth, meanDayLength, RXM)
USE PARM
USE MisceFun_Module
IMPLICIT NONE

! argument list 
CHARACTER(LEN = 80), INTENT(IN):: FSITE
INTEGER, INTENT(IN):: siteID
REAL, DIMENSION(6), INTENT(IN)::coeffChannelGeometry  ! coefficients for channel geometry: 
                                          ! X = Channel_geometry(N)*areaWshed** channel geometry(N+1) 
                                       ! X = depth for Channel_geometry(3) and (4); X = length for Channel_geometry(5) and (6) 
REAL, INTENT(OUT):: peakRateFactor, snowCont0, &  ! water con of snow on ground at start of simulation (mm)
                  siteAzimuth,           & ! azimuth orientation of land slope (degrees clockwise from north)
                  channelLen,            & ! Mainstream channel length (km, blank if unknown) 
                  channelDep,            & ! Channel depth (m)    
                  channelSlp,            & ! Mainstream channel slope (m/m, blank if unknown) 
                  surfNForChannel,       & ! Surface N for channel (blank if unknown)
                  nManningCoeff,         & ! Manning's roughness coefficient (APEX doc Eq. 37 and Eq. 76 )
                  meanDayLength, RXM
INTEGER,INTENT(OUT)::drainageDepth   
! local variables declaration
INTEGER:: I, II
REAL:: siteCO2, & ! CO2 concentration at the site (ppm)-- none zero value overrides CO2 input in EPICCONT. DAT 
       siteNO3, & ! NO3 concentration in irrigation water at this site (ppm) -- none zero value overrides CO2 input in EPICCONT. DAT
       rainNcon, CH, H, XX, X1, X2, X4, XM, SFL, TSF

! ///////////////////////////////////////////////////////////////////////////////
! ========== 1. Search the sitecom.dat to according to the Site_ID ================
! check if a list of site files exists
CALL OPENV(KR(23),FSITE,IDIR(2),KW(MSO) ) ! Catalog of site files  

II=-1  
DO WHILE(II/=siteID)                       ! search table   
      READ(KR(23),*,IOSTAT=NFL)II,SITEFILE    ! Site table: sitefile = umstead.sit   
      IF(NFL/=0)THEN
            IF(runMode==0)THEN                                               
                  WRITE(*,*)'SITE NO = ',siteID,' NOT IN SITE LIST FILE'
                  ERROR STOP                    
            ELSE
                  WRITE(KW(MSO),'(A,A8,A,I4,A)')' !!!!! ',siteName,&
                        'SITE NO = ',siteID,' NOT IN SITE LIST FILE'
                  ERROR STOP
            END IF
    END IF
END DO                                                                              
REWIND KR(23)    

! ================= 2. read data in .SIT file =================================
      !=========================== 13. read site file =======================================
      !  .sit file: site information      
CALL OPENV(KR(1), SITEFILE, IDIR(2), KW(MSO) ) ! check the site file: management practice                                      
!     TITLE=PROBLEM DESCRIPTION(3 LINES)                                             
!     LINES 1/3                                                                      
READ(KR(1),1)(TITLE(I),I=1,60)    
 !     LINE 4             11 variables                                                            
READ(KR(1),2)YLAT,XLOG,ELEV,peakRateFactor,siteCO2,siteNO3,rainNcon,X1,X2,snowCont0,siteAzimuth                                                                                              
!     LINE 5             10 variables                                                            
READ(KR(1),2)areaWshed,channelLen,channelSlp,channelDep,nManningCoeff,surfNForChannel, &
            uplandSlopLen,uplandSteep,conservFactor,timeIntervalGasDiff                             
  
! Default values applied if not provided in site file
IF(conservFactor<1.E-10)conservFactor=conservFactor0
IF(uplandSteep<1.E-10)uplandSteep=.001               
IF(nManningCoeff<1.E-10)nManningCoeff=.05     ! ManningN_Channel: Manning' roughness coefficient                                          
IF(surfNForChannel<1.E-10)surfNForChannel=.15 ! surfNForChannel:  Eq 131a in APEX-doc P37                            
IF(peakRateFactor<1.E-10)peakRateFactor=1.
IF(timeIntervalGasDiff<1.E-10)timeIntervalGasDiff=1.
snowWater=snowCont0                           ! snowCont0: water content of snow on groud at start of simulation (mm) 
!     READ MANAGEMENT INFORMATION                                                    
!     LINE 6                                                                         
READ(KR(1),3)irrCode,autoIrrInterval,fertApplyMinInterval,furrowFlag,drainageDepth,fertTimes, &
            autoManureTrigger,mowMinInterval,numAutoP  

! split the Irrigation code: IrrInput_choice = 0:   applies volume defined by ARMX(MaxIrr_sigle); 
!                                             =1:   applies input or ARMX                                                
CALL split_int(irrCode, irrChoice)

! update data according to the specific site 
! CO20, CNO30, rainNcon0: read from EPICCONT file
! CO2X, CNO3X, rainNcon: read from .sit file 
CO2=CO20                                                                            
IF(siteCO2>0.) CO2=siteCO2  
CNO3I=CNO30 
IF(siteNO3>0.) CNO3I=siteNO3 
CNO3I=CNO3I*.01  
rainNconX=rainNcon0                 ! from EPICCONT file  
IF(rainNcon>0.)rainNconX = rainNcon ! from site file 
rainNconX=rainNconX*.01 

!========================== 16. Potential Evaporation: Penman =================================
barometricPressure=101.3-ELEV*(.01152-5.44E-7*ELEV)   ! Eq. 93e in APEX doc    PB: barometric pressure in kPa
psychrometerConst=6.595E-4*barometricPressure         ! Eq. 93d in APEX doc    GMA: psychrometer constant in kPa/degree  

! ---------------------- Prepare data for Water Erosion ------------------------------
WSX=1.+ areaWshed                  ! why + 1 ?  Watershed Area(ha)  
! APEX-doc 112-113 P34
IF(channelLen<1.E-10)channelLen=coeffChannelGeometry(5)*areaWshed**coeffChannelGeometry(6)   
IF(channelSlp<1.E-10)channelSlp=uplandSteep*WSX**(-.3)  
IF(channelDep<1.E-10)channelDep=coeffChannelGeometry(3)*areaWshed**coeffChannelGeometry(4)   

! Compute SL for different water erosion models  Ref. APEX-doc 112- 113
! Eq. 2.94 in EPIC-doc on P26/ Eq. 112 and 112a APEX-doc P34                                              
XM=.3*uplandSteep/(uplandSteep+EXP(-1.47-61.09*uplandSteep))+.2   
uplandSlopLen=MIN(uplandSlopLen,SQRT(10000.*areaWshed))    
slopeLength=uplandSlopLen/22.127  
SL=slopeLength**XM*(uplandSteep*(65.41*uplandSteep+4.56)+.065) ! SL is calcuated for USLE model  
! calculation for RUSLE model: Eq. 113d - 113g in APEX-doc 
X1=3.*uplandSteep**.8+.56  
BETA=uplandSteep/(.0896*X1) 
RXM=BETA/(1.+BETA)  
RLF=slopeLength**RXM              ! SL = RSF*RLF is for RUSLE model

IF(uplandSlopLen>4.57)THEN        ! Eq. 113a-b   
      IF(uplandSteep>.09)THEN
            RSF=16.8*uplandSteep-.5  
      ELSE                                                                
            RSF=10.8*uplandSteep+.03  
      END IF
ELSE
      RSF=X1   
END IF  
! ---------------------------- cal peak rate --------------------------------
SX=SQRT(uplandSteep)                      ! uplandSteep: upperland slop steepness                                                                 
SX=SX/surfNForChannel                      ! SX = WSX in Eq. 131a APEX-doc P37 
IF(peakRateMethod > 0)THEN                ! peak rate estimated using TR55 method (Ref. APEX-doc P24)
      IF(channelLen>.1)THEN             
            SFL= 0.05                     ! SFL: shallow flow length in km     
      ELSE                             
            IF(channelLen>.05)THEN
                  SFL= channelLen-.05                                                             
            ELSE                                      
                  SFL=0.                                                                         
            END IF
      END IF
      TSF=SFL/MIN(2.16, 17.7*SX*surfNForChannel)      ! SFV < 2.19 km/h  Ref: APEX-doc Eq. 74                                       
      X1=MAX(channelLen-uplandSlopLen/1000.-SFL,0.) ! TCS is the time of concentration fro surface flow in h;                                
      TCC=X1*nManningCoeff/(3.6*channelDep**.66667*SQRT(channelSlp)) ! Eq. 72 TCC is the time of concentration for channel flow                                      
      TCS=.0913*(uplandSlopLen*surfNForChannel)**.8/uplandSteep**.4    ! questional: not match with Eq. 76 missing RFV: dailyrainfall                                   
      ! TSF: time of concentration for shallow channel flow                                                                                       
      ! TC: time of concentration                                                                                                                                                                        
      TC=TCC+TSF+TCS         ! Eq. 71 in APEX-doc  
ELSE                                                            
      ! in APEX-doc using rational eqs, using Manning'n for TCS; but ManningN_Surface for modified rational eqs
      TCS=.0216*(uplandSlopLen*surfNForChannel)**.75/uplandSteep**.375  ! Eq. 68 in APEX-doc on P23                                   
      TCC=1.75*channelLen*nManningCoeff**.75/(areaWshed**.125*channelSlp**.375)  ! Eq. 60  ManningN_Channel: Manning' roughness coefficient                                    
      X4=MIN(uplandSlopLen/360.,TCS)                                                         
      TC=X4+TCC                                                                      
END IF

!======== site determined day length========================================================
      ! Eq. 15c-h in APEX doc on page 14
      ! also see  .f90: calcuate day length and max radiation
! Here day length is calculated by latitude, meanDayLength
XX=YLAT/CLT                                                     
SIN1=SIN(XX)
COS1=COS(XX)
YLTS=SIN1
YLTC=COS1
YTN=TAN(XX)
YTN1=YTN

CH=.4349*ABS(YTN)         
IF(CH>=1.)THEN
      H=0.
ELSE
      H=ACOS(CH)
END IF
meanDayLength=7.72*H   ! day length in h  at the site
HR0=meanDayLength

1 FORMAT(20A4)
2 FORMAT(20F8.2)   
3 FORMAT(20I4) 


END SUBROUTINE Site_Info