SUBROUTINE Model_Setting(FPARM, FPRNT, FTR55, WINDFILE, otherCost, SCRX, coeffChannelGeometry)
USE PARM
USE MisceFun_Module
USE Initialize
USE Print_Module
IMPLICIT NONE 
! added by TXH
! this subroutine read model settings, parameters, and some defalut values
! read EPICFILE.DAT, EPICCON.DAT, and PRAM.DAT

! data dictionary
CHARACTER(LEN = 80):: FPARM, FPRNT, FTR55, WINDFILE
REAL, DIMENSION(2):: otherCost
REAL, DIMENSION(30, 2):: SCRX
REAL, DIMENSION(6):: coeffChannelGeometry
! declaration of local variables 
CHARACTER(LEN = 80):: fileName, AGSM
CHARACTER(LEN = 8)::AGIX
REAL, DIMENSION(6):: XTP
INTEGER, DIMENSION(150):: simuYears2
INTEGER:: i, j, yearIdx, N1
REAL:: X1, X2

! ======================= 1. read EPICCONT file ===================================
! Settings for the simulation
fileName='EPICCONT.DAT'                                            
CALL OPENV(KR(21),fileName,0,KW(MSO))  ! KW(MSO): the error file       
!     LINE 1/2                          40 variables                                                               
READ(KR(21),'(20I4)')numSimuYear, Year0, Month0, Day0, printChoice, weatherVarID, randomCycles, &
      weatherInputMode,leapYrFlag, ET_Method,curveNumFlag, peakRateMethod, erosionMode, &
      operationMode, rootRespirFlag, CN_Method, runoff_Method0, massPestFlag, soilP_runoff_Method,DOY_realtime,&
      numSeedInitialized, Enrich_Method, biomsConver_Method, limeFlag, erosionC_factor, FC_Method, &
      weatherVarForRuns,CO2_Method, N_vol_Method,correctionForRuns, deNitri_Method, NP_uptakeCurve, &
      O2_Method, dirChoice, soilSaturate_Method, LatChoice, P_apply_mode, PAR_Method, percolate_Method, NC_Miner_Method

!     LINE 3/6                           40 variables                                                  
READ(KR(21),1) rainNcon0, CO20,CNO30,saltConIrr,pestDamgFactor, periodHalfHrRain, wetDayCoeff, rainPara, &
      fieldLength,fieldWidth,fieldAngle0, deadCropRes0, windPowerPara, soilDiameter, windErosionFactor,&
      irrTrigger,irrRunoffRatio,irrAnnualMaxVol,irrSingleMinVol, irrSingleMaxVol, fertTrigger0, fertRate, maxAnnualN,&
      dayReduceStress,furrowDikeSafeFactor0,conservFactor0,lagoonVol0,lagoonInput,dayReduceLagoon,liqManureRatio,&
      grazeLimit,herdFrac,layerThickness,waterEroEq,baseStockRate,storLeachFracNO3

fieldAngle=fieldAngle0/CLT     ! convert from degree to rad
  ! default value for some variables
IF(layerThickness<1.E-10) layerThickness=.1   
IF(irrSingleMaxVol<1.E-10) irrSingleMaxVol=1000.                                                      
IF(maxAnnualN<1.E-10) maxAnnualN=200.                                                         
IF(grazeLimit<.01) grazeLimit=.01                                                           
IF(windPowerPara <1.E-10) windPowerPara =.5                                                           
IF(soilDiameter<1.E-10) soilDiameter=500.    
IF(deNitri_Method==0) deNitri_Method =2   
IF(fieldLength<1.E-10) fieldLength=.632                                                           
IF(fieldWidth<1.E-10) fieldWidth=.316                             
IF(periodHalfHrRain<1.E-10) periodHalfHrRain=10.                                                          
IF(wetDayCoeff <1.E-10) wetDayCoeff =.75                                                          
IF(rainPara<1.E-10) rainPara=1.3    
!     READ ECONOMIC DATA                                                             
!     LINE 7                                                                         
READ(KR(21),1)(XTP(i), i=1,6)                                          
CLOSE(KR(21)) ! ------------------ End reading the control file KR(21) -------------------

SELECT CASE(dirChoice+1)
CASE(1)
      IDIR=0                  
CASE(2)
      IDIR=1              
CASE(3)
      IDIR = [(i, i = 1,4)] 
END SELECT 

! ================================ 2. read PARM1102.DAT =============================    
CALL OPENV(KR(2), FPARM, 0, KW(MSO))        ! Equation parameters file  
      !     S_Curve = S CURVE SHAPE PARAMETERS (CONSTANTS EXCEPT FOR EXPERIMENTAL PURPOSES)                                                  
      !     LINE 1/30                                                                      
READ(KR(2),'(2F8.3)')((S_Curve(i, j),j=1,2),i=1,30) !MISCELLANEOUS PARAMETERS(CONSTANTS EXCEPT FOR EXPERIMENTAL PURPOSES                                       
!     LINE 31/40                                                                     
READ(KR(2),1)modelPara           
!     LINE 41                                                                          
READ(KR(2),1)irrCost,limeCost,fuelCost,laborCost,otherCost(1),otherCost(2)
! if economic data are specified in EPICCONT file, then use them to override the PARM file 
IF(XTP(1)>0.)irrCost=XTP(1)   
IF(XTP(2)>0.)limeCost=XTP(2) 
IF(XTP(3)>0.)fuelCost=XTP(3) 
IF(XTP(4)>0.)laborCost=XTP(4)  
IF(XTP(5)>0.)otherCost(1)=XTP(5) 
IF(XTP(6)>0.)otherCost(2)=XTP(6)  
      
!     LINE 42      
READ(KR(2),1)XKN50,XKN30,XKN10,CBVT0               ! What are they?                  
CLOSE(KR(2))                   
! ===================== close KR(2): PARM1102.DAT =============================

! copy S surve   
SCRX = S_Curve
! calcuate the coefficents of the S curve
DO I=1,29                                                                      
    IF(S_Curve(I,1)<1.E-10)CYCLE                                                    
    X1 = split(S_Curve(I,1))   ! !PAY ATTENTION, the split function changes values of S_Curve
    X2 = split(S_Curve(I,2))
    CALL Coeff_Scurve(S_Curve(I,1),S_Curve(I,2),X1,X2)          
END DO

!=========================== 3. read TR55COM.DAT ================================  
CALL OPENV(KR(10),FTR55,IDIR(1),KW(MSO) ) ! data for stochastic runoff estimation 
! coeffTR55  = COEFS OF 7TH DEG POLY IN TR55 PeakRunoff_rate EST                                    
READ(KR(10),'(8F10.6)') coeffTR55                               
READ(KR(10),1) coeffChannelGeometry                                                           
CLOSE(KR(10)) !---------------Close TR55COM.DAT ---------------------------------      

IF(coeffChannelGeometry(3)<1.E-10)coeffChannelGeometry(3)=.0208                                                
IF(coeffChannelGeometry(4)<1.E-10)coeffChannelGeometry(4)=.4                                                   
IF(coeffChannelGeometry(5)<1.E-10)coeffChannelGeometry(5)=.0803                                                
IF(coeffChannelGeometry(6)<1.E-10)coeffChannelGeometry(6)=.6   
                  
!======================= 4. read the AYEAR.DAT file ==============================
! AYEAR not used, but it is connected to the print settings via NFL
fileName='AYEAR.DAT'
CALL OPENV(KR(28),fileName,0,KW(MSO))  
READ(KR(28),*,IOSTAT=NFL)simuYears    ! character 
IF(NFL/=0)REWIND KR(28)               ! if the years are not stored as characters, then read them as integers
READ(KR(28),*,IOSTAT=NFL)simuYears2   ! integer 

MYR=numSimuYear 
!=========================== Allocated Variables ==================================
CALL Allocate_Parms ! Allocate matrix with different sizes 
!============================ 5. read the PRINT1102 file ==========================

CALL OPENV(KR(5),FPRNT,0,KW(MSO) )        ! output control file  
!     LINE 1/5                                                                       
IF(NFL/=0)READ(KR(5),2)(varID(I),I=1,numPrintVar)    ! NFL is questionable?  
!  varID can be selected from 100 parameters, but why the ID can be larger than 100?      
!  stateVarID   = OUTPUT VARIABLE ID NO (ACCUMULATED AND AVERAGE VALUES)    

! ========== collect ID of variables that are not 0 in  varID ===============
DO I=1,numPrintVar  
      IF(varID(I)<=0)EXIT  
      stateVarID(I)= varID(I)            ! stateVarID: flags of output for state variables                                                
END DO  
numPrintVar=I-1                                                                                          
!     LINE 6                                                                         
READ(KR(5),2) (varID(I),I=1,numConVar)  ! 4 concentration variables   

DO I=1,numConVar  
      IF(varID(I)<=0)EXIT  
      conVarID(I)= varID(I)  
END DO 
numConVar=I-1  
!     LINE 7/8                                                                       
READ(KR(5),2) (varID(I),I=1,numMonVar)  ! 40 monthly variables      
DO I=1,numMonVar  
      IF(varID(I)<=0)EXIT 
      monVarID(I)= varID(I) 
END DO  
numMonVar=I-1   
!     LINE 9/10                                                                      
READ(KR(5),2) (varID(I),I=1,numDayVar)  ! 40 daily variables   
DO I=1,numDayVar       
      IF(varID(I)<=0)EXIT  
      dailyVarID(I)= varID(I)  
END DO    
numDayVar=I-1  
!     LINE 11/12                                                                      
READ(KR(5),2) (varID(I),I=1,numAnnuVar) ! 40 annual variables  
DO I=1,numAnnuVar  
      IF(varID(I)<=0)EXIT  
      annualVarID(I)= varID(I)  
END DO    
numAnnuVar=I-1                         
!     LINE 13/14                                                                     
READ(KR(5),2) (varID(I),I=1,numEconVar) ! 40 variables for Flipsim economic analysis                                            
DO I=1,numEconVar  
      IF(varID(I)<=0)EXIT 
      econVarID(I)= varID(I) 
END DO  
numEconVar=I-1  
!  LINE 15/16                                                                     
READ(KR(5),2) KFL       ! KFL:  output file control   total MSO+5 = 38

CLOSE(KR(5)) ! ------------------ END of Print Settings ---------------------------

! ======================== 6. Prepare data for cal field capacity =================== 

IF(FC_Method==7.OR.FC_Method==8)THEN 
      fileName='SOIL35K.DAT'  
      CALL OPENV(KR(29), fileName, IDIR(1), KW(MSO)) 
      READ(KR(29),3)XAV,XDV,XRG,BRNG,NSX  
      NSNN=INT(.655*NSX**.493)  
      EXNN=.767*NSX**.049  
      DO I=1,NSX  
            READ(KR(29),3)(XSP(I,J),J=1,5) 
      END DO  
      CLOSE(KR(29)) 
END IF  

! ========================== 6. write RUN file ==========================      
IF(KFL(MSO+1)/=0)THEN
      fileName='RUN1102.SUM'
      CALL OPENV(KW(MSO+1), fileName, 0, KW(MSO) ) 
      WRITE(KW(MSO+1),4)IYER,IMON,IDAY,Hour,Mint,Sec   ! year, month, day, hour, Min, second  
       
      WRITE(KW(MSO+1),5)(varName(annualVarID(J)),J=1,numAnnuVar)  
END IF  

! ========================== 9. write data to soil file ===================== 
fileName='RTSOIL.DAT'
IF(DOY_realtime>0) CALL OPENV(KW(MSO+3),fileName,0,KW(MSO))  ! DOY_realtime: real time day of year  

! ========================== 10. GIS file output =============================
IF(KFL(NGF)>0)THEN  ! NGF : control ID for gis file   NGF = MSO+5 = 38
      DO I=1,150  
            IF(Year0==simuYears2(I)) THEN
                  yearIdx=I 
                  N1=NGF-1  
                  AGSM=simuYears(I)//"-"//simuYears(I+numSimuYear-1)//".TXT"   
                  OPEN(KW(N1),FILE=AGSM)  
                  IF(KFL(NGF)>=0)THEN
                        DO                                                                                  
                              READ(KW(N1),'(A80)',IOSTAT=NFL)WINDFILE  
                              IF(NFL/=0)EXIT  
                        END DO
                  END IF   
                  IGIS=NGF 
                  DO J=1,numSimuYear 
                        AGIX=simuYears(yearIdx)//".TXT" 
                        yearIdx=yearIdx+1                                                                         
                        OPEN(KW(IGIS),FILE=AGIX)                                                            
                        IF(KFL(NGF)<0)THEN
                              IGIS=IGIS+1                                                                       
                        ELSE
                              DO
                                    READ(KW(IGIS),'(A80)',IOSTAT=NFL)WINDFILE                                         
                                    IF(NFL/=0)EXIT                                                                
                              END DO
                        END IF
                  END DO
            END IF
      END DO 
END IF

1 FORMAT(10F8.2)  
2 FORMAT(20I4)  
3 FORMAT(10F10.0,I10)  
4 FORMAT(/T5,'EEPIC1102',2X,3I4,2X,2(I2,':'),I2/) 
5 FORMAT(13X,'Crop_Name',3X,'YLDG',26(4X,A4),5(4X,'Crop_Name',4X,'YLDG',4X,&               
            'YLDF',4X,'BIOM',4X,' Yield_N',4X,' YLP',4X,' YLC',4X,' FTN',4X,' FTP',&
            4X,'IRGA',4X,'IRDL',4X,'WUEF',4X,' Crop_AvailWater',4X,' CRF',4X,' CQV',4X,&
            ' THU',4X,' potentialHeatUnit',4X,'COST',4X,'COOP',4X,'RETN',4X,'PSTF',4X,'  Water_Stress',&
            4X,'  NS',4X,'  PS',4X,'  KS',4X,'  TS',4X,'  AS',4X,'  SS',4X,&
            '  bulkDensity',4X,'ALSA',4X,' SRT'))
END SUBROUTINE Model_Setting




