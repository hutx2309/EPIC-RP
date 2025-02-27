MODULE Print_Module
USE PARM
USE MisceFun_Module
IMPLICIT NONE
 
CONTAINS

! -------------------------- Timer ---------------------------------
SUBROUTINE Timer(ITR)
      !     EPIC1102
      !     THIS SUBPROGRAM SETS DATE AND TIME FOR OUTPUT AND CALCULATES
      !     ELAPSED TIME.
      IMPLICIT NONE
      ! local variables
      INTEGER:: ITR, IEH, IEM, IES, IEDT, IBDX, I1, IBT, IEX, II, ITS,&
                ITM, ITH 
      INTEGER, DIMENSION(8):: Values
            IF(ITR==0)THEN
                !CALL GETDAT(IYER,IMON,IDAY)
                !CALL GETTIM(Hour,Min,Sec,I100)
                CALL DATE_AND_TIME(VALUES = Values)
                IYER = Values(1)
                IMON = Values(2)
                IDAY = Values(3)
                Hour = Values(5)
                Mint = Values(6)
                Sec = Values(7) 
                RETURN
            END IF
            ! GETTIM and GETDAT is not standard intrinc procedure of fortran 
            !CALL GETTIM(IEH,IEM,IES,I100)  ! program end hour, minutes, sec, and milliseconds
            !CALL GETDAT(IY,MO,Date)
            CALL DATE_AND_TIME(VALUES = Values)
            IY = Values(1)
            MO = Values(2)
            Date = Values(3)
            IEH = Values(5)
            IEM = Values(6)
            IES = Values(7)
       
            CALL Leap_Yr_Check(IY,leapYr,leapYrFlag)
            IEDT = Cal_DOY(MO,Date, leapYr)
            CALL Leap_Yr_Check(IYER,leapYr,leapYrFlag)
            IBDX = Cal_DOY(IMON,IDAY, leapYr)
            I1=86400*((IY-IYER)*(366-leapYr)+IEDT-IBDX)
            IBT=Hour*3600+Mint*60+Sec
            IEX=IEH*3600+IEM*60+IES
            WRITE(KW(1),5000)Hour,Mint,Sec
            WRITE(KW(1),5100)IEH,IEM,IES
            WRITE(KW(1),5200)
            II=I1+IEX-IBT
            ITS=MOD(II,60)
            II=(II-ITS)/60
            ITM=MOD(II,60)
            ITH=(II-ITM)/60
            WRITE(KW(1),5320)ITH,ITM,ITS
            WRITE(*,5320)ITH,ITM,ITS
            RETURN
       5000 FORMAT(10X,'BEGINNING TIME: ',I2,2(':',I2),'.',I2)
       5100 FORMAT(10X,'ENDING    TIME: ',I2,2(':',I2),'.',I2)
       5200 FORMAT(10X,'----------------------------')
       5320 FORMAT(10X,'TOTAL RUN TIME: ',I2,2(':',I2),'.',I2)
       
END SUBROUTINE Timer

! ----------------------------------- 3. Open files for writing ----------------------------
SUBROUTINE Open_OutputFiles 
      !     EPIC1102
      !     THIS SUBPROGRAM OPENS FILES.
      IMPLICIT NONE
 
      ! local variables
      INTEGER:: I
      CHARACTER(LEN = 80):: outpath = "G:\EEPIC_PilotStudy\LM\output\EEPIC_"
      CHARACTER(4), DIMENSION(32)::AXT = [".OUT",".ACM",".SUM",".DHY",".DPS",".MFS",".MPS",".ANN",&
                             ".SOT",".DTP",".MCM",".DCS",".SCO",".ACN",".DCN",".SCN",".DGN",&
                             ".DWT",".ACY",".ACO",".DSL",".MWC",".ABR",".ATG",".MSW",".APS",&
                             ".DWC",".DHS",".DGZ",".DNC",".ASL",".DDN"]
 
      !     SELECT OUTPUT FILES--KFL=0(NO OUTPUT); KFL>0(GIVES OUTPUT FOR                 
      !     SELECTED FILE)                                                                 
      !     1  OUT = STANDARD OUTPUT                                                       
      !     2  ACM = ANNUAL CROPMAN                                                        
      !     3  SUM = AVE ANNUAL SUMMARY                                                    
      !     4  DHY = DAILY HYDROLOGY                                                       
      !     5  DPS = DAILY PESTICIDE                                                       
      !     6  MFS = MONTHLY FLIPSIM                                                       
      !     7  MPS = MONTHLY PESTICIDE                                                     
      !     8  ANN = ANNUAL                                                                
      !     9  SOT = ENDING SOIL TABLE                                                     
      !    10  DTP = DAILY SOIL TEMPERATURE                                                
      !    11  MCM = MONTHLY CROPMAN                                                       
      !    12  DCS = DAILY CROP STRESS                                                     
      !    13  SCO = SUMMARY OPERATION COST                                                
      !    14  ACN = ANNUAL SOIL ORGANIC C & N TABLE                                       
      !    15  DCN = DAILY SOIL ORGANIC C & N TABLE                                        
      !    16  SCN = SUMMARY SOIL ORGANIC C & N TABLE                                      
      !    17  DGN = DAILY GENERAL OUTPUT                                                  
      !    18  DWT = DAILY SOIL WATER IN CONTROL SECTION AND .5M SOIL T                    
      !    19  ACY = ANNUAL CROP YIELD                                                     
      !    20  ACO = ANNUAL COST                                                           
      !    21  DSL = DAILY SOIL TABLE.                                                     
      !    22  MWC = MONTHLY WATER CYCLE + N CYCLE                                         
      !    23  ABR = ANNUAL BIOMASS ROOT WEIGHT                                            
      !    24  ATG = ANNUAL TREE GROWTH.                                                   
      !    25  MSW = MONTHLY OUTPUT TO SWAT.                                               
      !    26  APS = ANNUAL PESTICIDE                                                      
      !    27  DWC = DAILY WATER CYCLE                                                     
      !    28  DHS = DAILY HYDROLOGY/SOIL                                                  
      !    29  DGZ = DAILY GRAZING FILE
      !    30  DNC = DAILY NITROGEN/CARBON CESAR IZAURRALDE
      !    31  ASL = ANNUAL SOIL TABLE
      !    32  DDN = DAILY DENITRIFICATION
      !   MSO  ERX = ERROR FILE                                                            
      !  MSO+1 RUN1102.SUM = AVE ANNUAL SUMMARY FILE FOR ALL SIMULATIONS IN A              
      !        BATCH.                                                                      
      !  MSO+2     = RTCROP.DAT                                                            
      !  MSO+3     = RTSOIL.DAT                                                            
      !  MSO+4 SGI = SUMMARY GIS FILE                                                      
      !  MSO+5 ANNUAL FILES FOR GIS   

    DO I=1,MSO-1
        IF(AXT(I)/="    ".AND.KFL(I)>0)OPEN(KW(I),FILE= TRIM(ADJUSTL(outpath))//TRIM(ADJUSTL(siteName))//AXT(I))
    END DO
END SUBROUTINE Open_OutputFiles

! --------------------------- 2. Write head info for each output file ---------------
SUBROUTINE Write_HeadInfo1(soilName, WPM1FILE, WINDFILE, XCC, soilHydGroup, fracFC, groundWaterResidenceTime0, &
                         splitedSoilLayer, soilCondCode, soilGroup, minThickSplit, fracSOC, fracHumus, YLAZ)
      IMPLICIT NONE
      ! arguments list
      CHARACTER(LEN = 80):: soilName, WPM1FILE, WINDFILE
      REAL:: XCC, soilHydGroup, fracFC, groundWaterResidenceTime0, splitedSoilLayer, soilCondCode, &
            soilGroup, minThickSplit, fracSOC, fracHumus, YLAZ
      ! local variables
      CHARACTER(LEN = 80):: A1, A2      
      INTEGER:: I, J
      REAL:: X1
      DO I=2,MSO           ! MSO = 33     number if files 
            IF(KFL(I)==0) CYCLE 
            ! .SOT file : the first three lines of a soil table after simulation 
            IF(I==9)THEN  
                  XCC=1. 
                  X1=0.  
                  WRITE(KW(9),'(T10,2A20,2X,3I4)')soilName,soilOrder,IYER,IMON,IDAY  
                  WRITE(KW(9),'(10F8.2)')soilAlbedo,soilHydGroup,fracFC,&
                        waterTableMinDep,waterTableMaxDep,waterTableHigt,&
                        groundWaterStor,groundWaterMaxStor,groundWaterResidenceTime0,returnFlowFrac,&  
                        splitedSoilLayer,soilCondCode,X1,soilGroup,minThickMaxLayer,minThickProfile,&
                        minThickSplit,fracSOC,fracHumus,XCC 
                  CYCLE
            END IF  
            WRITE(KW(I), 1)IYER,IMON,IDAY,Hour,Mint,Sec  
            WRITE(KW(I),'(T10,A8)')siteName  
            WRITE(KW(I),'(5X,3I4)')IRUN,IRO0,randomCycles  
            WRITE(KW(I),2)SITEFILE   
            WRITE(KW(I),2)WPM1FILE  
            WRITE(KW(I),2)WINDFILE  
            WRITE(KW(I),2)SOILFILE 
            WRITE(KW(I),2)OPSCFILE  
            ! .ACN file: Annual soil organic C and N table
            IF(I==14) CALL Print_OrgNCvar(Year0,1,1 )  
            ! .DCN file : Daily soil organic C and N table
            IF(I==15) CALL Print_OrgNC(1 )
      END DO

      IF(KFL(2)>0) WRITE(KW(2),9)varName(4),varName(10),varName(11),varName(14),varName(16),varName(17),&
      varName(29),varName(33),varName(42),varName(48),varName(47),varName(50),varName(51),&
      varName(52),varName(49),varName(43),varName(44),varName(45),varName(46),varName(56),&
      varName(54),varName(55),varName(57),varName(66)   

      ! Annual Crop Yield file   .ACY
      IF(KFL(19)>0)WRITE(KW(19),11) 

      ! **************** Add by TXH for testing **************************************     
      ! OPEN(KW(50), file = "I:\Diag_EEPIC\EEPIC-ph2\Yield_RP_Y6.txt")     
      ! WRITE(KW(50), "(1X,' YR  RT#',1X,'Crop_Name',3X,'YLDG',4X,'  HI', 4X,'tBio')")       
      ! ******************************************************************************

      IF(KFL(8)>0)WRITE(KW(8),10)(varName(annualVarID(J)),J=1,numAnnuVar)  
      IF(KFL(13)>0)WRITE(KW(13), 4)  
      IF(KFL(17)>0)WRITE(KW(17),12)(varName(dailyVarID(J)),J=1,numDayVar),'ZNH3','totNO3_N',&
            'NO31','PRK1','LN31',' Albedo',' HUI','AJHI',' LAI','BIOM','  Tot_RootWeight',&
            'AboveGround_Biomass','  HI','YLDX','YLDF','PotentialN_Demand',' NPP',' NEE',&
            (HEDS(I),I=23,28) 
      ! print .DWT :daily soil water in control section and 0.5m soil 
      IF(KFL(18)>0)THEN
            IF(weatherVarID>0)THEN
                  A1='INPUT WEATHER'
                  A2=FWTH
            ELSE
                  A1='GENERATED WEATHER'
                  A2=' '
            END IF
            WRITE(KW(18),13)siteName,YLAT,YLAZ,uplandSteep,ZCS,A1,A2  
      END IF

      IF(KFL(4)>0)WRITE(KW(4),3)       ! .DHY : daily hydrology   
      IF(KFL(20)>0)WRITE(KW(20),4)     ! .ACO : annual cost    
      IF(KFL(24)>0)WRITE(KW(24),5)     ! .ATG : annual tree growth   
      IF(KFL(26)>0)WRITE(KW(26),6)     ! .APS : annual pesticide  
      IF(KFL(29)>0)WRITE(KW(29),7)     ! .DGZ : daily grazing file  
      IF(KFL(23)>0)THEN                ! .ABR : annual biomass root weight  
            WRITE(KW(23),8)(SID(LORG(Layer_ID(J))),J=1,Actual_SoilLayers),SID(16)                               
            WRITE(KW(23),'(10X,A15,5X,16F8.2)')'DEPTH(m)',(Z(Layer_ID(I)),I=1,Actual_SoilLayers)                                    
      END IF     

    1 FORMAT(/T5,'EEPIC1102_Ph',2X,3I4,2X,2(I2,':'),I2/) 
    2 FORMAT(T10,A12)            
    3 FORMAT(/3X,'DATE',8X,'CN',7X,'RAIN',7X,'Q',9X,'TC',8X,'PeakRunoff_rate',7X,&                
            'Rainfall_Duration',6X,'ALTC',7X,'Alpha_05L'/4X,'Y M D',T25,'(mm)',6X,'(mm)',7X,'(H)',&            
            5X,'(mm/h)',6X,'(H)')    
    4 FORMAT(T61,'COTL',6X,'COOP',6X,'MTCO',6X,'MASS',6X,'FUEL'/4X,'Y',&             
            3X,'M',3X,'D',5X,'stableMineralP',14X,'CROP',2X,'Mon_Tillage#  HC',2X,'EQ  TR',2X,&                 
            '|----------($/ha)-----------|',2X,'(kg/ha)',4X,'(L/HA)') 
    5 FORMAT(3X,'Y',3X,'Y#',1X,'CROP',4X,' YLD',4X,'BIOM',4X,' RWT',4X,&             
            ' LAI',4X,' standCropResi') 
    6 FORMAT(2X,'YR   YR#',5X,'Q',5X,'subsurfLaterFlow',5X,'SEP_Flow',4X,'QDRN',7X,'Y',5X,&            
            'YOC',5X,'pestName',11X,'Crop_Name',4X,'PAPL',4X,'PSRO',4X,'Pest_Leach',4X,'PSSF',&
            4X,'PSED',4X,'PDGF',4X,'PDGS',4X,'PDRN',4X,'CMX4D'/12X,&
            '|------------(mm)------------|  (t/ha)',1X,'(kg/ha)',17X,&
            '|-----------------------------(g/ha)--------------------------------|',3X,'(ppb)') 
    7 FORMAT(4X,'Y',3X,'M',3X,'D',1X,'OPERATION',2X,'CROP',2X,'BIOMt/ha'&            
            ,3X,'RWTt/ha',6X,'LAI',4X,'STLt/ha',2X,'AGPMt/ha',5X,'overRideHI',4X,&               
            'YLDt/ha',2X,'YLSDt/ha',6X,'fracHeatUnit')  
    8 FORMAT(T55,'RWT SOIL LAYER #'/3X,'Y',3X,'Y#','   M   D',1X,'CROP',&
            4X,'BIOM',16(4X,A4))  
    9 FORMAT(1X,' YR  RT#',24(4X,A4),4X,'plawDepSOC',4X,' totSOC',4X,'APBC',4X,&              
            ' totLabileP',4X,'totNO3_N') 
   10 FORMAT(1X,'RUN   YR',4X,'AP15',4X,'PMTE',41(3X,A4,1X)) 
   11 FORMAT(1X,' YR  RT#',1X,'Crop_Name',3X,'YLDG',4X,'YLDF',4X,'WaterFrac_Yield',4X, &
            '  HI',4X,'BIOM',4X,'  Tot_RootWeight',4X,' Yield_N',4X,' YLP',4X,' YLC',4X,' FTN',&
            4X,' FTP',4X,' FTK',4X,'IRGA',4X,'IRDL',4X,'WUEF',4X,'GSET',4X,&
            ' Crop_AvailWater',4X,' CRF',4X,' CQV',4X,'COST',4X,'COOP',4X,'RYLG',4X,'RYLF',&
            4X,'PSTF',4X,'  Water_Stress',4X,'  NS',4X,'  PS',4X,'  KS',4X,'  TS',4X,&
            '  AS',4X,'  SS',4X,'PPOP',5X,'IPLD',5X,'IGMD',5X,'IHVD ') 
   12 FORMAT(4X,'Y   M   D',5X,'PDSW',50(6X,A4))              ! PDSW: plawDepthSoilWater  
   13 FORMAT(4X,'RUN #=',A8,2X,'LAT=',F7.2,' deg',2X,'EQUIV LAT=',F7.2,&
            ' deg',2X,'SLOPE=',F7.4,' m/m'/4X,'CONTROL SECTION DEPTHSm =',&
            3F6.3,2X,2A20/9X,'DATE',T23,'CONTROL SECTION'/2X,'Y#',4X,'Y   M   D',&
            5X,'SW1mm',5X,'SW2mm',2X,'PAW10m/m',2X,'PAW20m/m',2X,'PAW50m/m',1X,&
            'PAW100m/m',6X,'T10c',6X,'T20c',6X,'T50c',5X,'T100c',5X,'SW10%')             

END SUBROUTINE Write_HeadInfo1

! --------------------------- prepare daily/montly output files -------------------- 
SUBROUTINE Write_HeadInfo2
      IMPLICIT NONE
      ! arguments list
      ! local variables
      INTEGER:: I, J, K
      CHARACTER(LEN = 4), DIMENSION(15, 10):: HD28 = RESHAPE(&
            ['  Z1','  Z2','  Z3','  Z4','  Z5','  Z6','  Z7','  Z8','  Z9',' Z10',' Z11',' Z12',' Z13',' Z14',' Z15',&
             ' SW1',' SW2',' SW3',' SW4',' SW5',' SW6',' SW7',' SW8',' SW9','SW10','SW11','SW12','SW13','SW14','SW15',&
             ' WU1',' WU2',' WU3',' WU4',' WU5',' WU6',' WU7',' WU8',' WU9','WU10','WU11','WU12','WU13','WU14','WU15',&
             ' EV1',' EV2',' EV3',' EV4',' EV5',' EV6',' EV7',' EV8',' EV9','EV10','EV11','EV12','EV13','EV14','EV15',&
             ' PK1',' PK2',' PK3',' PK4',' PK5',' PK6',' PK7',' PK8',' PK9','PK10','PK11','PK12','PK13','PK14','PK15',&
             ' SF1',' SF2',' SF3',' SF4',' SF5',' SF6',' SF7',' SF8',' SF9','SF10','SF11','SF12','SF13','SF14','SF15',&
             ' N31',' N32',' N33',' N34',' N35',' N36',' N37',' N38',' N39','N310','N311','N312','N313','N314','N315',&
             ' UN1',' UN2',' UN3',' UN4',' UN5',' UN6',' UN7',' UN8',' UN9','UN10','UN11','UN12','UN13','UN14','UN15',&
             ' LN1',' LN2',' LN3',' LN4',' LN5',' LN6',' LN7',' LN8',' LN9','LN10','LN11','LN12','LN13','LN14','LN15',&
             '  T1','  T2','  T3','  T4','  T5','  T6','  T7','  T8','  T9',' T10',' T11',' T12',' T13',' T14',' T15'], [15, 10])
      CHARACTER(LEN = 7), DIMENSION(30):: HD30 = ['     ZC','    VWC','    AFP','   HKPO','   DPRO','   HKPC',&
             '   DPRC','   HKPN','   DPRN','   SMEA','   SMES','   WNH3','   WNO3','   WNO2','   WO2L','   WO2G', &
             'DO2CONS','  SSFO2','    VO2','  WCO2L','  WCO2G','DCO2GEN',' SSFCO2','   VCO2','  WN2OL','  WN2OG', &
             ' SSFN2O','   VN2O','  DN2OG','   DN2G']
      CHARACTER(LEN = 7), DIMENSION(15):: HDG = ['   GFO2','  GFCO2','  GFN2O','  DFO2S',' DBFO2B','  DFO2T',&
             '    QO2',' DFCO2S',' DFCO2B',' DFCO2T','   QCO2',' DFN2OS',' DFN2OB', ' DFN2OT','   QN2O']
      CHARACTER(LEN = 2), DIMENSION(30):: HDN = ['1 ','2 ','3 ','4 ','5 ','6 ','7 ','8 ','9 ','10','11',&
             '12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30']
             
      !////////////////////////////////////////////////////////////////////////////////
      IF(KFL(5)>0)WRITE(KW(5),1)HEDP       ! outputs of P related variables                                            
      IF(KFL(10)>0)WRITE(KW(10),2)         ! Temperature                                              
      IF(KFL(11)>0)WRITE(KW(11),3)HEDS(6),varName(4),varName(11),varName(14),varName(17),varName(16)! Hydrology  
      IF(KFL(12)>0)WRITE(KW(12),4)                                                                               
      IF(KFL(6)>0)WRITE(KW(6),5)(varName(econVarID(J)),J=1,numEconVar),(HEDS(I),I=4,6)                                                                         
      IF(KFL(22)>0)WRITE(KW(22),6)varName(4),varName(10),varName(11),varName(13),&                      
                   varName(14),(varName(I),I=16,20),(HEDS(I),I=6,8),'RNO3',(varName(I),I=43,46),&             
                   varName(49),varName(89),varName(52),varName(85),varName(50),(varName(I),I=59,61),&
                   'PotentialN_Demand',' Yield_N','Crop_Name',' YLD','TOTN'                                                    
      IF(KFL(27)>0)WRITE(KW(27),7)(varName(dailyVarID(K)),K=1,numDayVar),(HEDS(monVarID(K)),K=1,numMonVar)                                       
      IF(KFL(28)>0)WRITE(KW(28),8)varName(4),varName(10),varName(11),(varName(I),I=13,20),(HEDS(I),I=6,8),&
                   ((HD28(I,J),I=1,Actual_SoilLayers),J=1,10)
      IF(KFL(30)>0)THEN
          WRITE(KW(30),'(T15,3(A,E16.6))')'XKN1=',XKN1,' XKN3=',XKN3, ' XKN5=',XKN5                                            
          WRITE(KW(30),9)varName(4),((HD30(J),HDN(I),I=1,layersEqualThick),J=2,3),&
                ((HD30(J),HDN(I),I=1,layersEqualThick),J=12,17),(HDG(I),'  ',I=1,4),&
                ((HD30(J),HDN(I),I=1,layersEqualThick),J=18,22),(HDG(I),'  ',I=5,8),&
                ((HD30(J),HDN(I),I=1,layersEqualThick),J=23,26),(HDG(I),'  ',I=9,12),&
                ((HD30(J),HDN(I),I=1,layersEqualThick),J=27,30)
      END IF   
    1 FORMAT(4X,'Y M D  RT#  pestName',12X,10(4X,A4,2X),5X,'Q',8X,'SubSurf_LateralFlow',6X,&             
            'Percolation_Flow',4X,'ROCONC') 
    2 FORMAT(T27,'|------------------------------------TEMP(C) ----------------------------------|', &
            /6X,'DATE',8X,'DAMP',15X, &
            '|________________________________@ CENTER OF SOIL LAYERS___________________________|', &
            /T5,'Y   M   D',3X,'DEPTH(m)',5X,'SURF',7X,'1',9X,'2',9X,'3',9X,'4',9X,'5',9X,'6',9X,'7',9X,'8',9X,'9',8X,'10')  
    3 FORMAT(4X,'Y M RT#  Crop_Name',6X,'waterStress',6X,'NS',6X,'PS',6X,'KS',6X,'TS', 6X,'AS',6X,'SS',6(4X,A4))
    4 FORMAT(4X,'Y   M   D  RT#',4X,'Crop_Name',6X,' HUI',6X,'AJHI',6X,' LAI',6X,'BIOM',6X,'  Tot_RootWeight',6X,&
             ' AboveGround_Biomass',6X,'  HI',6X,'YLDX',6X,'YLDF',6X, 'PotentialN_Demand',6X,' NPP',6X,' NEE', 8X,&
             'waterStress',8X,'NS',8X,'PS',8X,'KS',8X,'TS',8X,'AS',8X,'SS',7X,'RW1',7X,'RW2',7X,'RW3',7X,'RW4',7X,&
             'RW5',7X,'RW6',7X,'RW7',7X,'RW8',7X,'RW9',6X,'RW10',6X,'RW11',6X,'RW12',6X,'RW13',6X,'RW14',6X,'RW15')     
    5 FORMAT(4X,'Y M RT#',43(4X,A4))  
    6 FORMAT(4X,'Y   M  ',32(A4,4X))
    7 FORMAT(4X,'Y   M   D  ',140(6X,A4,2X))  
    8 FORMAT(4X,'Y',3X,'M',3X,'D',8X,'SW15',8X,'SW30',7X,'NO315',7X, 'NO330',7X,'NH315',7X,'NH330',200(8X,A4))  
    9 FORMAT(4X,'Y',3X,'M',3X,'D',5X,A4,1X,1000(A7,A2,4X))            
END SUBROUTINE Write_HeadInfo2

! ---------------------------2. Print general informaton for simulation ------------------------
SUBROUTINE Print_GeneralInfo(channelLen, channelSlp, nManningCoeff, channelDep, siteAzimuth,&
                             surfNForChannel, RXM, peakRateFactor, snowCont0, otherCost, KRX)
      IMPLICIT NONE
      ! arguments list
      REAL, INTENT(IN):: channelLen,channelSlp, nManningCoeff, channelDep,siteAzimuth,surfNForChannel, &
                         RXM, peakRateFactor, snowCont0 
      REAL, DIMENSION(2), INTENT(IN):: otherCost
      INTEGER, INTENT(OUT):: KRX
      ! local variables
      CHARACTER(LEN = 80):: fileName
      CHARACTER(LEN = 2), DIMENSION(4)::RFPT= [' 1','1A',' 2',' 3']   
      ! //////////////////////////////////////////////////////////////////////////////////////
      WRITE(KW(1),'(//1X,A,A/)')'____________________GENERAL INFORMATION______________________' 
      KRX=KR(16)  
      IF(DOY_realtime>0)THEN 
            WRITE(KW(1),'(T10,A)')'REAL TIME SIMULATION MODE' ! What does "real time simulation mode" mean?  
            READ(KW(MSO+3),'(I4)')ISTP 
            IF(ISTP>0)THEN 
                  KRX=KR(22)
                  fileName='RTOPSC.DAT' 
                  CALL OPENV(KR(22),fileName,0,KW(MSO) )  
            END IF  
      ELSE
            WRITE(KW(1),'(T10,A)')'NORMAL SIMULATION MODE'  
      END IF  

      IF(modelPara(50)>0.)THEN  
            WRITE(KW(1),'(T10,A)')'DYNAMIC TECHNOLOGY'  
      ELSE  
            WRITE(KW(1),'(T10,A)')'STATIC TECHNOLOGY'  
      END IF  
      
      WRITE(KW(1),1)numSimuYear,Year0,Month0,Day0 

      IF(leapYrFlag>0)THEN  
            WRITE(KW(1),'(T10,A)')'LEAP YEAR IGNORED'                                         
      ELSE 
            WRITE(KW(1),'(T10,A)')'LEAP YEAR CONSIDERED'  
      END IF  

      WRITE(KW(1),2)areaWshed, YLAT, XLOG, ELEV, channelLen, channelSlp, nManningCoeff, channelDep, &
                  uplandSlopLen, uplandSteep, siteAzimuth, surfNForChannel   
      
      waterEroModel= INT(waterEroEq+1.1)   ! WaterEro_Eq: water erosion driving Eq  
      SELECT CASE(waterEroModel)                                                              
            CASE(1)                                                                      
                  NDVSS=34                    ! MUST   
            CASE(3)                                                                      
                  NDVSS=30                    ! USLE   
            CASE(4)                                                                      
                  NDVSS=33                    ! MUSS  
            CASE(5)                                                                      
                  NDVSS=31                    ! MUSL  
            CASE(7)                                                                      
                  NDVSS=36                    ! RUSL  
            CASE(8)
                  NDVSS=35                    ! RUS2
            CASE DEFAULT
                  NDVSS=31                    ! MUSL  
      END SELECT    
      
      WRITE(KW(1),'(T10,A,A4)')'WATER EROSION FACTORS--DRIVING EQ = ', varName(NDVSS) 
      WRITE(KW(1),3)SL,RXM,RLF,RSF,TC                                              
      WRITE(KW(1),'(T10,A)')'DAILY RUNOFF ESTIMATION'    
      
      ! FOR GREEN & AMPT ESTIMATE OF Q, RF EXP DST, PEAK RF RATE SIMULATED.
      runoff_Method = runoff_Method0+1   
      SELECT CASE(runoff_Method) 
            CASE(1)  
                  WRITE(KW(1),'(T15,A)')'NRCS CURVE NUMBER EQ'  
            CASE(2)                                                                      
                  WRITE(KW(1),'(T15,A/T15,A)')'GREEN & AMPT EQ','RF EXP DST--PEAK RF RATE SIM'  
            CASE(3)  
                  WRITE(KW(1),'(T15,A/T15,A)')'GREEN & AMPT EQ','RF EXP DST--PEAK RF RATE INPUT' 
            CASE(4) 
                  WRITE(KW(1),'(T15,A/T15,A)')'GREEN & AMPT EQ','RF UNIF DST--PEAK RF RATE printInterval' 
      END SELECT 

      SELECT CASE(CN_Method+1)  
            CASE(1)    
                  WRITE(KW(1),'(T15,A)')'VARIABLE CN DEPTH/SOIL-WATER WEIGHT' 
            CASE(2) 
                  WRITE(KW(1),'(T15,A)')'VARIABLE CN NO DP/totSoilWater WEIGHT'  
            CASE(3)  
                  WRITE(KW(1),'(T15,A)')'VARIABLE CN LINEAR NO DP/totSoilWater WEIGHT'  
            CASE(4)  
                  WRITE(KW(1),'(T15,A)')'CONSTANT CN'  
            CASE(5) 
                  WRITE(KW(1),'(T15,A)')'VARIABLE CN RETN PAR INDEX NO DP/totSoilWater WEIGHT'  
      END SELECT

      IF(curveNumFlag==0)THEN  
            WRITE(KW(1),'(T15,A)')'DAILY CN--STOCHASTIC'  
      ELSE                                                                                
            WRITE(KW(1),'(T15,A)')'DAILY CN--DETERMINISTIC' 
      END IF  

      IF(peakRateMethod >0)THEN  
            WRITE(KW(1),'(T10,A,A2)')'PEAK RATE EST WITH TR55--RF TYPE =', &  
                RFPT(peakRateMethod) 
      ELSE 
            WRITE(KW(1),'(T10,A)')'PEAK RATE EST WITH MOD RATIONAL EQ'  
      END IF

      IF(percolate_Method==0)THEN
            WRITE(KW(1),'(T10,A)')'SEP--HPERC'
      ELSE
            WRITE(KW(1),'(T10,A)')'SEP--HPERC1 (4mm SLUG FLOW)'
      END IF   

      IF(Enrich_Method>0)THEN                                                                 
            WRITE(KW(1),'(T10,A)')'GLEAMS ENRICHMENT RATIO'   
      ELSE  
            WRITE(KW(1),'(T10,A)')'EPIC ENRICHMENT RATIO'  
      END IF

      IF(biomsConver_Method>0)THEN  
            WRITE(KW(1),'(T10,A)')'WATER USE-BIOMASS CONVERSION'   
      ELSE 
            WRITE(KW(1),'(T10,A)')'RADIATION-BIOMASS CONVERSION'
      END IF

      IF(PAR_Method==0)THEN
            WRITE(KW(1),'(T10,A)')'PAR DRIVEN BY LAI'
      ELSE
            WRITE(KW(1),'(T10,A)')'PAR DRIVEN BY REMOTE SENSING EVI'
      END IF

      IF(N_vol_Method==0)THEN
            WRITE(KW(1),'(T10,A)')'ORIGINAL EPIC NITVOL EQS'
      ELSE
            WRITE(KW(1),'(T10,A)')'IZAURRALDE REVISED NITVOL EQS'
      END IF 
 
      SELECT CASE(deNitri_Method)
            CASE(1)
                  WRITE(KW(1),'(T10,A)')'EPIC DNIT'
            CASE(2)              
                  WRITE(KW(1),'(T10,A)')'KEMANIAN DNIT'
            CASE(3)
                  WRITE(KW(1),'(T10,A,F5.2,A)')'IZAURRALDE DNIT LayerThickness=',layerThickness,' m'
                  WRITE(KW(1),'(T15,3(A,E16.6),A)')'XKN1=',XKN1,' XKN3=',XKN3,&
                        ' XKN5=', XKN5,' ORIGINAL DW'   
                  IF(ABS(rootRespirFlag) < 1.0E-6)THEN
                        WRITE(KW(1),'(T10,A)')'ROOT RESPIRATION CALCULATED IN NDNITCI'
                  ELSE    
                        WRITE(KW(1),'(T10,A)')'ROOT RESPIRATION NOT CONSIDERED IN NDNITCI'  
                  END IF
            CASE(4)
                  WRITE(KW(1),'(T10,A,F5.2,A)')'IZAURRALDE DNIT LayerThickness=',layerThickness,' m'
                  WRITE(KW(1),'(T15,3(A,E16.6),A)')'XKN1=',XKN1,' XKN3=',XKN3,' XKN5=',XKN5,' NEW DW'
                  IF(rootRespirFlag == 0)THEN
                        WRITE(KW(1),'(T10,A)')'ROOT RESPIRATION CALCULATED IN NDNITCI'
                  ELSE  
                        WRITE(KW(1),'(T10,A)')'ROOT RESPIRATION NOT CONSIDERED IN NDNITCI'   
                  END IF
      END SELECT

      IF(NP_uptakeCurve>0)THEN  
            WRITE(KW(1),'(T10,A)')'N & P UPTAKE CONC S CURVE'  
      ELSE  
            WRITE(KW(1),'(T10,A)')'N & P UPTAKE CONC SMITH CURVE'   
      END IF

      IF(NC_Miner_Method==0)THEN
            WRITE(KW(1),'(T10,A)')'PHOENIX N-C MODEL'
      ELSE    
            WRITE(KW(1),'(T10,A)')'CENTURY N-C MODEL'
      END IF   

      SELECT CASE(O2_Method)   
            CASE(1)
                  WRITE(KW(1),'(T10,A)')'O2=F(C/clayFrac)'  
            CASE(2)    
                  WRITE(KW(1),'(T10,A)')'O2=gasO2Con(ISL)/gasO2Con(LD1)'  
            CASE DEFAULT                                                                                
                  WRITE(KW(1),'(T10,A)')'EPIC O2=F(Z)'  
      END SELECT   

      WRITE(KW(1),4)peakRateFactor,snowCont0,rainNcon0,CNO3I,saltConIrr  

      IF(massPestFlag>0)THEN                                                                 
            WRITE(KW(1),'(T10,A)')'NUTRIENT & PESTICIDE OUTPUT (MASS & CONC)'  
      ELSE                                                                                
            WRITE(KW(1),'(T10,A)')'NUTRIENT & PESTICIDE OUTPUT (MASS)'
      END IF      
 
      IF(soilP_runoff_Method>0)THEN    
            WRITE(KW(1),'(T10,A)')'MODIFIED NONLINEAR EQ SOL P RUNOFF' 
            ELSE                                                                                
            WRITE(KW(1),'(T10,A)')'GLEAMS PESTICIDE EQ SOL P RUNOFF'   
      END IF  

      IF(erosionC_factor==0)THEN                                                                 
            WRITE(KW(1),'(T10,A)')'RUSLE C FACTOR USED FOR ALL EROS EQS'   
      ELSE                                                                                
            WRITE(KW(1),'(T10,A)')'EPIC C FACTOR USED EXCEPT FOR RUSLE'   
      END IF 

      WRITE(KW(1),5)limeCost,irrCost,fuelCost,laborCost,otherCost(1),otherCost(2) 

    1 FORMAT(T10,'SIMULATION DURATION = ',I4,' Y'/T10,'BEGINNING DATE = ',3I4)  
    2 FORMAT(T10,'DRAINAGE AREA = ',F9.2,' ha'/T10,'LATITUDE = ',F7.2,' deg'/T10,'LONGITUDE = ', &
            F7.2,' deg'/T10,'ELEVATION = ',F7.1,' m'/T10,'CHANNEL'/T15,'LENGTH = ',F6.2,' km',5X,&
            'GRAD = ',F7.4,' m/m',5X,'MANNINGS N = ',F6.3,5X,'DEPTH = ',F6.3,' m'/T10, &
            'LAND SLOPE'/T15,'LENGTH = ',F6.1,' m',5X,'GRAD = ',F7.4,' m/m',5X,'AZIMUTH = ',&
            F4.0,' deg',5X,'MANNINGS N = ',F6.3) 
    3 FORMAT(T15,'LS= ',F6.3/T10,'RUSLE XM = ',F6.3/T10,'RUSLE SLP LG FAC = ',F6.3/T10,&
            'RUSLE SLP FAC = ',F6.3/T10,'TIME OF FLOW CONC = ',F5.2,' h')   
    4 FORMAT(T10,'PEAK RATE-Rainfall_EnergyFact ADJ FACTOR = ',F6.3/T10, &
            'INITIAL WATER CONTENT OF SNOW = ',F5.1,' mm'/T10,'AVE N CONC IN RAINFALL = ',F6.2,&            
            ' ppm'/T10,'CONC OF NO3 IN IRRIGATION WATER = ',F6.0,' ppm'/T10,&              
            'CONC OF SALT IN IRRIGATION WATER = ',F8.0,' ppm') 
    5 FORMAT(T10,'COSTS'/T15,'LIME = ',F5.2,' $/t'/T15,'IRR WATER = ',&              
            F6.2,' $/mm'/T15,'FUEL PRICE = ',F5.2,' $/l'/T15,'WAGE PRICE =',&              
            F5.2,' $/h'/T15,'MISCEL COST 1 = ',F7.2,' $/ha'/T15,'MISCEL COST 2 = ',F7.2,' $/ha')     
END SUBROUTINE Print_GeneralInfo

! ----------------------------------3. Print random numbers --------------------------------
SUBROUTINE Print_RandomNums
      IMPLICIT NONE
      ! arguments list

      ! local variables
      INTEGER:: I, II, KK, J
      REAL:: RN, XX
      ! RANDOM NUMBER GENERATOR ID NUMBERS  
      ! Random_ID = 1 DETERMINES WET AND DRY DAYS   
      !           = 2 RELATES WEATHER VARIABLES TO RAIN  
      !           = 3 RAINFALL AMOUNT     
      !           = 4 RAINFALL ENERGY (Rainfall_EnergyFact)- PEAK RUNOFF RATE (PeakRunoff_rate)   
      !           = 5 WIND SPEED   
      !           = 6 WIND DIRECTION  
      !           = 7 RELATIVE HUMIDITY  
      !           = 8 RUNOFF CURVE NUMBER  
      !           = 9 WITHIN DAY WIND SPEED DIST   
      IF(randomCycles>0)THEN
            DO J=1,20  
                  RN=Generate_Random(21)  
                  II=INT(100*randomCycles*RN)  
                  DO KK=1,II   
                        XX=Generate_Random(21)  
                  END DO  
                  IX(J)=IX(21)  
            END DO  
            CALL Shuffle 
      END IF  
                                                
      IX0=IX                       ! IX is 1:21  copy IX
      randID0=randID               ! copy randID  
      
      WRITE(KW(1),1)randomCycles,(IX(randID(I)),I=1,10),(randID(I),I=1,10) 

    1 FORMAT(///1X,'-----GENERATOR SEEDS AFTER',I9,' CYCLES'/(5X,10I12))  
END SUBROUTINE Print_RandomNums   
! ------------------------------------4. write crop variables -------------------------------
! ------------------------------ . Update Crop Variables -------------- ---------------------
SUBROUTINE Update_CropVar 
      !     EPIC1102
      !     THIS SUBPROGRAM ALLOWS REAL TIME UPDATES OF CROP VARIABLES
      IMPLICIT NONE 
	! local variables potential problems with local variables (SAVE attribute XRTC)
      INTEGER:: I
      REAL, DIMENSION(8,12)::XRTC
      REAL:: X1
      CHARACTER(80)::file_name
 
      WRITE(KW(1),170)IY,MO,DayOfMon
     
      ! CALL Pest_Damage  
      ! Below is Pest_Damage subroutine. Modified by TXH
      IF(pestGrowIndx>0.)THEN
            X1=pestDamgFactor*pestGrowIndx/GrowingSeason_Days
            PSTF(Crop_Num)=1.-(1.-pestDmagFactor(Crop_Num))*X1/(X1+EXP(S_Curve(9,1)-S_Curve(9,2)*X1))
      ELSE
            PSTF(Crop_Num)=1.
      END IF

      WRITE(KW(1),174)Crop_Name(Crop_Num),waterStress,N_StressFactor,P_StressFactor,SK,VAR(72),SSLT,PSTF(Crop_Num)
      WRITE(KW(1),191)totCropBio(Crop_Num),Current_LAI(Crop_Num),PPL0(Crop_Num),actualCropN(Crop_Num),&
                      actualCropP(Crop_Num),actualCropK(Crop_Num),seedYieldPrice(Crop_Num),forageYieldPrice(Crop_Num)
      
      IF(ISTP==2)THEN
          IF(IYS(1)/=0)THEN
              totCropBio(Crop_Num)=XRTC(1,Crop_Num)
              Current_LAI(Crop_Num)=XRTC(2,Crop_Num)
                PPL0(Crop_Num)=XRTC(3,Crop_Num)
                XLAI(Crop_Num)=maxLAI(Crop_Num)*PPL0(Crop_Num)/(PPL0(Crop_Num)+EXP(coeffPopuCurve(1,Crop_Num)- &
                             coeffPopuCurve(2,Crop_Num)*PPL0(Crop_Num)))
              actualCropN(Crop_Num)=XRTC(4,Crop_Num)
              actualCropP(Crop_Num)=XRTC(5,Crop_Num)
              !actualCropK(Crop_Num)=XRTC(6,Crop_Num)
              seedYieldPrice(Crop_Num)=XRTC(7,Crop_Num)
              forageYieldPrice(Crop_Num)=XRTC(8,Crop_Num)
          END IF
          WRITE(KW(1),191)totCropBio(Crop_Num),Current_LAI(Crop_Num),PPL0(Crop_Num),actualCropN(Crop_Num), &
                          actualCropP(Crop_Num),actualCropK(Crop_Num),seedYieldPrice(Crop_Num),forageYieldPrice(Crop_Num)      
          RETURN
      END IF
      
      IF(ISTP==0)THEN
          file_name='RTCROP.DAT'
          CALL OPENV(KW(MSO+2),file_name,0,KW(MSO))
          WRITE(KW(MSO+2),3)IYS(1)
          WRITE(KW(MSO+2),2)totCropBio(Crop_Num),Current_LAI(Crop_Num),PPL0(Crop_Num),actualCropN(Crop_Num), &
                            actualCropP(Crop_Num),actualCropK(Crop_Num),seedYieldPrice(Crop_Num),forageYieldPrice(Crop_Num)
          ISX=1
      ELSE
          file_name='RTCROP.DAT'       
          CALL OPENV(KW(MSO+2),file_name,0,KW(MSO))
          READ(KW(MSO+2),3)IYS(1)
          READ(KW(MSO+2),2)(XRTC(I,Crop_Num),I=1,8)
      END IF
      CLOSE(KW(MSO+2))   
      RETURN
    2 FORMAT(10F8.3)
    3 FORMAT(I4)
  170 FORMAT(///'*****UPDATE  YR = ',I2,'  MO = ',I2,'  DA = ',I2)
  174 FORMAT(5X,A4,10F7.2)
  191 FORMAT(5X,'BIOM=',F8.3,' t/ha   LAI=',F5.2,3X,'PPOP=',F6.1,'N_Rate_bySoil=',&
      F5.1,' kg/ha P_Rate_bySoil=',F5.1,' kg/ha',3X,'K_Rate_bySoil=',F5.1,' kg/ha',3X,'seedYieldPrice=',&
      F5.0,'$/t',3X,'forageYieldPrice=',F5.0,'$/t')
 
END SUBROUTINE Update_CropVar

! -------------------------------- 3. Updata Soil variables --------------------------------

SUBROUTINE Update_SoilVar 
      !     EPIC1102
      !     THIS SUBPROGRAM ALLOWS REAL TIME UPDATES OF SOIL VARIABLES
      IMPLICIT NONE
 			
      ! local variables potential problems with local variables (SAVE attribute XRTS)
      INTEGER:: IYCD = 1, L, I, J
      REAL:: XRTS(7,15), Z1, X1, X2, X3, X4
 
      WRITE(KW(1),120)(soilWater(Layer_ID(L)),L=1,Actual_SoilLayers)
      WRITE(KW(1),120)(soilTem(Layer_ID(L)),L=1,Actual_SoilLayers)
      WRITE(KW(1),120)(NO3_N_Soil(Layer_ID(L)),L=1,Actual_SoilLayers)
      WRITE(KW(1),120)(NH3_Weight(Layer_ID(L)),L=1,Actual_SoilLayers)
      WRITE(KW(1),120)(labileP(Layer_ID(L)),L=1,Actual_SoilLayers)
      WRITE(KW(1),120)(exchangeK(Layer_ID(L)),L=1,Actual_SoilLayers)
      WRITE(KW(1),120)(saltWeight(Layer_ID(L)),L=1,Actual_SoilLayers)
      IF(ISTP==2)GO TO 2
      IF(ISTP==0)GO TO 198
      READ(KW(MSO+3),9)(IYS(I),I=2,8)
      DO J=1,Actual_SoilLayers
          L=Layer_ID(J)
          READ(KW(MSO+3),6)(XRTS(I,L),I=1,7)
      END DO
      REWIND KW(MSO+3)
      ISTP=0
      WRITE(KW(MSO+3),9)ISTP
      ISTP=2
    2 DO J=1,Actual_SoilLayers
          L=Layer_ID(J)
          IF(IYS(2)==1)soilWater(L)=XRTS(1,L)*(fieldCapacity(L)-Wilt_Point(L))+Wilt_Point(L)
          IF(IYS(3)==1)soilTem(L)=XRTS(2,L)
          IF(IYS(4)==1)NO3_N_Soil(L)=XRTS(3,L)
          IF(IYS(5)==1)NH3_Weight(L)=XRTS(4,L)
          IF(IYS(6)==1)labileP(L)=XRTS(5,L)
          IF(IYS(7)==1)exchangeK(L)=XRTS(6,L)
          IF(IYS(8)==1)saltWeight(L)=XRTS(7,L)
      END DO
      GO TO 196
  198 REWIND KW(MSO+3)
      IF(ISX/=0)THEN
          ISTP=1
          WRITE(KW(MSO+3),9)ISTP
      END IF
      ISX=1
      IF(IYCD/=0)THEN
            DO I=2,8
                  IYS(I)=1
            END DO
      END IF
      WRITE(KW(MSO+3),9)(IYS(I),I=2,8)
      soilWaterRZ=0.
      potentialAvailWater=0.
      Z1=0.
      DO J=1,Actual_SoilLayers
            L=Layer_ID(J)
            X1=soilWater(L)-Wilt_Point(L)
            X4=.001*X1/(Z(L)-Z1)
            X2=fieldCapacity(L)-Wilt_Point(L)
            X3=X1/X2
            soilWaterRZ=soilWaterRZ+X1
            potentialAvailWater=potentialAvailWater+X2
            WRITE(KW(MSO+3),11)Z(L),X1,X4,X2,X3,soilTem(L),NO3_N_Soil(L),NH3_Weight(L),&
                  labileP(L),exchangeK(L),saltWeight(L)
      END DO
      WRITE(KW(MSO+3),11)RZ,soilWaterRZ,potentialAvailWater
  196 WRITE(KW(1),120)(soilWater(Layer_ID(L)),L=1,Actual_SoilLayers)
      WRITE(KW(1),120)(soilTem(Layer_ID(L)),L=1,Actual_SoilLayers)
      WRITE(KW(1),120)(NO3_N_Soil(Layer_ID(L)),L=1,Actual_SoilLayers)
      WRITE(KW(1),120)(NH3_Weight(Layer_ID(L)),L=1,Actual_SoilLayers)
      WRITE(KW(1),120)(labileP(Layer_ID(L)),L=1,Actual_SoilLayers)
      WRITE(KW(1),120)(exchangeK(Layer_ID(L)),L=1,Actual_SoilLayers)
      WRITE(KW(1),120)(saltWeight(Layer_ID(L)),L=1,Actual_SoilLayers)
 
    6 FORMAT(40X,F10.3,F10.2,5F10.0)
    9 FORMAT(20I4)
   11 FORMAT(F10.2,F10.1,F10.3,F10.1,2F10.2,5F10.0)
  120 FORMAT(1X,10F10.2)
 
END SUBROUTINE Update_SoilVar

!-----------------------------4. Print daily Organic C and N -----------------------
SUBROUTINE Print_OrgNC(KK) 
      !     EPIC1102
      !     THIS SUBPROGRAM OUTPUTS THE SOIL ORGANIC C AND N VARIABLES DAILY
      IMPLICIT NONE
 
      ! local variables
      INTEGER:: J, I, KK
 
      WRITE(KW(15),30)IYR,MO,KK
      WRITE(KW(15),2)(SID(LORG(Layer_ID(J))),J=1,Actual_SoilLayers),SID(16)
      WRITE(KW(15),3)(Z(Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(15),21)(SOIL(12,Layer_ID(I)),I=1,Actual_SoilLayers),totSoilWater
      WRITE(KW(15),27)(soilTem(Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(15),28)(cropResidu(Layer_ID(I)),I=1,Actual_SoilLayers),totCropResi ! Why crop residu is distributed in different soil layers?
      WRITE(KW(15),8)WBMX
      WRITE(KW(15),1)(soilResp(Layer_ID(I)),I=1,Actual_SoilLayers),Tot_SoilRes
      WRITE(KW(15),4)(XNetNMineralize(Layer_ID(I)),I=1,Actual_SoilLayers),mineralizedN
      !WRITE(KW(15),5)(DN2G(Layer_ID(I)),I=1,Actual_SoilLayers),SN2
      !WRITE(KW(15),6)(DN2OG(Layer_ID(I)),I=1,Actual_SoilLayers),totN2O
      WRITE(KW(15),7)(DN2G(I),I=1,layersEqualThick),SN2
      RETURN
    1 FORMAT(T5,'CO2 LOSS(kg/ha)',T20,16F8.3)
    2 FORMAT(T52,'SOIL LAYER NO'/T18,16(4X,A4))
    3 FORMAT(T5,'DEPTH(m)',T20,16F8.2)
    4 FORMAT(T5,'NET MN(kg/ha)',T20,16F8.2)
    !5 FORMAT(T5,'DN2G(kg/ha)',T20,16F8.4)
    !6 FORMAT(T5,'DN2OG(kg/ha)',T20,16F8.4) 
    7 FORMAT(T5,'DN2G(kg/ha)',T20,20F8.4)
    8 FORMAT(T5,'BIOMIX',T20,16F8.4) 
   21 FORMAT(T5,'totSoilWater(m/m)',T20,16F8.3)
   27 FORMAT(T5,'TEMP(c)',T20,16F8.2)
   28 FORMAT(T5,'cropResidu(t/ha)',T20,16F8.2)
   30 FORMAT(//T10,3I4)
END SUBROUTINE Print_OrgNC
 
! --------------------------- 5. Print Soil Organic C and N variables -----------------------------
SUBROUTINE Print_OrgNCvar(IYR1,MZ,KK)
      !     EPIC1102
      !     THIS SUBPROGRAM OUTPUTS THE SOIL ORGANIC C AND N VARIABLES
      IMPLICIT NONE
 
      ! local variables
      INTEGER:: IYR1,MZ,KK, I, J
      REAL:: X1

      WRITE(KW(14),30)IYR1,MZ,KK,CO2
      WRITE(KW(14),2)(SID(LORG(Layer_ID(J))),J=1,Actual_SoilLayers),SID(16)
      WRITE(KW(14),3)'DEPTH(m)',(Z(Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(14),3)'bulkDensity 33kpa(t/m3)',(bulkDensity(Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(14),1)'SAND(%)',(sandFrac(Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(14),1)'SILT(%)',(siltFrac(Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(14),1)'CLAY(%)',(clayFrac(Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(14),1)'ROCK(%)',(rockFrac(Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(14),1)'structLitt(kg/ha)',(structLitt(Layer_ID(I)),I=1,Actual_SoilLayers),totStructLitt
      WRITE(KW(14),1)'metabLitt(kg/ha)',(metabLitt(Layer_ID(I)),I=1,Actual_SoilLayers),totMetabLitt
      WRITE(KW(14),1)'lgStructLitt(kg/ha)',(lgStructLitt(Layer_ID(I)),I=1,Actual_SoilLayers),totLgStructLitt
      WRITE(KW(14),1)'CStructLitt(kg/ha)',(CStructLitt(Layer_ID(I)),I=1,Actual_SoilLayers),totCStructLitt
      WRITE(KW(14),1)'CMetabLitt(kg/ha)',(CMetabLitt(Layer_ID(I)),I=1,Actual_SoilLayers),totCMetabLitt
      WRITE(KW(14),1)'CLgStructLitt(kg/ha)',(CLgStructLitt(Layer_ID(I)),I=1,Actual_SoilLayers),totCLgStructLitt
      WRITE(KW(14),1)'NLgStructLitt(kg/ha)',(NLgStructLitt(Layer_ID(I)),I=1,Actual_SoilLayers),totNLgStructLitt
      WRITE(KW(14),1)'CBiomass(kg/ha)',(CBiomass(Layer_ID(I)),I=1,Actual_SoilLayers),totCBiomass
      WRITE(KW(14),4)'CSlowHumus(kg/ha)',(CSlowHumus(Layer_ID(I)),I=1,Actual_SoilLayers),totCSlowHumus
      WRITE(KW(14),4)'CPassiveHumus(kg/ha)',(CPassiveHumus(Layer_ID(I)),I=1,Actual_SoilLayers),totCPassiveHumus
      X1=.001*totSOC
      WRITE(KW(14),4)'SOC(kg/ha)',(SOC(Layer_ID(I)),I=1,Actual_SoilLayers),X1
      WRITE(KW(14),1)'NStructLitt(kg/ha)',(NStructLitt(Layer_ID(I)),I=1,Actual_SoilLayers),totNStructLitt
      WRITE(KW(14),1)'NMetabLitt(kg/ha)',(NMetabLitt(Layer_ID(I)),I=1,Actual_SoilLayers),totNMetabLitt
      WRITE(KW(14),1)'NBiomass(kg/ha)',(NBiomass(Layer_ID(I)),I=1,Actual_SoilLayers),totNBiomass
      WRITE(KW(14),4)'NSlowHumus(kg/ha)',(NSlowHumus(Layer_ID(I)),I=1,Actual_SoilLayers),totNSlowHumus
      WRITE(KW(14),4)'NPassiveHumus(kg/ha)',(NPassiveHumus(Layer_ID(I)),I=1,Actual_SoilLayers),totNPassiveHumus
      WRITE(KW(14),4)'SON(kg/ha)',(SON(Layer_ID(I)),I=1,Actual_SoilLayers),totSON
      WRITE(KW(14),4)'CFEM(kg/ha)',SMY(96)
      RETURN
    1 FORMAT(1X,A14,16F12.1)
    2 FORMAT(T52,'SOIL LAYER NO'/T18,16(4X,A4,4X))
    3 FORMAT(1X,A14,16F12.2)
    4 FORMAT(1X,A14,16F12.0)
   30 FORMAT(//T10,3I4,5X,'ATMOS CO2 =',F5.0,' ppm')
 
END SUBROUTINE Print_OrgNCvar
                           
! 6-------------------------  Print  DENITRIFICATION ----------------------------------------
SUBROUTINE Print_NO3Loss
      !     EPIC1102
      !     THIS SUBPROGRAM OUTPUTS DAILY OUTPUT FROM  Denitrification 
 
      IMPLICIT NONE
      ! local variables
      INTEGER:: I, J
 
      WRITE(KW(32),30)IYR,MO,DayOfMon
      WRITE(KW(32),2)(SID(LORG(Layer_ID(J))),J=1,Actual_SoilLayers),SID(16)
      WRITE(KW(32),3)(Z(Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(32),21)(SOIL(12,Layer_ID(I)),I=1,Actual_SoilLayers),totSoilWater
      WRITE(KW(32),27)(soilTem(Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(32),28)(cropResidu(Layer_ID(I)),I=1,Actual_SoilLayers),totCropResi
      WRITE(KW(32),8)WBMX
      WRITE(KW(32),1)(soilResp(Layer_ID(I)),I=1,Actual_SoilLayers),Tot_SoilRes
      WRITE(KW(32),4)(XNetNMineralize(Layer_ID(I)),I=1,Actual_SoilLayers),mineralizedN
      WRITE(KW(32),9)(NO3_N_Soil(Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(32),10)(H2OF(Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(32),11)(WNO3F(Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(32),12)(CBNF(Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(32),5)(weightDenitriN(Layer_ID(I)),I=1,Actual_SoilLayers),totDenitriN
 
      RETURN
    1 FORMAT(T5,'CO2 LOSS(kg/ha)',T20,16F8.3)
    2 FORMAT(T52,'SOIL LAYER NO'/T18,16(4X,A4))
    3 FORMAT(T5,'DEPTH(m)',T20,16F8.2)
    4 FORMAT(T5,'NET MN(kg/ha)',T20,16F8.2)
    5 FORMAT(T5,'DN(kg/ha)',T20,16F8.4)
    ! 6 FORMAT(T5,'DN2OG(kg/ha)',T20,16F8.4) 
    ! 7 FORMAT(T5,'DN2G(kg/ha)',T20,20F8.4)
    8 FORMAT(T5,'BIOMIX',T20,16F8.4) 
    9 FORMAT(T5,'NO3_N_Soil(kg/ha)',T20,16F8.4)      
   10 FORMAT(T5,'H2OF',T20,16F8.4)      
   11 FORMAT(T5,'WNO3F',T20,16F8.4)            
   12 FORMAT(T5,'CBNF',T20,16F8.4)      
   21 FORMAT(T5,'totSoilWater(m/m)',T20,16F8.3)
   27 FORMAT(T5,'TEMP(c)',T20,16F8.2)
   28 FORMAT(T5,'cropResidu(t/ha)',T20,16F8.2)
   30 FORMAT(//T10,3I4)
END SUBROUTINE Print_NO3Loss

! --------------------- 7. write the names of input weather variables ---------------
SUBROUTINE Print_WeatherVar 
      !     EPIC1102
      !     THIS SUBPROGRAM WRITES THE NAMES OF THE INPUT WEATHER VARIABLES.
      IMPLICIT NONE
 
      ! local variables 
      INTEGER:: K, N1, N2, J
 
      K=1
      N1=weatherVarID
      WRITE(KW(1),673)FWTH
      DO J=4,1,-1
          CALL split_int(N1,N2)
          IF(N1==0)GO TO 1
          KGN(N1)=1
          SELECT CASE(N1)
              CASE(1)
                  GO TO 42
              CASE(2)
                  KDT2(K)=67
              CASE(3)
                  KDT2(K)=3
              CASE(4) 
                  KDT2(K)=7
              CASE(5)
                  KDT2(K)=8
              CASE DEFAULT
          END SELECT
          K=K+1
   42     N1=N2
      END DO
    1 SELECT CASE(K)
          CASE(1)
              WRITE(KW(1),222)
          CASE(2)
              WRITE(KW(1),356)varName(KDT2(1))
          CASE(3)
              WRITE(KW(1),357)(varName(KDT2(J)),J=1,2)
          CASE(4)
              WRITE(KW(1),358)(varName(KDT2(J)),J=1,3)
          CASE(5)
              WRITE(KW(1),359)(varName(KDT2(J)),J=1,4)
      END SELECT
      RETURN
  222 FORMAT(/T10,'**********RAIN IS INPUT**********')
  356 FORMAT(/T10,'**********RAIN,',1X,A4,', ARE INPUT**********')
  357 FORMAT(/T10,'**********RAIN,',2(1X,A4,','),' ARE INPUT**********')
  358 FORMAT(/T10,'**********RAIN,',3(1X,A4,','),' ARE INPUT**********')
  359 FORMAT(/T10,'**********RAIN,',4(1X,A4,','),' ARE INPUT**********')
  673 FORMAT(T10,'DAILY WEATHER FILE = ',A80)

END SUBROUTINE Print_WeatherVar
 
! --------------------- Print date, time, title, page number -------------
SUBROUTINE Print_Page(IP)
      !     EPIC1102
      !     THIS SUBPROGRAM CHANGES PAGES, WRITES VERSION, DATE, & TITLE
      IMPLICIT NONE
      ! local variables
      INTEGER:: IP
      WRITE(KW(1),1)IYER,IMON,IDAY,Hour,Mint,Sec
      WRITE(KW(1),2)TITLE
      WRITE(KW(1),3)IRUN,IRO0,randomCycles
      IF(IP==0)RETURN
      WRITE(KW(1),4)SITEFILE
      WRITE(KW(1),4)SOILFILE
      WRITE(KW(1),4)OPSCFILE
      RETURN
    1 FORMAT('1'/T5,'EEPIC1102 - 2020 V1',2X,3I4,2X,2(I2,':'),I2)
    2 FORMAT(/(10X,20A4))
    3 FORMAT(5X,'RUN #=',I4,2X,'ROT #=',I4,2X,'GNSD #=',I4)
    4 FORMAT(5X,A80)
 
END SUBROUTINE Print_Page
! -------------------------Print Soil Information ---------------------------
SUBROUTINE Print_SoilVar1(soilName, ASG, adjustCN2, splitedSoilLayer, minThickSplit, &
                          fracSOC, fracHumus, cultivaYrs, soilCondCode, soilGroup)
      IMPLICIT NONE
      ! arguments list
      CHARACTER(LEN = 1), DIMENSION(4), INTENT(IN):: ASG
      CHARACTER(LEN = 80), INTENT(IN):: soilName
      REAL, INTENT(IN):: adjustCN2, splitedSoilLayer, minThickSplit, fracSOC, fracHumus, cultivaYrs, &
                         soilCondCode, soilGroup   
      ! local variables

      WRITE(KW(1),'(//1X,A/)')'____________________SOIL DATA____________________'  
      WRITE(KW(1),'(T10,A,A80)')'SOIL = ',SOILFILE 
      WRITE(KW(1),'(T10,A,A20)')'SOIL SERIES = ',soilName  
      WRITE(KW(1),'(T10,A,A20)')'SOIL ORDER = ',soilOrder 
      WRITE(KW(1),1)soilAlbedo,splitedSoilLayer,minThickMaxLayer,minThickProfile,minThickSplit,fracSOC,fracHumus,&
                      cultivaYrs,soilCondCode,soilGroup,plawDepth,plawDepSOC   
      WRITE(KW(1),2)waterTableMinDep,waterTableMaxDep,waterTableHigt,groundWaterStor,groundWaterMaxStor,groundWaterResidenceTime  
      
      IF(erosionMode>0)WRITE(KW(1),'(T10,A)')'STATIC SOIL PROFILE' 

      SELECT CASE(FC_Method+1)  
            CASE(1)  
                  WRITE(KW(1),'(T10,A)')'fieldCapacity/WP EST RAWLS METHOD DYNAMIC' 
            CASE(2)   
                  WRITE(KW(1),'(T10,A)')'fieldCapacity/WP EST BAUMER METHOD DYNAMIC'!? this method is not introduced by APEX-doc 
            CASE(3)  
                  WRITE(KW(1),'(T10,A)')'fieldCapacity/WP printInterval RAWLS METHOD DYNAMIC'  
            CASE(4)  
                  WRITE(KW(1),'(T10,A)')'fieldCapacity/WP printInterval BAUMER METHOD DYNAMIC' 
            CASE(5) 
                  WRITE(KW(1),'(T10,A)')'fieldCapacity/WP EST RAWLS METHOD STATIC'                         
            CASE(6)  
                  WRITE(KW(1),'(T10,A)')'fieldCapacity/WP EST BAUMER METHOD STATIC' 
            CASE(7)  
                  WRITE(KW(1),'(T10,A)')'fieldCapacity/WP printInterval STATIC' 
            CASE(8)  
                  WRITE(KW(1),'(T10,A)')'fieldCapacity/WP printInterval NEAREST NEIGHBOR METHOD &  
                                          &DYNAMIC'                                                               
            CASE(9) 
                  WRITE(KW(1),'(T10,A)')'fieldCapacity/WP printInterval NEAREST NEIGHBOR METHOD &  
                                          &STATIC'                                                            
            CASE(10) 
                  WRITE(KW(1),'(T10,A)')'fieldCapacity/WP EST BNW METHOD DYNAMIC'
            CASE(11)  
                  WRITE(KW(1),'(T10,A)')'fieldCapacity/WP EST BNW METHOD STATIC' 
      END SELECT
	
      IF(soilSaturate_Method>0)THEN
            WRITE(KW(1),'(T10,A)')'SAT COND ESTIMATED WITH RAWLS METHOD'
      ELSE
            WRITE(KW(1),'(T10,A)')'SAT COND INPUT'
      END IF  
      satuateCondv0=satuateCondv(Layer_ID(2))
      WRITE(KW(1),'(T10,A,F10.3)')'INITIAL STANDING DEAD = ',deadCropRes0                                                              
       
      WRITE(KW(1),3)ASG(hydroGroup),landUseCode,CN2,adjustCN2,S_Curve(30,1),S_Curve(30,2),S_Curve(4,1),&              
      S_Curve(4,2)  
    
    1 FORMAT(T10,'SOIL ALBEDO = ',F6.2/T10,'MAX NUMBER SOIL LAYERS = ',F3.0/T10,&
            'MIN THICKNESS FOR LAYER SPLITTING = ',F6.2,' m'/T10, &
            'MIN PROFILE THICKNESS--STOPS SIMULATION = ',F6.2,' m'/T10, &
            'SPLITTING PRIORITY THICKNESS = ',F6.2,' m'/T10,'BIOMASS/ORG C = ',F6.3/T10, &              
            'PASSIVE HUMUS/TOTAL HUMUS = ',F6.3/T10,'CULTIVATION HISTORY = ', F5.0,' y'/T10,&
            'WEATHERING CODE = ',F3.0/T10,'SOIL GROUP CODE =',1X,F3.0/T10,'ORG C IN TOP ',F6.2,' m = ',F7.1,' t/ha') 
    2 FORMAT(T10,'WATER TABLE DEPTH'/T15,'MIN = ',F6.2,' m'/T15,'MAX = ', &            
            F6.2,' m'/T15,'INITIAL = ',F6.2,' m'/T10,'GROUNDWATER STORAGE = ', &            
            F5.0,' mm'/T10,'MAX GROUNDWATER STORE = ',F5.0,' mm'/T10,'RETURN', &            
            ' FLOW TT = ',F7.2,' d')  
    3 FORMAT(T10,'HYDROLOGIC SOIL GROUP = ',A1/T10,'LAND USE NUMBER = ', I3/T10,'RUNOFF CN2 = ',F4.1/T10,&
            'SLP ADJ CN2 = ',F4.1/T10,'CN SCRVS_Curve(30)= ',2F6.0/T10,'CN S_Curve(4)= ',2E13.5)  
END SUBROUTINE Print_SoilVar1
! ------------------------ 10.2 Print Soil Variables ------------------------
SUBROUTINE Print_SoilVar2(YTP,L)
      !EPIC1102
      !THIS SUBPROGRAM OUTPUTS THE SOIL TABLE
      IMPLICIT NONE
      ! local variables
      INTEGER:: I, J, L
      REAL:: YTP(16), X1
      
      WRITE(KW(L),'(1X,A,1X,I4)')'YR=',IY
      WRITE(KW(L),2)(SID(LORG(Layer_ID(J))),J=1,Actual_SoilLayers),SID(16)
      WRITE(KW(L),1)'DEPTH(m)',(Z(Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(L),4)'POROSITY(m/m)',(SOIL(8,Layer_ID(I)),I=1,Actual_SoilLayers),YTP(4)
      WRITE(KW(L),4)'fieldCapacity totSoilWater(m/m)',(SOIL(9,Layer_ID(I)),I=1,Actual_SoilLayers),YTP(1)
      WRITE(KW(L),4)'WP totSoilWater(m/m)',(SOIL(20,Layer_ID(I)),I=1,Actual_SoilLayers),YTP(3)
      WRITE(KW(L),4)'totSoilWater(m/m)',(SOIL(12,Layer_ID(I)),I=1,Actual_SoilLayers),YTP(2)
      WRITE(KW(L),3)'SAT COND(mm/h)',(satuateCondv(Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(L),3)'H SC(mm/h)',(lateralFlowCondv(Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(L),1)'bulkDensity 33kpa(t/m3)',(bulkDensity(Layer_ID(I)),I=1,Actual_SoilLayers),YTP(5)
      WRITE(KW(L),1)'bulkDensity DRY(t/m3)',(SOIL(13,Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(L),5)'SAND(%)',(sandFrac(Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(L),5)'SILT(%)',(siltFrac(Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(L),5)'CLAY(%)',(clayFrac(Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(L),5)'ROCK(%)',(rockFrac(Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(L),5)'PH',(PH(Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(L),5)'SM BS(cmol/kg)',(totBase(Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(L),5)'CEC(cmol/kg)',(CEC(Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(L),5)'AL SAT(%)',(ALS(Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(L),5)'CaCO3 (%)',(CaCO3(Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(L),6)'LAB P(g/t)',(SOIL(1,Layer_ID(I)),I=1,Actual_SoilLayers),totLabileP
      WRITE(KW(L),3)'P SORP RTO',(sorbRatioP(Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(L),6)'MN P AC(g/t)',(SOIL(2,Layer_ID(I)),I=1,Actual_SoilLayers),totActiveP
      WRITE(KW(L),6)'MN P soilWater(g/t)',(SOIL(3,Layer_ID(I)),I=1,Actual_SoilLayers),TOP
      WRITE(KW(L),6)'ORG P(g/t)',(SOIL(4,Layer_ID(I)),I=1,Actual_SoilLayers),initialTotSOP
      WRITE(KW(L),6)'NO3(g/t)',(SOIL(5,Layer_ID(I)),I=1,Actual_SoilLayers),totNO3_N
      WRITE(KW(L),6)'solubleK(g/t)',(SOIL(14,Layer_ID(I)),I=1,Actual_SoilLayers),totSolubleK
      WRITE(KW(L),6)'exchangeK(g/t)',(SOIL(15,Layer_ID(I)),I=1,Actual_SoilLayers),totExchangeK
      WRITE(KW(L),6)'fixK(g/t)',(SOIL(16,Layer_ID(I)),I=1,Actual_SoilLayers),totFixK
      WRITE(KW(L),6)'ORG N(g/t)',(SOIL(6,Layer_ID(I)),I=1,Actual_SoilLayers),totSON
      X1=.001*totSOC
      WRITE(KW(L),3)'ORG C(%)',(SOIL(7,Layer_ID(I)),I=1,Actual_SoilLayers),X1
      WRITE(KW(L),3)'CROP cropResidu(t/ha)',(cropResidu(Layer_ID(I)),I=1,Actual_SoilLayers),totCropResi
      WRITE(KW(L),5)'structLitt(kg/ha)',(structLitt(Layer_ID(I)),I=1,Actual_SoilLayers),totStructLitt
      WRITE(KW(L),5)'metabLitt(kg/ha)',(metabLitt(Layer_ID(I)),I=1,Actual_SoilLayers),totMetabLitt
      WRITE(KW(L),5)'lgStructLitt(kg/ha)',(lgStructLitt(Layer_ID(I)),I=1,Actual_SoilLayers),totLgStructLitt
      WRITE(KW(L),5)'CStructLitt(kg/ha)',(CStructLitt(Layer_ID(I)),I=1,Actual_SoilLayers),totCStructLitt
      WRITE(KW(L),5)'CMetabLitt(kg/ha)',(CMetabLitt(Layer_ID(I)),I=1,Actual_SoilLayers),totCMetabLitt
      WRITE(KW(L),5)'CLgStructLitt(kg/ha)',(CLgStructLitt(Layer_ID(I)),I=1,Actual_SoilLayers),totCLgStructLitt
      WRITE(KW(L),5)'NLgStructLitt(kg/ha)',(NLgStructLitt(Layer_ID(I)),I=1,Actual_SoilLayers),totNLgStructLitt
      WRITE(KW(L),5)'CBiomass(kg/ha)',(CBiomass(Layer_ID(I)),I=1,Actual_SoilLayers),totCBiomass
      WRITE(KW(L),5)'CSlowHumus(kg/ha)',(CSlowHumus(Layer_ID(I)),I=1,Actual_SoilLayers),totCSlowHumus
      WRITE(KW(L),5)'CPassiveHumus(kg/ha)',(CPassiveHumus(Layer_ID(I)),I=1,Actual_SoilLayers),totCPassiveHumus
      X1=.001*totSOC
      WRITE(KW(L),5)'SOC(kg/ha)',(SOC(Layer_ID(I)),I=1,Actual_SoilLayers),X1
      WRITE(KW(L),5)'NStructLitt(kg/ha)',(NStructLitt(Layer_ID(I)),I=1,Actual_SoilLayers),totNStructLitt
      WRITE(KW(L),5)'NMetabLitt(kg/ha)',(NMetabLitt(Layer_ID(I)),I=1,Actual_SoilLayers),totNMetabLitt
      WRITE(KW(L),5)'NBiomass(kg/ha)',(NBiomass(Layer_ID(I)),I=1,Actual_SoilLayers),totNBiomass
      WRITE(KW(L),5)'NSlowHumus(kg/ha)',(NSlowHumus(Layer_ID(I)),I=1,Actual_SoilLayers),totNSlowHumus
      WRITE(KW(L),5)'NPassiveHumus(kg/ha)',(NPassiveHumus(Layer_ID(I)),I=1,Actual_SoilLayers),totNPassiveHumus
      WRITE(KW(L),6)'SON(kg/ha)',(SON(Layer_ID(I)),I=1,Actual_SoilLayers),totSON
      WRITE(KW(L),3)'elecCondv(mmho/cm)',(elecCondv(Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(L),6)'saltWeight(kg/ha)',(saltWeight(Layer_ID(I)),I=1,Actual_SoilLayers),totSalt
      WRITE(KW(L),4)'fracNO3Leach',(fracNO3Leach(Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(L),4)'gasO2Con(kg/ha)',(gasO2Con(Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(L),4)'gasCO2Con(kg/ha)',(gasCO2Con(Layer_ID(I)),I=1,Actual_SoilLayers)
      WRITE(KW(L),4)'gasN2OCon(kg/ha)',(gasN2OCon(Layer_ID(I)),I=1,Actual_SoilLayers)
      RETURN
    1 FORMAT(4X,A30,1X,16F12.5)  
    2 FORMAT(T52,'SOIL LAYER NO'/T42,16(A5,7X))
    3 FORMAT(4X,A30,1X,16F12.2)
    4 FORMAT(4X,A30,1X,16F12.3)
    5 FORMAT(4X,A30,1X,16F12.1)
    6 FORMAT(4X,A30,1X,16F12.0)
 
END SUBROUTINE Print_SoilVar2

! -------------------------- Print Management data -----------------------
SUBROUTINE Print_ManagementInfo(maxRainfall, yearTotPPT, drainageDepth, lagoonVolx, manureWeightX)
      IMPLICIT NONE
      ! arguments list
      REAL, INTENT(IN):: maxRainfall, yearTotPPT
      INTEGER, INTENT(IN):: drainageDepth
      REAL, INTENT(OUT):: lagoonVolx, manureWeightX
      ! local variables
      INTEGER:: NCOW, I, L
      REAL:: X1, X3, runoffToLagoon, QWW

      ! ///////////////////////////////////////////////////////////////////////
      WRITE(KW(1),'(//1X,A/)')'____________________MANAGEMENT DATA_____________________'   
      WRITE(KW(1),'(T10,A,A80)')'OPSC = ',OPSCFILE
      NCOW=INT(areaWshed/baseStockRate+.99) 
      ! ------------------------------ irrigation -------------------------
      IF(irrCode==0)THEN
            irrAnnualMaxVol=0.                                                  
            WRITE(KW(1),'(T10,A)')'DRYLAND AGRICULTURE'  
            GO TO 142
      ELSE
            IF(irrAnnualMaxVol < 1.E-10) irrAnnualMaxVol=2000.  
            IF(irrCode==4)THEN
                  WRITE(KW(1), 1)  
                  GO TO 137
            END IF
      END IF
     
      IF(irrTrigger<0.)THEN
            WRITE(KW(1),1)                                                               
            WRITE(KW(1),2)irrTrigger,autoIrrInterval                                                        
            GO TO 137
      END IF      
 
      IF(irrTrigger>0.)THEN
            WRITE(KW(1),1)                                                               
            IF(irrTrigger>1.)THEN                                                                 
                  WRITE(KW(1),3)irrTrigger,autoIrrInterval                                                      
            ELSE                                                                           
                  WRITE(KW(1),4)irrTrigger,autoIrrInterval                                                      
            END IF                                                                         
      ELSE
            WRITE(KW(1),'(T10,A)')'USER SPECIFIED IRRIGATION'
      END IF    
      
  137 WRITE(KW(1),5)irrAnnualMaxVol,irrSingleMinVol,irrSingleMaxVol   
      
      IF(irrChoice==0)THEN  
            WRITE(KW(1),'(T15,A)')'VARIABLE APPL VOLUMES'  
      ELSE                                                                           
            WRITE(KW(1),'(T15,A)')'FIXED APPL VOLUMES'   
      END IF
 
      SELECT CASE(irrCode)  
            CASE(1)   
                  WRITE(KW(1),'(T10,A)')'SPRINKLER IRRIGATION'   
            CASE(2)   
                  WRITE(KW(1),'(T10,A)')'FURROW IRRIGATION'  
            CASE(3)  
                  WRITE(KW(1),'(T10,A)')'IRRIGATION WITH FERT ADDED'  
                  WRITE(KW(1),'(T15,A,F6.3)')'RUNOFF RATIO = ',IrrRunoffRatio   
            CASE(5)  
                  WRITE(KW(1),'(T10,A)')'DRIP IRRIGATION'  
            CASE(4)
                  X3=lagoonInput*NCOW  
                  X1=maxRainfall-12.7 
                  runoffToLagoon=X1*X1/(maxRainfall+50.8)  
                  DALG=liqManureRatio*areaWshed  
                  X1=10.*DALG   
                  QWW=30.*X3/X1   
                  lagoonMaxVol=runoffToLagoon+QWW    
                  lagoonVol0=lagoonVol0*lagoonMaxVol   
                  WRITE(KW(1),6)DALG,lagoonVol0,lagoonMaxVol,dayReduceLagoon,lagoonInput  
                  lagoonInput=X3  
                  lagoonVolX=lagoonVol0   
                  lagoonVol0=X1*lagoonVol0  
                  lagoonMaxVol=X1*lagoonMaxVol 
                  lagoonVol=lagoonVol0  
                  fertCon=100.  
                  manureWeight=fertCon*lagoonVol    
                  manureWeightX=manureWeight           ! copy of manureWeight   
                  irrFromLagoon=(lagoonMaxVol-lagoonVol0)/(dayReduceLagoon+1.E-5) 
                  manureInputToLagoon=NCOW*fertRate*herdFrac  
                  AFLG=365.*manureInputToLagoon/areaWshed  
                  aveWaterFromLagoon=.1*(yearTotPPT*DALG+365.*X3)/areaWshed  
                  lagoonInput=X3  
                  GO TO 433
            CASE DEFAULT
                  irrCode=0                                                                              
      END SELECT    
      ! ----------------------------- fertilizer ---------------------------    
  142 IF(fertTrigger0>0.)THEN
            WRITE(KW(1),'(T10,A/T15,A,I3,A)')'AUTO SCHEDULED FERT','MIN APPL INTERVAL = ',fertApplyMinInterval,' D' 
            IAUF=1  
            IF(fertTrigger0>1.)THEN  
                  WRITE(KW(1),'(T15,A,F4.0,A)')'SOIL NO3 CONC TRIGGER = ',fertTrigger0,' g/t'                                    
            ELSE                                                                           
                  WRITE(KW(1),'(T15,A,F6.2)')'PLANT STRESS TRIGGER = ',fertTrigger0  
            END IF     
            IF(fertRate>1.)THEN   
                  WRITE(KW(1),'(T15,A,F8.1,A)')'FIXED RATE = ',fertRate,' kg/ha'   
            ELSE    
                  WRITE(KW(1),'(T15,A)')'VARIABLE RATE'  
            END IF  
      ELSE
            IAUF=0  
            IF(autoManureTrigger==0)WRITE(KW(1),'(T10,A)')'USER SCHEDULED FERT' 
      END IF  
      
  433 WRITE(KW(1),'(T15,A,F7.0,A)')'MAX N FERT/CROP = ',maxAnnualN,' kg/ha'
      ! P fertilizer
      IF(P_apply_mode>0)THEN
          WRITE(KW(1),'(T10,A)')'AUTO P FERT'
      ELSE
          WRITE(KW(1),'(T10,A)')'MANUAL P FERT'
      END IF  
      
      WRITE(KW(1),8)NCOW,grazeLimit,herdFrac 

      NII=autoIrrInterval                           ! N day application interval for automatic irrigation                                                                       
      furrowDikeSafeFactor=furrowDikeSafeFactor0    ! Furrow dike safety factor (0.1 - 1.0)                                                                    
      IF(furrowDikeSafeFactor0<1.E-10)furrowDikeSafeFactor=.9                                                         
      ! drainage
      IDR=drainageDepth  
      
      IF(drainageDepth>0)THEN
            X1=.001*drainageDepth                                                                   
            DO I=1,Actual_SoilLayers                                                                    
                  L=Layer_ID(I)                                                                     
                  IF(Z(L)>X1)EXIT                                                              
            END DO                                                                              
            IDR=L                                                                          
            lateralCondvHyN=lateralFlowCondv(L)                                                                    
            lateralFlowCondv(L)=MAX(10.*satuateCondv(L),(Porosity(L)-Wilt_Point(L))/(24.*dayReduceStress))                               
            lateralCondvHyD=lateralFlowCondv(L)                                                                         
            WRITE(KW(1),7) L, dayReduceStress, lateralFlowCondv(L)                                                   
            IF(dayReduceStress<groundWaterResidenceTime) groundWaterResidenceTime=dayReduceStress
      END IF   
 
      IF(furrowFlag>0) WRITE(KW(1),'(T10,A,F5.2)')'FURROW DIKE SYSTEM SAFETY FACTOR = ', furrowDikeSafeFactor
      IF(mowMinInterval>0) WRITE(KW(1),'(T10,A,I4,A)')'AUTO MOW INTERVAL = ',mowMinInterval,' d'  
      
      WRITE(KW(1),'(T10,A,F6.3)')'USLE P FACTOR = ',conservFactor 
      
      groundWaterResidenceTime=1.-EXP(-1./groundWaterResidenceTime)          ! APEX-doc Eq. 368 
      
      furrowDikeSafeFactor=furrowDikeSafeFactor*1000.   
      
      IF(limeFlag==0)THEN 
            WRITE(KW(1),'(T10,A)')'LIME APPLIED AS NEEDED'   
      ELSE                                                                           
            WRITE(KW(1),'(T10,A)')'NO LIME APPLICATIONS' 
      END IF    

    1 FORMAT(T10,'AUTOMATIC IRRIGATION') 
    2 FORMAT(T15,'potentialAvailWater DEFICIT TRIGGER = ',F5.0,' mm'/T15, &
            'MIN APPL INTERVAL = ',I3,' d')
    3 FORMAT(T15,'SOIL-WATER TENSION TRIGGER = ',F5.0,' kpa'/T15,'MIN APPL INTERVAL = ',I3,' d')  
    4 FORMAT(T15,'PLANT WATER STRESS TRIGGER = ',F6.2/T15,'MIN APPL INTERVAL = ',I3,' d')  
    5 FORMAT(T15,'MAX ANNUAL VOL APPL TO A CROP = ',F6.0,' mm'/T15,&
            'MIN SINGLE APPL VOL = ',F5.0,' mm'/T15,'MAX SINGLE APPL VOL = ',F5.0,' mm')  
    6 FORMAT(T15,'BIG GUN IRRIGATION FROM LAGOON'/T15,'CONFINEMENT AREA = ',F7.2,' ha'/T15, &
            'NORMAL LAGOON VOL = ',F5.0,' mm'/T15,'MAX LAGOON VOL = ',F5.0,' mm'/T15,&
            'DRAW DOWN TIME = ',F5.0,' d'/T15,'WASH WATER = ',F6.3,' m3/cow/d')
    7 FORMAT(T10,'TILE DRAIN IN SOIL LAYER',I3/T15,'DRAIN TIME = ',F5.2,&
            ' D'/T15,'FLOW RATE = ',F5.1,' mm/h')  
    8 FORMAT(T10,'NUMBER OF COWS = ',I8,' hd'/T10,'GRAZING LIMIT = ',&               
            F6.3,' t/ha'/T10,'FRACTION TIME COWS IN FEED AREA = ',F6.3)      
END SUBROUTINE Print_ManagementInfo

! --------------------------- Print crop paramters -----------------------
SUBROUTINE Print_CropParms1
      IMPLICIT NONE
      INTEGER:: I
      WRITE(KW(1),1)                                                               
      WRITE(KW(1),2)(Crop_Name(I),I=1,LC)                                               
      WRITE(KW(1),3)'RUE  ',(RUE(I),I=1,LC)        ! RUE factor                                   
      WRITE(KW(1),4)'WUE ',(WUE(I),I=1,LC)         ! WUE water use to biomass                                  
      WRITE(KW(1),4)'HI  ',(HI(I),I=1,LC)                                           
      WRITE(KW(1),3)'TOPT',(optimalTem(I),I=1,LC)                                         
      WRITE(KW(1),3)'TBAS',(plantMinT(I),I=1,LC)                                         
      WRITE(KW(1),7)'germinateHeatUnit',(germinateHeatUnit(I),I=1,LC)                                         
      WRITE(KW(1),4)'maxLAI',(maxLAI(I),I=1,LC)                                         
      WRITE(KW(1),4)'fracGrowSeasonLAIDecline',(fracGrowSeasonLAIDecline(I),I=1,LC)                                         
      WRITE(KW(1),5)'LAP1',(pointsLAIDevp(1,I),I=1,LC)                                       
      WRITE(KW(1),5)'LAP2',(pointsLAIDevp(2,I),I=1,LC)                                       
      WRITE(KW(1),5)'PPL1',(pointsPlantPopulation(1,I),I=1,LC)                                       
      WRITE(KW(1),5)'PPL2',(pointsPlantPopulation(2,I),I=1,LC)                                       
      WRITE(KW(1),5)'FRS1',(pointsFrostDamg(1,I),I=1,LC)                                       
      WRITE(KW(1),5)'FRS2',(pointsFrostDamg(2,I),I=1,LC)                                       
      WRITE(KW(1),4)'factorLAIDecline',(factorLAIDecline(I),I=1,LC)                                         
      WRITE(KW(1),4)'factorBioDecline',(factorBioDecline(I),I=1,LC)                                         
      WRITE(KW(1),3)'cropTolerantAl ',(cropTolerantAl(I),I=1,LC)                                          
      WRITE(KW(1),4)'areationThreshold ',(areationThreshold(I),I=1,LC)                                          
      WRITE(KW(1),6)'maxStomaCond ',(maxStomaCond(I),I=1,LC)                                          
      WRITE(KW(1),4)'CO2EffOnBio',(CO2EffOnBio(2,I),I=1,LC)                                       
      WRITE(KW(1),3)'VPDPara',(VPDPara(I),I=1,LC)                                         
      WRITE(KW(1),3)'VPDThreshold',(VPDThreshold(I),I=1,LC)                                         
      WRITE(KW(1),4)'VPDPara2',(VPDPara2(I),I=1,LC)                                         
      WRITE(KW(1),3)'seedRate ',(seedRate(I),I=1,LC)                                          
      WRITE(KW(1),4)'maxCropHeight ',(maxCropHeight(I),I=1,LC)                                          
      WRITE(KW(1),4)'maxRootDep',(maxRootDep(I),I=1,LC) 
      ! partitioning parameters to split biomss between aboveground and roots                                         
      WRITE(KW(1),5)'RWP1',(rootPartition(1,I),I=1,LC) ! fraction of root weight at emergence                                  
      WRITE(KW(1),5)'RWP2',(rootPartition(2,I),I=1,LC) ! fraction of root weight at maturity                                      
      WRITE(KW(1),6)'fracYieldN ',(fracYieldN(I),I=1,LC)                                          
      WRITE(KW(1),6)'fracYieldP ',(fracYieldP(I),I=1,LC)                                          
      WRITE(KW(1),6)'fracYieldK ',(fracYieldK(I),I=1,LC)                                          
      WRITE(KW(1),6)'lowerLimitHI',(lowerLimitHI(I),I=1,LC)                                         
      WRITE(KW(1),4)'pestDmagFactor ',(pestDmagFactor(I),I=1,LC)                                          
      WRITE(KW(1),4)'seedCost',(seedCost(I),I=1,LC)                                         
      WRITE(KW(1),4)'seedYieldPrice',(seedYieldPrice(I),I=1,LC)                                         
      WRITE(KW(1),4)'forageYieldPrice',(forageYieldPrice(I),I=1,LC)                                         
      WRITE(KW(1),4)'WCYS',(waterFracYield(I),I=1,LC)                                          
      WRITE(KW(1),6)'N_frac_Stage1 ',(uptakeParaN(1,I),I=1,LC)                                         
      WRITE(KW(1),6)'N_frac_Stage2 ',(uptakeParaN(2,I),I=1,LC)                                         
      WRITE(KW(1),6)'N_frac_Stage3 ',(uptakeParaN(3,I),I=1,LC)                                         
      WRITE(KW(1),6)'P_frac_Stage1 ',(uptakeParaP(1,I),I=1,LC)                                         
      WRITE(KW(1),6)'P_frac_Stage2 ',(uptakeParaP(2,I),I=1,LC)                                         
      WRITE(KW(1),6)'P_frac_Stage3 ',(uptakeParaP(3,I),I=1,LC)                                         
      WRITE(KW(1),6)'K_frac_Stage1 ',(uptakeParaK(1,I),I=1,LC)                                         
      WRITE(KW(1),6)'K_frac_Stage2 ',(uptakeParaK(2,I),I=1,LC)                                         
      WRITE(KW(1),6)'K_frac_Stage3 ',(uptakeParaK(3,I),I=1,LC)                                         
      WRITE(KW(1),5)'BW1 ',(windEroCrop(1,I),I=1,LC)                                        
      WRITE(KW(1),5)'BW2 ',(windEroCrop(2,I),I=1,LC)                                        
      WRITE(KW(1),5)'BW3 ',(windEroCrop(3,I),I=1,LC)                                        
      WRITE(KW(1),5)'Salinity_Yd1',(salinityThreshold(1,I),I=1,LC)                                        
      WRITE(KW(1),5)'Salinity_Yd2',(salinityThreshold(2,I),I=1,LC)                                        
      WRITE(KW(1),5)'LigninFrac_Plant1',(ligninFrac(1,I),I=1,LC)                                        
      WRITE(KW(1),5)'LigninFrac_Plant2',(ligninFrac(2,I),I=1,LC)                                        
      WRITE(KW(1),6)'turnoutFracCottn ',(turnoutFracCottn(I),I=1,LC)                                          
      WRITE(KW(1),6)'lintFracCottn ',(lintFracCottn(I),I=1,LC)
      WRITE(KW(1),6)'carbonEmissionSeedWeight',(carbonEmissionSeedWeight(I),I=1,LC)
      WRITE(KW(1),6)'leafWeightFrac',(leafWeightFrac(I),I=1,LC)
      WRITE(KW(1),8)(cropCode(I),I=1,LC)  
      
      1 FORMAT(//1X,'____________________CROP PARAMETERS_____________________'/)  
      2 FORMAT(10X,11(6X,A4))
      3 FORMAT(T7,A30,11F10.1)
      4 FORMAT(T7,A30,11F10.2)
      5 FORMAT(T7,A30,11F10.3)  
      6 FORMAT(T7,A30,11F10.4)
      7 FORMAT(T7,A30,11F10.0)
      8 FORMAT(T7,'cropCode ',11I10,'  cropCode')          
END SUBROUTINE Print_CropParms1

SUBROUTINE Print_CropParms2(LRG)
      IMPLICIT NONE
      ! arguments list
      INTEGER, INTENT(IN):: LRG
      ! local variables
      INTEGER:: I, J
      ! ------------------- print out variables -------------------------                                                             
      DO J=1,LRG                                                                          
            WRITE(KW(1),1)'POP ',(POP(I,J),I=1,LC)                   ! plant population                                
            WRITE(KW(1),'(T7,A30,11F10.2)')'MXLA',(PPLA(I,J),I=1,LC) ! LAI                               
      END DO     
      
      DO J=1,LRG                                                                     
            WRITE(KW(1),'(T7,A30,11F10.0)')'potentialHeatUnit ',(potentialHeatUnit(I,J),I=1,LC)  ! Potential Heat Unit                              
      END DO      
      
      ! ----------------------Print updated crop parameters --------------------
      CALL Print_Page(1)    
      
      WRITE(KW(1),2)                                                               
      WRITE(KW(1),'(10X,11(6X,A4))')(Crop_Name(I),I=1,LC)                                               
      WRITE(KW(1),1)'N_frac_Stage1 ',(uptakeParaN(1,I),I=1,LC)                                         
      WRITE(KW(1),1)'N_frac_Stage2 ',(uptakeParaN(2,I),I=1,LC)                                         
      WRITE(KW(1),1)'N_frac_Stage3 ',(uptakeParaN(3,I),I=1,LC)                                         
      WRITE(KW(1),1)'N_frac_Stage4 ',(uptakeParaN(4,I),I=1,LC)                                         
      WRITE(KW(1),1)'P_frac_Stage1 ',(uptakeParaP(1,I),I=1,LC)                                         
      WRITE(KW(1),1)'P_frac_Stage2 ',(uptakeParaP(2,I),I=1,LC)                                         
      WRITE(KW(1),1)'P_frac_Stage3 ',(uptakeParaP(3,I),I=1,LC)                                         
      WRITE(KW(1),1)'P_frac_Stage4 ',(uptakeParaP(4,I),I=1,LC)                                         
      WRITE(KW(1),1)'K_frac_Stage1 ',(uptakeParaK(1,I),I=1,LC)                                         
      WRITE(KW(1),1)'K_frac_Stage2 ',(uptakeParaK(2,I),I=1,LC)                                         
      WRITE(KW(1),1)'K_frac_Stage3 ',(uptakeParaK(3,I),I=1,LC)                                         
      WRITE(KW(1),3)'LAP1',(pointsLAIDevp(1,I),I=1,LC)                                       
      WRITE(KW(1),3)'LAP2',(pointsLAIDevp(2,I),I=1,LC)                                       
      WRITE(KW(1),3)'FRS1',(pointsFrostDamg(1,I),I=1,LC)                                       
      WRITE(KW(1),3)'FRS2',(pointsFrostDamg(2,I),I=1,LC)                                       
      WRITE(KW(1),3)'WAC1',(CO2EffOnBio(1,I),I=1,LC)                                       
      WRITE(KW(1),3)'CO2EffOnBio',(CO2EffOnBio(2,I),I=1,LC)                                       
      WRITE(KW(1),3)'PPC1',(coeffPopuCurve(1,I),I=1,LC)                                       
      WRITE(KW(1),3)'PPC2',(coeffPopuCurve(2,I),I=1,LC)    

    1 FORMAT(T7,A30,11F10.4) 
    2 FORMAT(//1X,'____________________UPDATED CROP PARAMETERS_____________________'/)
    3 FORMAT(T7,A30,11F10.3)   
END SUBROUTINE Print_CropParms2

SUBROUTINE Print_SoilTable(sulfCon, fracVPip, fracHPip)
      IMPLICIT NONE
      ! arguments list
      REAL, DIMENSION(15), INTENT(IN):: sulfCon, fracVPip, fracHPip
      ! local variables
      INTEGER:: I
      ! LINES 4/51                                                                     
      WRITE(KW(9),1)(Z(Layer_ID(I)),I=1,MSL),'  4 DEPTH(m)        '
      WRITE(KW(9),1)(bulkDensity(Layer_ID(I)),I=1,MSL),'  5 bulkDensity 33kpa(t/m3)  '                                   
      WRITE(KW(9),1)(SOIL(20,Layer_ID(I)),I=1,MSL),'  6 WP totSoilWater(m/m)       '                                       
      WRITE(KW(9),1)(SOIL(9,Layer_ID(I)),I=1,MSL),'  7 fieldCapacity totSoilWater(m/m)      '                                   
      WRITE(KW(9),1)(sandFrac(Layer_ID(I)),I=1,MSL),'  8 SAND(%)         '                                  
      WRITE(KW(9),1)(siltFrac(Layer_ID(I)),I=1,MSL),'  9 SILT(%)         '                                   
      WRITE(KW(9),2)(SOIL(6,Layer_ID(I)),I=1,MSL),' 10 ORG N(g/t)      '                                   
      WRITE(KW(9),1)(PH(Layer_ID(I)),I=1,MSL),' 11 PH              '                                   
      WRITE(KW(9),1)(totBase(Layer_ID(I)),I=1,MSL),' 12 SM BS(cmol/kg)  '                                   
      WRITE(KW(9),1)(SOIL(7,Layer_ID(I)),I=1,MSL),' 13 ORG C(%)        '                                   
      WRITE(KW(9),1)(CaCO3(Layer_ID(I)),I=1,MSL),' 14 CaCO3(%)          '                                   
      WRITE(KW(9),1)(CEC(Layer_ID(I)),I=1,MSL),' 15 CEC(cmol/kg)    '                                   
      WRITE(KW(9),1)(rockFrac(Layer_ID(I)),I=1,MSL),' 16 ROCK(%)         '                                   
      WRITE(KW(9),1)(SOIL(5,Layer_ID(I)),I=1,MSL),' 17 NO3(g/t)        '                                   
      WRITE(KW(9),1)(SOIL(1,Layer_ID(I)),I=1,MSL),' 18 LAB P(g/t)      '                                   
      WRITE(KW(9),1)(cropResidu(Layer_ID(I)),I=1,MSL),' 19 CROP cropResidu(t/ha)  '                                   
      WRITE(KW(9),1)(SOIL(13,Layer_ID(I)),I=1,MSL),' 20 bulkDensity DRY(t/m3)    '                                   
      WRITE(KW(9),1)(sorbRatioP(Layer_ID(I)),I=1,MSL),' 21 P SORP RTO       '                                   
      WRITE(KW(9),1)(satuateCondv(Layer_ID(I)),I=1,MSL),' 22 SAT COND(mm/h)   '                                   
      WRITE(KW(9),1)(lateralFlowCondv(Layer_ID(I)),I=1,MSL),' 23 H SC(mm/h)       '                                        
      WRITE(KW(9),1)(SOIL(4,Layer_ID(I)),I=1,MSL),' 24 ORG P(g/t)       '                                   
      WRITE(KW(9),1)(SOIL(15,Layer_ID(I)),I=1,MSL),' 25 exchangeK(g/t)        '                                   
      WRITE(KW(9),1)(elecCondv(Layer_ID(I)),I=1,MSL),' 26 elecCondv(mmho/cm)   '                                   
      WRITE(KW(9),1)(fracNO3Leach(Layer_ID(I)),I=1,MSL),' 27 fracNO3Leach            '                                        
      WRITE(KW(9),1)(SOIL(12,Layer_ID(I)),I=1,MSL),' 28 totSoilWater(m/m)         '                                   
      WRITE(KW(9),2)(fracVPip(Layer_ID(I)),I=1,MSL),' 29 FRACT V PIPE F  '                                        
      WRITE(KW(9),2)(fracHPip(Layer_ID(I)),I=1,MSL),' 30 FRACT H PIPE F  '                                        
      WRITE(KW(9),2)(structLitt(Layer_ID(I)),I=1,MSL),' 31 structLitt(kg/ha)      '                                   
      WRITE(KW(9),2)(metabLitt(Layer_ID(I)),I=1,MSL),' 32 metabLitt(kg/ha)      '                                   
      WRITE(KW(9),1)(lgStructLitt(Layer_ID(I)),I=1,MSL),' 33 lgStructLitt(kg/ha)     '                                   
      WRITE(KW(9),1)(CStructLitt(Layer_ID(I)),I=1,MSL),' 34 CStructLitt(kg/ha)     '                                   
      WRITE(KW(9),1)(CMetabLitt(Layer_ID(I)),I=1,MSL),' 35 CMetabLitt(kg/ha)     '                                   
      WRITE(KW(9),1)(CLgStructLitt(Layer_ID(I)),I=1,MSL),' 36 CLgStructLitt(kg/ha)    '                                   
      WRITE(KW(9),1)(NLgStructLitt(Layer_ID(I)),I=1,MSL),' 37 NLgStructLitt(kg/ha)   '                                   
      WRITE(KW(9),1)(CBiomass(Layer_ID(I)),I=1,MSL),' 38 CBiomass(kg/ha)     '                                   
      WRITE(KW(9),2)(CSlowHumus(Layer_ID(I)),I=1,MSL),' 39 CSlowHumus(kg/ha)     '                                   
      WRITE(KW(9),2)(CPassiveHumus(Layer_ID(I)),I=1,MSL),' 40 CPassiveHumus(kg/ha)     '                                   
      WRITE(KW(9),1)(NStructLitt(Layer_ID(I)),I=1,MSL),' 41 NStructLitt(kg/ha)     '                                   
      WRITE(KW(9),1)(NMetabLitt(Layer_ID(I)),I=1,MSL),' 42 NMetabLitt(kg/ha)     '                                   
      WRITE(KW(9),1)(NBiomass(Layer_ID(I)),I=1,MSL),' 43 NBiomass(kg/ha)     '                                   
      WRITE(KW(9),2)(NSlowHumus(Layer_ID(I)),I=1,MSL),' 44 NSlowHumus(kg/ha)     '                                   
      WRITE(KW(9),2)(NPassiveHumus(Layer_ID(I)),I=1,MSL),' 45 NPassiveHumus(kg/ha)     '
      WRITE(KW(9),1)(ironCon(Layer_ID(I)),I=1,MSL),' 46 IRON(%)         '
      WRITE(KW(9),1)(sulfCon(Layer_ID(I)),I=1,MSL),' 47 SULPHUR(%)      '
      WRITE(KW(9),'(15A8,A20)')(soilHorizon(Layer_ID(I)),I=1,MSL),' 48 SOIL HORIZON    '
      WRITE(KW(9),1)(gasO2Con(Layer_ID(I)),I=1,MSL),' 49 gasO2Con(kg/ha)     '
      WRITE(KW(9),1)(gasCO2Con(Layer_ID(I)),I=1,MSL),' 50 gasCO2Con(kg/ha)    '
      WRITE(KW(9),1)(gasN2OCon(Layer_ID(I)),I=1,MSL),' 51 gasN2OCon(kg/ha)    '
   
    1 FORMAT(15F8.2,A20) 
    2 FORMAT(15F8.0,A20)  
END SUBROUTINE Print_SoilTable

SUBROUTINE Print_SCN(XYR, XZP)
      IMPLICIT NONE
      ! arguments list
      INTEGER, INTENT(IN):: XYR
      REAL, DIMENSION(13,16), INTENT(IN):: XZP
      ! local variables
      INTEGER:: I, J, MS1
      REAL, DIMENSION(MSL + 1):: XTP, XYP, YTP

      ! /////////////////////////////////////////////////////////////////
      DO J=1,6                                                                       
            XTP(J)=0.                                                                    
            DO I=1,Actual_SoilLayers                                                                  
                  ISL=Layer_ID(I)                                                               
                  SMS(J,ISL)=SMS(J,ISL)/(SMS(11,ISL)+1.E-5)                                
                  XTP(J)=XTP(J)+SMS(J,ISL)                                                 
            END DO                                                                            
      END DO                                                                         
      DO J=7,10                                                                  
            XTP(J)=0.                                                                      
            DO I=1,Actual_SoilLayers                                                                    
                  ISL=Layer_ID(I)                                                                   
                  SMS(J,ISL)=SMS(J,ISL)/XYR                                                    
                  XTP(J)=XTP(J)+SMS(J,ISL)                                                     
            END DO
      END DO                                                                                   
      DO I=1,Actual_SoilLayers                                                                    
            ISL=Layer_ID(I)                                                                   
            XYP(ISL)=SOC(ISL)-XZP(6,ISL)                                                 
            YTP(ISL)=SON(ISL)-XZP(12,ISL)                                                
      END DO                                                                              
      XYP(16)=totSOC-XZP(6,16)                                                          
      YTP(16)=totSON-XZP(12,16)                                                         
      WRITE(KW(16),3)(SID(J),J=1,MSL),SID(MSL+1)                        
                                                                        
      DO J=1,6                                                                       
            XTP(J)=XTP(J)/Actual_SoilLayers                                                             
      END DO                                                                              
      MS1=MSL+1                                                                           
      WRITE(KW(16),1)'   Z',(Z(Layer_ID(I)),I=1,MSL),Z(Layer_ID(MSL))                        
      WRITE(KW(16),1)' SWF',(SMS(1,Layer_ID(I)),I=1,MSL),XTP(1)                         
      WRITE(KW(16),1)'TEMP',(SMS(2,Layer_ID(I)),I=1,MSL),XTP(2)                         
      WRITE(KW(16),1)'SWTF',(SMS(3,Layer_ID(I)),I=1,MSL),XTP(3)                         
      WRITE(KW(16),1)'TLEF',(SMS(4,Layer_ID(I)),I=1,MSL),XTP(4)                         
      WRITE(KW(16),1)'SPDM',(SMS(5,Layer_ID(I)),I=1,MSL),XTP(5)                         
      WRITE(KW(16),2)'RSDC',(SMS(7,Layer_ID(I)),I=1,MSL),XTP(7)                         
      WRITE(KW(16),2)'Soil_Respiration',(SMS(8,Layer_ID(I)),I=1,MSL),XTP(8)                         
      WRITE(KW(16),2)'XNetNMineralize',(SMS(9,Layer_ID(I)),I=1,MSL),XTP(9)                         
      WRITE(KW(16),2)'DNO3',(SMS(10,Layer_ID(I)),I=1,MSL),XTP(10)                       
      WRITE(KW(16),2)'HSC0',(XZP(1,Layer_ID(I)),I=1,MSL),XZP(1,MS1)                     
      WRITE(KW(16),2)'HSCF',(CSlowHumus(Layer_ID(I)),I=1,MSL),totCSlowHumus                            
      WRITE(KW(16),2)'HPC0',(XZP(2,Layer_ID(I)),I=1,MSL),XZP(2,MS1)                     
      WRITE(KW(16),2)'HPCF',(CPassiveHumus(Layer_ID(I)),I=1,MSL),totCPassiveHumus                            
      WRITE(KW(16),2)'LSC0',(XZP(3,Layer_ID(I)),I=1,MSL),XZP(3,MS1)                     
      WRITE(KW(16),2)'LSCF',(CStructLitt(Layer_ID(I)),I=1,MSL),totCStructLitt                            
      WRITE(KW(16),2)'LMC0',(XZP(4,Layer_ID(I)),I=1,MSL),XZP(4,MS1)                     
      WRITE(KW(16),2)'LMCF',(CMetabLitt(Layer_ID(I)),I=1,MSL),totCMetabLitt                            
      WRITE(KW(16),2)'BMC0',(XZP(5,Layer_ID(I)),I=1,MSL),XZP(5,MS1)                     
      WRITE(KW(16),2)'BMCF',(CBiomass(Layer_ID(I)),I=1,MSL),totCBiomass                            
      WRITE(KW(16),2)'Org_C_con0',(XZP(6,Layer_ID(I)),I=1,MSL),XZP(6,MS1)                     
      WRITE(KW(16),2)'WOCF',(SOC(Layer_ID(I)),I=1,MSL),totSOC                              
      WRITE(KW(16),2)'DWOC',(XYP(Layer_ID(I)),I=1,MSL),XYP(MS1)                         
      WRITE(KW(16),2)'HSN0',(XZP(7,Layer_ID(I)),I=1,MSL),XZP(7,MS1)                     
      WRITE(KW(16),2)'HSNF',(NSlowHumus(Layer_ID(I)),I=1,MSL),totNSlowHumus                            
      WRITE(KW(16),2)'HPN0',(XZP(8,Layer_ID(I)),I=1,MSL),XZP(8,MS1)                     
      WRITE(KW(16),2)'HPNF',(NPassiveHumus(Layer_ID(I)),I=1,MSL),totNPassiveHumus                            
      WRITE(KW(16),2)'LSN0',(XZP(9,Layer_ID(I)),I=1,MSL),XZP(9,MS1)                     
      WRITE(KW(16),2)'LSNF',(NStructLitt(Layer_ID(I)),I=1,MSL),totNStructLitt                            
      WRITE(KW(16),2)'LMN0',(XZP(10,Layer_ID(I)),I=1,MSL),XZP(10,MS1)                   
      WRITE(KW(16),2)'LMNF',(NMetabLitt(Layer_ID(I)),I=1,MSL),totNMetabLitt                            
      WRITE(KW(16),2)'BMN0',(XZP(11,Layer_ID(I)),I=1,MSL),XZP(11,MS1)                   
      WRITE(KW(16),2)'BMNF',(NBiomass(Layer_ID(I)),I=1,MSL),totNBiomass                            
      WRITE(KW(16),2)'Org_Nit_con0',(XZP(12,Layer_ID(I)),I=1,MSL),XZP(12,MS1)                   
      WRITE(KW(16),2)'Org_Nit_conF',(SON(Layer_ID(I)),I=1,MSL),totSON                              
      WRITE(KW(16),2)'DOrg_Nit_con',(YTP(Layer_ID(I)),I=1,MSL),YTP(MS1)                         
      WRITE(KW(16),1)'C/N0',(XZP(13,Layer_ID(I)),I=1,MSL),XZP(13,MS1)                   
      DO I=1,16                                                                    
            XTP(I)=0.                                                                    
      END DO                                                                              
      DO I=1,Actual_SoilLayers                                                                    
            ISL=Layer_ID(I)                                                                   
            XTP(ISL)=SOC(ISL)/SON(ISL)                                                   
      END DO                                                                              
      XTP(MS1)=totSOC/totSON                                                               
      WRITE(KW(16),1)'C/NF',(XTP(Layer_ID(I)),I=1,MSL),XTP(MSL+1)    

    1 FORMAT(1X,A4,20F10.3)
    2 FORMAT(1X,A4,20F10.0) 
    3 FORMAT(///T52,'SOIL LAYER NO'/T5,16(6X,A4))   
END SUBROUTINE Print_SCN

END MODULE Print_Module