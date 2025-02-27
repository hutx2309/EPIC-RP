MODULE Management
USE PARM
USE Soil_Module
USE MisceFun_Module
USE ReadVarTable_Module
USE Print_Module
IMPLICIT NONE
CONTAINS
SUBROUTINE Management_Info(IY1, KRX,  operationID, JRT, soilHydGroup, adjustCN2)
      IMPLICIT NONE
      ! Argument list
      INTEGER, INTENT(IN):: operationID, KRX
      INTEGER, INTENT(INOUT):: IY1, JRT
      REAL, INTENT(INOUT):: soilHydGroup, adjustCN2

      ! local variables declaration
      INTEGER:: I, II, IJ, I1, I3, J, J1, J2, JJ, K, K1, KK, L, M, NCRP
      REAL:: BASE, HU0, X1, X2, X4, XY, XZ
      ! ///////////////////////////////////////////
      ! ------------------------------- Select operation file --------------------------
      II=-1
      DO WHILE(II/=operationID)  
            READ(KR(15),*,IOSTAT=NFL)II,OPSCFILE  
            IF(NFL/=0)THEN
                  IF(runMode==0)THEN
                        WRITE(*,*)'OPS NO = ',operationID,' NOT IN OPSC LIST FILE'
                        ERROR STOP  
                  ELSE
                        WRITE(KW(MSO),'(A,A8,A,I4,A)')' !!!!! ',siteName,&
                              ' OPS NO = ',operationID,' NOT IN OPSC LIST FILE'
                        ERROR STOP
                  END IF
            END IF 
      END DO  
      REWIND KR(15)     
      
      ! ----------------------- read operation file ----------------------------
      CALL OPENV(KR(16),OPSCFILE,IDIR(4),KW(MSO))  
      !     LINE 1      
      READ(KR(16),'()')
      !     LINE 2  
      READ(KR(16),'(2I4)')landUseCode, autoIrrCode  ! landUseCode: Land use number  Ref EPIC0810  P36  
      ! autoIrrCode (1 - ...): apply irrigation from Till.dat. If auto irrigation is used, this irrigation will be used 
      !              to apply irrigation water. If none is specified, the default is operation #500 
      IF(autoIrrCode==0)autoIrrCode=500 
      hydroGroup=INT(soilHydGroup) 
      CALL CN2_Table                     ! read the SCS HYDROLOGIC SOIL GROUP-CURVE NUMBER TABLE   
      CALL Adjust_CN(CN2, adjustCN2)  
 
      CN2=adjustCN2  
      CN0=CN2  
      WRITE(KW(1),'(/1X,A)')'-----OPERATION SCHEDULE---------------'    

      HU(1)=0.                       ! heat unit?  
      IBGN=1                         ! Day of year?  
      Crop_Num=1      
      I=1                            ! I is the year now   
      K1=1  
  189 NCRP=IGO                       ! What are these variables?   How was IGO increased?                                                  
      I1=I-1  
      K=1   
      WRITE(KW(1),'(/T10,A,I4)') 'YR', I           

      numOfIrr=0                 ! counting times of irrigation  
      numOfFert=0                ! counting times of fertilization 
      numOfPest=0                ! counting times of pesticide 
      HU0=0.   
      IF(JDHU<366) HU(Crop_Num)=0. 
      IF(IGO>0)THEN              ! How is IGO increased to 1?
            DO J=1,MNC           ! MNC = 30  
                  IF(NHU(J)==0)CYCLE                                                           
                  LY(I,K)=NHU(J) 
                  K=K+1                                                                        
            END DO                                                                              
            J=0
      END IF  
      WRITE(KW(1),1)      
      !------------------------------ continue to read operation schedule ----------------------
      J=0   
      Operation_Schedule: DO       ! the DO ... END DO  loop read operation schedule year by year
            J=J+1                  ! J is the number of operations
            !     READ OPERATION SCHEDULE                                                        
            !  1  JX(1)= YR OF OPERATION                                                         
            !  2  JX(2)= MO OF OPERATION                                                         
            !  3  JX(3)= DAY OF OPERATION                                                        
            !  4  JX(4)= EQUIPMENT ID NO                                                         
            !  5  JX(5)= TRACTOR ID NO                                                           
            !  6  JX(6)= CROP ID NO                                                              
            !  7  JX(7)= XMTU--TIME FROM PLANTING TO MATURITY (Y)(FOR TREE CROPS AT PLANTING ONLY).                                                
            !          = TIME FROM PLANTING TO HARVEST (Y)(HARVEST ONLY)                         
            !          = PESTICIDE ID NO (FOR PESTICIDE APPLICATION ONLY)                        
            !          = FERTILIZER ID NO (FOR FERTILIZER APPLICATION ONLY)                      
            !  8  OPV1 = POTENTIAL HEAT UNITS FOR PLANTING (BIU)                    
            !          = APPLICATION VOLUME (mm)FOR IRRIGATION                                   
            !          = FERTILIZER APPLICATION RATE (kg/ha) = 0 FOR VARIABLE RATE
            !          = PESTICIDE APPLICATION RATE (g/ha)                                      
            !          = LIME APPLICATION RATE (t/ha)               
            !  9  OPV2 = LINE NUMBER FOR SCS HYDROLOGIC SOIL GROUP/RUNOFF CURVE NUMBER TABLE                                                            
            !          = SCS CURVE NUMBER (NEGATIVE)                                             
            !          = PEST CONTROL FACTOR FOR PEST APPLICATION (FRACTION OF PESTS CONTROLLED)                                                             
            !  10 OPV3 = PLANT WATER STRESS FACTOR(0-1); SOIL WATER TENSION(>1kpa);              
            !            OR PLANT AVAILABLE WATER DEFICIT IN ROOT ZONE(-mm)TO                    
            !            TRIGGER AUTO IRR. (0. OR BLANK DOES NOT CHANGE TRIGGER)                 
            !  11 OPV4 = RUNOFF VOL/VOL IRRIGATION WATER APPLIED                                 
            !  12 OPV5 = PLANT POPULATION (PLANTS/M**2 OR PLANTS/HA IF p/m2<1.) (FOR PLANTING ONLY)                                                     
            !  13 OPV6 = MAX ANNUAL N FERTILIZER APPLIED TO A CROP (0. OR BLANK                  
            !            DOES NOT CHANGE FMX; > 0 SETS NEW FMX)(FOR PLANTING ONLY                
            !  14 OPV7 = TIME OF OPERATION AS FRACTION OF GROWING SEASON (ENTER                  
            !            EARLIEST POSSIBLE MO & DAY -- JX(2) & JX(3))
            !  15 OPV8 = MINIMUM USLE C FACTOR
            !  16 OPV9 = MOISTURE CONTENT OF GRAIN REQUIRED FOR HARVEST                            
            !     LINE 3/L                                                                       
            IF(I==1.OR.J>1)READ(KRX,2,IOSTAT=NFL)JX,(OPV(L),L=1,9)  
            IF(NFL/=0)EXIT
            IF(I==1.AND.J==1)JJK0=JX(6)       ! Crop ID  
            IF(JX(1)/=IY1)EXIT                 
            CALL Tillage_Table                ! read tillage from the tillage operation table, return NDT and irrigation operation  
            LT(I,J)=NDT                       ! Initialized in AINLZ, modified in Tillage_Table  
            JH(I,J)=JX(6)    
            IJ=LT(I,J)                        ! Is it necessary to use IJ?  
            tillDOY(I,J)= Cal_DOY(JX(2),JX(3), 1)! calculate DOY given the month and the day of the month of operation   
            X4=tillDepth(IJ)*1000.            ! tillDepth: Tillage Depth (mm), given in Tillage_Table  
            I3=currentOps(IJ)                 ! currentOps: operation code  
            IF(IBGN<tillDOY(I,J))GO TO 419   
            IF(IBGN==tillDOY(I,J))GO TO 422  
            IF(IGO>0) HU(Crop_Num)=HU(Crop_Num)+ Cal_Heat(IBGN,365,BASE,0)/(potentialHeatUnit(Crop_Num,IHU(Crop_Num))+1.)                       
                  !???? questionable Base has no value? How is IGO changed?                                                                  
            IBGN=1                                                                         
        419 IF(IGO>0) HU(Crop_Num)=HU(Crop_Num)+ Cal_Heat(IBGN,tillDOY(I,J),BASE,0)/(potentialHeatUnit(Crop_Num,IHU(Crop_Num))+1.) 
          
            HU0=HU0 + Cal_Heat(IBGN,tillDOY(I,J),0.,1)/yearTotHeat   ! cal accumulated heat when the operation is conducted                                      
            IBGN=tillDOY(I,J)  
        422 IF(OPV(7)>0.)GO TO 420 
             ! Operation_Flag is read from EPICCONT.DAT   0: normal operation, 1: automatic heat unit schedule                                                          
            IF(operationMode==0)GO TO 420  
            IF(IGO==0)GO TO 441  
            IF(cropCode(Crop_Num)==plantCategoryCode(1).OR.cropCode(Crop_Num)==plantCategoryCode(2).OR.&
            cropCode(Crop_Num)==plantCategoryCode(4).OR.cropCode(Crop_Num)==plantCategoryCode(5).OR.&
            cropCode(Crop_Num)==plantCategoryCode(9))GO TO 423  
        441 fracHeatUnit(I,J)=HU0 
            GO TO 421  
        423 fracHeatUnit(I,J)=HU(Crop_Num)  
            GO TO 421                                                                      
        420 fracHeatUnit(I,J)=OPV(7)                                                               
        421 CALL Operation_Inputs(I3,I,J,JRT)                                                         
            X1=MAX(0.,COOP(IJ))              ! for economic analysis                                                     
            X2=MAX(0.,COTL(IJ))     
          
            ! PRINTOUT OPERATION SCHEDULE  
          
            WRITE(KW(1),3)I,JX(2),JX(3),equipmentName(IJ),I3,JX(4),JX(5),X1,X2,mixEfficiency(IJ), &            
                        tillRoughness(IJ),X4,soilCompactFrac(IJ),ridgeHeight(IJ),ridgeInterval(IJ), &
                        furrowHeight(IJ),furrowInterval(IJ),harvestEffi(IJ),overRideHI(IJ),CN2, &
                        irrTrigger,IrrRunoffRatio,fracHeatUnit(I,J)                                                     
            X1=0.                                                                          
            IF(tillDepth(IJ)>BIG)BIG=tillDepth(IJ)                                                     
            IF(I3/=operationCode(5).AND.I3/=operationCode(6))GO TO 180 ! 5: plant in rows 6: plant with drill ?                              
          
            NCRP=NCRP+1          ! number of crops?                                                             
            IGO=IGO+1            ! growing times
          
            CALL Read_CropTable                  ! read plant/crop parameters  
          
            NBCX(I,Crop_Num)=NBCX(I,Crop_Num)+1  ! What is NBCX?                                                      
            BASE=plantMinT(Crop_Num)                                                               
            IPL=tillDOY(I,J)+365*I1              ! IPL: operation date of the next year                                                      
            IPLD(1,Crop_Num)=IPL                                                                
            LY(I,K)=Crop_Num                                                                    
            NHU(Crop_Num)=Crop_Num                                                                   
            K=K+1                                                                          
            X1=seedRate(Crop_Num)*seedCost(Crop_Num)                                                          
            X2=X1+X2  
          
            WRITE(KW(1),4)X1,Crop_Name(Crop_Num) ! X1 is cost here.
          
            GO TO 584                                                                      
        180 IF(I3/=operationCode(1).AND.I3/=operationCode(2).AND.I3/=operationCode(3))GO TO 584                          
            IF(varID(JX(6))==0)GO TO 584                                                    
            Crop_Num= varID(JX(6))              ! confusing                                                       
            IF(I3==operationCode(1))THEN
                  IHV=tillDOY(I,J)+365*I1                                                            
                  IHVD(1,Crop_Num)=IHV                                                                
                  NHU(Crop_Num)=0                                                                     
                  IGO=MAX(0,IGO-1)              ! IGO was modified here                                                        
                  IF(NBCX(I,Crop_Num)==0) NBCX(I,Crop_Num)=1
            END IF                                                     
            HU(Crop_Num)=0.                                                                     
            LYR(I,J)=MAX(1,JX(7))                                                          
            WRITE(KW(1),4)HU(Crop_Num),Crop_Name(Crop_Num)                                              
        584 IF(KFL(13)>0)WRITE(KW(13),5)I,JX(2),JX(3),equipmentName(IJ),JX(6),I3,JX(4),JX(5),X2,COOP(IJ),X1
      END DO Operation_Schedule     
      
      tillDOY(I,J)=367                                                                   
      fracHeatUnit(I,J)=10.   ! It should be with in 1~1.2? Why the initial is 10?                                                                
      NBC(I)=NCRP                                                                 
      J1=J-1                                                                              
      numOfTill(I)=J1                                                                      
      numOfPestArr(I)=numOfPest                                                                     
      LT(I,J)=1                                                                      
      JH(I,J)=0                                                                      
      CND(I,J)=CN2                                                                   
      QIR(I,J)=IrrRunoffRatio                                                                   
      irrTriggerArr(I,J)=irrTrigger                                                                   
      fertCFactirArr(I,J)=minCFactor
      HWC(I,J)=0. 
      baseStockRateArr(I,J)=baseStockRate                  
      
      ! ------------------- print Irrigation schedule ---------------------
      MO=1        
      IF(numOfIrr==0) GO TO 185    
      
      WRITE(KW(1),6)    
      
      DO L=1,J1                                                                  
            J2=LT(I,L)                                                                     
            IF(currentOps(J2)/=operationCode(8))CYCLE
            DayOfYear=tillDOY(I,L)                                                                   
            IF(leapYr==0)DayOfYear=DayOfYear+1                                                            
            CALL Cal_Mon(DayOfYear,MO)                                                             
            DayOfMon = Cal_DOM(DayOfYear, MO, leapYr)                                                                    
            XZ=irrCost*irrVolumeArr(I,L)                                                              
            XY=MAX(0.,COTL(J2))+XZ                                                         
            !     PRINTOUT IRRIGATION SCHEDULE                                                   
            WRITE(KW(1),7)I,MO,DayOfMon,irrVolumeArr(I,L),XY,fracHeatUnit(I,L)                                
            IF(KFL(13)>0)WRITE(KW(13),8)I,MO,DayOfMon,JH(I,L),currentOps(J2),XY,XZ,&                 
                             irrVolumeArr(I,L)
      END DO  
      
      ! ---------------------------print fertilizer schedule ------------------------
      MO=1                                                                           
  185 IF(numOfFert==0)GO TO 187  
      WRITE(KW(1),9)  
      JJ=367   
      KK=0  
      
      DO L=1,J1   
            J2=LT(I,L) 
            IF(currentOps(J2)/=operationCode(9))CYCLE
            X1=MAX(0.,COTL(J2))   
            DayOfYear=tillDOY(I,L)  
            IF(leapYr==0)DayOfYear=DayOfYear+1   
            IF(DayOfYear==JJ.AND.NBT(J2)==0.AND.NBE(J2)==KK)X1=0.
            JJ=DayOfYear  
            KK=NBE(J2)    
            CALL Cal_Mon(DayOfYear,MO)   
            DayOfMon = Cal_DOM(DayOfYear, MO, leapYr)  
            M=LFT(I,L)  
            XZ=FCST(M)*WFA(I,L)  
            XY=X1+XZ   
            ! PRINTOUT FERTILIZER SCHEDULE                                                   
            WRITE(KW(1),10)I, MO, DayOfMon, fertName(M), KDF(M), NBE(J2), NBT(J2),XY,WFA(I,L),tillDepth(J2), &
                        FN(M), FNH3(M), FNO(M), FP(M), FPO(M), FK(M), fracHeatUnit(I,L)                
            IF(KFL(13)>0)WRITE(KW(13),11)I, MO, DayOfMon, fertName(M), JH(I,L), KDF(M), currentOps(J2), &
                              NBE(J2), NBT(J2), XY, XZ, WFA(I,L)                                         
      END DO
      ! --------------------------------print pesticide schedule------------------------------ 
      MO=1   
  187 JJ=367 
      KK=0  
      IF(numOfPest>0)THEN
            WRITE(KW(1),12)  
            DO L=1,J1                                                                  
                  J2=LT(I,L)  
                  IF(currentOps(J2)/=operationCode(7))CYCLE
                  X1=MAX(0.,COTL(J2))  
                  DayOfYear=tillDOY(I,L)  
                  IF(leapYr==0)DayOfYear=DayOfYear+1  
                  IF(DayOfYear==JJ.AND.NBT(J2)==0.AND.NBE(J2)==KK)X1=0.
                  JJ=DayOfYear  
                  KK=NBE(J2)   
                  CALL Cal_Mon(DayOfYear,MO)  
                  DayOfMon = Cal_DOM(DayOfYear, MO, leapYr)  
                  M=LPC(I,L)                                                                     
                  XZ=Pest_Cost(M)*pestRate(I,L)                                                           
                  XY=X1+XZ                                                                       
                  !     PRINTOUT PESTICIDE SCHEDULE                                                    
                  WRITE(KW(1),13)I,MO,DayOfMon,pestName(M),KDP(M),NBE(J2),NBT(J2),XY,pestRate(I,L),pestKillEffi(I,L), &
                              fracHeatUnit(I,L)                                                      
                  IF(KFL(13)>0)WRITE(KW(13),14)I, MO, DayOfMon, pestName(M), JH(I,L), KDP(M), currentOps(J2), &
                              NBE(J2),NBT(J2),XY,XZ,pestRate(I,L)                                        
            END DO         
      END IF
      
      IF(NFL==0.AND.JX(1)>0)THEN
            I=JX(1)  
            IY1=I    
            GO TO 189 
      END IF

      REWIND KR(16) !=================== End of reading .OPC ============================================  

      ! =============== print auto irrigation, lime, fertilizer information if they were scheduled =======  

      cropRotationYrs=IY1             ! crop rotation duration   yr (1-10)  
      IGSD=0    
      
      IF(DOY_realtime>0.AND.DOY_realtime<366)IGSD=cropRotationYrs 
      JX(4)=autoIrrCode 
      JX(5)=0
      
      CALL Tillage_Table                    ! return NDT and irrigation operation                                                
      autoIrrCode=NDT 
      WRITE(KW(1), 15)equipmentName(NDT),tillDepth(NDT)  
      ! IF(IAUF==0) GO TO 689                                                           
      JX(4)=261 
      JX(5)=12  
      CALL Tillage_Table                    ! return NDT and fertilizer operation                                                            
      IAUF=NDT     
      
      IF(limeFlag==0)THEN                   ! if lime operation is not specified, then read from the TILL table.
            JX(4)=267  
            JX(5)=12  
            CALL Tillage_Table  
            IAUL=NDT                                                                       
            WRITE(KW(1), 16)equipmentName(NDT),tillDepth(NDT)
      END IF
      
      L=1    
      
      IF(fertTrigger0>0.)THEN 
            IDFT(1)=fertTimes    
            IF(IDFT(1)==0)THEN  
                  IDFT(1)=52                                                                   
                  IDFT(2)=52                                                                   
            ELSE                                                                           
                  IDFT(2)=fertTimes                                                                 
            END IF
            DO K=1,2                                                                       
                  JX(7)=IDFT(K)                                                                
                  CALL Fert_Table(L)                                                              
                  IDFT(K)=L                                                                    
            END DO                                                                         
          WRITE(KW(1),17)equipmentName(IAUF),tillDepth(IAUF),fertName(IDFT(1))
      END IF
      
      IF(P_apply_mode>0)THEN
            IF(numAutoP==0) numAutoP=53
            JX(7)=numAutoP
            CALL Fert_Table(L)
            numAutoP=L
      END IF  

    1 FORMAT(/T10,'TILLAGE OPERATIONS'/6X,'DATE',T21,'stableMineralP',4X,'EQ',4X,'TR',3X, &
            'OCST',3X,'TCST',2X,'SCST',5X,'MX',4X,'RR',4X,'DP',4X,'FR',3X,'RHT',3X,'RIN',3X,&
            'furrowHeight',3X,'furrowInterval',4X,'HV',4X,'HV',8X,'NRCS', 4X,'IRR',5X,'Q/', /6X, &
            'Y M D',3X,'NAME',2X,'CD',4X,'NO',4X,'NO',1X,'________$/ha________',4X,'EF',4X,'mm',4X,'mm',2X,'COMP',&
            4X,'mm',5X, 'm',4X,'mm',5X,'m',4X,'EF',3X,'IDX',2X,'CROP',4X,'CN',3X,'TRGR',3X,'irrVolumeArr',2X,'fracHeatUnit')
    2 FORMAT(3I3,4I5,10F8.0)
    3 FORMAT(4X,I3,2I2,1X,A8,I2,2I6,2F7.2,7X,F6.2,2F6.0,F6.3,F6.0,F6.2,F6.0, &            
            2F6.2,F6.3,5X,F7.1,2F7.2,2F6.3)
    4 FORMAT('+',T50,F7.2,T118,A4) 
    5 FORMAT(1X,3I2,2X,A8,8X,I6,6X,3I4,3F7.2)
    6 FORMAT(/T10,'IRRIGATION WATER APPLIED'/T7,'DATE',T14,'VOL',3X,'COST',3X,'HU'/T7,'Y M D   mm   $/ha   SCD')
    7 FORMAT(5X,3I2,F6.0,2F6.2) 
    8 FORMAT(1X,3I2,2X,'IRGA',12X,I6,6X,I4,8X,F7.2,7X,F7.2,F7.0)     
    9 FORMAT(/T10,'FERTILIZER APPLIED'/T7,'DATE',T25,'FT',4X,'EQ',4X, 'TR',3X,'COST',4X,'WT',2X,'DPTH',1X, &
            '-------------FRACTION OF WT--------------'/T7,'Y M D',3X,'NAME',6X,'NO',4X,'NO',4X,'NO',3X, &                
            '$/ha',1X,'kg/ha',5X,'m', 4X,'MN',3X,'NH3',4X,'ON',4X,'MP',4X,'stableMineralP',&             
            4X,'MK',2X,'fracHeatUnit')
   10 FORMAT(5X,3I2,1X,A8,3I6,F7.2,F6.0,8F6.3) 
   11 FORMAT(1X,3I2,2X,A8,8X,I6,2X,4I4,F7.2,7X,F7.2,F7.0) 
   12 FORMAT(/T10,'PESTICIDES APPLIED'/T7,'DATE',T33,'PS',4X,'EQ',4X,&               
            'TR',3X,'COST', 2X,'RATE',2X,'KILL'/T7,'Y M D   NAME',14X,'NO',4X,&  
            'NO',4X,'NO',3X, '$/ha',1X,'kg/ha',4X,'EF',2X,'fracHeatUnit')
   13 FORMAT(5X,3I2,1X,A16,3I6,F7.2,3F6.2) 
   14 FORMAT(1X,3I2,2X,A16,I6,2X,4I4,F7.2,7X,F7.2,F7.2)                     
   15 FORMAT(/T10,'AUTO IRR EQUIP  = ',A8,2X,'DEPTH = ',F6.3,' M')
   16 FORMAT(T10,'AUTO LIME EQUIP = ',A8,2X,'DEPTH = ',F6.3,' M')                    
   17 FORMAT(T10,'AUTO FERT EQUIP = ',A8,2X,'DEPTH = ',F6.3,' M',2X,&                
            'FERT = ',A8/)                            
END SUBROUTINE Management_Info

! ************************** I Residual Operation *********************************  
!--------------------I.1. Burn all standing crop residue ----------------------------
SUBROUTINE Burn_Residue                           
      IMPLICIT NONE
      !     EPIC1102
      !     THIS SUBPROGRAM BURNS ALL STANDING AND FLAT CROP RESIDUE.
	 
	! local variables
      INTEGER:: J
      REAL:: X1, X2,X3, X4, X5, ADD, SUM, RTO
	  
      X2=1.- modelPara(49)
      ADD=0.
      SUM=0.
      DO J=1,LC
          RTO=MIN(.99,abvGroundBiom(J)/(totCropBio(J)+1.E-10))
          X1=modelPara(49)*abvGroundBiom(J)
          X4=420.*X1
          VAR(99)=VAR(99)-X4
          SMM(99,MO)=SMM(99,MO)-X4
          totCropBio(J)=totCropBio(J)-X1
          abvGroundBiom(J)=abvGroundBiom(J)-X1
          X3=modelPara(49)*actualCropN(J)*RTO
          actualCropN(J)=actualCropN(J)-X3
          X5=standCropResi(J)*X2
          standCropResi(J)=standCropResi(J)-X5
          ADD=ADD+420.*X5
          X1=modelPara(49)*standDeadResiN(J)
          standDeadResiN(J)=standDeadResiN(J)-X1
          SUM=SUM+X1+X3
      END DO
      structLitt(LD1)=structLitt(LD1)*X2
      metabLitt(LD1)=metabLitt(LD1)*X2
      X1=modelPara(49)*NStructLitt(LD1)
      NStructLitt(LD1)=NStructLitt(LD1)-X1
      X3=modelPara(49)*NMetabLitt(LD1)
      NMetabLitt(LD1)=NMetabLitt(LD1)-X3
      SUM=SUM+X1+X3
      lgStructLitt(LD1)=lgStructLitt(LD1)*X2
      X1=modelPara(49)*CStructLitt(LD1)
      CStructLitt(LD1)=CStructLitt(LD1)-X1
      X3=modelPara(49)*CMetabLitt(LD1)
      CMetabLitt(LD1)=CMetabLitt(LD1)-X3
      CLgStructLitt(LD1)=CLgStructLitt(LD1)*X2
      NLgStructLitt(LD1)=CStructLitt(LD1)-CLgStructLitt(LD1)
      ADD=ADD+X1+X3
      SMM(97,MO)=SMM(97,MO)+ADD
      VAR(97)=ADD
      SMM(98,MO)=SMM(98,MO)+SUM
      VAR(98)=SUM
      cropResidu(LD1)=.001*(structLitt(LD1)+metabLitt(LD1))
      standDeadResOrganism=0.
      standDeadResiOrgN=0.
      X1=standDeadResiP*modelPara(49)
      FreshOrgP_Residu(LD1)=FreshOrgP_Residu(LD1)+X1+StandDeadRes_OrgP
      standDeadResiP=standDeadResiP-X1
      StandDeadRes_OrgP=0.
 
END SUBROUTINE Burn_Residue
                         
! ------------------ I.2. CN from residue to soil -------------------------------------
SUBROUTINE CN_TopSoil(X1, X2, LCD)
      IMPLICIT NONE
      !     THIS SUBPROGRAM REMOVES C AND N FROM STANDING RESIDUE AND ADDS IT TO THE TOP SOIL 
      !      LAYER AS A RESULT OF A TILLAGE OPERATION.
      ! local variables
      INTEGER:: LCD, JSL
      REAL:: BX1, CNRTO, Y1, Y2, X1, X2, X3, X6, X7, X8, XX, XY, XZ, C7, RLN, RLR
	  
      IF(LCD>0)THEN
            JSL=ISL
            BX1=.1
      ELSE
            JSL=LD1
            BX1=.05
      END IF
   
      CNRTO=(CStructLitt(JSL)+CMetabLitt(JSL))/(NStructLitt(JSL)+NMetabLitt(JSL)+1.E-5)
      ! Organic residues are added to soil, a fraction of soil mineral N is sorbed on to the litter N component
      IF(CNRTO>=10.)THEN
            Y1=BX1*NO3_N_Soil(JSL)     ! APEX-doc Eq. 164
            Y2=BX1*NH3_Weight(JSL)
            X2=X2+Y1+Y2                ! X2: Litter N component 
            NO3_N_Soil(JSL)=NO3_N_Soil(JSL)-Y1
            NH3_Weight(JSL)=NH3_Weight(JSL)-Y2
      END IF    
      
      RLN=1000.*standDeadLignin/(standDeadResiN(Crop_Num)+1.E-5)
      RLR=standDeadLignin/(standCropResi(Crop_Num)+1.E-5)
      
      IF(RLR>.8)THEN
            RLR=.8
      ELSE
            IF(RLR<.1)RLR=.1
      END IF        
      
      X7=1000.*X1                 ! X1: organic materials in residue
      C7=.42*X7                   ! 0.42 is carbon fraction of organic materials
      SMS(7,JSL)=SMS(7,JSL)+C7   
      SMM(73,MO)=SMM(73,MO)+C7    ! SMM(73,MO), VAR(73): Carbon in organic materials of residues
      VAR(73)=VAR(73)+C7
      CN_Ratio=C7/(X2+1.E-5)           ! CN ratio 
      ! APEX-doc Eq. 165
      X8=.85-.018*RLN
      IF(X8<.01)THEN
          X8=.01
      ELSE
          IF(X8>.7)X8=.7
      END IF
      XX=X7*X8
      metabLitt(JSL)=metabLitt(JSL)+XX
      XZ=X7-XX
      structLitt(JSL)=structLitt(JSL)+XZ
      lgStructLitt(JSL)=lgStructLitt(JSL)+XZ*RLR
      X6=X2
      SMM(86,MO)=SMM(86,MO)+X6       ! SMM(86,MO):  N in litter component/ organic materials
      XY=.42*XZ                      ! Carbon in structural litter
      CStructLitt(JSL)=CStructLitt(JSL)+XY
      CLgStructLitt(JSL)=CLgStructLitt(JSL)+XY*RLR
      NLgStructLitt(JSL)=CStructLitt(JSL)-CLgStructLitt(JSL)
      X3=MIN(X6,XY/150.)
      NStructLitt(JSL)=NStructLitt(JSL)+X3
      CMetabLitt(JSL)=CMetabLitt(JSL)+.42*XX
      NMetabLitt(JSL)=NMetabLitt(JSL)+X6-X3
      cropResidu(JSL)=.001*(structLitt(JSL)+metabLitt(JSL))
 
END SUBROUTINE CN_TopSoil
 
! ******************************* II. Irrigation ***********************************
! --------------------------------II.1 irrgation ------------------------------

REAL FUNCTION Cal_WTN( )
      IMPLICIT NONE

      ! local variables
      INTEGER:: J
      REAL:: XX, X1
      DO J=1,Actual_SoilLayers
            ISL=Layer_ID(J)
            IF(Z(ISL)>=.15)GO TO 1
      END DO
      ISL=Layer_ID(Actual_SoilLayers)
    1 XX=LOG10(Wilt_Point(ISL))
      X1=3.1761-1.6576*(LOG10(soilWater(ISL))-XX)/(LOG10(fieldCapacity(ISL))-XX)
      IF(X1<4.)THEN
            Cal_WTN=MAX(5.,10.**X1)
      ELSE
            Cal_WTN=10000.
      END IF    
 
END FUNCTION Cal_WTN
! --------- -------- II.2 cal water needed for irrigation ----------------------------

SUBROUTINE Cal_IrrWater(IrrWater,EFD,ZX,JRT,IRX,IRY)                
      IMPLICIT NONE
      !     EPIC1102
      !     THIS SUBPROGRAM IS USED TO SIMULATE AUTOMATIC OR USER SPECIFIED IRRIGATION APPLICATIONS. 
      !     COMPUTES THE AMOUNT OF IRRIGATION WATER NEEDED TO BRING THE ROOT ZONE WATER CONTENT TO FIELD CAPACITY, FOR
      !     AUTOMATIC OPTION.  USER SPECIFIED AMOUNT IS APPLIED FOR MANUAL OPTION.  EROSION AND RUNOFF ARE ESTIMATED.

      ! local variables
      INTEGER:: N1, I, IRY, J, JRT, IRX 
      REAL:: IrrWater, X3, X4, EFD, YX,XX, X1, X2, QXM, ZX, QPX, DX, PX, AX, VX, CY, XEF
  
      IF(VIRT>=irrAnnualMaxVol.OR.NII<autoIrrInterval)THEN
            IrrWater=0.
            JRT=1                     ! no irrigation
            RETURN
      END IF
      X3=irrCost*IrrWater
      N1=MAX(1,NCP(Crop_Num))
      I=1
      X4=IrrWater      
      IrrWater=IrrWater*EFD
      YX=0.
      IF(IRY==1)THEN
            IF(irrChoice>0)GO TO 5
            X1=IrrWater
            XX=(potentialAvailWater-soilWaterRZ)/(1.-IrrRunoffRatio)
      ELSE
            X1=10000.
            IF(irrChoice>0)THEN
                  XX=10000.
            ELSE
                  XX=(potentialAvailWater-soilWaterRZ)/(1.-IrrRunoffRatio)
            END IF
      END IF
      X4=MIN(X1, irrAnnualMaxVol-VIRT,XX, irrSingleMaxVol)
      IrrWater=X4*EFD
      X3=irrCost*X4
      IF(IrrWater> irrSingleMinVol)GO TO 5
      IrrWater=0.
      GO TO 10
    5 NII=0
      QXM=IrrWater*IrrRunoffRatio  
      IF(irrCode==5)THEN
            DO J=1,Actual_SoilLayers
                  I=Layer_ID(J)
                  IF(ZX<Z(I))EXIT
            END DO
            soilWater(I)=soilWater(I)+IrrWater
      ELSE
            IF(QXM>0.)THEN    ! Ref: APEX-doc  Eq. 131 a -g 
                  QPX=QXM/24.   ! QXM: application rate in    mm/d    QPX: flow rate in  mm/h
                  IF(irrCode/=1.AND.RidgeHeight2>0..AND.Ridge_IntervalX>0.)THEN
                        X1=1000.*Ridge_IntervalX/RidgeHeight2                ! RidgeHeight2: mm   ; Ridge_IntervalX: m
                        QPX=2.778E-6*QPX*Ridge_IntervalX*areaWshed/fieldWidth! QPX: flow rate in m3/s
                        DX=(2.*QPX/(SX*X1*(1./(4.+16./(X1*X1)))**.3333))**.375 ! APEX-doc Eq. 131g
                        X2=DX*X1
                        PX=2.*SQRT(DX*DX+.25*X2*X2)
                        AX=.5*DX*X2
                        VX=(AX/PX)**.6667*SX
                        CY=modelPara(36)*VX**modelPara(31)
                        YX=10.*QXM*CY*EK         ! Furrow erosion  
                  ELSE
                        YX=2.5*CropManage_Factor*USL*SQRT(QPX*QXM)
                  END IF
            END IF
      END IF    
      Runoff=Runoff+QXM  
      totSoilErosion=totSoilErosion+YX                          ! Accumulate soil erosion
      Sediment(waterEroModel)=Sediment(waterEroModel)+YX
      X1=soilWaterRZ-potentialAvailWater
      verticleFlowSalt=.01*IrrWater*saltConIrr
      SMM(69,MO)=SMM(69,MO)+verticleFlowSalt              ! SMM(69,MO), VAR(69) : Salt in irrigation water   
      VAR(69)=verticleFlowSalt
      IF(KG(Crop_Num)>0.OR.JPL(Crop_Num)>0)XHSM=HU(Crop_Num)/potentialHeatUnit(Crop_Num,IHU(Crop_Num))
      IF(NOP>0)WRITE(KW(1),9)IYR,MO,DayOfMon,IrrWater,waterStress,WTN,X1,irrTrigger,QXM,YX,XHSM,X3
      VIRT=VIRT+X4   
      COST=COST+X3
      IF(VIRT>0.)THEN
            IF(IRY==0)THEN
                  X1=COTL(IRX)
                  X2=X1-COOP(IRX)
                  COST=COST+X1
                  CSFX=CSFX+X2
            END IF
            IF(KFL(20)>0)THEN
                  WRITE(KW(20),14)IYR,MO,DayOfMon,equipmentName(IRX),cropID(Crop_Num),currentOps(IRX),NBE(IRX),X3,&
                  X3,IrrWater
                  IF(IRY==0)WRITE(KW(20),50)IYR,MO,DayOfMon,equipmentName(IRX),cropID(Crop_Num),currentOps(IRX),&
                  NBE(IRX),NBT(IRX),X1,X2,Fuel_Use(IRX)
            END IF
      END IF
      SMM(19,MO)=SMM(19,MO)+X4                        ! SMM(19, MO), VAR(19) : irrigation water applied
      VAR(19)=X4
      SMM(92,MO)=SMM(92,MO)+Fuel_Use(IRX)             ! SMM(19, MO) : fuel use
      XEF=X4-IrrWater
      SMM(84,MO)=SMM(84,MO)+XEF                       ! SMM(84, MO), VAR(84) : Irrigation distribution loss
      VAR(84)=XEF
      IrrWater=IrrWater-QXM
      VIR(N1,Crop_Num)=VIR(N1,Crop_Num)+X4
      VIL(N1,Crop_Num)=VIL(N1,Crop_Num)+XEF
      NO3_IrrWater=IrrWater*CNO3I
      ANA(Crop_Num)=ANA(Crop_Num)+NO3_IrrWater        ! ANA : NO3 in water
      NO3_N_Soil(I)=NO3_N_Soil(I)+NO3_IrrWater
      SMM(60,MO)=SMM(60,MO)+NO3_IrrWater               ! SMM(60, MO), VAR(60) : NO3 
      VAR(60)=VAR(60)+NO3_IrrWater
      IF(NOP>0.AND.NO3_IrrWater>0.)WRITE(KW(1),133)IYR,MO,DayOfMon,NO3_IrrWater,IrrWater,XHSM
   10 JRT=0
      IF(irrCode==5)THEN
            AWC(Crop_Num)=AWC(Crop_Num)+IrrWater
            IrrWater=0.
      END IF
 
    9 FORMAT(1X,3I4,2X,'IRRIGATE',10X,'VOL=',F5.0,'MM',2X,'waterStress=',F6.2,&
      2X,'WTN=',E12.4,'kpa',2X,'PWDF=',F6.0,'MM',2X,'TRGR=',F7.2,&
      2X,'Q=',F5.0,'MM',2X,'Y=',F5.1,'t/ha',2X,'fracHeatUnit=',F6.2/2X,'COST=',&
      F7.2,'$/ha')
   14 FORMAT(1X,3I4,2X,A8,8X,I6,6X,2I4,4X,F10.2,10X,2F10.2)
   50 FORMAT(1X,3I4,2X,A8,8X,I6,6X,3I4,2F10.2,20X,F10.2)
  133 FORMAT(1X,3I4,2X,'NO3 FERT = ',F5.0,'kg/ha',2X,'IRR VOL = ',F5.0,'MM',2X,'fracHeatUnit = ',F6.2)
  
END SUBROUTINE Cal_IrrWater

! ******************************* III. Fertilier ***********************************
! ---------------------------------III.1 Apply Fertilizer ----------------------------
SUBROUTINE Apply_Fert(IRC,JFT)
      !     EPIC1102
      !     THIS SUBPROGRAM APPLIES N AND P FERTILIZER AT SPECIFIED DATES, RATES, AND DEPTH.
      IMPLICIT NONE 
      !   local variables
      INTEGER:: I, J, JFT, IRC, N1
      REAL:: ADD, X, X1,X2, X3, X4, X5,X6, X7, X8,X9, X10, X11, ZZ,XN, W1, RLN,XX, YY, XZ, YZ, &
            Y1, Y2, MinThick_profT, SOM, DZB, BC, BCF, CECM, OM, CECA ! XHSM is not used, which may cause problem
                
      IF(ANA(Crop_Num)>=FNMX(Crop_Num))RETURN
      X1=0.
      I=LD1
      MinThick_profT=tillDepth(JFT)
      SELECT CASE(IRC)
            CASE(1)
                  numOfFert=LFT(IRO,KT)
                  X1=WFA(IRO,KT)
            CASE(2)
                  numOfFert=IDFT(1)
                  X1=fertRate
            CASE(3)
                  IF(irrCode/=4)THEN
                        numOfFert=IDFT(1)
                        X1=fertRate
                  ELSE
                        numOfFert=IDFT(2)
                  END IF
            CASE(4)
                  numOfFert=IDFT(2)
                  X3=NStressTrigger-TNOR
                  GO TO 1
            CASE(5)
                  numOfFert=numAutoP
                  X1=APMU
                  X3=APMU*FN(numOfFert)
                  GO TO 13              
      END SELECT
      DO J=1,Actual_SoilLayers
            I=Layer_ID(J)
            IF(MinThick_profT<Z(I))EXIT
      END DO
      IF(X1>0.)THEN
            IF(FN(numOfFert)>0.)THEN
                  X3=X1*FN(numOfFert)
            ELSE
                  X3=0.
                  GO TO 13
            END IF
      ELSE
            X1=MAX(modelPara(28)*UNA(Crop_Num)-TNOR,0.)
            IF(IAP(Crop_Num)==0)THEN
                  IAP(Crop_Num)=1
            ELSE
                  X1=X1*(1.-HUI(Crop_Num))
            END IF    
            X3=X1
            IF(FN(numOfFert)<1.E-10)RETURN
      END IF
    1 X3=MIN(X3,FNMX(Crop_Num)-ANA(Crop_Num))
      X1=X3/FN(numOfFert)
   13 X2=X1*FP(numOfFert)        ! fration of N, P, K
      X4=X3*FNH3(numOfFert)
      X5=X1*FNO(numOfFert)
      X6=X1*FPO(numOfFert)
      X7=X3-X4
      X8=X1*FOC(numOfFert)
      SMM(65,MO)=SMM(65,MO)+X8        ! SMM(65,MO), VAR(65): organic carbon fraction in fertilizer    
      VAR(65)=X8
      
      IF(FOC(numOfFert)>.5)THEN
            ZZ=0.
            DO I=1,Actual_SoilLayers
                  ISL=Layer_ID(I)
                  IF(Z(ISL)<=.15)THEN
                        SOM=.172*SOC(ISL)/WT(ISL)
                        DZB=(Z(ISL)-ZZ)/.15
                        BC=.1*X1*DZB/WT(ISL) ! % BIOCHAR IN LAYER ISL
                        bulkDensity(ISL)=100./(SOM/.244+(100.-SOM-BC)/mineralBulkDensity(ISL)+BC/.55)
                        BCF=.01*BC           ! FRACTION BIOCHAR IN LAYER ISL
                        CECM=(CEC(ISL)+BCF*187.)/(1.+BCF)
                        ADD=CECM-CEC(ISL)
                        OM=(SOM+BC)*WT(ISL)*1.E4
                        CECA=10.*ADD/OM
                        CEC(ISL)=CECM 
                        X=6.6-(LOG(3.8056/(PH(ISL)-3.495))-1.)/1.08
                        XN=X+CECA
                        PH(ISL)=3.8056/(1.+EXP(1.08*(6.6-XN)))+3.495
                        W1=X8*DZB
                        metabLitt(ISL)=metabLitt(ISL)+.02*X1*DZB
                        CMetabLitt(ISL)=CMetabLitt(ISL)+.02*W1
                        CSlowHumus(ISL)=CSlowHumus(ISL)+.6*W1
                        CPassiveHumus(ISL)=CPassiveHumus(ISL)+.38*W1                                      
                  ELSE
                        EXIT          
                  END IF
                  ZZ=Z(ISL)  
            END DO
          
            IF(I<Actual_SoilLayers)THEN
                  SOM=.172*SOC(ISL)/WT(ISL)
                  DZB=(.15-ZZ)/.15
                  BC=.1*X1*DZB/WT(ISL) ! % BIOCHAR IN LAYER ISL
                  bulkDensity(ISL)=100./(SOM/.244+(100.-SOM-BC)/mineralBulkDensity(ISL)+BC/.55)
                  BCF=.01*BC           ! FRACTION BIOCHAR IN LAYER ISL
                  CECM=(CEC(ISL)+BCF*187.)/(1.+BCF)
                  ADD=CECM-CEC(ISL)
                  OM=(SOM+BC)*WT(ISL)*1.E4
                  CECA=10.*ADD/OM
                  CEC(ISL)=CECM
                  X=6.6-(LOG(3.8056/(PH(ISL)-3.495))-1.)/1.08
                  XN=X+CECA
                  PH(ISL)=3.8056/(1.+EXP(1.08*(6.6-XN)))+3.495
                  W1=X8*DZB
                  metabLitt(ISL)=metabLitt(ISL)+.02*X1*DZB
                  CMetabLitt(ISL)=CMetabLitt(ISL)+.02*W1
                  CSlowHumus(ISL)=CSlowHumus(ISL)+.6*W1
                  CPassiveHumus(ISL)=CPassiveHumus(ISL)+.38*W1                                                     
            END IF      
      ELSE
            IF(X8>.1)THEN
                  RLN=.175*X8/(X3+X5)
                  X10=.85-.018*RLN
                  IF(X10<.01)THEN
                        X10=.01
                  ELSE
                        IF(X10>.7)X10=.7
                  END IF
                  XX=X8*X10
                  CMetabLitt(I)=CMetabLitt(I)+XX
                  YY=X1*X10
                  metabLitt(I)=metabLitt(I)+YY
                  ZZ=X5*X10
                  NMetabLitt(I)=NMetabLitt(I)+ZZ
                  NStructLitt(I)=NStructLitt(I)+X5-ZZ
                  XZ=X8-XX
                  CStructLitt(I)=CStructLitt(I)+XZ
                  CLgStructLitt(I)=CLgStructLitt(I)+XZ*.175
                  NLgStructLitt(I)=CStructLitt(I)-CLgStructLitt(I)
                  YZ=X1-YY
                  structLitt(I)=structLitt(I)+YZ
                  lgStructLitt(I)=lgStructLitt(I)+YZ*.175
            END IF
      END IF
      X9=X1*FK(numOfFert)
      X11=X1*FSLT(numOfFert)
      NH3_Weight(I)=NH3_Weight(I)+X4
      labileP(I)=labileP(I)+X2
      SOP(I)=SOP(I)+X6
      SMM(59,MO)=SMM(59,MO)+X5         ! SMM(59,MO), VAR(59) : Organic Nitrogen Fertilizer (manure) kg/ha
      VAR(59)=X5
      SMM(62,MO)=SMM(62,MO)+X6         ! SMM(62,MO), VAR(62) : Organic P  Fertilizer                kg/ha
      VAR(62)=X6
      SMM(61,MO)=SMM(61,MO)+X4         ! SMM(61,MO), VAR(61) : NH3   Fertilizer                     kg/ha
      VAR(61)=X4
      SMM(63,MO)=SMM(63,MO)+X2         ! SMM(63,MO), VAR(63) : Minearl P Fertilizer                 kg/ha   
      VAR(63)=X2
      NO3_N_Soil(I)=NO3_N_Soil(I)+X7
      ANA(Crop_Num)=ANA(Crop_Num)+X3
      SMM(60,MO)=SMM(60,MO)+X7          ! SMM(60,MO), VAR(60) : NO3   Fertilizer                     kg/ha
      VAR(60)=VAR(60)+X7
      solubleK(I)=solubleK(I)+X9
      SMM(64,MO)=SMM(64,MO)+X9          ! SMM(64,MO), VAR(64) : K   Fertilizer                     kg/ha
      VAR(64)=X9
      saltWeight(I)=saltWeight(I)+X11
      SMM(72,MO)=SMM(72,MO)+X11         ! SMM(72,MO), VAR(72) : Salt   Fertilizer                     kg/ha
      VAR(72)=X11
      N1=MAX(1,NCP(Crop_Num))
      TotN_Frac(N1,Crop_Num)=TotN_Frac(N1,Crop_Num)+X3+X5
      TotP_Frac(N1,Crop_Num)=TotP_Frac(N1,Crop_Num)+X2+X6
      TotK_Frac(N1,Crop_Num)=TotK_Frac(N1,Crop_Num)+X9
      XX=X1*FCST(numOfFert)
      COST=COST+XX
      SMM(96,MO)=SMM(96,MO)+FCEM(numOfFert)*X1   ! SMM(96,MO), VAR(96) : carbon emission                    kg/ha
      IF(IRC==3)THEN
            Y1=COTL(JFT)
            Y2=Y1-COOP(JFT)
            COST=COST+Y1
            CSFX=CSFX+Y2
      END IF
      IF(KFL(20)>0)THEN
            WRITE(KW(20),18)IYR,MO,DayOfMon,fertName(numOfFert),cropID(Crop_Num),KDF(numOfFert),&
                          currentOps(JFT),NBE(JFT),NBT(JFT),XX,XX,X1
            IF(IRC==3)THEN
                  WRITE(KW(20),50)IYR,MO,DayOfMon,equipmentName(JFT),cropID(Crop_Num),currentOps(JFT),&
                              NBE(JFT),NBT(JFT),Y1,Y2,Fuel_Use(JFT)
            END IF
      END IF
      IF(NOP>0)WRITE(KW(1),90)IYR,MO,DayOfMon,fertName(numOfFert),X1,MinThick_profT,X3,X4,X5,X2,&
      X6,X9,XHSM,XX
      NFA=0
      RETURN
   18 FORMAT(1X,3I4,2X,A8,8X,I6,2X,4I4,F10.2,10X,2F10.2,10X,F10.2)
   50 FORMAT(1X,3I4,2X,A8,8X,I6,6X,3I4,2F10.2,20X,F10.2)
   90 FORMAT(1X,3I4,2X,A8,2X,'RATE=',F6.0,'kg/ha',1X,'DPTH=',F5.2,'m', &
             1X,'MN=',F5.0,1X,'NH3=',F5.0,1X,'ON=',F5.0,1X,'MP=',F5.0,1X,'stableMineralP=', &
             F5.0,1X,'MK=',F5.0,1X,'fracHeatUnit=',F5.2,2X,'COST=',F8.2,'$/ha')
END SUBROUTINE Apply_Fert
	
! --------------------- III.2 Estimate aluminum saturation ---------------------------
	
SUBROUTINE Est_Aluminum(SB,DSB,C1,soilPH, soilALS, OC, BSA)
      IMPLICIT NONE
      !     EPIC1102
      !     THIS SUBPROGRAM ESTIMATES ALUMINUM SATURATION USING BASE
      !     SATURATION, ORGANIC C, AND PH.
      ! Ref: APEX-doc Eq 309 
      ! global variables
      REAL:: soilPH, soilALS, BSA
	! local variables
      REAL:: SB, DSB, C1, OC

      SB=SB-DSB
      SB=MAX(.02,SB)
      BSA=C1*SB
      IF(soilPH<=5.6)THEN
            soilALS=154.2-1.017*BSA-3.173*OC-14.23*soilPH
            IF(soilALS<.01)THEN
                  soilALS=0.
                  RETURN
            ELSE
                  IF(soilALS>95.)soilALS=95.
            END IF
      ELSE
            soilALS=0.
      END IF        
END SUBROUTINE Est_Aluminum
 
!------------------------------- III.3 Apply Lime ------------------------------
	
SUBROUTINE Apply_Lime(TLX)
      IMPLICIT NONE
      !     EPIC1102
      !     THIS SUBPROGRAM APPLIES LIME WHEN THE SUM OF THE SOIL LIME REQUIRE
      !     MENT AND ACCUMULATED LIME REQUIREMENT CAUSED BY N FERTILIZER EXCEED
      !     4 t/ha.
 
      ! local variables
      INTEGER:: J, L1,K
      REAL:: SMFN = 0.0, OC, TOT,XY, XZ, ZZ, XX, W3, W2, RTO, X1, BS, DSB, TLX, EAL, TLZ, &
            DBS, ALN, PHN, W1, ALSX, BSA

      OC=0.
      TOT=0.
      XZ=0.
      ZZ=0.
      XX=0.
      DO J=1,Actual_SoilLayers
            ISL=Layer_ID(J)
            IF(Z(ISL)>BIG)EXIT
            XY=WT(ISL)
            OC=OC+.1*SOC(ISL)   
            XZ=XZ+CEC(ISL)*XY
            TOT=TOT+PH(ISL)*XY
            ZZ=ZZ+totBase(ISL)*XY
            XX=XX+XY
      END DO
      IF(J>Actual_SoilLayers)THEN
            J=Actual_SoilLayers
            ISL=Layer_ID(Actual_SoilLayers)
      ELSE
            L1=Layer_ID(J-1)
            W3=Z(ISL)-Z(L1)
            W2=BIG-Z(L1)
            RTO=W2*WT(ISL)/W3
            X1=.1*SOC(ISL)/WT(ISL)
            OC=OC+RTO*X1
            TOT=TOT+RTO*PH(ISL)
            ZZ=ZZ+RTO*totBase(ISL)
            XZ=XZ+RTO*CEC(ISL)
            XX=XX+RTO
      END IF
      XZ=XZ/XX
      OC=OC/XX
      TOT=TOT/XX
      ZZ=ZZ/XX
      XY=.001*XX
      X1=SMY(60)+SMY(61)
      DSB=.036*(SMY(50)+X1-SMFN)/XX
      SMFN=X1
      BS=100./XZ
      TOT=TOT-.05*DSB*BS
      CALL Est_Aluminum(ZZ,DSB,BS,TOT,ALSX, OC, BSA)
      IF(TOT>6.5.AND.TLX<1.E-5)GO TO 6
      IF(IDSP==4)THEN
            EAL=.01*ALSX*XZ
            IF(TLX>0.)THEN
                  TLZ=TLX
            ELSE
                  TLZ=EAL*XY
                  IF(TLZ<1.)GO TO 6
            END IF
            TOT=5.4
            CALL Est_Aluminum(ZZ,-EAL,BS,TOT,ALSX, OC, BSA)
            GO TO 7
      END IF
      DBS=MIN((6.5-TOT)/.023,90.-BSA)
      ALN=0.
      RTO=1.
      IF(TLX>0.)THEN
            TLZ=TLX
      ELSE          
            TLZ=DBS*XY/BS
            IF(TLZ>2.)THEN
                  RTO=2./TLZ
                  TLZ=2.
            ELSE
                  IF(TOT>5.)GO TO 6
            END IF
      END IF
      PHN=(6.5-TOT)*RTO+TOT
      DBS=MIN((PHN-TOT)/.023,90.-BSA)
      BSA=(BSA+DBS)/BS
      GO TO 8
    6 TLZ=0.
    7 ALN=ALSX
      PHN=TOT
      BSA=ZZ
    8 limeRate=TLZ  
      DO K=1,J
            ISL=Layer_ID(K)
            TOT=totBase(ISL)
            XZ=PH(ISL)
            ALSX=ALS(ISL)
            totBase(ISL)=BSA
            PH(ISL)=PHN
            ALS(ISL)=ALN
      END DO
      IF(J==Actual_SoilLayers)RETURN
      ISL=Layer_ID(J)
      W1=Z(ISL)-BIG
      totBase(ISL)=(W1*TOT+W2*totBase(ISL))/W3
      PH(ISL)=(W1*XZ+W2*PH(ISL))/W3
      ALS(ISL)=MAX(.001,(W1*ALSX+W2*ALS(ISL))/W3)
 
END SUBROUTINE Apply_Lime  
      
! ******************************* IV.  Pesticide **********************************  
! -------------------------------IV. 1 apply pesticide ----------------------------
SUBROUTINE Apply_Pesticide 
      !     EPIC1102
      !     THIS SUBPROGRAM APPLIES PESTICIDES TO CROP CANOPY & SOIL
      IMPLICIT NONE
 
      ! local variables
      INTEGER:: K
      REAL:: XX, X1
 
      numOfPest=LPC(IRO,KT)
      XX=pestRate(IRO,KT)*harvestEffi(JT1)
      SMMP(1,numOfPest,MO)=SMMP(1,numOfPest,MO)+XX       ! Pest controlled
      totPestControl(Crop_Num,numOfPest)=totPestControl(Crop_Num,numOfPest)+XX
      pestVar(1,numOfPest)=XX
      X1=Pest_Cost(numOfPest)*pestRate(IRO,KT)
      COST=COST+X1
      SMM(96,MO)=SMM(96,MO)+Pest_C_emission(numOfPest)       ! SMM(96,MO): C emission from pestiside
      
      IF(KFL(20)>0)WRITE(KW(20),1)IYR,MO,DayOfMon,pestName(numOfPest),cropID(Crop_Num),&
                      KDP(numOfPest),currentOps(JT1),NBE(JT1),NBT(JT1),X1,X1,pestRate(IRO,KT)
	  
      pestGrowIndx=pestGrowIndx-pestKillEffi(IRO,KT)*modelPara(37)
      
      IF(NOP>0)WRITE(KW(1),6)IYR,MO,DayOfMon,pestName(numOfPest),pestRate(IRO,KT),harvestEffi(JT1),&
                             pestKillEffi(IRO,KT),pestGrowIndx,X1
      IF(tillDepth(JT1)<1.E-10)THEN
            X1=XX*GroundCover_Frac
            pestPlantIntercep(numOfPest)=pestPlantIntercep(numOfPest)+X1
            pestInSoil(numOfPest,LD1)=pestInSoil(numOfPest,LD1)+XX-X1
            RETURN
      ELSE
            DO K=1,Actual_SoilLayers
                  ISL=Layer_ID(K)
                  IF(tillDepth(JT1)<=Z(ISL))EXIT
            END DO
            pestInSoil(numOfPest,ISL)=pestInSoil(numOfPest,ISL)+XX
      END IF
    1 FORMAT(1X,3I4,2X,A16,I6,2X,4I4,F10.2,10X,3F10.2)
    6 FORMAT(1X,3I4,2X,A8,2X,'APPL RATE = ',F5.1,'kg/ha',2X,'APPL EFF = ',F6.2,2X,'KILL EFF = ',F6.2,2X,&
            'PST IDX = ',E12.4,2X,'COST=',F7.0,'$/ha')
END SUBROUTINE Apply_Pesticide    

!  ------------------------- IV.2 PESTICIDES transport and degradation -------------------------
SUBROUTINE Pest_Trans 
      !     EPIC1102
      !     THIS SUBPROGRAM SIMULATES PESTICIDE TRANSPORT & DEGRADATION
      IMPLICIT NONE
 
      ! local variables
      INTEGER:: II, I, I1, K, L1
      INTEGER, DIMENSION(5):: NY= [1,4,21,60,90]
      INTEGER, DIMENSION(90):: NXP = (/(I, I = 1, 90)/)
      REAL:: SQB(5), SYB(5), QQ, Y1, Y2, ADD, SUM, TOT, X1, X2, X3, X4, PQ, PY, WO, DGF, DGS, DK, XX, &
            V, VPST, CO, Pest_WaterFlow, Bio_contrl_fact       
 
      II=NXP(90)
      QQ=Runoff+lateralFlow
      Y2=Sediment(waterEroModel) ! selected a way of soil loss by water erosion
      SQB(1)=QQ                          ! totoal water flow
      SYB(1)=Y2                          ! soil loss by water erosion (sediment)
      DO I=2,5
            SQB(I)=SQB(I)+QQ-VQ(NXP(NY(I)))
            SYB(I)=SYB(I)+Y2-VY(NXP(NY(I)))
      END DO
      DO 9 K=1,NDP
      ADD=0.
      SUM=0.
      TOT=0.
      X3=0.
      PQ=0.
      PY=0.
      Y1=pestInSoil(K,LD1)         ! Pesticide in soil layers
      IF(IGO>0)THEN
            IF(pestPlantIntercep(K)>.01)THEN
                  IF(totRunoff>2.54)THEN
                  ! COMPUTE PESTICIDE WASH OFF FROM FOLIAGE
                  ! WO: the amount of pesticide washed off the plants by a rainstorm      
                  WO=WashOff_Frac(K)*pestPlantIntercep(K)          
                  pestPlantIntercep(K)=pestPlantIntercep(K)-WO
                  Y1=Y1+WO
                  END IF
                  !  COMPUTE PESTICIDE DEGRADATION FROM FOLIAGE
                  DGF=pestPlantIntercep(K)*pestHfLifePlant(K)
                  pestPlantIntercep(K)=pestPlantIntercep(K)-DGF
                  SMMP(6,K,MO)=SMMP(6,K,MO)+DGF              ! SMMP(6,K,MO), pestVar(6): pestcide degraded on plant
                  pestVar(6,K)=DGF  
            ELSE
                  pestPlantIntercep(K)=0.
            END IF
      END IF
      !     COMPUTE PESTICIDE LOSS FROM TOP SOIL LAYER IN RUNOFF,
      !     LATERAL SUBSURFACE FLOW, & PERCOLATION
      IF(Y1>.01)THEN                  
            DK=.0001*Coeff_OrgCarbon(K)*SOC(LD1)       ! APEX-doc  Eq. 257
            X1=Porosity(LD1)-Wilt_Point(LD1)
            XX=X1+DK
            V=Runoff+subsurfLaterFlow(LD1)+percolateFlow(LD1)  ! total runoff
            IF(V>0.)THEN
                  VPST=Y1*(1.-EXP(-V/XX))
                  CO=MIN(pestSolubity(K),VPST/(percolateFlow(LD1)+modelPara(18)*(Runoff+subsurfLaterFlow(LD1))))
                  Bio_contrl_fact=modelPara(18)*CO
                  X3=CO*percolateFlow(LD1)        ! leach from percolation
                  PQ=Bio_contrl_fact*Runoff          ! leach from runoff
                  SMMP(2,K,MO)=SMMP(2,K,MO)+PQ           ! SMMP(2,K,MO), pestVar(2): pesticide in surface runoff
                  pestVar(2,K)=PQ
                  SUM=Bio_contrl_fact*subsurfLaterFlow(LD1) ! leach from subsurface lateral flow
                  Y1=Y1-X3-PQ-SUM
                  ! WRITE(KW(1),3)IYR,MO,DayOfMon,Runoff,subsurfLaterFlow(LD1),Y1,VPST,CO,SUM,PQ
                  ! COMPUTE PESTICIDE LOSS WITH SEDIMENT
                  IF(enrichRatio>0.)THEN                
                        Bio_contrl_fact=DK*Y1/XX
                        PY=enrichRatio*Bio_contrl_fact
                        SMMP(5,K,MO)=SMMP(5,K,MO)+PY        ! SMMP(5,K,MO), pestVar(5): pesticide in sediment
                        pestVar(5,K)=PY  
                        Y1=Y1-PY
                  END IF
            END IF
            ! COMPUTE PESTICIDE DEGRADATION IN TOP SOIL LAYER
            DGS=Y1*pestHfLifeSoil(K)
            Y1=Y1-DGS
            TOT=DGS
            ADD=Y1
      ELSE
            Y1=0.
      END IF
      pestInSoil(K,LD1)=Y1       
      !     COMPUTE PESTICIDE MOVEMENT THRU SOIL LAYERS BY LATERAL
      !     SUBSURFACE FLOW & PERCOLATION
      !     COMPUTE PESTICIDE MOVEMENT THRU SOIL LAYERS BY LATERAL
      X2=0.
      DO L1=2,Actual_SoilLayers
            ISL=Layer_ID(L1)
            Y1=pestInSoil(K,ISL)
            Y1=Y1+X3
            X3=0.
            IF(Y1>.01)THEN
                  V=percolateFlow(ISL)+subsurfLaterFlow(ISL)
                  IF(V>0.)THEN
                        VPST=Y1*(1.-EXP(-V/(Porosity(ISL)-Wilt_Point(ISL)+.0001*Coeff_OrgCarbon(K)*SOC(ISL))))
                        CO=MIN(pestSolubity(K),VPST/(percolateFlow(ISL)+modelPara(18)*subsurfLaterFlow(ISL)))
                        Bio_contrl_fact=modelPara(18)*CO
                        X4=Bio_contrl_fact*subsurfLaterFlow(ISL)
                        IF(ISL==IDR)THEN                        ! if in drainage system
                              SMMP(10,K,MO)=SMMP(10,K,MO)+X4      ! SMMP(10,K,MO) : pesticide drainaged
                              pestVar(10,K)=X4
                        END IF
                        SUM=SUM+X4
                        X3=CO*percolateFlow(ISL)
                        IF(L1==Actual_SoilLayers)X2=X3
                        Y1=Y1-X4-X3
                  ELSE
                        ! COMPUTE PESTICIDE DEGRADATION IN SOIL LAYERS
                        DGS=Y1*pestHfLifeSoil(K)
                        Y1=Y1-DGS
                        TOT=TOT+DGS
                        ADD=ADD+Y1
                  END IF
            ELSE
                  Y1=0.
            END IF
            pestInSoil(K,ISL)=Y1
      END DO
      SMMP(3,K,MO)=SMMP(3,K,MO)+X2               ! SMMP(3,K,MO), pestVar(3) : pesticide  percolated
      pestVar(3,K)=X2
      SMMP(4,K,MO)=SMMP(4,K,MO)+SUM              ! SMMP(4,K,MO), pestVar(4) : pesticide  drainged
      pestVar(4,K)=SUM
      SMMP(7,K,MO)=SMMP(7,K,MO)+TOT              ! SMMP(7,K,MO), pestVar(7) : pesticide degraded in soil
      pestVar(7,K)=TOT
      SMMP(9,K,MO)=ADD                           ! SMMP(9,K,MO), pestVar(9) : total pesticide in soil layers
      pestVar(9,K)=ADD
      
      Pest_Leach(K)=X2
      Pest_Drainged(K)=SUM
      Pest_WaterFlow=PQ+SUM
      !  ------???  below do not understand --------
      SPQ(1,K)=Pest_WaterFlow
      Pest_Sediment(1,K)=PY
      DO I=2,5                      ! Why is it 5 ?
            SPQ(I,K)=SPQ(I,K)+Pest_WaterFlow-PVQ(K,NXP(NY(I)))   ! Why add Pest_WaterFlow ?
            IF(SPQ(I,K)<1.E-3.OR.SQB(I)<1.E-3)THEN
                  SPQC(I,K)=0.
            ELSE
                SPQC(I,K)=100.*SPQ(I,K)/SQB(I)
            END IF
            Pest_Sediment(I,K)=Pest_Sediment(I,K)+PY-PVY(K,NXP(NY(I)))
      END DO
      DO I=1,5
            IF(APQ(I,K,IY)<=SPQ(I,K))THEN
                  APQ(I,K,IY)=SPQ(I,K)
                  AQB(I,K,IY)=SQB(I)
            END IF
            IF(APQC(I,K,IY)<SPQC(I,K))APQC(I,K,IY)=SPQC(I,K)
            IF(APY(I,K,IY)>Pest_Sediment(I,K))CYCLE
            APY(I,K,IY)=Pest_Sediment(I,K)
            AYB(I,K,IY)=SYB(I)
      END DO
      PVQ(K,II)=Pest_WaterFlow
      PVY(K,II)=PY
      VQ(II)=QQ
      VY(II)=Y2
    9 CONTINUE
      DO I=90,2,-1
            I1=I-1
            NXP(I)=NXP(I1)
      END DO
      NXP(1)=II
END SUBROUTINE Pest_Trans
!  ---------------------IV.3 Pest factor and damage --------------------------------
SUBROUTINE Pest_Damage 
      IMPLICIT NONE
      !     EPIC1102
      !     THIS SUBPROGRAM CALCULATES THE PEST FACTOR BASED ON THE MINIMUN
      !     PEST FACTOR (PST(Crop_Num)) FOR CROP Crop_Num AND THE SUM OF DAILY PEST
      !     DAMAGE(pestGrowIndx)

	! local variables
      REAL:: X1
	  
      IF(pestGrowIndx>0.)THEN
            X1=pestDamgFactor*pestGrowIndx/GrowingSeason_Days
            PSTF(Crop_Num)=1.-(1.-pestDmagFactor(Crop_Num))*X1/(X1+EXP(S_Curve(9,1)-S_Curve(9,2)*X1))
      ELSE
            PSTF(Crop_Num)=1.
      END IF
END SUBROUTINE Pest_Damage
! ******************************* V. Harvest ****************************************** 
! --------------------- Harvest Root Crops --------------------------------------
SUBROUTINE Har_RootCrop(YY,X3,X1,X6,X7,N1)
      IMPLICIT NONE
 
	! local variables
      INTEGER:: N1, J
      REAL:: AD1, AD2, X1, X3, X6, X7, X8, X9, X10, X11,  X12, X13, XX, XZ, YY, U2
	  
      JD=Crop_Num
      X1=overRideHI(JT1)
      X3=totRootWeight(Crop_Num)
      XX=totCropBio(Crop_Num)
      XZ=X1*X3
      CALL Pest_Damage 
      X6=PSTF(Crop_Num)
      TPSF(N1,Crop_Num)=TPSF(N1,Crop_Num)+X6
      NPSF(N1,Crop_Num)=NPSF(N1,Crop_Num)+1
      YY=XZ*harvestEffi(JT1)*X6
      X12=actualCropN(Crop_Num)/XX
      X13=actualCropP(Crop_Num)/XX
      X8=actualCropK(Crop_Num)/XX
      X9=.9*abvGroundBiom(Crop_Num)
      X10=X12*X9
      AD1=X10
      CALL CN_TopSoil(X9,X10,0)
      U2=X13*X9 
      AD2=U2
      FreshOrgP_Residu(LD1)=FreshOrgP_Residu(LD1)+U2
      solubleK(LD1)=solubleK(LD1)+actualCropK(Crop_Num)
      Yield_N=YY*X12
      YLP=YY*X13
      YLC=.42*YY
      YLK=YY*X8
      XX=X1*(1.-harvestEffi(JT1))
      totRootWeight(Crop_Num)=0.
      DO J=1,Layer_RD_Reach
            ISL=Layer_ID(J)
            X11=RWT(ISL,Crop_Num)*XX
            RWT(ISL,Crop_Num)=RWT(ISL,Crop_Num)*(1.-X1)
            X10=X11*X12
            AD1=AD1+X10
            CALL CN_TopSoil(X11,X10,1)
            U2=X11*X13
            AD2=AD2+U2 
            FreshOrgP_Residu(ISL)=FreshOrgP_Residu(ISL)+U2
            totRootWeight(Crop_Num)=totRootWeight(Crop_Num)+RWT(ISL,Crop_Num)
      END DO
      YLD(Crop_Num)=YY
      YLD1(N1,Crop_Num)=YLD1(N1,Crop_Num)+YY
      YLNF(N1,Crop_Num)=YLNF(N1,Crop_Num)+Yield_N
      YLPF(N1,Crop_Num)=YLPF(N1,Crop_Num)+YLP
      YLCF(N1,Crop_Num)=YLCF(N1,Crop_Num)+YLC
      YLKF(N1,Crop_Num)=YLKF(N1,Crop_Num)+YLK
      X7=.001*X12
      abvGroundBiom(Crop_Num)=abvGroundBiom(Crop_Num)-X9
      totCropBio(Crop_Num)=totRootWeight(Crop_Num)+abvGroundBiom(Crop_Num)
      actualCropN(Crop_Num)=actualCropN(Crop_Num)-Yield_N-AD1
      actualCropP(Crop_Num)=actualCropP(Crop_Num)-YLP-AD2
      !actualCropK(Crop_Num)=actualCropK(Crop_Num)-YLK
      HU(Crop_Num)=HU(Crop_Num)*.1
      Current_LAI(Crop_Num)=.001
  
END SUBROUTINE Har_RootCrop
	  
! -------------------------------------- 2. Tillage ---------------------------------
SUBROUTINE Tillage_Mix(Mix_Efficiency, Tillage_Dep, NMIX, IBMX)
      IMPLICIT NONE
      !     EPIC1102
      !     THIS SUBPROGRAM MIXES N,P, AND CROP RESIDUE WITHIN THE PLOW DEPTH ACCORDING TO THE MIXING harvestEffi OF THE IMPLEMENT, CALCULATES
      !     THE CHANGE IN BULK DENSITY, CONVERTS STANDING RESIDUE TO FLAT RESIDUE, AND ESTIMATES THE IMPLEMENT'S EFFECT ON RIDGE HEIGHT AND
      !     INTERVAL.
      ! Ref: APEX-doc Chapter 2.9 Tillage   Eqs. 304 - 308

	! local variables potential problems with local variables (SAVE attribute TST XTP)
      INTEGER:: I, ISM = 58, IBMX, NMIX, J, K, J1, LL, I1, LD2
      REAL:: TST(58),DUM(15),XTMP(15),XTP(8),YTP(8), AD1, AD2, ZLSC1, ZLMC1, ZBMC1, ZHPC1, ZHSC1, & 
            Tillage_Dep, XX, X1, X2, X3, X4, RTO, W1, ZZ, ZO, ZL,  ZON, ZN, DX, ZP, ZOP, ZOK, ZK, F, &
            Mix_Efficiency, SUM, TOT, RE , RX, PMA 
 
      AD1=0.
      ZLSC1=0.
      ZLMC1=0.
      ZBMC1=0.
      ZHSC1=0.
      ZHPC1=0.
      DO I=1,Actual_SoilLayers
          ISL=Layer_ID(I)
          ZLSC1=ZLSC1+CStructLitt(ISL)      ! Total C content of stuctural litter(Kg/ha)                                                        
          ZLMC1=ZLMC1+CMetabLitt(ISL)      ! Total C content of metabolic litter(Kg/ha)                                                           
          ZBMC1=ZBMC1+CBiomass(ISL)           ! Total C content of biomass (Kg/ha)                                                      
          ZHSC1=ZHSC1+CSlowHumus(ISL)        ! Total C CONTENT OF SLOW HUMUS(kg/ha)                                                       
          ZHPC1=ZHPC1+CPassiveHumus(ISL)        ! Total C CONTENT OF passive HUMUS(kg/ha)                                                          
          AD1=AD1+CBiomass(ISL)+CPassiveHumus(ISL)+CSlowHumus(ISL)+CMetabLitt(ISL)+CStructLitt(ISL)
      END DO
      ! ------------- if 1 --------------
      IF(IBMX==0)THEN         ! what is IBMX?
            ! ---------- if 2 -------------
            IF(NMIX==0)THEN     ! what is NMIX?  times of mixing ?
                  ! ------- if 3 -------------
                  IF(Tillage_Dep<0.)THEN  ! Tillage_Dep == ZMAX, how could it be less than zero?
                  !     MOW, SHRED, ETC
                        DO J=1,LC    ! crop loop
                              IF(cropCode(J)==plantCategoryCode(7).OR.cropCode(J)==plantCategoryCode(8).OR. &
                              cropCode(J)==plantCategoryCode(10)) CYCLE                     
                              ! 7, 8 , 10 are trees
                              IF(CPHT(J)+Tillage_Dep<0..OR.abvGroundBiom(J)<.001)CYCLE
                              XX=(CPHT(J)+Tillage_Dep)/CPHT(J)
                              ZZ=XX*standCropResi(J)
                              ZO=XX*standDeadResOrganism
                              ZL=XX*abvGroundBiom(J)
                              standCropResi(J)=standCropResi(J)-ZZ
                              standDeadResOrganism=standDeadResOrganism-ZO
                              abvGroundBiom(J)=abvGroundBiom(J)-ZL
                              X1=1.-XX
                              Current_LAI(J)=Current_LAI(J)*X1
                              HU(J)=HU(J)*X1
                              standDeadLignin=standDeadLignin*X1
                              DX=MIN(.99,ZL/(totCropBio(J)+1.E-10))
                              X1=ZZ+ZL+ZO
                              ZZ=XX*standDeadResiN(J)
                              ZON=XX*standDeadResiOrgN
                              ZN=DX*actualCropN(J)
                              standDeadResiN(J)=standDeadResiN(J)-ZZ
                              standDeadResiOrgN=standDeadResiOrgN-ZON
                              X2=ZZ+ZN+ZON
                              CALL CN_TopSoil(X1, X2, 0)
                              ZZ=XX*standDeadResiP
                              standDeadResiP=standDeadResiP-ZZ
                              ZOP=XX*StandDeadRes_OrgP
                              StandDeadRes_OrgP=StandDeadRes_OrgP-ZOP
                              ZP=DX*actualCropP(J)
                              FreshOrgP_Residu(LD1)=FreshOrgP_Residu(LD1)+ZZ+ZP+ZOP
                              ZZ=XX*standDeadResiK
                              standDeadResiK=standDeadResiK-ZZ
                              ZOK=XX*StandDeadRes_OrgK
                              StandDeadRes_OrgK=StandDeadRes_OrgK-ZOK
                              ZK=DX*actualCropK(J)
                              solubleK(LD1)=solubleK(LD1)+ZZ+ZK+ZOK
                              totCropBio(J)=totCropBio(J)-ZL
                              actualCropN(J)=MAX(1.E-5,actualCropN(J)-ZN)
                              actualCropP(J)=actualCropP(J)-ZP
                              !actualCropK(J)=MAX(1.E-5,actualCropK(J)-ZK)
                              CPHT(J)=-Tillage_Dep
                        END DO
                        RETURN
                  END IF  ! ----------------- end if 3
            END IF  ! --------------------- end if 2 ---------------------   
          
            IF(Z(LD1)>=Tillage_Dep)RETURN
            Tot_Rainfall=0.
            RCF=1.
          
            IF(ridgeHeight(JT1)<ridgeHeight(JT2))THEN
                  RidgeHeight2=ridgeHeight(JT1)+(ridgeHeight(JT2)-ridgeHeight(JT1))*EXP(-Tillage_Dep/tillDepth(JT2))
            ELSE
                  RidgeHeight2=ridgeHeight(JT1)
                  Ridge_IntervalX=ridgeInterval(JT1)
            END IF
            F=1.-EXP(-56.9*Tillage_Dep*Mix_Efficiency)         ! APEX-doc  Eq. 307 
            SUM=0.
            TOT=0.
            DO K=1,LC
                  X1=standCropResi(K)*F                            ! APEX-doc  Eq. 307   ! X1: residues ----> soil
                  SUM=SUM+X1
                  standCropResi(K)=MAX(1.E-10,standCropResi(K)-X1) ! remaining standing residue weights after tillage
                  XX=F*standDeadResiN(K)                           ! XX: amount of N in residue ---> soil
                  TOT=TOT+XX
                  standDeadResiN(K)=MAX(1.E-10,standDeadResiN(K)-XX) ! remaining N in residue
                  XX=F*standDeadResiP                              ! XX: amount of P in residue ---> soil (fresh organic P)
                  standDeadResiP=MAX(1.E-10,standDeadResiP-XX)     ! remaining N in residue
                  FreshOrgP_Residu(LD1)=FreshOrgP_Residu(LD1)+XX
                  XX=F*standDeadResiK                              ! XX: amount of K in residue ---> soil
                  standDeadResiK=MAX(1.E-10,standDeadResiK-XX)
                  solubleK(LD1)=solubleK(LD1)+XX
                  XX=F*standDeadLignin                             ! XX: amount of lignin in residue ---> soil
                  standDeadLignin=MAX(.1*standCropResi(K),standDeadLignin-XX)
            END DO
            XX=standDeadResOrganism*F
            standDeadResOrganism=MAX(1.E-5,standDeadResOrganism-XX)
            X1=SUM+XX                      ! total organism
            ZON=F*standDeadResiOrgN
            standDeadResiOrgN=MAX(1.E-5,standDeadResiOrgN-ZON)
            X2=TOT+ZON                     ! total organic N
            AD1=AD1+420.*X1
          
            CALL CN_TopSoil(X1, X2, 0 )
          
            XX=F*StandDeadRes_OrgP
            StandDeadRes_OrgP=MAX(1.E-5,StandDeadRes_OrgP-XX)
            FreshOrgP_Residu(LD1)=FreshOrgP_Residu(LD1)+XX
            XX=F*StandDeadRes_OrgK
            StandDeadRes_OrgK=MAX(1.E-5,StandDeadRes_OrgK-XX)
            solubleK(LD1)=solubleK(LD1)+XX
            RandRough=MAX(1.,tillRoughness(JT1))
            TLMF=0.
      END IF  ! --------------- end if 1 -----------
      
      DO I=1,ISM
            TST(I)=0.
      END DO
      XX=0.
      XTP(1)=structLitt(LD1)
      XTP(2)=metabLitt(LD1)
      XTP(3)=lgStructLitt(LD1)
      XTP(4)=CStructLitt(LD1)
      XTP(5)=CMetabLitt(LD1)
      XTP(6)=CLgStructLitt(LD1)
      XTP(7)=NStructLitt(LD1)
      XTP(8)=NMetabLitt(LD1)
      
      DO J=1,Actual_SoilLayers
            ISL=Layer_ID(J)
            IF(IBMX>0)Mix_Efficiency=WBMX
            XTMP(ISL)=rockFrac(ISL)
            ZZ=Z(ISL)-XX
            IF(Z(ISL)>=Tillage_Dep)EXIT
            IF(NMIX<=0)THEN                              ! APEX-docEq. 305
                  postTillBulkDensity(ISL)=postTillBulkDensity(ISL)-(postTillBulkDensity(ISL)-.6667*bulkDensity(ISL))*Mix_Efficiency
                  clayFrac(ISL)=clayFrac(ISL)*ZZ
                  siltFrac(ISL)=siltFrac(ISL)*ZZ
                  rockFrac(ISL)=rockFrac(ISL)*ZZ
            END IF
            PMA=activeMineralP(ISL)+labileP(ISL)
            DUM(ISL)=sorbRatioP(ISL)*PMA
            soilSupplyRateP(ISL)=PMA-DUM(ISL)
          
            ! EXTRACT THE FRACTION OF MATERIAL TO BE MIXED AND PLACE IN TST STORAGE
            TST(1)=Soil_Move(NO3_N_Soil(ISL),Mix_Efficiency)+TST(1)
            TST(4)=Soil_Move(NBiomass(ISL),Mix_Efficiency)+TST(4)     ! When ISL = 3 there is a small problem
            TST(5)=Soil_Move(NStructLitt(ISL),Mix_Efficiency)+TST(5)
            TST(6)=Soil_Move(NMetabLitt(ISL),Mix_Efficiency)+TST(6)
            TST(9)=Soil_Move(CBiomass(ISL),Mix_Efficiency)+TST(9)
            TST(14)=Soil_Move(structLitt(ISL),Mix_Efficiency)+TST(14)
            TST(15)=Soil_Move(metabLitt(ISL),Mix_Efficiency)+TST(15)
            TST(16)=Soil_Move(lgStructLitt(ISL),Mix_Efficiency)+TST(16)
            TST(10)=Soil_Move(CStructLitt(ISL),Mix_Efficiency)+TST(10)
            TST(11)=Soil_Move(CMetabLitt(ISL),Mix_Efficiency)+TST(11)
            TST(12)=Soil_Move(CLgStructLitt(ISL),Mix_Efficiency)+TST(12)
            IF(J==1)THEN     ! A copy of the first layer?
                  YTP(1)=structLitt(LD1)
                  YTP(2)=metabLitt(LD1)
                  YTP(3)=lgStructLitt(LD1)
                  YTP(4)=CStructLitt(LD1)
                  YTP(5)=CMetabLitt(LD1)
                  YTP(6)=CLgStructLitt(LD1)
                  YTP(7)=NStructLitt(LD1)
                  YTP(8)=NMetabLitt(LD1)
            END IF
            TST(17)=Soil_Move(SOP(ISL),Mix_Efficiency)+TST(17)
            TST(19)=Soil_Move(labileP(ISL),Mix_Efficiency)+TST(19)
            TST(20)=Soil_Move(activeMineralP(ISL),Mix_Efficiency)+TST(20)
            TST(21)=Soil_Move(FreshOrgP_Residu(ISL),Mix_Efficiency)+TST(21)
            TST(22)=Soil_Move(stableMineralP(ISL),Mix_Efficiency)+TST(22)
            IF(NMIX==0)THEN
                  TST(23)=Soil_Move(clayFrac(ISL),Mix_Efficiency)+TST(23)
                  TST(24)=Soil_Move(siltFrac(ISL),Mix_Efficiency)+TST(24)
                  TST(27)=Soil_Move(rockFrac(ISL),Mix_Efficiency)+TST(27)
            END IF
            TST(25)=Soil_Move(DUM(ISL),Mix_Efficiency)+TST(25)
            TST(26)=Soil_Move(soilSupplyRateP(ISL),Mix_Efficiency)+TST(26)
            TST(28)=Soil_Move(NH3_Weight(ISL),Mix_Efficiency)+TST(28)
            I1=29
            DO I=1,NDP
                  TST(I1)=Soil_Move(pestInSoil(I,ISL),Mix_Efficiency)+TST(I1)
                  I1=I1+1
            END DO 
            XX=Z(ISL)
      END DO   
      
      IF(J<=Actual_SoilLayers)THEN          ! This is the layer that till depth reach
            RTO=(Tillage_Dep-XX)/ZZ
            IF(IBMX>0)Mix_Efficiency=WBMX
            RE=RTO*Mix_Efficiency
            IF(NMIX==0)THEN
                  postTillBulkDensity(ISL)=postTillBulkDensity(ISL)-(postTillBulkDensity(ISL)-.6667*bulkDensity(ISL))*RE
                  clayFrac(ISL)=clayFrac(ISL)*ZZ
                  siltFrac(ISL)=siltFrac(ISL)*ZZ
                  rockFrac(ISL)=rockFrac(ISL)*ZZ
            END IF
            PMA=activeMineralP(ISL)+labileP(ISL)
            DUM(ISL)=sorbRatioP(ISL)*PMA
            soilSupplyRateP(ISL)=PMA-DUM(ISL)
            TST(1)=Soil_Move(NO3_N_Soil(ISL),RE)+TST(1)
            TST(4)=Soil_Move(NBiomass(ISL),RE)+TST(4)
            TST(5)=Soil_Move(NStructLitt(ISL),RE)+TST(5)
            TST(6)=Soil_Move(NMetabLitt(ISL),RE)+TST(6)
            TST(9)=Soil_Move(CBiomass(ISL),RE)+TST(9)
            TST(10)=Soil_Move(CStructLitt(ISL),RE)+TST(10)
            TST(11)=Soil_Move(CMetabLitt(ISL),RE)+TST(11)
            TST(12)=Soil_Move(CLgStructLitt(ISL),RE)+TST(12)
            TST(14)=Soil_Move(structLitt(ISL),RE)+TST(14)
            TST(15)=Soil_Move(metabLitt(ISL),RE)+TST(15)
            TST(16)=Soil_Move(lgStructLitt(ISL),RE)+TST(16)
            TST(17)=Soil_Move(SOP(ISL),RE)+TST(17)
            TST(19)=Soil_Move(labileP(ISL),RE)+TST(19)
            TST(20)=Soil_Move(activeMineralP(ISL),RE)+TST(20)
            TST(21)=Soil_Move(FreshOrgP_Residu(ISL),RE)+TST(21)
            TST(22)=Soil_Move(stableMineralP(ISL),RE)+TST(22)
            IF(NMIX==0)THEN
                  TST(23)=Soil_Move(clayFrac(ISL),RE)+TST(23)
                  TST(24)=Soil_Move(siltFrac(ISL),RE)+TST(24)
                  TST(27)=Soil_Move(rockFrac(ISL),RE)+TST(27)
            END IF
            TST(25)=Soil_Move(DUM(ISL),RE)+TST(25)
            TST(26)=Soil_Move(soilSupplyRateP(ISL),RE)+TST(26)
            TST(28)=Soil_Move(NH3_Weight(ISL),RE)+TST(28)
            I1=29
            DO I=1,NDP
                  TST(I1)=Soil_Move(pestInSoil(I,ISL),RE)+TST(I1)
                  I1=I1+1
            END DO
      ELSE
            J=Actual_SoilLayers
            Tillage_Dep=Z(Layer_ID(Actual_SoilLayers))
      END IF
      
      J1=J-1
      ! COMPUTE MATERIAL PER DEPTH (kg/ha/m)
      DO I=1,ISM
            TST(I)=TST(I)/Tillage_Dep
      END DO
      XX=0.
      DO J=1,J1     
            LL=Layer_ID(J)
            ZZ=Z(LL)-XX
            ! DISTRIBUTE MIXED MATERIAL UNIFORMLY THRU PLOW DEPTH
            NO3_N_Soil(LL)=TST(1)*ZZ+NO3_N_Soil(LL)
            NBiomass(LL)=TST(4)*ZZ+NBiomass(LL)
            NStructLitt(LL)=TST(5)*ZZ+NStructLitt(LL)
            NMetabLitt(LL)=TST(6)*ZZ+NMetabLitt(LL)
            CBiomass(LL)=TST(9)*ZZ+CBiomass(LL)
            CStructLitt(LL)=TST(10)*ZZ+CStructLitt(LL)
            CMetabLitt(LL)=TST(11)*ZZ+CMetabLitt(LL)
            CLgStructLitt(LL)=TST(12)*ZZ+CLgStructLitt(LL)
            structLitt(LL)=TST(14)*ZZ+structLitt(LL)
            metabLitt(LL)=TST(15)*ZZ+metabLitt(LL) ! Result will be a little different from origional EPIC when LL =2 
            lgStructLitt(LL)=TST(16)*ZZ+lgStructLitt(LL)
            NLgStructLitt(LL)=CStructLitt(LL)-CLgStructLitt(LL)
            cropResidu(LL)=.001*(structLitt(LL)+metabLitt(LL))
            SOP(LL)=TST(17)*ZZ+SOP(LL)
            labileP(LL)=TST(19)*ZZ+labileP(LL)
            activeMineralP(LL)=TST(20)*ZZ+activeMineralP(LL)
            FreshOrgP_Residu(LL)=TST(21)*ZZ+FreshOrgP_Residu(LL)
            stableMineralP(LL)=TST(22)*ZZ+stableMineralP(LL)
            DUM(LL)=TST(25)*ZZ+DUM(LL)
            soilSupplyRateP(LL)=TST(26)*ZZ+soilSupplyRateP(LL)
            IF(NMIX==0)THEN
                  rockFrac(LL)=TST(27)+rockFrac(LL)/ZZ
                  clayFrac(LL)=TST(23)+clayFrac(LL)/ZZ
                  siltFrac(LL)=TST(24)+siltFrac(LL)/ZZ
            END IF
            NH3_Weight(LL)=TST(28)*ZZ+NH3_Weight(LL)
            I1=29
            DO I=1,NDP
                  pestInSoil(I,LL)=TST(I1)*ZZ+pestInSoil(I,LL)
                  I1=I1+1
            END DO
            sorbRatioP(LL)=DUM(LL)/(soilSupplyRateP(LL)+DUM(LL))
            RX=MIN(1.,(100.-rockFrac(LL))/(100.-XTMP(LL)))
            fieldCapacity(LL)=fieldCapacity(LL)*RX
            Wilt_Point(LL)=Wilt_Point(LL)*RX
            Porosity(LL)=Porosity(LL)*RX
            CALL WiltPoint_Sat(LL)
            sandFrac(LL)=100.-clayFrac(LL)-siltFrac(LL)
            WT(LL)=bulkDensity(LL)*ZZ*1.E4
            XX=Z(LL)
      END DO
      
      XX=Tillage_Dep-Z(Layer_ID(J1))
      NO3_N_Soil(ISL)=NO3_N_Soil(ISL)+TST(1)*XX
      NBiomass(ISL)=NBiomass(ISL)+TST(4)*XX
      NStructLitt(ISL)=NStructLitt(ISL)+TST(5)*XX
      NMetabLitt(ISL)=NMetabLitt(ISL)+TST(6)*XX
      CBiomass(ISL)=CBiomass(ISL)+TST(9)*XX
      CStructLitt(ISL)=CStructLitt(ISL)+TST(10)*XX
      CMetabLitt(ISL)=CMetabLitt(ISL)+TST(11)*XX
      CLgStructLitt(ISL)=CLgStructLitt(ISL)+TST(12)*XX
      structLitt(ISL)=structLitt(ISL)+TST(14)*XX
      metabLitt(ISL)=metabLitt(ISL)+TST(15)*XX
      lgStructLitt(ISL)=lgStructLitt(ISL)+TST(16)*XX
      NLgStructLitt(ISL)=CStructLitt(ISL)-CLgStructLitt(ISL)
      cropResidu(ISL)=.001*(structLitt(ISL)+metabLitt(ISL))
      SOP(ISL)=SOP(ISL)+TST(17)*XX
      labileP(ISL)=labileP(ISL)+TST(19)*XX
      activeMineralP(ISL)=activeMineralP(ISL)+TST(20)*XX
      FreshOrgP_Residu(ISL)=FreshOrgP_Residu(ISL)+TST(21)*XX
      stableMineralP(ISL)=stableMineralP(ISL)+TST(22)*XX
      DUM(ISL)=DUM(ISL)+TST(25)*XX
      soilSupplyRateP(ISL)=soilSupplyRateP(ISL)+TST(26)*XX
      IF(NMIX==0)THEN
            rockFrac(ISL)=rockFrac(ISL)+TST(27)*XX
            clayFrac(ISL)=clayFrac(ISL)+TST(23)*XX
            siltFrac(ISL)=siltFrac(ISL)+TST(24)*XX
      END IF
      NH3_Weight(ISL)=NH3_Weight(ISL)+TST(28)*XX
      I1=29
      DO I=1,NDP
            pestInSoil(I,ISL)=pestInSoil(I,ISL)+TST(I1)*XX
            I1=I1+1
      END DO
      sorbRatioP(ISL)=DUM(ISL)/(soilSupplyRateP(ISL)+DUM(ISL))
      ZZ=Z(ISL)-Z(Layer_ID(J1))
      IF(NMIX==0)THEN
            rockFrac(ISL)=rockFrac(ISL)/ZZ
            clayFrac(ISL)=clayFrac(ISL)/ZZ
            siltFrac(ISL)=siltFrac(ISL)/ZZ
      END IF
      IF(XTMP(ISL)>0.)THEN
            RX=MIN(1.,(100.-rockFrac(ISL))/(100.-XTMP(ISL)))
            fieldCapacity(ISL)=fieldCapacity(ISL)*RX
            Wilt_Point(ISL)=Wilt_Point(ISL)*RX
            Porosity(ISL)=Porosity(ISL)*RX
            CALL WiltPoint_Sat(ISL)
      END IF
      sandFrac(ISL)=100.-clayFrac(ISL)-siltFrac(ISL)
      WT(ISL)=bulkDensity(ISL)*ZZ*1.E4
      AD2=0.
      totCStructLitt=0.
      totCMetabLitt=0.
      totCBiomass=0.
      totCSlowHumus=0.
      totCPassiveHumus=0.
      
      DO I=1,Actual_SoilLayers
            ISL=Layer_ID(I)
            totCStructLitt=totCStructLitt+CStructLitt(ISL)                                                              
            totCMetabLitt=totCMetabLitt+CMetabLitt(ISL)                                                              
            totCBiomass=totCBiomass+CBiomass(ISL)                                                              
            totCSlowHumus=totCSlowHumus+CSlowHumus(ISL)                                                              
            totCPassiveHumus=totCPassiveHumus+CPassiveHumus(ISL)                                                              
            AD2=AD2+CBiomass(ISL)+CPassiveHumus(ISL)+CSlowHumus(ISL)+CMetabLitt(ISL)+CStructLitt(ISL)
      END DO
 
      IF(Mix_Efficiency<1.)RETURN
      LD2=Layer_ID(2)
      X1=abvGroundBiom(Crop_Num)+standCropResi(Crop_Num)+standDeadResOrganism
      totCropBio(Crop_Num)=totCropBio(Crop_Num)-abvGroundBiom(Crop_Num)
      XX=abvGroundBiom(Crop_Num)/(totCropBio(Crop_Num)+1.E-10)
      X2=XX*actualCropN(Crop_Num)
      X3=XX*actualCropP(Crop_Num)
      W1=XX*actualCropK(Crop_Num)
      X4=standDeadResiN(Crop_Num)+standDeadResiOrgN+X2
      
      CALL CN_TopSoil(X1, X4, 0 )
      
      NMetabLitt(LD2)=NMetabLitt(LD2)+NMetabLitt(LD1)
      NMetabLitt(LD1)=0.
      NStructLitt(LD2)=NStructLitt(LD2)+NStructLitt(LD1)
      NStructLitt(LD1)=0.
      structLitt(LD2)=structLitt(LD2)+structLitt(LD1)
      structLitt(LD1)=0.
      metabLitt(LD2)=metabLitt(LD2)+metabLitt(LD1)
      metabLitt(LD1)=0.
      lgStructLitt(LD2)=lgStructLitt(LD2)+lgStructLitt(LD1)
      lgStructLitt(LD1)=0.
      CStructLitt(LD2)=CStructLitt(LD2)+CStructLitt(LD1)
      CStructLitt(LD1)=0.
      CMetabLitt(LD2)=CMetabLitt(LD2)+CMetabLitt(LD1)
      CMetabLitt(LD1)=0.
      CLgStructLitt(LD2)=CLgStructLitt(LD2)+CLgStructLitt(LD1)
      CLgStructLitt(LD1)=0.
      NLgStructLitt(LD2)=NLgStructLitt(LD2)+NLgStructLitt(LD1)
      NLgStructLitt(LD1)=0.
      FreshOrgP_Residu(LD2)=FreshOrgP_Residu(LD2)+FreshOrgP_Residu(LD1)+standDeadResiP+StandDeadRes_OrgP+X3
      FreshOrgP_Residu(LD1)=0.
      labileP(LD2)=labileP(LD2)+labileP(LD1)
      labileP(LD1)=0.
      solubleK(LD2)=solubleK(LD2)+standDeadResiK+StandDeadRes_OrgK+W1
      cropResidu(LD1)=0.
      abvGroundBiom(Crop_Num)=0.
      standDeadResOrganism=0.
      standCropResi(Crop_Num)=0.
      standDeadResiN(Crop_Num)=0.
      standDeadResiOrgN=0.
      standDeadLignin=0.
      standDeadResiP=0.
      StandDeadRes_OrgP=0.
      standDeadResiK=0.
      StandDeadRes_OrgK=0.
      actualCropN(Crop_Num)=MAX(1.E-5,actualCropN(Crop_Num)-X2)
      actualCropP(Crop_Num)=actualCropP(Crop_Num)-X3
      !actualCropK(Crop_Num)=MAX(1.E-5,actualCropK(Crop_Num)-W1)
END SUBROUTINE Tillage_Mix 

! ------------------------- RUE convert to residue when crops are killed -----------------
SUBROUTINE Bio2_Res
      IMPLICIT NONE
      !     EPIC1102
      !     THIS SUBPROGRAM CONVERTS LIVE BIOMASS TO RESIDUE WHEN A CROP IS KILLED.
	! local variables
      REAL:: X1, X2, X3, XX, W1 
      INTEGER:: J
      
      abvGroundBiom(Crop_Num)=totCropBio(Crop_Num)-totRootWeight(Crop_Num)
      standCropResi(Crop_Num)=standCropResi(Crop_Num)+abvGroundBiom(Crop_Num)
      X1=abvGroundBiom(Crop_Num)+totRootWeight(Crop_Num)
      XX=actualCropN(Crop_Num)/X1
      X3=actualCropP(Crop_Num)/X1
      W1=actualCropK(Crop_Num)/X1
      standDeadResiN(Crop_Num)=standDeadResiN(Crop_Num)+XX*abvGroundBiom(Crop_Num)
      standDeadResiP=standDeadResiP+X3*abvGroundBiom(Crop_Num)
      standDeadResiK=standDeadResiK+W1*abvGroundBiom(Crop_Num)
      standDeadLignin=standDeadLignin+CLG*abvGroundBiom(Crop_Num)
      DO J=1,Layer_RD_Reach
            ISL=Layer_ID(J)
            X1=RWT(ISL,Crop_Num)
            X2=X1*XX
            CALL CN_TopSoil(X1, X2, 1 )
            FreshOrgP_Residu(ISL)=FreshOrgP_Residu(ISL)+X1*X3
            solubleK(ISL)=solubleK(ISL)+X1*W1
            RWT(ISL,Crop_Num)=0.
      END DO
      YLD(Crop_Num)=0.
      totCropBio(Crop_Num)=0.
      abvGroundBiom(Crop_Num)=0.
      actualCropN(Crop_Num)=0.
      UNMX(Crop_Num)=0.
      actualCropP(Crop_Num)=0.
      UPMX(Crop_Num)=0.
      actualCropK(Crop_Num)=0.
      totRootWeight(Crop_Num)=0.
      RD(Crop_Num)=0.
      RETURN
END SUBROUTINE Bio2_Res
 
! ---------- Crop operations --------------------------
SUBROUTINE Till_OPS(CSTX,COX,JRT)             
      !     EPIC1102
      !     THIS SUBPROGRAM CONTROLS ALL TILLAGE OPERATIONS INCLUDING PLANTING, HARVESTING, 
      !     AND AUTOMATIC FERTILIZER APPLICATIONS AT PLANTING.
      IMPLICIT NONE
      ! local variables potential problems with local variables (SAVE attribute YTP)
      REAL:: YTP(16), CSTX, COX, Mix_Efficiency, Tillage_Dep, RWL, HVWC, RNR, RPR, RKR, RLR,&
            F, FT, GCOW, YLSD, X1,X2,X3, X4, X5, X6, X7, X8, X9, X10, X11,XX, XZ, YZ, ZZ, Y4, Y5, Y6, YY, &
            Z2, Z3, Z4, RTO
      INTEGER:: JRT, JHV, II, NN, N1, J, K, I2
 
	! define a function or formula
      !FNPP(X)=maxLAI(Crop_Num)*X/(X+EXP(coeffPopuCurve(1,Crop_Num)-coeffPopuCurve(2,Crop_Num)*X))
      
      JRT=0
      II=currentOps(JT1)
      NN=NBC(IRO)
      N1=MAX(1,NCP(Crop_Num))
      X1=CND(IRO,KT)
      IF(ABS(X1-CN0)>0.)THEN
            X2=SMX
            CALL Adjust_CN(X1, X3)
            CN0=X1
            CN2=X1
            SCI=SMX*SCI/X2
      END IF
      IF(II==operationCode(1).OR.II==operationCode(2).OR.II==operationCode(3).OR.II==operationCode(19)&
      .OR.II==operationCode(22))GO TO 10
      IF(II==operationCode(5))GO TO 61
      IF(II==operationCode(6))GO TO 53
      IF(II==operationCode(7).OR.II==operationCode(8))GO TO 57
      IF(II==operationCode(10))GO TO 51
      IF(II==operationCode(11))GO TO 52
      IF(II==operationCode(12).OR.II==operationCode(13))GO TO 56
      IF(II==operationCode(14))GO TO 59
      IF(II==operationCode(23))GO TO 63
      IF(II==operationCode(24))GO TO 64
      GO TO 6
   51 CSTX=-CSTX*YLD1(N1,Crop_Num)/(1.-waterFracYield(Crop_Num))
      COX=CSTX
      GO TO 57
   52 CSTX=-CSTX*YLD(Crop_Num)/(1.-waterFracYield(Crop_Num))
      COX=CSTX
      GO TO 57
   56 IF(ICUS(JT1)==0)GO TO 57
      CSTX=-CSTX*YLD(Crop_Num)/(1.-waterFracYield(Crop_Num))
      COX=CSTX
      GO TO 57
   53 IDRL=1
   61 ISL=Layer_ID(2)            ! Why the second layer?
      DO K=1,NN
            I2=LY(IRO,K)
            IF(JH(IRO,KT)==cropID(I2))GO TO 27
      END DO
      GO TO 26
   27 IF(KG(I2)>0)GO TO 26
      IF(soilTem(ISL)<plantMinT(I2)+2..AND.MO<12)THEN
            KOMP(KT)=0
            JRT=1
            RETURN
      END IF
      AQV=0.
      ARF=0.
      IGO=IGO+1
      DO KC=1,NN
            IF(JE(KC)>MNC)EXIT
      END DO
      JE(KC)=I2
      Crop_Num=I2
      AWC(Crop_Num)=soilWaterRZ
      KG(Crop_Num)=Crop_Num
      JP(Crop_Num)=0
      IYH(Crop_Num)=1
      GSEP=0.
      GSVP=0.
      SRA=0.
      sumActuralEP(Crop_Num)=0.
      sumPotentialEP(Crop_Num)=0.
      ACET(Crop_Num)=0.
      XDLAI(Crop_Num)=fracGrowSeasonLAIDecline(Crop_Num)
      XDLA0(Crop_Num)=0.
      WaterFrac_Yield=.3
      standDeadResOrganism=standDeadResOrganism+standCropResi(Crop_Num)
      StandDeadRes_OrgK=StandDeadRes_OrgK+standDeadResiK
      standDeadResiOrgN=standDeadResiOrgN+standDeadResiN(Crop_Num)
      StandDeadRes_OrgP=StandDeadRes_OrgP+standDeadResiP
      standCropResi(Crop_Num)=0.
      standDeadResiK=0.
      standDeadLignin=0.
      standDeadResiN(Crop_Num)=0.
      standDeadResiP=0.
      RD(Crop_Num)=tillDepth(JT1)
      ROSP=ridgeInterval(JT1)
      HU(Crop_Num)=0.
      totCropBio(Crop_Num)=seedRate(Crop_Num)*5.E-4
      totCropBioX(Crop_Num)=totCropBio(Crop_Num)
      SM99=420.*totCropBio(Crop_Num)
      VAR(99)=VAR(99)+SM99
      SMM(99,MO)=SMM(99,MO)+SM99                                                                             
      totRootWeight(Crop_Num)=.4*totCropBio(Crop_Num)
      RWT(ISL,Crop_Num)=totRootWeight(Crop_Num)
      CPHT(Crop_Num)=0.
      PPL0(Crop_Num)=POP(Crop_Num,IHU(Crop_Num))
      XLAI(Crop_Num)=maxLAI(Crop_Num)*PPL0(Crop_Num)/(PPL0(Crop_Num)+EXP(coeffPopuCurve(1,Crop_Num)-coeffPopuCurve(2,Crop_Num)&
                                                                          *PPL0(Crop_Num)))
      DMLX(Crop_Num)=XLAI(Crop_Num)
      X1=seedRate(Crop_Num)*seedCost(Crop_Num)
      COST=COST+X1
      SMM(96,MO)=SMM(96,MO)+carbonEmissionSeedWeight(Crop_Num)*seedRate(Crop_Num)
      Layer_RD_Reach=2
      JPL(Crop_Num)=1
      IF(NCP(Crop_Num)==0)NCP(Crop_Num)=1
      N1=MAX(1,NCP(Crop_Num))
      IPLD(N1,Crop_Num)=IYR*10000+MO*100+DayOfMon
      IHVD(N1,Crop_Num)=0
      IF(NOP>0)WRITE(KW(1),32)IYR,MO,DayOfMon,Crop_Name(Crop_Num),groundCover,X1
      IF(KFL(20)>0)WRITE(KW(20),49)IYR,MO,DayOfMon,equipmentName(JT1),cropID(Crop_Num),II,&
      NBE(JT1),NBT(JT1),X1,X1,seedRate(Crop_Num)
    6 Mix_Efficiency=mixEfficiency(JT1)
      PPL0(Crop_Num)=(1.-PlantRedu_frac(JT1))*PPL0(Crop_Num)
      XLAI(Crop_Num)=maxLAI(Crop_Num)*PPL0(Crop_Num)/(PPL0(Crop_Num)+EXP(coeffPopuCurve(1,Crop_Num)-coeffPopuCurve(2,Crop_Num)&
                                                                          *PPL0(Crop_Num)))
      DMLX(Crop_Num)=XLAI(Crop_Num)
      Tillage_Dep=tillDepth(JT1)
      IF(II/=operationCode(19).AND.II/=operationCode(2)) CALL Tillage_Mix(Mix_Efficiency,Tillage_Dep, 0, 0)
      
      IF(Tillage_Dep>BIG)tillDepth(JT1)=BIG
      IF(II==operationCode(15))THEN
            satuateCondv(Layer_ID(2))=modelPara(33)
      ELSE
            IF(II==operationCode(16))satuateCondv(Layer_ID(2))=satuateCondv0
      END IF
   57 IF(IDR>0)THEN
            IF(II==operationCode(25))THEN
                  lateralFlowCondv(IDR)=lateralCondvHyN
            ELSE
                  IF(II==operationCode(26))lateralFlowCondv(IDR)=lateralCondvHyD
            END IF
      END IF
      IF(KFL(20)>0)WRITE(KW(20),50)IYR,MO,DayOfMon,equipmentName(JT1),cropID(Crop_Num),II,&
      NBE(JT1),NBT(JT1),CSTX,COX,Fuel_Use(JT1)
      SMM(92,MO)=SMM(92,MO)+Fuel_Use(JT1) 
      IF(II==operationCode(2).OR.II==operationCode(3).OR.II==operationCode(19))GO TO 26
    7 XX=tillDepth(JT1)*1000.
      IF(NOP>0)WRITE(KW(1),28)IYR,MO,DayOfMon,equipmentName(JT1),XX,XHSM,CSTX
      IF(II/=operationCode(17).AND.II/=operationCode(18))GO TO 26
      IF(II/=operationCode(18))THEN
            dikeHeight=furrowHeight(JT1)
            DKHL=dikeHeight
            dikeInterval=furrowInterval(JT1)
            furrowFlag=1
      ELSE
            furrowFlag=0
      END IF
      IF(NOP>0)WRITE(KW(1),30)dikeHeight,dikeInterval,XHSM
      GO TO 26
      59 CALL Burn_Residue 
      GO TO 7
   63 ICV=1
      GO TO 57
   64 ICV=0
      GO TO 57
   10 DO K=1,NN
            IF(JE(K)>MNC)CYCLE
            IF(JH(IRO,KT)==cropID(JE(K)))EXIT
      END DO
      IF(K>NN)THEN
            CSTX=0.
            COX=0.
            JRT=1
            RETURN
      END IF
      Crop_Num=JE(K)
      RWL=totRootWeight(Crop_Num)
      IF(cropCode(Crop_Num)==plantCategoryCode(7).OR.cropCode(Crop_Num)==plantCategoryCode(8).OR. &
                                                       cropCode(Crop_Num)==plantCategoryCode(10))THEN
            IF(IYH(Crop_Num)/=LYR(IRO,KT).AND.LYR(IRO,KT)/=1)THEN
                  KOMP(KT)=1
                  RETURN
            END IF
      END IF
      N1=MAX(1,NCP(Crop_Num))
      JHV=K
      KHV=1
      IF(II==operationCode(1))GO TO 22
      IF(II/=operationCode(2).AND.II/=operationCode(19).AND.II/=operationCode(22))THEN
      IF(IHT(JT1)>0)GO TO 26
            IHT(JT1)=1
      END IF
      HVWC=HWC(IRO,KT)                                                                 
      IF(HVWC>0..AND.WaterFrac_Yield>HVWC.AND.MO<12)THEN
            JRT=1
            KOMP(KT)=0
            RETURN
      END IF
      IF(JP(Crop_Num)==0)THEN
            JP(Crop_Num)=1
            IF(II/=operationCode(3))NCR(Crop_Num)=NCR(Crop_Num)+1
      END IF
      HUF(Crop_Num)=MAX(HUF(Crop_Num),HU(Crop_Num))
      DMF(N1,Crop_Num)=totCropBioX(Crop_Num)
      TRA(Crop_Num)=SRA
      IF(RD(Crop_Num)>RDF(Crop_Num))RDF(Crop_Num)=RD(Crop_Num)
      X9=totCropBio(Crop_Num)+.001
      X2=MAX(uptakeParaN(3,Crop_Num),actualCropN(Crop_Num)/X9)
      X7=.001*X2
      X3=actualCropP(Crop_Num)/X9
      X8=actualCropK(Crop_Num)/X9
      XX=standCropResi(Crop_Num)+1.E-10
      RNR=standDeadResiN(Crop_Num)/XX
      RPR=standDeadResiP/XX
      RKR=standDeadResiK/XX
      standDeadLignin=CLG*XX
      RLR=MIN(.8,standDeadLignin/(standCropResi(Crop_Num)+1.E-5))
      IF(overRideHI(JT1)<1.E-10)THEN
            IF(cropCode(Crop_Num)==plantCategoryCode(8))THEN
                  F=1.
            ELSE
                  XX=100.*sumActuralEP(Crop_Num)/(sumPotentialEP(Crop_Num)+1.E-5)
                  F=XX/(XX+EXP(S_Curve(10,1)-S_Curve(10,2)*XX))
            END IF
            XX=MAX(AJHI(Crop_Num)-lowerLimitHI(Crop_Num),0.)
            FT=MAX(.1,1.+modelPara(50)*(IYR-2000))
            X1=(F*XX+lowerLimitHI(Crop_Num))*FT
            IF(cropCode(Crop_Num)==plantCategoryCode(8))THEN
                  X2=modelPara(76)/AWC(Crop_Num)
                  X1=MIN(HI(Crop_Num),X1*X2)
            END IF
            X2=1000.*fracYieldN(Crop_Num)*(X7/uptakeParaN(3,Crop_Num))**.1
            X3=1000.*fracYieldP(Crop_Num)*(.001*X3/uptakeParaP(3,Crop_Num))**.1
            GO TO 17
      END IF
      IF(II/=operationCode(19).AND.II/=operationCode(22))GO TO 16
      IF(cropCode(Crop_Num)==plantCategoryCode(7).OR.cropCode(Crop_Num)==plantCategoryCode(8).OR. &
                                                      cropCode(Crop_Num)==plantCategoryCode(10))THEN
            KTT=0
            GO TO 26
      END IF
      KOMP(KT)=0
      KTT=KT
      IF(II==operationCode(22))THEN
            XX=CPHT(Crop_Num)-HMO(JT1)
            IF(XX<.001.OR.Mow_DayLap<mowMinInterval)RETURN
            X1=XX/CPHT(Crop_Num)
            ZZ=HMO(JT1)/CPHT(Crop_Num)
            YZ=ZZ*standCropResi(Crop_Num)
            Mow_DayLap=0
            GO TO 45
      END IF
      IF(standLiveDeadPlant<grazeLimit)RETURN
      GCOW=areaWshed/baseStockRateArr(IRO,KT)
      XX=GCOW*overRideHI(JT1)/(areaWshed*harvestEffi(JT1))
      X1=MIN(XX/standLiveDeadPlant,.9)
      GO TO 17
   16 X1=overRideHI(JT1)
      IF(tillDepth(JT1)<=0.)GO TO 17
      YLSD=0.
      CALL Har_RootCrop(YY,X3,X1,X6,X7,N1)
      GO TO 25
   17 ZZ=MAX(.01,1.-X1)
      YZ=X1*standCropResi(Crop_Num)
   45 XZ=X1*abvGroundBiom(Crop_Num)
      HIF(N1,Crop_Num)=X1
      AJHI(Crop_Num)=0.
      CPHT(Crop_Num)=CPHT(Crop_Num)*ZZ
      IF(II==operationCode(19))THEN
            HU(Crop_Num)=HU(Crop_Num)*ZZ
      ELSE
            HU(Crop_Num)=HU(Crop_Num)*modelPara(69)
      END IF
      Current_LAI(Crop_Num)=Current_LAI(Crop_Num)*ZZ
      abvGroundBiom(Crop_Num)=abvGroundBiom(Crop_Num)*ZZ
      standCropResi(Crop_Num)=MAX(1.E-10,standCropResi(Crop_Num)-YZ)
      standDeadLignin=standDeadLignin*ZZ
      CALL Pest_Damage 
	  
      TPSF(N1,Crop_Num)=TPSF(N1,Crop_Num)+PSTF(Crop_Num)
      NPSF(N1,Crop_Num)=NPSF(N1,Crop_Num)+1
      !================================ Crop Yield ======================================
      YLD(Crop_Num)=XZ*harvestEffi(JT1)*PSTF(Crop_Num)    ! Yield is calcuated! 
      YLSD=YZ*harvestEffi(JT1)
      Y4=YZ*RNR
      Y5=YZ*RPR
      Y6=YZ*RKR
      standDeadResiN(Crop_Num)=MAX(1.E-5,standDeadResiN(Crop_Num)-Y4)
      standDeadResiP=MAX(1.E-5,standDeadResiP-Y5)
      standDeadResiK=MAX(1.E-5,standDeadResiK-Y6)
      standDeadLignin=MAX(standDeadLignin-YZ*RLR,.1*standCropResi(Crop_Num))
      X6=PSTF(Crop_Num)
      X4=MIN(XZ*X2,actualCropN(Crop_Num))
      X5=MIN(XZ*X3,actualCropP(Crop_Num))
      X11=XZ-YLD(Crop_Num)+YZ-YLSD
      Z2=YLSD*RNR
      Z3=YLSD*RPR
      Z4=YLSD*RKR
      Yield_N=MIN(.9*(actualCropN(Crop_Num)+standDeadResiN(Crop_Num)),YLD(Crop_Num)*X2+Z2)
      YLP=MIN(.9*(actualCropP(Crop_Num)+standDeadResiP),YLD(Crop_Num)*X3+Z3)
      X10=X4-Yield_N+Y4
      CALL CN_TopSoil(X11,X10,0)
      FreshOrgP_Residu(LD1)=FreshOrgP_Residu(LD1)+X5-YLP+Y5
      YY=YLD(Crop_Num)+YLSD
      YLD(Crop_Num)=YY
      YLC=.42*YY
      IF(overRideHI(JT1)>0.)THEN
            YLD2(N1,Crop_Num)=YLD2(N1,Crop_Num)+YY
            X11=YY*forageYieldPrice(Crop_Num)
      ELSE
            IF(cropCode(Crop_Num)==plantCategoryCode(9))THEN
                  YLD1(N1,Crop_Num)=YLD1(N1,Crop_Num)+turnoutFracCottn(Crop_Num)*YY
                  YLD2(N1,Crop_Num)=YLD2(N1,Crop_Num)+YLD1(N1,Crop_Num)*(1./lintFracCottn(Crop_Num)-1.)
            ELSE
                  YLD1(N1,Crop_Num)=YLD1(N1,Crop_Num)+YY
                  X10=YY*seedYieldPrice(Crop_Num)
            END IF            
      END IF
      JD=Crop_Num
      GS_Rad=SRA
      SRA=0.
      actualCropN(Crop_Num)=MAX(1.E-5,actualCropN(Crop_Num)-X4)
      actualCropP(Crop_Num)=actualCropP(Crop_Num)-X5
      totCropBio(Crop_Num)=totCropBio(Crop_Num)-XZ
      totRootWeight(Crop_Num)=totCropBio(Crop_Num)-abvGroundBiom(Crop_Num)
      RTO=totRootWeight(Crop_Num)/RWL
      DO J=1,Layer_RD_Reach
            ISL=Layer_ID(J)
            RWT(ISL,Crop_Num)=RWT(ISL,Crop_Num)*RTO
      END DO
      X3=MAX(0.,totRootWeight(Crop_Num))
      RWF(N1,Crop_Num)=X3
      YLNF(N1,Crop_Num)=YLNF(N1,Crop_Num)+Yield_N
      YLPF(N1,Crop_Num)=YLPF(N1,Crop_Num)+YLP
      YLKF(N1,Crop_Num)=YLKF(N1,Crop_Num)+YLK
      YLCF(N1,Crop_Num)=YLCF(N1,Crop_Num)+YLC
      GO TO 25
   22 IF(printChoice==5)THEN
            CALL Prep_SoilVar(YTP)
            WRITE(KW(1),'(T5,A)')'SOIL DATA'
            CALL Print_SoilVar2(YTP,1)
      END IF
      
      CALL Bio2_Res 
      NCP(Crop_Num)=MIN(NBCX(IRO,Crop_Num),NCP(Crop_Num)+1)
      IF(YLD1(N1,Crop_Num)+YLD2(N1,Crop_Num)<1.E-10)NCR(Crop_Num)=NCR(Crop_Num)+1
      JE(JHV)=MNC+1
      GS_PlantTranspiration=GSEP
      Var_GS(4)=Var_GS(4)+GS_PlantTranspiration
      GSEP=0.
      GS_VPD=GSVP                ! GS_VPD: growing season VPD (kPa)
      Var_GS(6)=Var_GS(6)+GS_VPD
      GSVP=0.
      IGO=MAX(0,IGO-1)
      KG(Crop_Num)=0
      IYH(Crop_Num)=0
      JPL(Crop_Num)=0
      HU(Crop_Num)=0.
      HUI(Crop_Num)=0.
      HSM=0.
      Current_LAI(Crop_Num)=0.
      WLV(Crop_Num)=0.
      ANA(Crop_Num)=0.
      NFA=0
      NII=autoIrrInterval
      CSTF(N1,Crop_Num)=COST
      CSOF(N1,Crop_Num)=COST-CSFX
      COST=0.
      CSFX=0.
      IHU(Crop_Num)=IHU(Crop_Num)+1
      IF(IHU(Crop_Num)>NHU(Crop_Num))IHU(Crop_Num)=1
      availWaterCrop(N1,Crop_Num)=AWC(Crop_Num)
      CQV(N1,Crop_Num)=AQV
      CRF(N1,Crop_Num)=ARF 
      ET_GrowingSeason(N1,Crop_Num)=ACET(Crop_Num)+ET_GrowingSeason(N1,Crop_Num)
      pestGrowIndx=0.
      GrowingSeason_Days=0
      waterStress=1.
      GroundCover_Frac=0.
      GroundCover_StandLiveBio=0.
      wiltPoint=0.
      soilSupplyRateN=0.
      cropVar=0.
      GO TO 6
   25 TYN=TYN+Yield_N
      TYP=TYP+YLP
      TYK=TYK+YLK
      TYC=TYC+YLC
      IHVD(N1,Crop_Num)=IYR*10000+MO*100+DayOfMon
      IF(ICUS(JT1)/=0.AND.CSTX<=0.)THEN
            CSTX=-CSTX*YLD(Crop_Num)
            COX=CSTX
      END IF
      IF(NOP>0.AND.II/=operationCode(19))WRITE(KW(1),29)IYR,MO,DayOfMon,equipmentName(JT1),&
      Crop_Name(JD),YY,YLSD,standLiveDeadPlant,totCropBioX(Crop_Num),X3,X1,X6,X7,WaterFrac_Yield,XHSM,CSTX
      GO TO 6
   26 IF(NOP>0.AND.II==operationCode(19))WRITE(KW(1),62)IYR,MO,DayOfMon,equipmentName(JT1),&
      Crop_Name(Crop_Num),YLD(Crop_Num),YLSD,standLiveDeadPlant,abvGroundBiom(Crop_Num),X3,X1,XHSM
      IF(KFL(29)>0.AND.II==operationCode(19))WRITE(KW(29),31)IYR,MO,DayOfMon,equipmentName(JT1),&
      Crop_Name(Crop_Num),totCropBio(Crop_Num),X3,Current_LAI(Crop_Num),abvGroundBiom(Crop_Num),standLiveDeadPlant,X1,&
      YLD(Crop_Num),YLSD,XHSM 
 
   28 FORMAT(1X,3I4,2X,A8,2X,'DPTH = ',F5.0,'mm',2X,'fracHeatUnit = ',F6.2,2X&
      ,'COST=',F7.0,'$/ha')
   29 FORMAT(1X,3I4,2X,A8,2X,A4,2X,'YLD=',F7.2,'t/ha',2X,'YLSD=',F7.2,&
      't/ha', 2X,'standLiveDeadPlant=',F7.2,'t/ha',2X,'BIOM=',F7.2,'t/ha',2X,'totRootWeight=',F7.2,&
      't/ha',2X,'HI=',F6.2,2X,'PSTF=',F5.2,2X,'NCN=',F6.3,'G/G',2X,'waterFracYield=',&
      F5.2,2X,'fracHeatUnit=',F5.2,2X,'COST=',F7.0,'$/ha')
   30 FORMAT('+',T45,'furrowHeight=',F5.0,'MM',2X,'furrowInterval=',F6.2,'M',2X,'fracHeatUnit=',F6.2)
   31 FORMAT(1X,3I4,2X,A8,2X,A4,9F10.3)
   32 FORMAT(1X,3I4,2X,A4,2X,'cropResidu = ',F5.1,'T',2X,'COST=',F7.0,'$/ha')
   49 FORMAT(1X,3I4,2X,A8,8X,I6,6X,3I4,F10.2,10X,3F10.2)
   50 FORMAT(1X,3I4,2X,A8,8X,I6,6X,3I4,2F10.2,20X,F10.2)
   62 FORMAT(1X,3I4,2X,A8,2X,A4,2X,'YLD=',F7.4,'t/ha',2X,'Sediment=',F7.2,&
            &'t/ha',2X,'standLiveDeadPlant=',F7.2,'t/ha',2X,'abvGroundBiom=',F7.2,'t/ha',2X,'RWT=',F7.2,&
            &'t/ha',2X,'HI=',F7.3,2X,'fracHeatUnit=',F6.2)
END SUBROUTINE Till_OPS

END MODULE Management