MODULE Weather_Generator
USE PARM
USE MisceFun_Module

IMPLICIT NONE
CONTAINS
! *********************** Parameters for simulator ***********************    
SUBROUTINE Para_SimuWeather    
      !     EPIC1102
      !     THIS SUBPROGRAM SIMULATES DAILY PRECIPITATION, MAXIMUM AND MINUMUM
      !     IrrWater TEMPERATURE, SOLAR RADIATION AND RELATIVE HUMIDITY.  ALSO
      !     ALSO PROVIDES AND OPTIONS TO SIMULATE VARIOUS COMBINATIONS
      !     GIVEN DAILY PRECIPITATION.
 IMPLICIT NONE
     
      ! local variables potential problems with local variables (SAVE attribute XX E)
      INTEGER:: I, J 
      REAL, DIMENSION(3,3):: A(3,3) = RESHAPE([.594,.454,-.004,.076,.261,-.037,-.018,-.129,.222], [3,3]), &
                             B(3,3) = RESHAPE([.767,.304,.274,0.,.692,-.33,0.,0.,.873], [3,3])
      REAL:: XX(3), E(3), XXX, III, Z2, YY, V2
 
      XXX=.5*(TMXM-TMNM)
      III=1
      Z2=wetDayFreq(NWI,MO)
      YY=.9*Z2
      TXXM=TMXM+XXX*Z2
      RHM=(RH(NWI,MO)-YY)/(1.-YY)       ! APEX-dco Eq. 22 : RH: long-term  average RH
      IF(RHM<.05)RHM=.5*RH(NWI,MO)
      RM=SRAM/(1.-.25*Z2)
      IF(Rainfall>2.)THEN
          TXXM=TXXM-XXX
          RM=.5*RM
          RHM=RHM*.1+.9
      END IF
      
      DO I=1,3
          V2=Generate_Random(randID(2))
          E(I)= SQRT(-2.*LOG(V1))*COS(6.283185*V2)
          V1=V2
      END DO
      DO I=1,3
          WX(I)=0.
          XX(I)=0.
              DO J=1,3
                  WX(I)=WX(I)+B(I,J)*E(J)
                  XX(I)=XX(I)+A(I,J)*XIM(J)
              END DO
      END DO
      DO I=1,3
          WX(I)=WX(I)+XX(I)
          XIM(I)=WX(I)
      END DO
   
END SUBROUTINE Para_SimuWeather    
!****************************** 1. Rainfall *****************************
!------------ 1.simulate rainfall using skewed distribution -------------
REAL FUNCTION Sim_Rain(R6,X)
      !     EPIC1102
      !     THIS SUBPROGRAM COMPUTES DAILY PRECIPITATION AMOUNT FROM SKEWED
      !     NORMAL DISTRIBUTION.
      IMPLICIT NONE
      ! local variables
      REAL:: XLV, R6, X 
      
      XLV=(X-R6)*R6+1.
      XLV=(XLV**3-1.)*2./rainfallSkew
      Sim_Rain=XLV*rainfallStd+rainfallVolum
      IF(Sim_Rain<.01)Sim_Rain=.01
      RETURN
END FUNCTION Sim_Rain
    
! Generate rainfall
SUBROUTINE Generate_Rainfall(IEX)
      IMPLICIT NONE
      INTEGER, INTENT(IN):: IEX
      ! local variables
      REAL:: RN, ZZ, V4, R6
    
      IF(IEX==0)THEN
            RN=1.-Generate_Random(randID(1))
            IF(RN>PBWM+.001)THEN
                  Rainfall=0.
                  LW=1
                  RETURN
            END IF
      END IF
      V4=Generate_Random(randID(3))
      IF(rainDistriFlag==0)THEN
            R6=rainfallSkew/6.
            ZZ=SQRT(-2.*LOG(V3))*COS(6.283185*V4)
            Rainfall=Sim_Rain(R6,ZZ)*PCF(NWI,MO)
            V3=V4
      ELSE
            Rainfall=rainfallVolum*(-LOG(V4))**rainPara
      END IF
      LW=2
 
END SUBROUTINE Generate_Rainfall

! **************************** 2. Temperature *****************************
! generate temperature
SUBROUTINE Generate_Tem 
      IMPLICIT NONE

      IF(TMX<100.)THEN
            TMN=MIN(TMNM+TNSD*WX(2),TMX-.2*ABS(TMX))
      ELSE
            TMX=MAX(TXXM+TXSD*WX(1),TMN+.2*ABS(TMN))
      END IF

END SUBROUTINE Generate_Tem
    
!  generate air tempeature 
    
SUBROUTINE Generate_AirTem 
      IMPLICIT NONE 
      
      TMX=TXXM+TXSD*WX(1)
      TMN=TMNM+TNSD*WX(2)
      IF(TMN>TMX) TMN=TMX-.2*ABS(TMX)
      
END SUBROUTINE Generate_AirTem
    
!***************************  Solar Radiation ***************************
!--------------- 4. generate solar radiation --------------------------
    
SUBROUTINE Generate_SolarRad    
      IMPLICIT NONE
 
      ! local variables
      REAL:: RX 
      
      RX=maxRad-RM
      SRAD=RM+WX(3)*RX/4.
      IF(SRAD<=0.)SRAD=.05*maxRad
END SUBROUTINE Generate_SolarRad
    
! ********************************** Wind *****************************
! generate wind
SUBROUTINE Generate_WindSpeed 
      IMPLICIT NONE
      ! local variables
      REAL:: V6 
    
      V6=Generate_Random(randID(5))
      U10=UAVM(MO)*(-LOG(V6))**windPowerPara 
     
END SUBROUTINE Generate_WindSpeed

! generate Wind direction 
    
REAL FUNCTION Wind_Dir() 
      IMPLICIT NONE 
 
      ! local variables
      REAL::FX, G
      INTEGER:: J, J1, XJ1
   
      FX=Generate_Random(randID(6))
      DO J=1,16
            J1=J-1
            IF(windDirection(MO,J)>FX)GO TO 3
      END DO
      J=16
    3 IF(J==1)THEN
            G=FX/windDirection(MO,J)
      ELSE
            G=(FX-windDirection(MO,J1))/(windDirection(MO,J)-windDirection(MO,J1))
      END IF
      XJ1=J1
      Wind_Dir=PI2*(G+XJ1-.5)/16.
      IF(Wind_Dir<0.)Wind_Dir=PI2+Wind_Dir
 
END FUNCTION Wind_Dir
    
! ************************************* RH ******************************
! generate RH
    
SUBROUTINE Generate_RH 
      ! Ref: APEX-doc Eq 23 and 24
      IMPLICIT NONE
 
      ! local variables
      REAL:: Q1, BLM  
      ! RHP: the peak of the triangular distribution (BLM or RHD)
      Q1=RHM-1.
      maxRH=RHM-Q1*EXP(Q1)       ! Eq. 23: maxRH, largest RH can be generated on the day  
      BLM=RHM*(1.-EXP(-RHM))      ! BLM, lowest RH can be generated on the day
      RHD=Triangle_Num(BLM,RHM,maxRH,7)
      
END SUBROUTINE Generate_RH
    
END MODULE Weather_Generator