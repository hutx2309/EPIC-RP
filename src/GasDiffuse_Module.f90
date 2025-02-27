MODULE GasDiffuse_Module
USE PARM
USE Initialize

! Gas diffusion in soil
IMPLICIT NONE
CONTAINS 
! ----------------------- Thomas algorithm -------------------
SUBROUTINE Thomas(B,D,A,C,N)
!     SUBPROGRAM TRIDIAG USES THE THOMAS ALHORITHM
!     D()...SOLUTION RETURNED AS C()
!     B()...BELOW DIAGONAL ELEMENTS
!     D()...DIAGONAL ELEMENTS
!     A()...ABOVE DIAGONAL ELEMENTS
!     C()...RIGHT HAND SIDE

      ! local variables
      INTEGER:: N, I, J
      REAL:: R
      REAL:: A(N),B(N),C(N),D(N)

      !     FORWARD ELIMINATION
      DO I=2,N
            R=B(I)/D(I-1)
            D(I)=D(I)-R*A(I-1)
            C(I)=C(I)-R*C(I-1)
      END DO
      !     BACK SUBSTITUTION
      C(N)=C(N)/D(N)
      DO I=2,N
            J=N-I+1
            C(J)=(C(J)-A(J)*C(J+1))/D(J)
      END DO           
END SUBROUTINE Thomas

! ----------------------- Gas transport  ----------------------
SUBROUTINE Solve_GasTrans(CONC,DPRM,CUP,DSURF,FLUXOUT)
!     SUBPROGRAM SOLVES THE GAS TRANSPORT EQUATION
IMPLICIT NONE
 
      ! local variables
      INTEGER:: ID
      ! arrays for storing tridiagonal matrix and Humus_To_Slow
      REAL:: A(MSC),B(MSC),C(MSC),D(MSC),CONC(MSC),CONCTMP(MSC),DPRM(MSC), & 
             ALX = 0.5, DTMIN = 0.01, DTSUB, T, FLUXOUT, R, R1, R2, DSURF, CUP, F1, F2 
      LOGICAL::CONVG        ! convergence flag
      
      ! ALX is the explicit fraction, 1-ALX the implicit 

      C=0.
      DTSUB=timeIntervalGasDiff
      CONVG = .FALSE.   
      DO WHILE(.NOT.CONVG)
            !! set P_Rate_bySoil the integration
            CONCTMP=CONC        ! set concentration to initial values
            T=0.                ! set time to top of the hour
            FLUXOUT=0.          ! zero the hourly flux
            CONVG=.TRUE.        ! will be unset when and if we observe a failure
            DO WHILE(T<timeIntervalGasDiff)
                  !! These must be initialized here because we may take a
                  !! short step at the end of the integration interval
                  R=DTSUB/(layerThickness**2)
                  R1=ALX*R
                  R2=(1.-ALX)*R
                  ! UPPER BOUNDARY CONDITION
                  A(1)=-R2*DPRM(1)
                  B(1)=0.0
                  C(1)=AFP(1)*CONCTMP(1)+R1*(DPRM(1)*(CONCTMP(2)-CONCTMP(1))-&
                  DSURF*(CONCTMP(1)-CUP))+R2*DSURF*CUP
                  D(1)=AFP(1)+R2*(DSURF+DPRM(1))
                  ! Calculate surface flux using the old concentration values
                  F1=-DSURF*(CONCTMP(1)-CUP)/layerThickness
                  !      print *,'start conctmp(1), cup, F1,alx : ', conctmp(1), cup, F1, alx !XXX
                  ! MAIN COMPUTATIONS
                  DO ID=2,IUN
                        A(ID)=-R2*DPRM(ID)
                        B(ID)=-R2*DPRM(ID-1)
                        C(ID)=AFP(ID)*CONCTMP(ID)+R1*(DPRM(ID)*(CONCTMP(ID+1)-&
                        CONCTMP(ID))-DPRM(ID-1)*(CONCTMP(ID)-CONCTMP(ID-1)))
                        D(ID)=AFP(ID)+R2*(DPRM(ID-1)+DPRM(ID))
                  END DO
                  ! LOWER BOUNDARY CONDITION
                  A(layersEqualThick)=0.0
                  B(layersEqualThick)=-R2*DPRM(layersEqualThick-1)
                  C(layersEqualThick)=AFP(layersEqualThick)*CONCTMP(layersEqualThick)-R1*DPRM(layersEqualThick-1)*&
                                      (CONCTMP(layersEqualThick)-CONCTMP(layersEqualThick-1))
                  D(layersEqualThick)=AFP(layersEqualThick)+R2*DPRM(layersEqualThick-1)
                  ! SOLVE TRIADIAGONAL SYSTEM
                  CALL Thomas(B,D,A,C,layersEqualThick)
                  !DO ID=1,layersEqualThick
                  !IF(ABS(C(ID)-CONCTMP(ID))/ABS(CONCTMP(ID)+1.E-4)>.5.AND.DTSUB>DTMIN)THEN
                  IF(ANY(ABS(C-CONCTMP)/(ABS(CONCTMP)+1.0E-4)>0.5).AND.DTSUB>DTMIN)THEN
                        !! excessive change.  reduce time step and try again
                        !! (but if the timestep is already too small, then we
                        !! just live with it)
                        DTSUB=DTSUB*.5
                        CONVG=.FALSE.
                        EXIT
                  ELSE 
                        CONCTMP=MAX(C,1.E-10)
                  END IF
                  !END DO
                  !IF(CONVG==.FALSE.)EXIT
                  ! Calculate surface flux using the new concentration values
                  F2=-DSURF*(C(1)-CUP)/layerThickness ! use C instead of conctmp, so we don't see the effect of MAX above
                  ! Combine the fluxes in the same ratio used in the Crank-Nicholson scheme
                  FLUXOUT=FLUXOUT+(ALX*F1+(1.-ALX)*F2)*DTSUB
                  ! update the time
                  T=T+DTSUB
                  IF(timeIntervalGasDiff-T<DTSUB)THEN
                        DTSUB=timeIntervalGasDiff-T    ! ensure that we don't overshoot our intended finish time
                  END IF 
            END DO                 ! while(t < timeIntervalGasDiff
      END DO                     ! while(.not. convg)
      CONC=CONCTMP
END SUBROUTINE Solve_GasTrans

! -------------------------- Hourly Gas CONCENTRATIONS ----------------------
SUBROUTINE Gas_Cyc 
! THIS LOOP CALCULATES THE PRODUCTION AND CONSUMPTION OF O2, CO2, AND N2O BY SOIL LAYER WITHIN ONE HOUR. 
! IT ALSO UPDATES POOLS OF NO3 AND NO2 FACTOR TO CONVERT MASS (kg/ha FOR A GIVEN SOIL LAYER) TO 
! GAS CONCENTRATION (G/M3 SOIL)
 
      IMPLICIT NONE
              					 
      ! local variables
      INTEGER:: J
      REAL:: AD1, AD2, AD3, XFCM, XMP, XSTN, EA, EAO2, EAO2R, ESRR, WN2G, ESMR,X1,X2, X3, X4, VWCE, &
        DW, SoilPartl_dia2, XFC, BX, A2, CONSTA, CONSTB, PVC, VEDW, HX, SA, &
        DAO2TC, XKT, QTB, QTC, O2M, RRTC, RRC, XKTR, QTBR, QTCR, O2R, SUM, ESD,&
        CNO3, CNO2, WN5, WN3, WN1, EAN5, EAN3, EAN1, EAD, GENN2O, O2CONS, CO2GEN, DF
        		
 
      REAL:: B1, B2, B3, COX, DAO2, O2MW, XN, FD, RGF, RMF, VU, VL, XL, SoilPartl_dia1, STN, &
        A1, T1, FCMP, VWCR, WPMP, A, B, C 
 
      ! XKN5,XKN3,XKN1/2.8,0.5,0.14/ *affinity coeff in gasdn0607                                                                        
      ! XKN5,XKN3,XKN1/10.0,0.5,0.0005/ *affinity coeff in epic0810v5                                                                        
      ! XKN5,XKN3,XKN1/10.0,2.5,2.5/ *as per Grant et al. (1993)                                                                   
      ! XKN5,XKN3,XKN1/10.0,2.5,2.5/ *latest approximation (2011/4/5)                                                        
      ! XKN5,XKN3,XKN1/0.077,0.126,0.057/ *latest approximation (2011/4/26)
      ! XKN5,XKN3,XKN1/10.,.5,.0005/,                                                                                                                               
      DATA B1,B2,B3/70.,140.,720./,COX/.032/,DAO2/7.2E-6/,O2MW/32./,&
      XN/2.58368E15/,FD/.19/,RGF,RMF/.547,.01785/,VU/1./,VL/-.15/,&
      SoilPartl_dia1/3.27E-7/,STN/.0728/,A1/.61821/,T1/3.27536/,FCMP/3.06/,&
      VWCR/.03/,WPMP/153.06/,A,B,C/.002,1.4,.5/
      ! THIS LOOP CALCULATES THE PRODUCTION AND CONSUMPTION OF O2, CO2, AND
      ! N2O BY SOIL LAYER WITHIN ONE HOUR. IT ALSO UPDATES POOLS OF NO3 AND
      ! NO2 FACTOR TO CONVERT MASS (kg/ha FOR A GIVEN SOIL LAYER) TO 
      ! GAS CONCENTRATION (G/M3 SOIL)
      AD1=0.
      AD2=0.
      AD3=0.
      XFCM=LOG10(FCMP)
      XMP=XFCM-LOG10(WPMP)
      XSTN=2.*STN
      DO J=1,layersEqualThick
            AD1=AD1+NO3_N_Soil(J)+NO2Weight(J)+WN2O(J)
            ! EA=TOTAL ELECTRONS ACCEPTED BY O2 AND N OXIDES
            EA=0.
            EAO2=0.
            EAO2R=0.
            ESRR=0.
            WN2G=0.
            ESMR=soilResp(J)/B3
            X4=SOT(J)+273.15
            IF(deNitri_Method==3)THEN
                  VWCE=MIN(.99,(VWC(J)-VWCR)/(TPOR(J)-VWCR))
                  VWCE=MAX(.01,VWCE)
                  X3=.001*((VWCE**(-1./C)-1.)**(1./B))/A
                  DW=MAX(1.01E-6,1.E-6+8.E-6*X3**(-0.945703126))
            ELSE
                  SoilPartl_dia2=1.E-6*(T1-(1./A1)*LOG((VU-VL)/(CBVT-VL)-1.))
                  XFC=LOG10(VFC(J))
                  BX=XMP/(LOG10(VWP(J))-XFC)
                  A2=10.**(BX*XFC+XFCM)
                  CONSTA=9.82*A2/XSTN
                  CONSTB=1./BX
                  PVC=(CONSTA*SoilPartl_dia2)**CONSTB
                  VEDW=VWC(J)-PVC
                  HX=.9549*PVC/((.5*SoilPartl_dia1)**2+(.5*SoilPartl_dia2)**2+SoilPartl_dia1*SoilPartl_dia2/4.)
                  XL=SQRT(HX**2+(.5*(SoilPartl_dia2-SoilPartl_dia1))**2)
                  SA=1.5708*(SoilPartl_dia1+SoilPartl_dia2)*XL
                  DW=MAX(1.1E-6,VEDW/SA)
            END IF
            DAO2TC=DAO2*(X4/293.15)**6
            XKT=1.5708E-10*XN*CBiomass(J)*DAO2TC*DW/(DW-1.E-6)
            IF(XKT>0.)THEN
                  QTB=XKT*(CLO2(J)-COX)-ESMR
                  QTC=XKT*COX*CLO2(J)
                  IF(QTB>1.E10)THEN
                        EAO2=ESMR
                  ELSE
                        O2M=(QTB+SQRT(QTB*QTB+4.*XKT*QTC))/(2.*XKT)
                        EAO2=ESMR*O2M/(O2M+COX)
                  END IF
            ELSE
                  EAO2=0.
            END IF
            IF(RWTZ(J)>1.E-5.AND. rootRespirFlag==0)THEN
                  ! NEW DERIVATION FOR ELECTRON SUPPLY DUE TO ROOT RESPIRATION
                  ! SEE MODEL DOCUMENTATION
                  ! ROOT RESPIRED C (KG C HA-1 D-1)
                  X1=(RGF/(1.-RGF))*MAX(0.,DRWX(J))+RMF*RWTZ(J)
                  RRTC=.42*X1
                  ! ESRR=MOLE E- M-2 H-1 FROM ROOT RESPIRATION - USE ESRR
                  ESRR=5.833E-4*X1
                  RRC=ESRR*B3 
                  !XKTR=54.573*DAO2TC*RWTZ(J)
                  X2=LOG(DW/.001)
                  IF(X2<=0.)X2=1.
                  XKTR=125.664*DAO2TC*RWTZ(J)/X2
                  QTBR=XKTR*(CLO2(J)-COX)-ESRR
                  IF(QTBR>1.E10)THEN
                        EAO2R=ESRR
                  ELSE
                        QTCR=XKTR*COX*CLO2(J)            
                        ! SOLVE QUADRATIC EQN FOR O2M AND O2MR
                        O2R=(QTBR+SQRT(QTBR*QTBR+4.*XKTR*QTCR))/(2.*XKTR)
                        ! ELECTRONS FROM MICROBE AND ROOT RESPIRATION ACCEPTED BY O2
                        EAO2R=ESRR*O2R/(O2R+COX)
                  END IF
            END IF 
            SUM=EAO2+EAO2R
            ESD=ESMR+ESRR-SUM
            ! ELECTRONS AVAILABLE FOR DENITRIFICATION
            ESD=FD*ESD
            ! COMPETITION FOR ELECTRONS AMONG OXIDES OF N
            ! CALCULATE WEIGHING FACTORS FIRST
            X1=DZ10*VWC(J)
            CNO3=MAX(1.E-5,NO3_N_Soil(J)/X1)
            CNO2=MAX(1.E-5,NO2Weight(J)/X1)
            WN5=5.*CNO3/(XKN5+CNO3)
            WN3=3.*CNO2/(XKN3*(1.+CNO3/XKN5)+CNO2)
            WN1=CLN2O(J)/(XKN1*(1.+CNO2/XKN3)+CLN2O(J))
            ! CALCULATE THE RATES OF REDUCTION OF OXIDES OF N
            X2=ESD/(WN1+WN3+WN5)
            X1=MAX(1.E-10,NO3_N_Soil(J)/B1)
            EAN5=MIN(X1,X2*WN5)
            X1=MAX(1.E-10,NO2Weight(J)/B1)
            EAN3=MIN(X1,X2*WN3)
            IF(WN2O(J)>0.)THEN
                  X1=WN2O(J)/B2
                  EAN1=MIN(X1,X2*WN1)
            ELSE
                  EAN1=0.
            END IF
            ! THESE ARE THE RESULTS BY LAYER AT THE END OF ONE HOUR
            ! IF NOT ALL ELECTRONS CAN BE ACCEPTED BY O2 (ESD>0.)
            ! TOTAL ELECTRONS ACCEPTED AND TRANSFORMATIONS OF N OXIDES
            EA=EA+EAO2+EAO2R+EAN5+EAN3+EAN1
            SMEA(J)=SMEA(J)+EA
            SMES(J)=SMES(J)+ESMR+ESRR
            ! EAD=TOTAL DEFICIT OF ELECTRONS
            EAD=EAN5+EAN3+EAN1
            ! LIQUID POOLS
            NO3_N_Soil(J)=NO3_N_Soil(J)-EAN5*B1
            NO2Weight(J)=NO2Weight(J)-(EAN3-EAN5)*B1
            ! GAS POOLS
            ! NITROUS OXIDE AND DINITROGEN
            ! GENN2O CALCULATES HOW MUCH N2O IS GENERATED (kg/ha/h per layer)
            GENN2O=EAN3*B1-EAN1*B2
            ! DN2OG(J) ACCUMULATES N2O GENERATED DURING A DAY (kg/ha)
            DN2OG(J)=DN2OG(J)+GENN2O
            ! WN2O(J) UPDATES THE N2O POOL (kg/ha). *WN2O(J)=denitriedN2O(J)*
            WN2O(J)=WN2O(J)+GENN2O
            ! AN2OC(J) CONVERTS MASS OF N2O INTO CONCENTRATION (G/M3)
            AN2OC(J)=WN2O(J)/DZ10
            ! WN2G CALCULATES THE MASS OF N2 GENERATED (kg/ha)
            WN2G=EAN1*B2
            !IF(WN2O(J)<0.)WRITE(KW(1),'(1X,A,4I4,10E16.6)')'~~~~~',IY,MO,&
            !DayOfMon,J,WN2O(J),GENN2O,WN2G
            ! DN2G(J) ACCUMULATES N2 GENERATED DURING A DAY (kg/ha)
            DN2G(J)=DN2G(J)+WN2G
            ! O2CONS= O2 CONSUMED (kg/ha)
            ! (MOL E/M2*1/4*MOL O2/MOL E*32 G O2/MOL O2*10.)
            O2CONS=2.5*SUM*O2MW
            ! DO2CONS ACCUMULATES O2 CONSUMED DAILY 
            ! CARBON DIOXIDE
            DO2CONS(J)=DO2CONS(J)+O2CONS
            ! AO2C(J) RECALCULATES O2 CONCENTRATION IN LAYER (G/M3)
            AO2C(J)=MAX(0.,AO2C(J)-O2CONS/DZ10)
            ! CO2GEN IS CO2 GENERATED (kg/ha)
            !(MOL E/M2*1/4*MOL C/MOL E*12 G C/MOL C*10.)
            ! FACTOR 10. ABOVE IS TO CONVERT G/M2 TO kg/ha
            CO2GEN=EA*30.
            !	  DCO2GEN(J) ACCUMULATES CO2 GENERATED DAILY 
            DCO2GEN(J)=DCO2GEN(J)+CO2GEN
            ! ACO2C(J) RECALCULATES CO2 CONCENTRATION IN LAYER (G/M3)
            ACO2C(J)=MAX(0.,ACO2C(J)+CO2GEN/DZ10)
            AD2=AD2+NO3_N_Soil(J)+NO2Weight(J)+WN2O(J)
            AD3=AD3+WN2G
      END DO
      DF=AD2+AD3-AD1
      IF(ABS(DF/AD1)>1.E-5)WRITE(KW(1),1)IY,MO,DayOfMon,AD1,AD2,AD3,DF
    1 FORMAT(1X,'NDNITCI',3I4,4E16.6)
 
END SUBROUTINE Gas_Cyc

! ------------------------- Ice fraction in soil ------------------------------

SUBROUTINE Cal_Ice(IFW,SOILTEMP,VOLWC,TOTPOR)
      ! CORRECT VOLWC FOR FREEZING WHEN TEMPERATURE IS BELOW 0C.  ICE
      ! FRACTION IS ASSUMED TO BE A LINEAR RAMP BETWEEN 0C AND -8C.  WE
      ! SET A FLOOR OF 1.0E-8 ON THE RESIDUAL WATER FRACTION IN ORDER TO
      ! PREVENT DIVIDE-BY-ZERO PROBLEMS AND TO ENSURE THAT THIS FUNCTION
      ! IS INVERTIBLE.

      IMPLICIT NONE
 
      ! local variables potential problems with local variables (SAVE attribute)
      INTEGER:: IFW, J
      REAL:: FAC, SOILTEMP(layersEqualThick),VOLWC(layersEqualThick),TOTPOR(layersEqualThick)

      DO J=1,layersEqualThick
            IF(SOILTEMP(J)<0.)THEN
                  ! SET A FLOOR ON FAC TO AVOID POSSIBLE DIVIDE-BY-ZERO AND TO ENSURE RELATION IS INVERTIBLE
                  FAC=MAX(1.+.125*SOILTEMP(J),1.0E-8)
            ELSE
                  FAC=1.
            END IF
            IF(IFW==0)THEN
                  TOTPOR(J)=TOTPOR(J)-(1.-FAC)*VOLWC(J) ! REDUCE TOTAL POROSITY BY THE VOLUME OF THE ICE COMPONENT
                  VOLWC(J)=VOLWC(J)*FAC
            ELSE
                  VOLWC(J)=VOLWC(J)/FAC
                  TOTPOR(J)=TOTPOR(J)+(1.0-FAC)*VOLWC(J) ! ADD THE ICE BACK TO THE TOTAL
            END IF
      END DO
 
END SUBROUTINE Cal_Ice

! -------------------------- interpolate concentrations of soil ----------------
SUBROUTINE Interp_Concentration(X,Y,N1,N2)
      ! EPIC1102   
      ! THIS SUBPROGRAM INTERPOLATES CONCENTRATIONS FROM LAYERS WITH N_Rate_bySoil
      ! EQUAL THICKNESS TO LAYERS OF EQUAL THICKNESS USED IN DIFFERENTIAL
      ! EQUATIONS OF GAS DIFFUSION.
      IMPLICIT NONE
 
      ! local variables    potential problems with variables (may with SAVE attribute)
      INTEGER:: N1, N2, K, J, L, I
      REAL:: X(15), Y, ZZ, Z1, TOT
      DIMENSION:: Y(N2)

      ZZ=0.
      Z1=0.
      TOT=0.
      J=1
      DO K=1,N1
            L=Layer_ID(K)
            DO 
                  IF(ZC(J)>Z(L))EXIT
                  Y(J)=TOT+X(L)*(ZC(J)-ZZ)
                  ZZ=ZC(J)
                  J=J+1
                  IF(J>N2)RETURN
                  TOT=0.
            END DO
            IF(J>N2)RETURN
            TOT=TOT+X(L)*(Z(L)-ZZ)
            Z1=Z(L) 
            ZZ=Z(L)
      END DO
      DO I=1,N2-1
            Y(I)=Y(I)/layerThickness
      END DO 
      Y(N2)=MAX(X(Layer_ID(N1)),TOT/layerThickness)
END SUBROUTINE Interp_Concentration

! -------------------- interpolate soil properties from equal to unequal layers -----
SUBROUTINE Interp_Soil3(X,Y,N1,N2)
      ! EPIC1102
      ! THIS SUBPROGRAM INTERPOLATES SOIL PROPERTIES FROM LAYERS WITH   
      ! EQUAL THICKNESS (OUTPUT FROM DIF EQ SOLN OF GAS DIFFUSION EQS) TO
      ! LAYERS OF UNEQUAL THICKNESS (INPUT SOIL LAYERS).
      IMPLICIT NONE
 
      ! local variables  potential problems with local variables (SAVE attribute)
      INTEGER:: J, K, L, N1, N2
      REAL:: X(30), Y(30), Z1, TOT 
        
      Z1=0.
      TOT=0.
      Y=0.
      J=1
      DO K=1,N2
            DO WHILE(J<=N1)
                  L=Layer_ID(J)
                  IF(Z(L)>ZC(K))EXIT
                  Y(L)=TOT+X(K)*(Z(L)-Z1)
                  Z1=Z(L)
                  J=J+1
                  TOT=0.
            END DO
            IF(J<=N1)THEN
                  TOT=TOT+X(K)*(ZC(K)-Z1)
                  Z1=ZC(K)
            ELSE
                  EXIT  
            END IF
      END DO
      Z1=0.
      DO J=1,N1
            L=Layer_ID(J)     
            Y(L)=Y(L)/(Z(L)-Z1)
            Z1=Z(L)
      END DO
END SUBROUTINE Interp_Soil3

! ------- interpolate soil properties from equal to unequal layers --- ----- 
SUBROUTINE Interp_Soil2(X,Y,N1,N2)
      ! EPIC1102
      ! THIS SUBPROGRAM INTERPOLATES SOIL PROPERTIES FROM LAYERS WITH   
      ! EQUAL THICKNESS (OUTPUT FROM DIF EQ SOLN OF GAS DIFFUSION EQS) TO
      ! LAYERS OF UNEQUAL THICKNESS (INPUT SOIL LAYERS).
      IMPLICIT NONE
 
      ! local variables   potential problems with local variables (SAVE attribute)
      INTEGER:: J, K, L, N1, N2   
      REAL:: X(30), Y(30), Z1, TOT, AD1, AD2

      Z1=0.
      TOT=0.
      AD1=0.
      AD2=0.
      Y=0.
      J=1
      DO K=1,N2
            AD1=AD1+X(K)
            DO WHILE(J<=N1)
                  L=Layer_ID(J)
                  IF(Z(L)>ZC(K))EXIT
                  Y(L)=TOT+X(K)*(Z(L)-Z1)/layerThickness
                  AD2=AD2+Y(L)
                  Z1=Z(L)
                  J=J+1
                  TOT=0.
            END DO
            IF(J<=N1)THEN
                  TOT=TOT+X(K)*(ZC(K)-Z1)/layerThickness
                  Z1=ZC(K)
            ELSE
                  EXIT  
            END IF
      END DO
      L=Layer_ID(N1)     
      Y(L)=Y(L)+AD1-AD2
END SUBROUTINE Interp_Soil2

! ---------------------------- interpolate soil properties -----------------

SUBROUTINE Interp_Soil(X, Y, N1, N2)
      ! EPIC1102
      ! THIS SUBPROGRAM INTERPOLATES SOIL PROPERTIES FROM LAYERS WITH  
      ! EQUAL THICKNESS TO LAYERS OF EQUAL THICKNESS USED IN DIFFERENTIAL
      ! EQUATIONS OF GAS DIFFUSION.
      IMPLICIT NONE
 
      ! local variables   potential problems with local variables (SAVE attribute)
      INTEGER:: N1, N2, K, J, L, I
      REAL:: X(15), Y(N2), ZZ, Z1, TOT, AD1, AD2
 
      ZZ=0.
      Z1=0.
      TOT=0.
      Y=0.
      J=1
      DO K=1,N1
            L=Layer_ID(K)
            DO 
                  IF(ZC(J)>Z(L))EXIT
                  Y(J)=TOT+X(L)*(ZC(J)-ZZ)/(Z(L)-Z1)
                  ZZ=ZC(J)
                  J=J+1
                  IF(J>N2)RETURN
                  TOT=0.
            END DO
            TOT=TOT+X(L)*(Z(L)-ZZ)/(Z(L)-Z1)
            IF(J>N2)EXIT
            Z1=Z(L) 
            ZZ=Z(L)
      END DO
      I=MIN(J,N2)
      Y(I)=MAX(Y(I),TOT)
      AD1=0.
      DO I=1,N1
            AD1=AD1+X(I)
      END DO
      AD2=0.
      DO I=1,N2
            AD2=AD2+Y(I)
      END DO
      IF(ABS((AD1-AD2)/(AD1+1.E-20))>.01)THEN
            WRITE(KW(1),'(A,3I4,A,2E13.5)')'!!!!!',IYR,MO,DayOfMon,'Interp_Soil ERR',AD1,AD2
      END IF    
END SUBROUTINE Interp_Soil

! ------------------------------ Solve Gas Diffusion Equation ---------------
SUBROUTINE Solve_GasDiff 
      ! why this subroutine can not be referenced? potential problem with local variables 
      !     EPIC1102
      !     THIS PROGRAM SOLVES THE GAS DIFFUSION EQUATION
      IMPLICIT NONE
      
      ! local variables should have SAVE attribute 
      REAL, DIMENSION(30):: YTP, XTP1,XTP2, XTP3, XTP4, XTP5, XTP6, XTP7, XTP8, XTP9, XTP10, &
                     XTP11, XTP12, XTP13, XTP14, XTP15, XTP16, XTP17, XTP18, XTP19
      !     WBM & RCI REVISED GAS CONC. IN AIR (g/m3) (8/25/2011)      
      !     PARTIAL PRESSURE AT SEA LEVEL (0 m)
      !     O2=299.2957143;
      !     CO2=.766071429; CO2-C=.208928571
      !     N2O=.000616786; N2O-N=.0003925
      !     CLOC,CUPN,CLON/279.,240.,.18,10.,.00018,1./
      REAL:: GASC = 0.08205783, DCAO = .06132, DCAC = .05556, DCAN = .05148, CUPO = 299., &
            CUPC = .2089, CUPN = .0003925 
      REAL::FLXO2,FLXCO2,FLXN2O, DZX, t_inter_gasdiff10,SM1,SM2, &
            AD1, AD2, AD3, AD4, AD5, AD6, AD7, AD8, AD9, AD10, AD11, AD12, &
            TOT1, TOT2, TOT3, TOT4, ABST, HKF, X1, X2, X3, X4, XXO, XXC, XXN, XTPOR, XVWC, XVFC,&
            XAFP, EPS, DF, DCSO, DCSC, DCSN, XHKPO, XHKPC, XHKPN, DFO2B, DFCO2B, DFN2OB, &
            TIME, WO2GB, WCO2GB, WN2OGB, WO2GA, WCO2GA, WN2OGA, ZZ
 	  
      INTEGER:: J , L, I, IT, K, numTimeIntervals
 
      DZX=1000.*layerThickness
      t_inter_gasdiff10=10.*timeIntervalGasDiff
      SM1=0.
      AD1=0.
      AD2=0.
      AD3=0.
      TOT1=0.
      XTP5=0.
      DO J=1,Actual_SoilLayers
            L=Layer_ID(J)
            YTP(L)=soilResp(L)
            XTP1(L)=NO3_N_Soil(L)
            XTP2(L)=NO2Weight(L)
            XTP3(L)=WN2O(L)
            XTP4(L)=CBiomass(L)
            XTP5(L)=RWT(L,Crop_Num)
            XTP6(L)=MAX(1.E-10,DRWX(L))
            XTP7(L)=Weight_CO2L(L)
            XTP8(L)=WCO2G(L)
            XTP9(L)=Weight_N2OL(L)
            XTP10(L)=WN2OG(L)
            XTP11(L)=Weight_O2L(L)
            XTP12(L)=WO2G(L)
            XTP13(L)=SubSurfFlow_CO2L(L)
            XTP14(L)=SubSurfFlow_N2O(L)
            XTP15(L)=SubSurfFlow_O2L(L)
            XTP16(L)=VerticalFlow_CO2(L)
            XTP17(L)=verticalFlowN2O(L)
            XTP18(L)=VerticalFlow_O2(L)
            XTP19(L)=NH3_Weight(L)
            SM1=SM1+NO3_N_Soil(L)+NO2Weight(L)+WN2O(L)
            AD1=AD1+WN2O(L)
            AD2=AD2+WN2OG(L)
            AD3=AD3+Weight_N2OL(L)
            TOT1=TOT1+CBiomass(L)
      END DO
      CALL Interp_Soil(Porosity,TPOR,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil(soilWater,VWC,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil(Wilt_Point,VWP,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Concentration(soilTem,SOT,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil(YTP,soilResp,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil(XTP1,NO3_N_Soil,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil(XTP2,NO2Weight,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil(XTP3,WN2O,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil(XTP4,CBiomass,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil(XTP5,RWTZ,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil(XTP6,DRWX,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil(XTP7,Weight_CO2L,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil(XTP8,WCO2G,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil(XTP9,Weight_N2OL,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil(XTP10,WN2OG,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil(XTP11,Weight_O2L,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil(XTP12,WO2G,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil(XTP13,SubSurfFlow_CO2L,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil(XTP14,SubSurfFlow_N2O,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil(XTP15,SubSurfFlow_O2L,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil(XTP16,VerticalFlow_CO2,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil(XTP17,verticalFlowN2O,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil(XTP18,VerticalFlow_O2,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil(XTP19,NH3_Weight,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil(fieldCapacity,VFC,Actual_SoilLayers,layersEqualThick)
      CALL Cal_Ice(0,SOT,VWC,TPOR)
      AD4=0.
      AD5=0.
      AD6=0.
      TOT2=0.
      DO I=1,layersEqualThick
            AD4=AD4+WN2O(I)
            AD5=AD5+WN2OG(I)
            AD6=AD6+Weight_N2OL(I)
            TOT2=TOT2+CBiomass(I)
      END DO
      ZZ=DZX
      !Z1=1000.*(Z(Layer_ID(Actual_SoilLayers))-ZC(layersEqualThick-1))
      !IF(Z1<1.E-5)layersEqualThick=layersEqualThick-1
      DO I=1,layersEqualThick
            IF(I==layersEqualThick)ZZ=1000.*(Z(Layer_ID(Actual_SoilLayers))-ZC(layersEqualThick-1))
            TPOR(I)=TPOR(I)/ZZ
            VWC(I)=VWC(I)/ZZ
            VFC(I)=VFC(I)/ZZ
            VWP(I)=VWP(I)/ZZ
            ABST=SOT(I)+273.15
            AFP(I)=MAX(1.E-5,TPOR(I)-VWC(I))
            HKF=.018/(GASC*ABST)
            X1=.01*ABST
            XXO=EXP(-66.7354+87.4755/X1+24.4526*LOG(X1))
            HKPO(I)=HKF*(1./XXO-1.)
            ! HERE WE OBTAIN IN HENRYS CONST AS A FUNCT OF SOIL T DIRECTLY
            ! UNITS ARE kpa M3 MOL-1
            XXC=EXP(-8.286*LOG(ABST)+46.742)*1.E-03
            HKPC(I)=HKF*(1./XXC-1.)
            XXN=EXP(-60.7467+88.828/X1+21.2531*LOG(X1))
            HKPN(I)=HKF*(1./XXN-1.)
      END DO
      DO I=1,layersEqualThick
            IF(I<layersEqualThick)THEN
                  ! Estimate soil properties at grid cell boundaries using
                  ! linear interpolation (NB: someday we might want to do
                  ! a spline interpolation here)
                  ABST=.5*(SOT(I)+SOT(I+1))+273.15
                  XTPOR=.5*(TPOR(I)+TPOR(I+1))
                  XVWC=.5*(VWC(I)+VWC( I+1))
                  XVFC=.5*(VFC(I)+VFC(I+1))
                  XAFP=.5*(AFP(I)+AFP(I+1))
            ELSE
                  ! for the last grid cell we don't have a value to
                  ! interpolate with, so just use the cell-center value.
                  ABST=SOT(I)+273.15
                  XTPOR=TPOR(I)
                  XVWC=VWC(I)
                  XVFC=VFC(I)
                  XAFP=AFP(I)
            END IF
            HKF=.018/(GASC*ABST)
            ! THIS IS THE MILLINGTON-QUIRCK COEFF
            EPS=XAFP**3.333/XTPOR**2
            ! DIFFUSION COEFFICIENT IN SOIL
            DCSO=DCAO*EPS
            X1=.01*ABST
            XXO=EXP(-66.7354+87.4755/X1+24.4526*LOG(X1))
            XHKPO=HKF*(1./XXO-1.)
            X3=ABST/273.15
            X4=X3**1.75
            DPRO(I)=DCSO*XAFP/(XTPOR+XVWC*(1./XHKPO-1.))*X3*X3
            DCSC=DCAC*EPS
            ! HERE WE OBTAIN IN HENRYS CONST AS A FUNCT OF SOIL T DIRECTLY
            ! UNITS ARE kpa M3 MOL-1
            XXC=EXP(-8.286*LOG(ABST)+46.742)*1.E-03
            XHKPC=HKF*(1./XXC-1.)
            DPRC(I)=DCSC*XAFP/(XTPOR+XVWC*(1./XHKPC-1.))*X4
            DCSN=DCAN*EPS
            XXN=EXP(-60.7467+88.828/X1+21.2531*LOG(X1))
            XHKPN=HKF*(1./XXN-1.)
            DPRN(I)=DCSN*XAFP/(XTPOR+XVWC*(1./XHKPN-1.))*X4
      END DO
      IF(IY==1.AND.Date==IBD) THEN
            !     THIS SUBPROGRAM DOES AN INITIAL PARTITION EACH GAS CONC. IN IrrWater INTO GAS CONCENTRATION IN AIR 
            !     AND LIQUID PHASE AND THEN CALCULATES THE TOTAL GAS CONCENTRATION (G/M3 SOIL)
            !     IT IS AN APPROXIMATION AND SHOULD BE CALCULATED ONLY ONCE AT THE BEGINNING OF EACH SIMULATION 
            !     - RCI - 7/19/09
            DO J=1,layersEqualThick
                  CLO2(J)=gasO2Con(J)/HKPO(J) 
                  CLCO2(J)=gasCO2Con(J)/HKPC(J)
                  CLN2O(J)=gasN2OCon(J)/HKPN(J) 
                  AO2C(J)=gasO2Con(J)*AFP(J)+CLO2(J)*VWC(J)
                  ACO2C(J)=gasCO2Con(J)*AFP(J)+CLCO2(J)*VWC(J)
                  AN2OC(J)=gasN2OCon(J)*AFP(J)+CLN2O(J)*VWC(J)
            END DO
	  
      END IF
      ! SET TO 0. LAYER ARRAYS OF DAILY PRODUCTION AND CONSUMPTION
      ! OF GASES (O2, CO2, N2O, N2 AND N2O+N2 [DDENIT])
       
      DO I=1,layersEqualThick
            DO2CONS(I)=0.
            DCO2GEN(I)=0.
            DN2OG(I)=0.
            DN2G(I)=0.
      END DO
      
      ! INITIALIZE DAILY GAS FLUXES (WBM & RCI, 8/25/11)      
      SMEA=0.
      SMES=0.
      DFO2S=0.
      DFCO2S=0.
      DFN2OS=0.
      DFO2B=0.
      DFCO2B=0.
      DFN2OB=0.
      DFO2T=0.
      DFCO2T=0.
      DFN2OT=0.      
      TIME=0.
      numTimeIntervals=INT(24./timeIntervalGasDiff+.9) 
      ! TIME LOOP TO CALCULATE GENERATION AND CONSUMPTION OF
      ! GASES, GAS TRANSPORT, AND FLUX AT THE SURFACE
      DO IT=1,numTimeIntervals
            TIME=TIME+timeIntervalGasDiff
            ! CALCULATE GENERATION AND CONSUMPTION OF GASES
            CALL Gas_Cyc 
            ! RE-CALCULATE GAS CONCENTRATIONS IN LIQUID AND AIR PHASES
            ! PRODUCTION AND CONSUMPTION OF GASES
            DO J=1,layersEqualThick
	            ! CLO2=CONC GAS IN LIQ PHASE (G/M3 WATER)
                  CLO2(J)=AO2C(J)/(AFP(J)*HKPO(J)+VWC(J))   
                  ! gasO2Con=CONC GAS IN GAS PHASE (G/M3 AIR)
                  gasO2Con(J)=AO2C(J)/(TPOR(J)+VWC(J)*(1./HKPO(J)-1.))
                  ! CLCO2=CONC GAS IN LIQ PHASE
                  CLCO2(J)=ACO2C(J)/(AFP(J)*HKPC(J)+VWC(J))   
                  ! gasCO2Con=CONC GAS IN GAS PHASE
                  gasCO2Con(J)=ACO2C(J)/(TPOR(J)+VWC(J)*(1./HKPC(J)-1.))
                  ! CLN2O=CONC GAS IN LIQ PHASE
                  CLN2O(J)=AN2OC(J)/(AFP(J)*HKPN(J)+VWC(J)) 
                  ! gasN2OCon=CONC GAS IN GAS PHASE
                  gasN2OCon(J)=AN2OC(J)/(TPOR(J)+VWC(J)*(1./HKPN(J)-1.))
            END DO
          
            WO2GB=0.
            WCO2GB=0.
            WN2OGB=0.
            WO2GA=0.
            WCO2GA=0.
            WN2OGA=0.
            DO J=1,layersEqualThick
                  X1=AFP(J)*DZ10
                  WO2GB=WO2GB+gasO2Con(J)*X1
                  WCO2GB=WCO2GB+gasCO2Con(J)*X1
                  WN2OGB=WN2OGB+gasN2OCon(J)*X1
                  !IF(J==layersEqualThick-1)THEN
                        !WO2GBT=WO2GB
                        !WCO2GBT=WCO2GB
                        !WN2OGBT=WN2OGB
                  !END IF
            END DO
            ! MOVE O2 AND STORE VALUES OF GAS CONC. IN 2-DIM ARRAY
            CALL Solve_GasTrans(gasO2Con,DPRO,CUPO,DCAO,FLXO2 ) 
            ! MOVE CO2 AND STORE VALUES OF GAS CONC. IN 2-DIM ARRAY
            CALL Solve_GasTrans(gasCO2Con,DPRC,CUPC,DCAC,FLXCO2 )
            ! MOVE N2O AND STORE VALUES OF GAS CONC. IN 2-DIM ARRAY
            CALL Solve_GasTrans(gasN2OCon,DPRN,CUPN,DCAN,FLXN2O )
            ! RECALCULATE GAS CONCENTRATIONS IN LIQUID AND GAS PHASES
            DO J=1,layersEqualThick
                  AO2C(J)=gasO2Con(J)*AFP(J)+CLO2(J)*VWC(J)
                  ACO2C(J)=gasCO2Con(J)*AFP(J)+CLCO2(J)*VWC(J)
                  AN2OC(J)=gasN2OCon(J)*AFP(J)+CLN2O(J)*VWC(J)
                  X2=AN2OC(J)*DZ10
                  WN2O(J)=X2
                  X1=AFP(J)*DZ10
                  WO2GA=WO2GA+gasO2Con(J)*X1
                  WCO2GA=WCO2GA+gasCO2Con(J)*X1
                  WN2OGA=WN2OGA+gasN2OCon(J)*X1
                  !IF(J==layersEqualThick-1)THEN
                  !WO2GAT=WO2GA
                  !WCO2GAT=WCO2GA
                  !WN2OGAT=WN2OGA
                  !END IF    
            END DO
            DFO2T=DFO2T+(WO2GA-WO2GB)*timeIntervalGasDiff
            DFCO2T=DFCO2T+(WCO2GA-WCO2GB)*timeIntervalGasDiff
            DFN2OT=DFN2OT+(WN2OGA-WN2OGB)*timeIntervalGasDiff
            DFO2S=DFO2S+FLXO2*t_inter_gasdiff10 
            DFCO2S=DFCO2S+FLXCO2*t_inter_gasdiff10
            DFN2OS=DFN2OS+FLXN2O*t_inter_gasdiff10
            DFO2B=DFO2T-DFO2S
            DFCO2B=DFCO2T-DFCO2S
            DFN2OB=DFN2OT-DFN2OS
            DO J=1,layersEqualThick
	            ! CLO2=CONC GAS IN LIQ PHASE (G/M3 WATER)
                  CLO2(J)=AO2C(J)/(AFP(J)*HKPO(J)+VWC(J))   
                  ! gasO2Con=CONC GAS IN GAS PHASE (G/M3 AIR)
                  gasO2Con(J)=AO2C(J)/(TPOR(J)+VWC(J)*(1./HKPO(J)-1.))
                  ! CLCO2=CONC GAS IN LIQ PHASE
                  CLCO2(J)=ACO2C(J)/(AFP(J)*HKPC(J)+VWC(J))   
                  ! gasCO2Con=CONC GAS IN GAS PHASE
                  gasCO2Con(J)=ACO2C(J)/(TPOR(J)+VWC(J)*(1./HKPC(J)-1.))
                  ! CLN2O=CONC GAS IN LIQ PHASE
                  CLN2O(J)=AN2OC(J)/(AFP(J)*HKPN(J)+VWC(J)) 
                  ! gasN2OCon=CONC GAS IN GAS PHASE
                  gasN2OCon(J)=AN2OC(J)/(TPOR(J)+VWC(J)*(1./HKPN(J)-1.))
            END DO
      END DO
      totN2O=0.
      SN2=0.
      totDenitriN=0.
      AD7=0.
      AD8=0.
      AD9=0.
      TOT3=0.
      DO J=1,layersEqualThick
            IF(SMES(J)>0.)THEN
                  EAR(J)=SMEA(J)/SMES(J)
            ELSE
                  EAR(J)=1.
            END IF
            X1=AFP(J)*DZ10
            WO2G(J)=gasO2Con(J)*X1
            WCO2G(J)=gasCO2Con(J)*X1
            WN2OG(J)=gasN2OCon(J)*X1
            X1=VWC(J)*DZ10
            Weight_O2L(J)=CLO2(J)*X1
            Weight_CO2L(J)=CLCO2(J)*X1
            Weight_N2OL(J)=CLN2O(J)*X1
            !IF(VWC(J)>VFC(J).AND.SOT(J)>0.)THEN
            SN2=SN2+DN2G(J)
            totN2O=totN2O+DN2OG(J)
            totDenitriN=totDenitriN+DN2OG(J)+DN2G(J)
            !END IF
            XTP1(J)=NO3_N_Soil(J)
            XTP2(J)=NO2Weight(J)
            XTP3(J)=WN2O(J)
            XTP4(J)=CBiomass(J)
            XTP5(J)=EAR(J)
            XTP6(J)=WN2OG(J)
            XTP7(J)=Weight_N2OL(J)
            XTP8(J)=DN2G(J)
            XTP9(J)=DN2OG(J)
            XTP10(J)=soilResp(J)
            XTP11(J)=WO2G(J)
            XTP12(J)=Weight_O2L(J)
            XTP13(J)=WCO2G(J)
            XTP15(J)=Weight_CO2L(J)
            XTP14(J)=SubSurfFlow_N2O(J)
            XTP17(J)=verticalFlowN2O(J)
            XTP19(J)=NH3_Weight(J)
            AD7=AD7+WN2O(J)
            AD8=AD8+WN2OG(J)
            AD9=AD9+Weight_N2OL(J)
            TOT3=TOT3+CBiomass(J)
      END DO
      IF(KFL(30)>0) WRITE(KW(30),19)IYR,MO,DayOfMon,VAR(4),(VWC(K),K=1,layersEqualThick),&
      (AFP(K),K=1,layersEqualThick),(NH3_Weight(K),K=1,layersEqualThick),(NO3_N_Soil(K),K=1,layersEqualThick),&
      (NO2Weight(K),K=1,layersEqualThick),(Weight_O2L(K),K=1,layersEqualThick),(WO2G(K),K=1,layersEqualThick),&
      (DO2CONS(K),K=1,layersEqualThick),DFO2S,DFO2B,DFO2T,QO2,(SubSurfFlow_O2L(K),K=1,layersEqualThick),&
      (VerticalFlow_O2(K),K=1,layersEqualThick),(Weight_CO2L(K),K=1,layersEqualThick),&
      (WCO2G(K),K=1,layersEqualThick),(DCO2GEN(K),K=1,layersEqualThick),DFCO2S,DFCO2B,DFCO2T,QCO2,&
      (SubSurfFlow_CO2L(K),K=1,layersEqualThick),(VerticalFlow_CO2(K),K=1,layersEqualThick),&
      (Weight_N2OL(K),K=1,layersEqualThick),(WN2OG(K),K=1,layersEqualThick),DFN2OS,DFN2OB,DFN2OT,QN2O,&
      (SubSurfFlow_N2O(K),K=1,layersEqualThick),(verticalFlowN2O(K),K=1,layersEqualThick),&
      (DN2OG(K),K=1,layersEqualThick),(DN2G(K),K=1,layersEqualThick),(SMEA(K),K=1,layersEqualThick),&
      (SMES(K),K=1,layersEqualThick),(EAR(K),K=1,layersEqualThick)
      CALL Cal_Ice(1,SOT, VWC,TPOR)
      CALL Interp_Soil2(XTP1,NO3_N_Soil,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil2(XTP2,NO2Weight,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil2(XTP3,WN2O,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil2(XTP4,CBiomass,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil3(XTP5,EAR,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil2(XTP6,WN2OG,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil2(XTP7,Weight_N2OL,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil2(XTP11,WO2G,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil2(XTP12,Weight_O2L,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil2(XTP13,WCO2G,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil2(XTP15,Weight_CO2L,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil2(XTP14,SubSurfFlow_N2O,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil2(XTP17,verticalFlowN2O,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil2(XTP19,NH3_Weight,Actual_SoilLayers,layersEqualThick)
      CALL Interp_Soil2(XTP10,soilResp,Actual_SoilLayers,layersEqualThick)
      SM2=0.
      AD10=0.
      AD11=0.
      AD12=0.
      TOT4=0.
      DO J=1,Actual_SoilLayers
            L=Layer_ID(J)
            gasO2Con(L)=100.*WO2G(L)/(Porosity(L)-soilWater(L))
            SM2=SM2+NO3_N_Soil(L)+NO2Weight(L)+WN2OG(L)+Weight_N2OL(L)
            AD10=AD10+WN2O(L)
            AD11=AD11+WN2OG(L)
            AD12=AD12+Weight_N2OL(L)
            TOT4=TOT4+CBiomass(L)
      END DO
      IF(ABS(AD1-AD4)>1.E-5)WRITE(KW(1),'(1X,3I4,A,E13.6,A,E13.6)')IY,MO,DayOfMon,'^^^^^ A1=',AD1,'  A4=',AD4
      IF(ABS(AD2-AD5)>1.E-5)WRITE(KW(1),'(1X,3I4,A,E13.6,A,E13.6)')IY,MO,DayOfMon,'^^^^^ A2=',AD2,'  A5=',AD5
      IF(ABS(AD3-AD6)>1.E-5)WRITE(KW(1),'(1X,3I4,A,E13.6,A,E13.6)')IY,MO,DayOfMon,'^^^^^ A3=',AD3,'  A6=',AD6
      IF(ABS(AD7-AD10)>1.E-5)WRITE(KW(1),'(1X,3I4,A,E13.6,A,E13.6)')IY,MO,DayOfMon,'^^^^^ A7=',AD7,'  A10=',AD10
      IF(ABS(AD8-AD11)>1.E-5)WRITE(KW(1),'(1X,3I4,A,E13.6,A,E13.6)')IY,MO,DayOfMon,'^^^^^ A8=',AD8,'  A11=',AD11
      IF(ABS(AD9-AD12)>1.E-5)WRITE(KW(1),'(1X,3I4,A,E13.6,A,E13.6)')IY,MO,DayOfMon,'^^^^^ A9=',AD9,'  A12=',AD12
      IF(ABS((TOT1-TOT2)/TOT1)>1.E-5)WRITE(KW(1),4)IY,MO,DayOfMon,'^^^^^ TOT1=',TOT1,' TOT2=',TOT2
      IF(ABS((TOT3-TOT4)/TOT3)>1.E-5)WRITE(KW(1),4)IY,MO,DayOfMon,'^^^^^ TOT3=',TOT3,' TOT4=',TOT4
      DF=SM1-SM2-SN2+DFN2OT
      IF(ABS(DF/SM1)>.1)THEN
            WRITE(KW(1),3)IY,MO1,DayOfMon,SM1,totDenitriN,SM2,totN2O,SN2,DFN2OT,DF
      END IF
 
 
    3 FORMAT(4X,'GAS XBAL ',3I4,2X,'BTOT=',F11.6,2X,'DNIT=',F11.6,2X,'FTOT=',&
      F11.6,2X,'totN2O=',F11.6,2X,'SN2=',F11.6,2X,'DFN2OT=',F11.6,2X,'DF=',F11.6)
    4 FORMAT(1X,3I4,A12,E13.6,A6,E13.6)          
   19 FORMAT(1X,3I4,F9.1,1000E13.5) 
 
END SUBROUTINE Solve_GasDiff

END MODULE GasDiffuse_Module