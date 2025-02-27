MODULE Nutrient_Module

USE PARM
USE Management
USE GasDiffuse_Module, ONLY: Interp_Concentration
IMPLICIT NONE

CONTAINS

SUBROUTINE Nutrient_Process(NCC, XX, XZP)
      IMPLICIT NONE
      ! Argument lists
      INTEGER, INTENT(IN):: NCC
      REAL, INTENT(IN):: XX  
      REAL, DIMENSION(13,16), INTENT(OUT):: XZP
      ! local variables declaration
      INTEGER:: I, J, L, L1, MXP
      REAL, DIMENSION(30)::ACO2, XTP1, XTP2
      REAL:: RTO, TOT, X1, Z1, ZZ, SUM, WPMX 

      ! ////////////////////////////////////////////////
      Z1=0.                                                                         
      Actual_Soil_Layer: DO J=1,Actual_SoilLayers  ! each soil layer, P_Rate_bySoil to 10 layers   
            L=Layer_ID(J)
            ZZ=Z(L)-Z1                             ! depth 
            ACO2(L)=gasO2Con(L)                    ! maximum 30 layers
            XTP1(L)=gasCO2Con(L) 
            XTP2(L)=gasN2OCon(L)
            IF(NCC>0)THEN                                                                  
                  SOC(L)=CStructLitt(L)+CMetabLitt(L)+CBiomass(L)+CSlowHumus(L)+CPassiveHumus(L) 
                  SON(L)=NStructLitt(L)+NMetabLitt(L)+NBiomass(L)+NSlowHumus(L)+NPassiveHumus(L)  
            ELSE 
                  CStructLitt(L)=.42*structLitt(L)  
                  CMetabLitt(L)=.42*metabLitt(L) 
                  CLgStructLitt(L)=.42*lgStructLitt(L) 
                  NLgStructLitt(L)=CStructLitt(L)-CLgStructLitt(L)  
            END IF   
            ! ---------------------- only C and N ----------------------------
            XZP(1,L)=CSlowHumus(L)          ! C content  
            XZP(2,L)=CPassiveHumus(L) 
            XZP(3,L)=CStructLitt(L) 
            XZP(4,L)=CMetabLitt(L) 
            XZP(5,L)=CBiomass(L)  
            XZP(6,L)=SOC(L)  
            XZP(7,L)=NSlowHumus(L) 
            XZP(8,L)=NPassiveHumus(L) 
            XZP(9,L)=NStructLitt(L)  
            XZP(10,L)=NMetabLitt(L)    
            XZP(11,L)=NBiomass(L) 
            XZP(12,L)=SON(L)  
            XZP(13,L)=SOC(L)/SON(L)  
            ! ------------------------- including C, N, P, K -------------------
            SOL(1,L)=CSlowHumus(L)  
            SOL(2,L)=CPassiveHumus(L)   
            SOL(3,L)=CStructLitt(L) 
            SOL(4,L)=CMetabLitt(L) 
            SOL(5,L)=CBiomass(L)  
            SOL(6,L)=SOC(L) 
            SOL(7,L)=NSlowHumus(L) 
            SOL(8,L)=NPassiveHumus(L)  
            SOL(9,L)=NStructLitt(L)  
            SOL(10,L)=NMetabLitt(L)  
            SOL(11,L)=NBiomass(L)   
            SOL(12,L)=SON(L)  
            SOL(13,L)=activeMineralP(L)  
            SOL(14,L)=SOP(L) 
            SOL(15,L)=stableMineralP(L)  
            SOL(16,L)=exchangeK(L)  
            SOL(17,L)=fixK(L)  
            SOL(18,L)=soilWater(L)  
            SOL(19,L)=structLitt(L)  
            SOL(20,L)=metabLitt(L)    
            SOL(21,L)=lgStructLitt(L)  
            SOL(22,L)=CLgStructLitt(L) 
            SOL(23,L)=NLgStructLitt(L)  

            ! ----------------------------- variables within plaw depth ---------------------          
            ! The IF...END IF is to calculate the some parameters that within the plaw depth.         
            IF(Z(L)<=plawDepth)THEN                                                              
                  SUM=SUM+WT(L)                                                                
                  plawDepLabileP=plawDepLabileP+labileP(L)  ! labileP: Labile P concentration (1) Estimated current labile P concentration at start of simulation.-- Zhao                                                  
                  plawDepSOC=plawDepSOC+SOC(L)              ! plawDepSOC: Organic C concentration within the plaw depth                                                  
                  plawDepSoilWater=plawDepSoilWater+soilWater(L)-Wilt_Point(L)                                                       
                  plawDepFieldCapa=plawDepFieldCapa+fieldCapacity(L)-Wilt_Point(L)                                                       
                  MXP=J                                     ! The max layer that does not exceed the plaw depth                                               
            END IF
            NO2Weight(L)=0.
            WN2O(L)=0.
            Z1=Z(L)  
      END DO  Actual_Soil_Layer  ! end of 10 soil layer loops
      
      IF(MXP==Actual_SoilLayers)THEN  
            plawDepth=Z(Layer_ID(Actual_SoilLayers))  
      ELSE  
            L1=Layer_ID(MXP+1) 
            X1=0.    
            IF(MXP>0) X1=Z(Layer_ID(MXP))     ! Looks like if plaw depth is between the MXP and MXP+1 layer   
            RTO=(plawDepth-X1)/(Z(L1)-X1)    ! the ratio between plaw depth and the first layer, but what does it mean?   
            SUM=SUM+WT(L1)*RTO  
            plawDepLabileP=plawDepLabileP+labileP(L1)*RTO   
            plawDepSOC=plawDepSOC+SOC(L1)*RTO   
            plawDepSoilWater=plawDepSoilWater+RTO*(soilWater(L1)-Wilt_Point(L1))  
            plawDepFieldCapa=plawDepFieldCapa+RTO*(fieldCapacity(L1)-Wilt_Point(L1))    
      END IF
      
      WPMX=.001*SUM                                          ! SUM is the total weight of soil within plaw depth
      ZMIX=MIN(modelPara(24),Z(Layer_ID(Actual_SoilLayers))) ! maximum depth for biologcal mixing(m): 0.1-0.3  
      ! plawDepSOC=.1*plawDepSOC/SUM   
      APBC=plawDepLabileP/WPMX   
      plawDepSOC=.001*plawDepSOC  
      aveBulkDensity=1.E-4*soilElement(3)/XX     ! Average bulk density    XX now is the 10th depth Z(10)   
      totSON=totNStructLitt+totNMetabLitt+totNBiomass+totNSlowHumus+totNPassiveHumus  ! Z... : accumulate of N, C in each layer  
      totSOC=totCStructLitt+totCMetabLitt+totCBiomass+totCSlowHumus+totCPassiveHumus  
      TWN0=totSON   
      ! the last column of XZP
      XZP(1,16)=totCSlowHumus  
      XZP(2,16)=totCPassiveHumus  
      XZP(3,16)=totCStructLitt  
      XZP(4,16)=totCMetabLitt   
      XZP(5,16)=totCBiomass  
      XZP(6,16)=totSOC  
      XZP(7,16)=totNSlowHumus  
      XZP(8,16)=totNPassiveHumus  
      XZP(9,16)=totNStructLitt  
      XZP(10,16)=totNMetabLitt  
      XZP(11,16)=totNBiomass   
      XZP(12,16)=totSON  
      XZP(13,16)=totSOC/totSON            ! ratio of C:N
       
      IF(deNitri_Method>2)THEN
            ! customerized soil layers with equal thickness 
            layersEqualThick=INT(Z(Layer_ID(Actual_SoilLayers))/layerThickness+.999)

            IF(layersEqualThick>30)THEN     ! maximum 30 layers
                  layersEqualThick=30
                  layerThickness=Z(Layer_ID(Actual_SoilLayers))/30.
            ELSE
                  IF(layersEqualThick<Actual_SoilLayers)THEN
                        layersEqualThick=Actual_SoilLayers
                        X1=layersEqualThick
                        layerThickness=Z(Layer_ID(Actual_SoilLayers))/X1
                  END IF              
            END IF
            DZ10=10.*layerThickness      ! Layer thickness for differential Eq. soln to gas diff Eqs (m)
            TOT=0.
            DO I=1,layersEqualThick
                  TOT=TOT+layerThickness
                  ZC(I)=TOT
            END DO
            CALL Interp_Concentration(ACO2,gasO2Con,Actual_SoilLayers,layersEqualThick )
            CALL Interp_Concentration(XTP1,gasCO2Con,Actual_SoilLayers,layersEqualThick)
            CALL Interp_Concentration(XTP2,gasN2OCon,Actual_SoilLayers,layersEqualThick)
            IUN=layersEqualThick-1      ! What is IUN? If deNitri_Method <= 2 then IUN will be undefined.!!!!
      END IF     


END SUBROUTINE Nutrient_Process

! ----------------------------------1.1 P flux -----------------------------------
SUBROUTINE P_Flux 
      IMPLICIT NONE
      ! EPIC1102
      ! THIS SUBPROGRAM COMPUTES P FLUX BETWEEN THE LABILE, ACTIVE MINERAL
      ! AND STABLE MINERAL P POOLS.

      ! local variables
      REAL:: RTO, RMN, X1, ROC

      RTO=MIN(.8,sorbRatioP(ISL)/(1.-sorbRatioP(ISL)))
      RMN=modelPara(77)*(labileP(ISL)-activeMineralP(ISL)*RTO)!RMN: mineral P flow rate kg/ha/d 
      X1=4.*activeMineralP(ISL)-stableMineralP(ISL)  ! if X1 < 0 , stable p ---> active P (very slow)   
      IF(X1>500.)THEN                                  ! if X1 >0 ,  active P ---> stable P 
            ROC=10.**(LOG10(flowCoeffP(ISL))+LOG10(X1))
      ELSE
            ROC=flowCoeffP(ISL)*X1
      END IF
      ROC=modelPara(78)*ROC    ! para(78): coefficient regulating p flux between active and stable pool
      stableMineralP(ISL)=stableMineralP(ISL)+ROC  ! ROC: active P ---> stable P
      activeMineralP(ISL)=activeMineralP(ISL)-ROC+RMN
      labileP(ISL)=labileP(ISL)-RMN    ! RMN: labile (soluable) P ---> active P
END SUBROUTINE P_Flux

! ---------------------------------- 1.2 P mineralization using Papran Eqs --------------------
SUBROUTINE P_mineral(Bio_contrl_fact)
      ! EPIC1102
      ! THIS SUBPROGRAM SIMULATES MINERALIZATION P USING PAPRAN EQS
      ! Ref: APEX-doc   Eqs 236 -245
      IMPLICIT NONE
 
      ! local variables
      REAL:: HMP, TKG, R4, X1, CPR, CPRF, DECR, RMP, Bio_contrl_fact

      HMP=HumusMine_RateCons(ISL)*SOP(ISL)/(NSlowHumus(ISL)+NPassiveHumus(ISL))  ! HMP: humus P minearalization rate in kg/ha/d
      TKG=cropResidu(ISL)*1000.                      ! t/ha convert to Kg/ha
      R4=.58*TKG
      X1=MAX(1.,FreshOrgP_Residu(ISL)+labileP(ISL))
      CPR=MIN(2000.,R4/X1)                            ! Eq 236c
      CPRF=EXP(-.693*(CPR-200.)/200.)                 ! Eq. 236b in APEX-doc Pg 52   CPRF: 
      DECR=MAX(.01,.05*CPRF*Bio_contrl_fact)          ! DECR: decay rate constant for fresh organic P in /d 
      RMP=DECR*FreshOrgP_Residu(ISL)                  ! Eq. 236 in APEX-doc Pg 52  RMP: is the mineralizatin rate of fresh organic P in Kg/ha/d
      IF(FreshOrgP_Residu(ISL)>1.E-5)THEN             ! FreshOrgP_Residu: fresh organic P in crop residue in kg/ha: but why it has a layer ISL?
            FreshOrgP_Residu(ISL)=FreshOrgP_Residu(ISL)-RMP
      ELSE
            RMP=FreshOrgP_Residu(ISL)
            FreshOrgP_Residu(ISL)=0.
      END IF
      LaibleP_Layer=.8*RMP+HMP      ! APEX-doc P53: 80% of RMP is added to labile P content
      X1=labileP(ISL)+LaibleP_Layer  
      IF(X1<0.)THEN
            LaibleP_Layer=labileP(ISL)
            HMP=LaibleP_Layer-.8*RMP
            labileP(ISL)=1.E-5
      ELSE
            labileP(ISL)=labileP(ISL)+LaibleP_Layer
      END IF
      SOP(ISL)=SOP(ISL)-HMP+.2*RMP   ! 20% of RMP is added to organic P pool (WPO)
END SUBROUTINE P_mineral
! -----------------------------------1.3 Daily P Loss --------------------------------------

SUBROUTINE PK_Loss 
      ! EPIC1102
      ! THIS SUBPROGRAM PREDICTS DAILY P and K LOSS, GIVEN SOIL LOSS AND
      ! ENRICHMENT RATIO
      IMPLICIT NONE
 
      ! local variables
      REAL:: X1, X2, X3, YAP, V, RTO, YMP 

      ! ------------------------- P loss in sediment -----------------------
      ! enrichRatioFinal: concentration of Organic nitrogen in the sediment devided by that in surface soil 
      sedimentLossP=enrichRatioFinal*SOP(LD1)
      SOP(LD1)=SOP(LD1)-sedimentLossP
     
      X2=labileP(LD1)
      YAP=enrichRatioFinal*X2                           ! Labile P in soil -----> sediment phase P
      X2=X2-YAP
      
      YMP=activeMineralP(LD1)*enrichRatioFinal         ! Active minearal P -----> sediment phase P
      activeMineralP(LD1)=activeMineralP(LD1)-YMP
      sedimentLossP=sedimentLossP+YMP+YAP
      ! ------------------------- P loss in runoff -------------------------
      V=Runoff+percolateFlow(LD1)+subsurfLaterFlow(LD1)
      
      IF(V>0.)THEN
            X1=MAX(5.*V,WT(LD1)*modelPara(8)) ! para(8): p con in sediment divided by that of the water (range: 10-20)
            IF(Runoff>0.)THEN
                  IF(soilP_runoff_Method>0)THEN
                        RTO=(10.*SOP(LD1)/WT(LD1))**modelPara(34)
                        runoffLossLabileP=MIN(.5*X2,X2*Runoff*RTO/X1)
                  ELSE
                        runoffLossLabileP=X2*Runoff/X1
                  END IF
                  X2=X2-runoffLossLabileP
            END IF
            X3=MIN(.5*X2,X2*(percolateFlow(LD1)+subsurfLaterFlow(LD1))/X1)
            labileP(LD1)=X2-X3
            VerticleFlow_LaibleP=X3*percolateFlow(LD1)/V
            SubSurfFlow_LaibleP=X3*subsurfLaterFlow(LD1)/V
      END IF    
      
      ! --------------------- K loss in sediment -----------------------
      X1=exchangeK(LD1)*enrichRatioFinal
      exchangeK(LD1)=exchangeK(LD1)-X1
      X3=fixK(LD1)*enrichRatioFinal
      fixK(LD1)=fixK(LD1)-X3                   ! SMM(78,MO), VAR(78): K in sediment 
      VAR(78)=X1+X3
      SMM(78,MO)=SMM(78,MO)+X1+X3
      
END SUBROUTINE PK_Loss
 
!-------------------------1.4 K Flux --------------------------
SUBROUTINE K_Flux 
      IMPLICIT NONE
      !     EPIC1102
      !     THIS SUBPROGRAM COMPUTES K FLUX BETWEEN THE WATER SOLUBLE,
      !     EXCHANGEABLE, & FIXED POOLS
 
      ! local variables
      REAL:: RSE, REF
      RSE=(solubleK(ISL)-exchangeK(ISL)*EQKS(ISL))*modelPara(29) ! Para(29) coefficients regulate flow between soluble and exchangeable K pools
      REF=(exchangeK(ISL)-fixK(ISL)*EQKE(ISL))*modelPara(22)     ! Para(22) regulates flow between exchangeable and fixed K pools
      solubleK(ISL)=MAX(.0001,solubleK(ISL)-RSE)
      exchangeK(ISL)=MAX(.0001,exchangeK(ISL)+RSE-REF)
      fixK(ISL)=fixK(ISL)+REF
END SUBROUTINE K_Flux

! -----------------------------1.5 NO3 Leaching ----------------------------

SUBROUTINE NPK_Leach(L1 )
      !     EPIC1102
      !     THIS SUBPROGRAM ESTIMATES DAILY NO3 LEACHING BY PERCOLATION AND
      !     LATERAL SUBSURFACE FLOW FOR ALL LAYERS EXCEPT THE SURFACE LAYER.
      IMPLICIT NONE
 				
      ! local variables
      INTEGER:: L1
      REAL:: V, VV, VNX, X1, X2, X3, XX 

      labileP(ISL)=labileP(ISL)+VerticleFlow_LaibleP    
      NO3_N_Soil(ISL)=NO3_N_Soil(ISL)+verticalFlowNO3(L1)
      NO2Weight(ISL)=NO2Weight(ISL)+VerticalFlow_NO2(L1)
      Weight_N2OL(ISL)=Weight_N2OL(ISL)+verticalFlowN2O(L1)
      Weight_CO2L(ISL)=Weight_CO2L(ISL)+VerticalFlow_CO2(L1)
      Weight_O2L(ISL)=Weight_O2L(ISL)+VerticalFlow_O2(L1)
      solubleK(ISL)=solubleK(ISL)+VerticleFlow_SolubleK
      saltWeight(ISL)=saltWeight(ISL)+verticleFlowSalt
      VerticleFlow_LaibleP=0.
      VerticleFlow_SolubleK=0.
      verticleFlowSalt=0.
      SubSurfFlow_LaibleP=0.
      verticalFlowNO3(ISL)=0.
      VerticalFlow_NO2(ISL)=0.
      verticalFlowN2O(ISL)=0.
      VerticalFlow_CO2(ISL)=0.
      VerticalFlow_O2(ISL)=0.
      
      V=percolateFlow(ISL)+subsurfLaterFlow(ISL)+1.E-10 ! percolation and subsurface lateral flow
      
      IF(V<1.E-5)RETURN      ! if there are no percolation and subsurface lateral flow, there is no leaching
      
      X2=1.-EXP(-V/(fracNO3Leach(ISL)*Porosity(ISL))) ! APEX-doc Eq. 148
      X1=NO3_N_Soil(ISL)-.001*WT(ISL)*modelPara(27) ! para(27) lower limit nitrate concentrations 
      
      IF(X1>0.)THEN
            X3=X1*X2           ! X3 is concentration of NO3-N after leaching
            VV=X3/V            ! Eq. 150 
            NO3_N_Soil(ISL)=NO3_N_Soil(ISL)-X3  ! Eq. 149
            verticalFlowNO3(ISL)=VV*percolateFlow(ISL)
            SubSurfFlow_NO3=VV*subsurfLaterFlow(ISL)
      END IF
      
      IF(NO2Weight(ISL)>0.)THEN
            VNX=X2*NO2Weight(ISL)
            VV=VNX/V
            NO2Weight(ISL)=NO2Weight(ISL)-VNX        
            VerticalFlow_NO2(ISL)=VV*percolateFlow(ISL)
            SubSurfFlow_NO2=VV*subsurfLaterFlow(ISL)
      END IF
      
      IF(Weight_N2OL(ISL)>0.)THEN
            VNX=X2*Weight_N2OL(ISL)
            VV=VNX/V
            Weight_N2OL(ISL)=Weight_N2OL(ISL)-VNX
            !WN2O(ISL)=WN2O(ISL)-VNX
            verticalFlowN2O(ISL)=VV*percolateFlow(ISL)
            SubSurfFlow_N2O(ISL)=VV*subsurfLaterFlow(ISL)
      END IF
      
      IF(Weight_CO2L(ISL)>0.)THEN
            VNX=X2*Weight_CO2L(ISL)
            VV=VNX/V
            Weight_CO2L(ISL)=Weight_CO2L(ISL)-VNX
            VerticalFlow_CO2(ISL)=VV*percolateFlow(ISL)
            SubSurfFlow_CO2L(ISL)=VV*subsurfLaterFlow(ISL)
      END IF
      
      IF(Weight_O2L(ISL)>0.)THEN
            VNX=X2*Weight_O2L(ISL)
            VV=VNX/V
            Weight_O2L(ISL)=Weight_O2L(ISL)-VNX
            VerticalFlow_O2(ISL)=VV*percolateFlow(ISL)
            SubSurfFlow_O2L(ISL)=VV*subsurfLaterFlow(ISL)
      END IF
      
      IF(solubleK(ISL)>0.)THEN
            X3=solubleK(ISL)*X2
            solubleK(ISL)=solubleK(ISL)-X3
            VerticleFlow_SolubleK=X3*percolateFlow(ISL)/V
            SubSurfFlow_SolubleK=X3-VerticleFlow_SolubleK
      END IF
      
      IF(saltWeight(ISL)>1.E-5)THEN
            X3=saltWeight(ISL)*X2
            saltWeight(ISL)=saltWeight(ISL)-X3
            verticleFlowSalt=X3*percolateFlow(ISL)/V
            SubSurfFlow_Salt=X3-verticleFlowSalt
      END IF
      
      IF(labileP(ISL)>1.E-5)THEN
            XX=MIN(.75,(percolateFlow(ISL)+subsurfLaterFlow(ISL))/WT(ISL))
            X3=XX*labileP(ISL)
            labileP(ISL)=labileP(ISL)-X3
            VerticleFlow_LaibleP=X3*percolateFlow(ISL)/V
            SubSurfFlow_LaibleP=X3*subsurfLaterFlow(ISL)/V
      END IF    
      
END SUBROUTINE NPK_Leach

! ------------------------ 1.6 Daily NO3 Leaching for surface layer ------------------
SUBROUTINE NO3_Leach_Surf 
      !     EPIC1102
      !     THIS SUBPROGRAM ESTIMATES DAILY NO3 LEACHING BY PERCOLATION AND
      !     LATERAL SUBSURFACE FLOW AND NO3 LOSS IN RUNOFF FOR THE SURFACE LAYER

      ! Ref: APEX-doc Eqs 145-153 
      IMPLICIT NONE
 
      ! local variables
      REAL:: V, X1, X2, X3, X4, X5, XX, Vertical_Con,COSL, CSSL, Vertical_Kcon, Horizon_Kcon, QSK, Horizon_Con 

      NO3_N_Soil(LD1)=NO3_N_Soil(LD1)+verticalFlowNO3(LD1)
      saltWeight(LD1)=saltWeight(LD1)+verticleFlowSalt
      
      verticalFlowNO3(LD1)=0.
      VerticleFlow_SolubleK=0.
      verticleFlowSalt=0.
      
      V=Runoff+subsurfLaterFlow(LD1)+percolateFlow(LD1)
      
      IF(V<=0.)RETURN
	  
      X1=NO3_N_Soil(LD1)-.001*WT(LD1)*modelPara(27)    ! NO3-N for leaching
	  
      X2=1.-EXP(-V/(fracNO3Leach(LD1)*Porosity(LD1)))
	 
      X4=percolateFlow(LD1)+modelPara(14)*(Runoff+subsurfLaterFlow(LD1)) ! Eq. 153 
 
      IF(X1>0.)THEN
            X3=X1*X2
            Vertical_Con =MIN(.01*modelPara(63),X3/X4)  ! para(63) : upper limit of N con in percolating water (ppm)
            Horizon_Con=modelPara(14)*Vertical_Con      ! para(14) : nitrate con in surface runoff / in percolate
            QNO3=Runoff*Horizon_Con
            SubSurfFlow_NO3=Horizon_Con*subsurfLaterFlow(LD1)
            verticalFlowNO3(LD1)=Vertical_Con*percolateFlow(LD1)
            NO3_N_Soil(LD1)=NO3_N_Soil(LD1)-QNO3-SubSurfFlow_NO3-verticalFlowNO3(LD1)
      END IF
      
      IF(NO2Weight(LD1)>0.)THEN
            X3=X2*NO2Weight(LD1)
            Vertical_Con=X3/X4
            Horizon_Con=modelPara(14)*Vertical_Con
            QNO2=Runoff*Horizon_Con
            SubSurfFlow_NO2=Horizon_Con*subsurfLaterFlow(LD1)
            VerticalFlow_NO2(LD1)=Vertical_Con*percolateFlow(LD1)
            NO2Weight(LD1)=NO2Weight(LD1)-QNO2-SubSurfFlow_NO2-VerticalFlow_NO2(LD1)
      END IF
      IF(Weight_N2OL(LD1)>0.)THEN
            X3=X2*Weight_N2OL(LD1)
            Vertical_Con=X3/X4
            Horizon_Con=modelPara(14)*Vertical_Con
            QN2O=Runoff*Horizon_Con
            SubSurfFlow_N2O(LD1)=Horizon_Con*subsurfLaterFlow(LD1)
            verticalFlowN2O(LD1)=Vertical_Con*percolateFlow(LD1)
            XX=QN2O+SubSurfFlow_N2O(LD1)+verticalFlowN2O(LD1)
            Weight_N2OL(LD1)=Weight_N2OL(LD1)-XX
         
      END IF
      IF(Weight_CO2L(LD1)>0.)THEN
            X3=X2*Weight_CO2L(LD1)
            Vertical_Con=X3/X4
            Horizon_Con=modelPara(14)*Vertical_Con
            QCO2=Runoff*Horizon_Con
            SubSurfFlow_CO2L(LD1)=Horizon_Con*subsurfLaterFlow(LD1)
            VerticalFlow_CO2(LD1)=Vertical_Con*percolateFlow(LD1)
            Weight_CO2L(LD1)=Weight_CO2L(LD1)-QCO2-SubSurfFlow_CO2L(LD1)-VerticalFlow_CO2(LD1)
      END IF
      IF(Weight_O2L(LD1)>0.)THEN
            X3=X2*Weight_O2L(LD1)
            Vertical_Con=X3/X4
            Horizon_Con=modelPara(14)*Vertical_Con
            QO2=Runoff*Horizon_Con
            SubSurfFlow_O2L(LD1)=Horizon_Con*subsurfLaterFlow(LD1)
            VerticalFlow_O2(LD1)=Vertical_Con*percolateFlow(LD1)
            Weight_O2L(LD1)=Weight_O2L(LD1)-QO2-SubSurfFlow_O2L(LD1)-VerticalFlow_O2(LD1)
      END IF
      IF(solubleK(LD1)>0.)THEN
            X3=solubleK(LD1)*X2
            solubleK(LD1)=solubleK(LD1)-X3
            Vertical_Kcon=X3/X4
            Horizon_Kcon=modelPara(14)*Vertical_Kcon
            VerticleFlow_SolubleK=Vertical_Kcon*percolateFlow(LD1)
            SubSurfFlow_SolubleK=Horizon_Kcon*subsurfLaterFlow(LD1)
            QSK=Horizon_Kcon*Runoff
            SMM(79,MO)=SMM(79,MO)+QSK                      ! SMM(79,MO), VAR(79): Soluble K in runoff
            VAR(79)=QSK
      END IF
      IF(saltWeight(LD1)>0.)THEN
            X5=(saltWeight(LD1)+1.E-5)*X2
            COSL=X5/X4
            CSSL=modelPara(14)*COSL
            SubSurfFlow_Salt=CSSL*subsurfLaterFlow(LD1)
            verticleFlowSalt=COSL*percolateFlow(LD1)
            XX=X5-verticleFlowSalt-SubSurfFlow_Salt
            SMM(70,MO)=SMM(70,MO)+XX                        ! SMM(70,MO), VAR(70): Soluble Salt in runoff 
            VAR(70)=XX
            saltWeight(LD1)=saltWeight(LD1)-X5
      END IF
END SUBROUTINE NO3_Leach_Surf

! ------------------------1.7 NH3 convert to NO3 or volatilization --------------------

SUBROUTINE NH3_Nir(Z5 )
      !     EPIC1102
      !     THIS SUBPROGRAM SIMULATES THE TRANSFORMATION FROM NH3 TO NO3, AND
      !     THE VOLATILIZATION OF NH3 USING MODIFIED METHODS OF REDDY AND OF THE CERES MODEL
 
      IMPLICIT NONE
 
      ! local variables
      REAL:: X1, XX, F, FAF, AKAV, AKAN, FCEC, FZ, FPH, Z5, AVOL, RNIT, RNV, GNO2, GNO3

      X1=.41*(soilTem(ISL)-5.)                 ! APEX-doc Eq.230b   X1: temperature factor
      IF(X1<=0.)RETURN
      
      IF(ISL==LD1)THEN
            FAF=.335+.16*LOG(U10+.2)              ! APEX-doc Eq.231  
            AKAV=X1*FAF     ! for top layer
      ELSE                                      
            FCEC=MAX(.3,1.-.038*CEC(ISL))
            FZ=1.-Z5/(Z5+EXP(S_Curve(12,1)-S_Curve(12,2)*Z5))
            AKAV=X1*FCEC*FZ ! other layers
      END IF
      
      IF(PH(ISL)>7.)THEN
            IF(PH(ISL)>7.4)THEN
                  FPH=5.367-.599*PH(ISL)
            ELSE
                  FPH=1.
            END IF
      ELSE
            FPH=.307*PH(ISL)-1.269
      END IF
      
      AKAN=X1*SoilWater_Fact(ISL)*FPH
      AKAV=AKAV*SoilWater_Fact(ISL)
      XX=AKAV+AKAN
      IF(XX<1.E-5)RETURN
      F=MIN(modelPara(64),1.-EXP(-XX))   ! F : coefficient of combined nitrification and volatilization   
      RNV=F*NH3_Weight(ISL)               ! RNV: combined nitrification and volization  (kg/ha/d) NH3_Weight: kg/ha
      AVOL=RNV*modelPara(57)             ! para(57): Volatilization/Nitrification partitioning coefficient (0.05-0.5) 
      Vol_NH3=Vol_NH3+AVOL                ! AVOL caused Nirt_NH3 difference when ISL = 3 at 01-07
     
      IF(N_vol_Method==0)THEN
            RNIT=RNV-AVOL                   ! RNIT: nitrification                      kg/ha/d
            NH3_Weight(ISL)=NH3_Weight(ISL)-AVOL-RNIT
            NO3_N_Soil(ISL)=NO3_N_Soil(ISL)+RNIT
            Nitri_NH3=Nitri_NH3+RNIT
      ELSE
            NH3_Weight(ISL)=NH3_Weight(ISL)-AVOL
            F=MIN(1.,modelPara(64)*AKAN)    ! para(64): upper limit of nitrification-volatilization as a fraction of NH3 percent
            GNO2=F*NH3_Weight(ISL)
            NH3_Weight(ISL)=NH3_Weight(ISL)-GNO2
            NO2Weight(ISL)=NO2Weight(ISL)+GNO2
            IF(PH(ISL)>5.5)THEN
                  IF(PH(ISL)>7.2)THEN
                        FPH=4.367-.5324*PH(ISL)
                  ELSE
                        FPH=1.
                  END IF
            ELSE
                  FPH=.307*PH(ISL)-1.269
            END IF
            AKAN=X1*SoilWater_Fact(ISL)*FPH
            F=MIN(1.,modelPara(64)*AKAN)
            GNO3=F*NO2Weight(ISL)
            NO2Weight(ISL)=NO2Weight(ISL)-GNO3
            NO3_N_Soil(ISL)=NO3_N_Soil(ISL)+GNO3
            Nitri_NH3=Nitri_NH3+GNO2
      END IF        
END SUBROUTINE NH3_Nir
! ---------------------------------1.8 Crop_Fall --------------------------------------
SUBROUTINE Crop_Fall(ZZ )
      !     EPIC1102
      !     THIS SUBPROGRAM SIMULATES THE CONVERSION OF STANDING DEAD CROP
      !     RESIDUE TO FLAT RESIDUE.
      IMPLICIT NONE
 
      ! local variables
      INTEGER:: K
      REAL:: SUM, TOT, X1, X2, ZZ, ZS

      SUM=0.
      TOT=0.
      DO K=1,LC
            IF(standCropResi(K)<.001)CYCLE
          
            X1=ZZ*standCropResi(K)
            standCropResi(K)=standCropResi(K)-X1
            SUM=SUM+X1
          
            X2=ZZ*standDeadResiN(K)
            standDeadResiN(K)=standDeadResiN(K)-X2
            TOT=TOT+X2
          
            ZS=ZZ*standDeadResiP
            FreshOrgP_Residu(LD1)=FreshOrgP_Residu(LD1)+ZS
            standDeadResiP=standDeadResiP-ZS
          
            ZS=ZZ*standDeadResiK
            solubleK(LD1)=solubleK(LD1)+ZS
            standDeadResiK=MAX(1.E-5,standDeadResiK-ZS)
          
            ZS=ZZ*standDeadLignin
            standDeadLignin=standDeadLignin-ZS
      END DO
      
      ZZ=MIN(1.,ZZ*10.)      
      ZS=ZZ*standDeadResOrganism
      standDeadResOrganism=MAX(1.E-5,standDeadResOrganism-ZS)
      SUM=SUM+ZS     
      ZS=ZZ*standDeadResiOrgN
      TOT=TOT+ZS
      
      CALL CN_TopSoil(SUM, TOT, 0 )
      
      standDeadResiOrgN=MAX(1.E-5,standDeadResiOrgN-ZS)
      ZS=ZZ*StandDeadRes_OrgP
      FreshOrgP_Residu(LD1)=FreshOrgP_Residu(LD1)+ZS
      StandDeadRes_OrgP=MAX(1.E-5,StandDeadRes_OrgP-ZS)
      ZS=ZZ*StandDeadRes_OrgK
      solubleK(LD1)=solubleK(LD1)+ZS
      StandDeadRes_OrgK=MAX(1.E-5,StandDeadRes_OrgK-ZS)
END SUBROUTINE Crop_Fall

! ------------------------------1.9 simulate C respiration by Century ------------------
SUBROUTINE Sim_SoilResp(Bio_contrl_fact)
      !     EPIC1102
      !     THIS SUBPROGRAM SIMULATES C RESPIRATION USING EQUATIONS TAKEN FROM CENTURY.
      ! Ref: APEX-doc Eqs. 180-190
      IMPLICIT NONE
      ! local variables
      real:: XZ, RLR, APCO2, ASCO2, ABCO2, A1CO2, A1, ABP, ASP,ASX,APX,  RBM, RLM, RLS, HPNC, XBM, BMNC, HSNC, &
            XLMNC, XLSNC, CRX, CR, X1, X2, TLSLNCA, TLSLCA, TLSNA, TLSCA, TLMCA, TLMNA, TBMCA, TBMNA, THSCA, &
            THSNA, THPCA, THPNA,Bio_contrl_fact
 
      XZ=NO3_N_Soil(ISL)+NO2Weight(ISL)+NH3_Weight(ISL)
      ! AD1=NStructLitt(ISL)+NMetabLitt(ISL)+NBiomass(ISL)+NSlowHumus(ISL)+NPassiveHumus(ISL)+XZ
      RLR=lgStructLitt(ISL)/(structLitt(ISL)+1.E-5)
      IF(RLR>.8)THEN
            RLR=.8
      ELSE
            IF(RLR<.1)RLR=.1
      END IF    
      APCO2=.55      ! APEX-doc : allocation from passive humus to CO2 (0.55)
      ASCO2=.60      ! APEX-doc : allocation from slow humus to CO2 (0.55)
      IF(ISL==LD1)THEN
            Bio_contrl_fact=Bio_contrl_fact*modelPara(51) ! para(51): O2 modifying microbial activity (0.8-0.95)
            ABCO2=.55       ! Allocation from biomass to CO2 (surface: 0.6)
            A1CO2=.55       ! ALMCO2 + ALSCO2 ? allocation from (Metabolic + Structural) CO2
            RBM=.0164       ! Transformation rate of microbial biomass under optimal cond (surface: 0.0164/d)
            RLM=.0405       ! Metabolic litter transformation rate under optimal cond (surface: 0.0405 /d)
            RLS=.0107       ! Structural litter transformation rate under optimal conditions (surface: 0.0107 /d)
            HPNC=.1
            XBM=1.          ! transformation rate of microbial biomass control by soil texture and structure (surface: 1)
      ELSE
            ABCO2=.17+.0068*sandFrac(ISL) ! other layers
            A1CO2=.55
            RBM=.02        ! other layers: 0.02 /d
            RLM=.0507      ! other layers: 0.0507 /d
            RLS=.0132      ! other layers: 0.0132 /d   
            XBM=.25+.0075*sandFrac(ISL) ! other layers: 
            HPNC=NPassiveHumus(ISL)/CPassiveHumus(ISL)  ! N/C ratio of the passive humus
      END IF
      BMNC=NBiomass(ISL)/CBiomass(ISL)                  ! APEX-doc NCBM: N/C ratio of biomass formed from surface litter 
      HSNC=NSlowHumus(ISL)/CSlowHumus(ISL)              ! N/C ratio of the slow humus
      CMetabLitt(ISL)=MAX(.01,CMetabLitt(ISL))
      XLMNC=NMetabLitt(ISL)/CMetabLitt(ISL)
      CStructLitt(ISL)=MAX(.01,CStructLitt(ISL))
      XLSNC=NStructLitt(ISL)/CStructLitt(ISL)
      ABP=.003+.00032*clayFrac(ISL)                      ! ABP: Allocation from biomass to passive humus: 0 surface layer, others: all other layers
      ASP=MAX(.001,modelPara(45)-.00009*clayFrac(ISL))   ! ASP: Allocation from slow humus to passive; 0.001: surface layer, others: all other layers
      A1=1.-A1CO2
      ASX=1.-ASCO2-ASP
      APX=1.-APCO2
      ! TRANSFORMATION STRUCTURAL LITTER WBM 2012-11-22
      ! CR = scaling factor for decomposition (0 - 1) basewd BMNC and &
      ! ratio of NC-substrate/Yield. 
      CRX=1.-(runoffAdj-BMNC)/(runoffAdj-dampDepAdj)
      CR=0.
      IF(BMNC>dampDepAdj.AND.BMNC<runoffAdj)CR=CRX
      IF(NStructLitt(ISL)/(CStructLitt(ISL)*(1-RLR))/A1>BMNC.OR.BMNC>runoffAdj)CR=1.
      X1=RLS*Bio_contrl_fact*CR*EXP(-3.*RLR)
      X2=1.-RLR
      TLSLNCA=X1*CStructLitt(ISL)*X2       ! APEX-doc  Eq.181
      TLSLCA=TLSLNCA*RLR/X2                   ! APEX-doc  Eq.180    Structual Litter ---> microbes (Active)
      TLSCA=TLSLNCA/X2                        ! APEX-doc  Eq.177
      TLSNA=TLSCA*XLSNC                       ! APEX-doc  Eq.182
      
      ! TRANSFORMATIONS METABOLIC LITTER; SIX LINES BELOW ADDED WBM 2012-11-14
      CR=0.
      IF(BMNC>dampDepAdj.AND.BMNC<runoffAdj)CR=CRX
      IF(XLMNC/A1>BMNC.OR.BMNC>runoffAdj)CR=1.
      X1=RLM*Bio_contrl_fact*CR
      TLMCA=CMetabLitt(ISL)*X1           ! APEX-doc  Eq.183   C:  Metabolic Litter ---> microbes (Active)
      TLMNA=TLMCA*XLMNC                      ! APEX-doc  Eq.184   N: 
      ! TRANSFORMATIONS MICROBIAL BIOMASS; NEXT SIX LINES ADDED WBM 2012-11-22
      ! CBiomass: C in soil microbial biomass and associate productes (kg/ha)
      CR=0.
      IF(BMNC>dampDepAdj.AND.BMNC<runoffAdj)CR=CRX
      IF(BMNC/(1.-ABCO2)>BMNC.OR.BMNC>runoffAdj)CR=1.
      TBMCA=CBiomass(ISL)*RBM*Bio_contrl_fact*CR*XBM   ! APEX-doc  Eq.185   transformation of microbial biomass C
      TBMNA=TBMCA*BMNC                                  ! APEX-doc  Eq.186   transformation of microbial biomass N 
      ! TRANSFORMATION OF SLOW HUMUS; NEXT SIX LINES ADDED WBM 2012-11-22
      CR=0.
      IF(BMNC>dampDepAdj.AND.BMNC<runoffAdj)CR=CRX
      IF(HSNC/(1.-ASCO2)>BMNC.OR.BMNC>runoffAdj)CR=1.
      THSCA=CSlowHumus(ISL)*slowHumusTranRate*Bio_contrl_fact*CR ! APEX-doc Eq.187 transformation of slow humus C
      THSNA=THSCA*HSNC                                         ! APEX-doc Eq.188 transformation of slow humus N
      ! TRANSFORMATIONS OF PASSIVE HUMUS; NEXT SIX LINES ADDED WBM 2012-11-22
      CR=0.
      IF(BMNC>dampDepAdj.AND.BMNC<runoffAdj)CR=CRX
      IF(HPNC/APX>BMNC.OR.BMNC>runoffAdj)CR=1.
      THPCA=CPassiveHumus(ISL)*Bio_contrl_fact*CR*passivHumusTranRate ! APEX-doc Eq.189 transformation of passive humus C
      THPNA=THPCA*HPNC                                         ! APEX-doc Eq.188 transformation of passive humus N
      IF(ISL==LD1) THEN
            soilResp(ISL)=MAX(0.,.3*TLSLCA+A1CO2*(TLSLNCA+TLMCA)+ABCO2*TBMCA+ASCO2*THSCA)   ! RSPC
      ELSE
            soilResp(ISL)=MAX(0.,.3*TLSLCA+A1CO2*(TLSLNCA+TLMCA)+ABCO2*TBMCA+ASCO2*THSCA+APCO2*THPCA)
      END IF    
END SUBROUTINE Sim_SoilResp

! ------------------------------------1 Master Nutrient Cycling ----------------------

SUBROUTINE Nutrient_Cyc 
      !     EPIC1102
      !     THIS SUBPROGRAM IS THE MASTER NUTRIENT CYCLING SUBPROGRAM. CALLS NPMIN, NYNIT, NLCH, NCNMI, AND NDNIT FOR EACH SOIL
      !     LAYER.
      IMPLICIT NONE
      ! local variables
      INTEGER:: K, IRTO, ISL0, L1,J
      REAL:: STDX,  QSK, TSFCO2, TSFO2,TSFP, XX, AD1, AD2, VNO31, Layer_Thickness, X1, X2, X3, X4, Z5, F, &
            RTO, ZZ, O2_Factor,Bio_contrl_fact
      
      STDX=0.     
      DO K=1,LC
            STDX=STDX+standCropResi(K)
      END DO 
      SGMN=0.
      totDenitriN=0.
      SN2=0.
      totN2O=0.
      SMP=0.
      SIP=0.
      QNO3=0.
      QNO2=0.
      QN2O=0.
      QCO2=0.
      QO2=0.
      QSK=0.
      VerticleFlow_LaibleP=0.
      VerticalFlow_NO2=0.
      verticalFlowN2O=0.
      VerticalFlow_CO2=0.
      VerticalFlow_O2=0.
      SubSurfFlow_N2O=0.
      SubSurfFlow_CO2L=0.
      SubSurfFlow_O2L=0.
      SubSurfFlow_SolubleK=0.
      SubSurfFlow_Salt=0.
      TSFNO3=0.
      TSFNO2=0.
      TSFN2O=0.
      TSFCO2=0.
      TSFO2=0.
      TSFK=0.
      TSFS=0.
      Vol_NH3=0.
      Nitri_NH3=0.
      Tot_SoilRes=0.
      XX=0.
      AD1=0.
      AD2=0.
      VNO31=verticalFlowNO3(LD1)
      IRTO=0
      WBMX=0.
      
      DO J=1,Actual_SoilLayers   ! soil layer loop
            ISL=Layer_ID(J)
            ISL0=ISL
            AD1=AD1+NO3_N_Soil(ISL)+NO2Weight(ISL)+Weight_N2OL(ISL)
            soilResp(ISL)=0.
            XNetNMineralize(ISL)=0.
            SubSurfFlow_NO3=0.
            SubSurfFlow_NO2=0.
            Layer_Thickness=Z(ISL)-XX      ! XX is Z(ISL-1)
            X1=soilWater(ISL)-Wilt_Point(ISL)
            IF(X1<0.)THEN
                  SoilWater_Fact(ISL)=.1*(soilWater(ISL)/Wilt_Point(ISL))**2    ! Eq. 163 and 163a in APEX-doc P44
            ELSE
                  SoilWater_Fact(ISL)=MIN(1.,.1+.9*SQRT(X1/(fieldCapacity(ISL)-Wilt_Point(ISL))))
            END IF
          
            CALL P_Flux 
          
            CALL K_Flux  
          
            IF(ISL/=LD1)THEN     ! Not top layer: NPK leach from the top layer to other layers
                  L1=Layer_ID(J-1) ! The last layer 
                  CALL NPK_Leach(L1)
            ELSE               ! Top layer: orgainc materials decay to soil
                  ZZ=MIN(.5,modelPara(66)*(1.+.1*Rainfall)) ! para(66): convert standing res to flat res
                  IF(STDX+standDeadResOrganism>.001)&
                  CALL Crop_Fall(ZZ)
                  IF(landUseCode/=35.AND.Rainfall>0.)THEN     ! 35: pavement and urban area
                        CALL NO3_Leach_Surf 
                        CALL PK_Loss 
                        SMM(54,MO)=SMM(54,MO)+sedimentLossP   ! SMM(54,MO), VAR(54): Sediment phase P lost (kg/ha)
                        VAR(54)=sedimentLossP
                  END IF
            END IF
          
            TSFNO3=TSFNO3+SubSurfFlow_NO3
            TSFNO2=TSFNO2+SubSurfFlow_NO2
            TSFN2O=TSFN2O+SubSurfFlow_N2O(ISL)
            TSFP=TSFP+SubSurfFlow_LaibleP
            TSFCO2=TSFCO2+SubSurfFlow_CO2L(ISL)
            TSFO2=TSFO2+SubSurfFlow_O2L(ISL)
            TSFK=TSFK+SubSurfFlow_SolubleK
            TSFS=TSFS+SubSurfFlow_Salt
            AD2=AD2+NO3_N_Soil(ISL)+NO2Weight(ISL)+Weight_N2OL(ISL)  
          
            IF(ISL==IDR)THEN
                  VAR(53)=SubSurfFlow_NO3+SubSurfFlow_NO2+SubSurfFlow_N2O(ISL)
                  SMM(53,MO)=SMM(53,MO)+VAR(53)                   ! SMM(53,MO), VAR(53): Soluble Nitrogen in drainage outflow
                  SMM(114,MO)=SMM(114,MO)+SubSurfFlow_LaibleP     ! SMM(114,MO), VAR(114): Soluble P in drainage outflow
                  VAR(114)=SubSurfFlow_LaibleP
            END IF
            ! ---------------------- NH3 Nitrification and Volatilization -----------------------
          
            Z5=.5*(Z(ISL)+XX)*1000    ! XX: the last layer Z(ISL-1) depth to the middle of a soil layer in mm
          
            IF(NH3_Weight(ISL)>.01)CALL NH3_Nir(Z5)  
		  
            ! ------------------------ C Transformation  ---------------------------
          
            IF(soilTem(ISL)>0.)THEN             
                  Tem_Fact(ISL)=soilTem(ISL)/(soilTem(ISL)+EXP(S_Curve(14,1)-S_Curve(14,2)*soilTem(ISL))) ! soil temperature factor used in regulating microbal processes
                  X1 =MIN(20.,EXP(modelPara(52)*(bulkDensity(ISL)-postTillBulkDensity(ISL))))            ! bulkDensity: bulk density, postTillBulkDensity: post tillage bulk density
                  O2_Factor=1.-modelPara(53)*Z5/(Z5+EXP(S_Curve(20,1)-S_Curve(20,2)*Z5))! O2_Factor: oxygen factor controlling biological processes as a function of depth
              
                  IF(deNitri_Method>2.AND.gasO2Con(LD1)>0..AND.O2_Method==2)O2_Factor=gasO2Con(ISL)/gasO2Con(LD1)
              
                  X2=SQRT(Tem_Fact(ISL)*SoilWater_Fact(ISL))*O2_Factor    ! Eq. 179  combination factor of T, SW, O2
              
                  X3=1.
              
                  IF(soilWater(ISL)>1.1*fieldCapacity(ISL))X3=0.
              
                  IF(ISL==LD1)THEN       ! top layer 
                        X4=cropResidu(LD1)
                        F=1.+10.*X4/(X4+EXP(S_Curve(27,1)-S_Curve(27,2)*X4))
                  ELSE                   ! other layers
                        F=1.
                  END IF         
              
                  IF(Z(ISL)<ZMIX)THEN
                        WBMX=WBMX+Layer_Thickness*X2*X3*F
                  ELSE
                        IF(IRTO==0)THEN
                              X1=ZMIX-XX
                              RTO=X1/Layer_Thickness
                              WBMX=WBMX+X1*X2*RTO*X3
                              IRTO=1
                        END IF
                  END IF    
                  Bio_contrl_fact=MIN(10.,X2*modelPara(20)*X1) 
                  CALL Sim_SoilResp(Bio_contrl_fact) 
                  CALL P_mineral(Bio_contrl_fact )
                  SMP=SMP+ LaibleP_Layer    ! accumulate soil laible P in each layer
            END IF          
            XX=Z(ISL)
      END DO   ! Soil layer loop
      ! WBMX: weight of biomass mixing?
      WBMX=WBMX*modelPara(25)/ZMIX   ! para(25): biological mixing efficiency simulates mixing in top soil by earth worms : 0.1-0.5
      solubleNConc=solubleNConc+verticalFlowNO3(ISL)
      SMM(57,MO)=SMM(57,MO)+VerticleFlow_LaibleP         ! SMM(57,MO), VAR(57) : labile P in verticle flow 
      VAR(57)=VerticleFlow_LaibleP
      SMM(81,MO)=SMM(81,MO)+VerticleFlow_SolubleK        ! SMM(81,MO), VAR(81) : soluble K in verticle flow
      VAR(81)=VerticleFlow_SolubleK
      SMM(82,MO)=SMM(82,MO)+verticleFlowSalt            ! SMM(82,MO), VAR(82) : salt in verticle flow  
      VAR(82)=verticleFlowSalt
END SUBROUTINE Nutrient_Cyc
!======================================================================================================
!----------------------------------------- 2. Phoenix ----------------------------------
SUBROUTINE Phoenix( Z5, RTOE )
      !     EPIC1102
      !     THIS SUBPROGRAM SIMULATES MINERALIZATION AND IMMOBILIZATION OF N AND C USING POOLS 
      !     FOLLOWING CENTURY (IZAURRALDE ET AL. 2006) 
      !     & C/N OF MICROBIAL BIOMASS FOLLOWING PHOENIX (MCGILL ET AL. 1981)
      IMPLICIT NONE
                          
      ! local variables
      REAL:: SDF = 0.0, XZ, RLR,AD1, AD2, XHN,X1,X2, X3, O2_Factor,Z5, Y1,Y2,Y3,RTOE,APCO2, ASCO2,RLM, RLS, &
      ABCO2, A1CO2,RBM, HPNC, XBM, BMNC, HSNC, XLMNC, XLSNC, ABP, ASP, A1, ASX, APX, &
      CRX, CR, CR0,TLSCA, TLSLNCA, TLSLCA, TLSNA, TLMCA, TLMNA, TBMCA, TBMNA, THPCA,THPNA, &
      THSNA, THSCA, XSF, CA, CU, XU, XDENOM, DMDNMAX, XWNH3, XWNO2, XWNO3, Bio_contrl_fact,&
      DMDNH3MX, DMDNO2MX, DMDNO3MX, DMDN, XX, TUPH3, TUPO2, TUPO3, XNIMO, RTO, DF, &
      WKA, WNCMAX, WNCMIN, VMU, WKMNH3, WKMNO2, WKMNO3

      ! NO2Weight (ISL) added by WBM 2012-11-05
      XZ=NO3_N_Soil(ISL)+NO2Weight(ISL)+NH3_Weight(ISL)
      !AD1=NStructLitt(ISL)+NMetabLitt(ISL)+NBiomass(ISL)+NSlowHumus(ISL)+NPassiveHumus(ISL)+XZ
      RLR=lgStructLitt(ISL)/(structLitt(ISL)+1.E-5)
      IF(RLR>.8)THEN
            RLR=.8
      ELSE
            IF(RLR<.1)RLR=.1
      END IF    
      AD1=CBiomass(ISL)+CPassiveHumus(ISL)+CSlowHumus(ISL)+CMetabLitt(ISL)+CStructLitt(ISL)
      XHN=NSlowHumus(ISL)+NPassiveHumus(ISL)
      IF(bulkDensity(ISL)>postTillBulkDensity(ISL))THEN
            X1=MIN(20.,EXP(modelPara(52)*(bulkDensity(ISL)-postTillBulkDensity(ISL)))) ! para(52)coefficient in equation expressing tillage effect on residue decay rate
      ELSE    
            X1=1.
      END IF
      SMS(4,ISL)=SMS(4,ISL)+X1
      O2_Factor=1.-modelPara(53)*Z5/(Z5+EXP(S_Curve(20,1)-S_Curve(20,2)*Z5))
      SELECT CASE(O2_Method)
            CASE(1)
                  Y1=.1*SOC(ISL)/WT(ISL)
                  Y2=2.11+.0375*clayFrac(ISL)
                  Y3=100.*Y1/Y2
                  O2_Factor=Y3/(Y3+EXP(S_Curve(24,1)-S_Curve(24,2)*Y3))
            CASE(2)
                  IF(gasO2Con(LD1)>0.)O2_Factor=gasO2Con(ISL)/gasO2Con(LD1)
      END SELECT
      IF(deNitri_Method>2)O2_Factor=O2_Factor*RTOE    ! PARA(20):no introduction of para(20) in EPIC0810 manual
      Bio_contrl_fact=MIN(10.,SQRT(Tem_Fact(ISL)*SoilWater_Fact(ISL))*modelPara(20)*O2_Factor*X1) ! Eq. 179 in APEX-doc P46
      SMS(11,ISL)=SMS(11,ISL)+1.
      SMS(1,ISL)=SMS(1,ISL)+SoilWater_Fact(ISL)
      APCO2=.55
      ASCO2=.60
      IF(ISL==LD1)THEN
            IF(cropCode(Crop_Num)==plantCategoryCode(1))THEN
                  Bio_contrl_fact=1.
                  RLM=.2
                  RLS=.2
            ELSE    
                  Bio_contrl_fact=Bio_contrl_fact*modelPara(51)
                  RLM=.0405
                  RLS=.0107
            END IF
            ABCO2=.55
            A1CO2=.55
            RBM=.0164
            HPNC=.1
            XBM=1.
      ELSE
            ABCO2=.17+.0068*sandFrac(ISL)
            A1CO2=.55
            RBM=.02
            RLM=.0507
            RLS=.0132
            XBM=.25+.0075*sandFrac(ISL)
            HPNC=NPassiveHumus(ISL)/CPassiveHumus(ISL) 
      END IF
      BMNC=NBiomass(ISL)/CBiomass(ISL)
      HSNC=NSlowHumus(ISL)/CSlowHumus(ISL)
      XLMNC=NMetabLitt(ISL)/CMetabLitt(ISL)
      XLSNC=NStructLitt(ISL)/CStructLitt(ISL)
      ABP=.003+.00032*clayFrac(ISL)
      SMS(3,ISL)=SMS(3,ISL)+Bio_contrl_fact
      ASP=MAX(.001,modelPara(45)-.00009*clayFrac(ISL)) ! Para(45): coefficient allocating slow to passive humus 0.001-0.05
      A1=1.-A1CO2
      ASX=1.-ASCO2-ASP
      APX=1.-APCO2
      ! TRANSFORMATION OF STRUCTURAL LITTER WBM 2012-11-14
      ! CR = scaling factor for decomposition (0 - 1) based BMNC and &
      ! ratio of NC-substrate/Yield.
      CRX=1.-(runoffAdj-BMNC)/(runoffAdj-dampDepAdj)
      CR=0.
      IF(BMNC>dampDepAdj.AND.BMNC<runoffAdj)CR=CRX
      CR0=CR
      X2=1.-RLR
      IF((NStructLitt(ISL)/(CStructLitt(ISL)*X2))/A1>BMNC.OR.BMNC>runoffAdj)CR=1.
      TLSCA=RLS*Bio_contrl_fact*CR*EXP(-3.*RLR)*CStructLitt(ISL)
      TLSLNCA=TLSCA*X2
      TLSLCA=TLSCA*RLR
      TLSNA=TLSCA*XLSNC
      ! TRANSFORMATIONS METABOLIC LITTER; SIX LINES BELOW ADDED WBM 2012-11-14
      CR=CR0
      IF(XLMNC/A1>BMNC.OR.BMNC>runoffAdj)CR=1.
      TLMCA=CMetabLitt(ISL)*RLM*Bio_contrl_fact*CR
      TLMNA=TLMCA*XLMNC
      ! TRANSFORMATIONS MICROBIAL BIOMASS; NEXT SIX LINES ADDED WBM 2012-11-14
      CR=CR0
      IF(BMNC/(1.-ABCO2)>BMNC.OR.BMNC>runoffAdj)CR=1.
      TBMCA=CBiomass(ISL)*RBM*Bio_contrl_fact*CR*XBM
      TBMNA=TBMCA*BMNC
      ! TRANSFORMATION OF SLOW HUMUS; NEXT SIX LINES ADDED WBM 2012-11-14
      CR=CR0
      IF(HSNC/(1.-ASCO2)>BMNC.OR.BMNC>runoffAdj)CR=1.
      THSCA=CSlowHumus(ISL)*slowHumusTranRate*Bio_contrl_fact*CR
      THSNA=THSCA*HSNC
      ! TRANSFORMATIONS OF PASSIVE HUMUS; NEXT SIX LINES ADDED WBM 2012-11-14
      IF(ISL/=LD1)THEN
            CR=CR0
            IF(HPNC/APX>BMNC.OR.BMNC>runoffAdj)CR=1.
            THPCA=CPassiveHumus(ISL)*Bio_contrl_fact*CR*passivHumusTranRate
            THPNA=THPCA*HPNC
            ! MINERALIZATION AND IMMOBILIZATION (32 LINES) ADDED WBM 2012-11-14
            ! CA = AMMONIFICATION RATE SCALING FACTOR BASED ON BMNC (0 - 1)
            ! CU = IMMOPBILIZATION (UPTAKE) RATE SCALING FACTOR BASED ON BMNC (0 - 1)
            ! WKA = SPECIFIC BASE RATE FOR AMMONIFICATION (d-1)
            ! WNCMIN = BMNC AT WHICH IMMOBILIZATION IS A MAXIMUM; BMNC AT WHICH AMMONIFICATION CEASES
            ! WNCMAX = BMNC AT WHICH IMMOBILIZATION CEASES; BMNC AT WHICH AMMONIFICATION IS A MAXIMUM
            ! VMU = MAXIMUM RATE OF UPTAKE OF N DURING IMMOBILIZATION (gN (gC-1) d-1)
            ! WKMNH3 = HALF SATURATION CONSTANT FOR AMMONIA IMMOBILIZATION (mg N L-1)
            ! WKMNO2 = HALF SATURATION CONSTANT FOR NITRITE IMMOBILIZATION (mg N L-1)
            ! WKMNO3 = HALF SATURATION CONSTANT FOR NITRATE IMMOBILIZATION (mg N L-1)
      END IF
      ! modified by TXH --- 2020-04-14
      WKA=modelPara(84)
      WNCMIN=modelPara(85)
      WNCMAX=modelPara(86)  
      VMU=modelPara(87)
      WKMNH3=modelPara(88)
      WKMNO2=modelPara(89)
      WKMNO3=modelPara(90)  
      XSF=(BMNC-WNCMIN)/(WNCMAX-WNCMIN)
      CA=0.
      IF(BMNC>WNCMIN.AND.BMNC<WNCMAX)CA=XSF
      IF(BMNC>WNCMAX)CA=1.
      CU=.001
      IF(BMNC>WNCMIN.AND.BMNC<WNCMAX)CU=1.-XSF
      IF(BMNC<WNCMIN)CU=1.
      ! AMMONIFICATION (N MINERLIZATION) PRODUCRT OF CA*Bio_contrl_fact*WKA*BMN
      SUP=MAX(0.,CA*Bio_contrl_fact*WKA*NBiomass(ISL))
      ! MICROBIAL UPTAKE OF N (N IMMOBILIZATION) MICHAELIS-MENTEN EXPRESION &
      ! DRIVEN BY MICROBIAL CARBON AND REGUALTED BY CU
      XU=Bio_contrl_fact*CU*VMU
      XDENOM=.01*soilWater(ISL)
      DMDNMAX=XU*CBiomass(ISL)
      XWNH3=NH3_Weight(ISL)/XDENOM
      XWNO2=NO2Weight(ISL)/XDENOM
      XWNO3=NO3_N_Soil(ISL)/XDENOM
      DMDNH3MX=MIN(NH3_Weight(ISL),XWNH3*DMDNMAX/(WKMNH3+XWNH3))
      DMDNO2MX=MIN(NO2Weight(ISL),XWNO2*DMDNMAX/(WKMNO2+XWNO2))
      DMDNO3MX=MIN(NO3_N_Soil(ISL),XWNO3*DMDNMAX/(WKMNO3+XWNO3))
      DMDN=DMDNH3MX+DMDNO2MX+DMDNO3MX
      XX=1.
      IF(DMDN>DMDNMAX)XX=DMDNMAX/DMDN
      TUPH3=DMDNH3MX*XX
      TUPO2=DMDNO2MX*XX
      TUPO3=DMDNO3MX*XX
      SGMN=SGMN+SUP
      XNIMO=-TUPH3-TUPO2-TUPO3
      XNetNMineralize(ISL)=SUP+XNIMO
      SMNIM=SMNIM+XNIMO        ! ???
      NO3_N_Soil(ISL)=MAX(1.E-10,NO3_N_Soil(ISL)-TUPO3)
      NO2Weight(ISL)=MAX(1.E-10,NO2Weight(ISL)-TUPO2)
      NH3_Weight(ISL)=MAX(1.E-10,NH3_Weight(ISL)-TUPH3+SUP)
      mineralizedN=mineralizedN+XNetNMineralize(ISL)
      SMS(9,ISL)=SMS(9,ISL)+XNetNMineralize(ISL)
      IF(TLSCA>CStructLitt(ISL))TLSCA=CStructLitt(ISL)
      CStructLitt(ISL)=CStructLitt(ISL)-TLSCA
      IF(TLSLCA>CLgStructLitt(ISL))TLSLCA=CLgStructLitt(ISL)
      CLgStructLitt(ISL)=CLgStructLitt(ISL)-TLSLCA
      NLgStructLitt(ISL)=CStructLitt(ISL)-CLgStructLitt(ISL)
      TLMCA=MIN(CMetabLitt(ISL),TLMCA)
      IF(metabLitt(ISL)>0.)THEN
            RTO=MAX(.42,CMetabLitt(ISL)/metabLitt(ISL))
            metabLitt(ISL)=metabLitt(ISL)-TLMCA/RTO
            CMetabLitt(ISL)=MAX(.01,CMetabLitt(ISL)-TLMCA)
      END IF
      lgStructLitt(ISL)=lgStructLitt(ISL)-TLSLCA/.42
      structLitt(ISL)=CStructLitt(ISL)/.42
      IF(ISL==LD1)THEN 
            X3=ASX*THSCA+A1*(TLMCA+TLSLNCA)
            X1=.7*TLSLCA+TBMCA*(1.-ABCO2)
            NBiomass(ISL)=NBiomass(ISL)-TBMNA+THSNA*(1.-ASP)+TLMNA+TLSNA-XNetNMineralize(ISL)
            !Soil_Respiration(ISL)=.3*TLSLCA+A1CO2*(TLSLNCA+TLMCA)+ABCO2*TBMCA+ASCO2*THSCA
      ELSE
            X3=APX*THPCA+ASX*THSCA+A1*(TLMCA+TLSLNCA)
            X1=.7*TLSLCA+TBMCA*(1.-ABP-ABCO2)
            X2=THSCA*ASP+TBMCA*ABP
            CPassiveHumus(ISL)=CPassiveHumus(ISL)-THPCA+X2
            NBiomass(ISL)=NBiomass(ISL)-TBMNA+THPNA+THSNA*(1.-ASP)+TLMNA+TLSNA-XNetNMineralize(ISL)
            NPassiveHumus(ISL)=NPassiveHumus(ISL)-THPNA+TBMNA*ABP+THSNA*ASP
            !Soil_Respiration(ISL)=.3*TLSLCA+A1CO2*(TLSLNCA+TLMCA)+ABCO2*TBMCA+ASCO2*THSCA+&
            !APCO2*THPCA
      END IF        
      CBiomass(ISL)=CBiomass(ISL)-TBMCA+X3
      CSlowHumus(ISL)=CSlowHumus(ISL)-THSCA+X1
      AD2=CBiomass(ISL)+CPassiveHumus(ISL)+CSlowHumus(ISL)+CMetabLitt(ISL)+CStructLitt(ISL)
      soilResp(ISL)=AD1-AD2
      NSlowHumus(ISL)=NSlowHumus(ISL)-THSNA+TBMNA*(1.-ABP)
      NStructLitt(ISL)=NStructLitt(ISL)-TLSNA
      NMetabLitt(ISL)=NMetabLitt(ISL)-TLMNA
      SMM(74,MO)=SMM(74,MO)+soilResp(ISL)
      VAR(74)=VAR(74)+soilResp(ISL)
      SMS(8,ISL)=SMS(8,ISL)+soilResp(ISL)
      Tot_SoilRes=Tot_SoilRes+soilResp(ISL)      ! Total respiration of CO2
      cropResidu(ISL)=.001*(structLitt(ISL)+metabLitt(ISL))
      HumusMine_RateCons(ISL)=XHN-NSlowHumus(ISL)-NPassiveHumus(ISL)
      DF=AD1-AD2-soilResp(ISL)
      SDF=SDF+DF
      IF(ABS(DF)>.01)WRITE(KW(1),3)IY,MO,DayOfMon,ISL,AD1,AD2,soilResp(ISL),DF,SDF
    3 FORMAT(1X,'NCNMI',4I4,10E16.6)      
END SUBROUTINE Phoenix

! 3. ---------------------------------- CENTURY model -----------------------------------------------

SUBROUTINE CENTURY(Z5)
      !     EPIC0810
      !     THIS SUBPROGRAM SIMULATES MINERALIZATION AND IMMOBILIZATION OF N AND C USING EQUATIONS TAKEN FROM CENTURY.
      !     Ref: APEX-doc P47
      IMPLICIT NONE
  
      ! local variables

      REAL:: XZ, RLR,AD1,XHN,X1, X3, O2_Factor,Z5, Y1,Y2,Y3,OXZ, OXK,APCO2, ASCO2, RBM, HPNC, XBM, BMNC, &
      HSNC, ABCO2, A1CO2, RLM, RLS, ABP, ASP, TLSCP, TLSLCP, TLSLNCP, TLSNP, TLMCP, TLMNP, &
      TBMCP, TBMNP, THSCP, THSNP, THPCP, THPNP, A1, ASX, APX, PN1, PN2, PN3, PN5, PN6, PN7,&
      PN8, PN9, SUM, CPN1, CPN2, CPN3, CPN4, CPN5, WMIN, DMDN, THSCA, THSNA, THPCA, THPNA, &
      TLSCA, TLSLNCA, TLSLCA, TLSNA, TLMCA, TLMNA, TBMCA, TBMNA, DF1, DF2, DF3, DF4, Bio_contrl_fact, &
      DF5, DF6, ADD, ADF1, ADF2, ADF3, ADF4, ADF5, TOT, XX, AD2 

      XZ=NO3_N_Soil(ISL)+NH3_Weight(ISL)
      ! AD1=NStructLitt(ISL)+NMetabLitt(ISL)+NBiomass(ISL)+NSlowHumus(ISL)+NPassiveHumus(ISL)+XZ
      RLR=MIN(.8,lgStructLitt(ISL)/(structLitt(ISL)+1.E-5))   ! ligin /sturctural litter 
      XHN=NSlowHumus(ISL)+NPassiveHumus(ISL)
      AD1=CStructLitt(ISL)+CMetabLitt(ISL)+CBiomass(ISL)+CSlowHumus(ISL)+CPassiveHumus(ISL)
	! ------------- O2 factor for microbial process -------------
      !IF(O2_meth_flag>0)THEN
      Y1=.1*SOC(ISL)/WT(ISL)
      Y2=2.11+.0375*clayFrac(ISL)
      Y3=100.*Y1/Y2
      OXK=Y3/(Y3+EXP(S_Curve(24,1)-S_Curve(24,2)*Y3))                  ! Oxygen content simulated by C and clay
	!ELSE            
      OXZ=1.-modelPara(53)*Z5/(Z5+EXP(S_Curve(20,1)-S_Curve(20,2)*Z5))! Oxygen content simulated by soil depth (Eq. 179a)
      !END IF
      O2_Factor=modelPara(30)*OXK+(1.-modelPara(30))*OXZ
      ! -------------- Tillage factor for microbial process -----------
      IF(bulkDensity(ISL)>postTillBulkDensity(ISL))THEN
            X1=MIN(20.,EXP(modelPara(52)*(bulkDensity(ISL)-postTillBulkDensity(ISL))))! tillage effect on residue decay rate
      ELSE
            X1=1.
      END IF
      ! ----------------------- Overall factor ------------------------ 
      Bio_contrl_fact=MIN(10.,SQRT(Tem_Fact(ISL)*SoilWater_Fact(ISL))*modelPara(20)*O2_Factor*X1)  ! APEX-doc  Eq.179
      
      SMS(4,ISL)=SMS(4,ISL)+X1                       ! SMS(4, ISL) : Tillage factor
      SMS(11,ISL)=SMS(11,ISL)+1.                     ! SMS(11,ISL) : times of running CENTURY ?
      SMS(1,ISL)=SMS(1,ISL)+SoilWater_Fact(ISL)      ! SMS(1, ISL) : Soil water factor
     
      APCO2=.55                                      ! Allocate rate: Passive humus  ----> CO2    
      ASCO2=.60                                      ! Allocate rate: Slow humus  ----> CO2 
      IF(ISL==LD1)THEN           ! Top layer
            IF(cropCode(JD)==plantCategoryCode(1))THEN  ! warm legimu crop 
                  Bio_contrl_fact=1.
                  RLM=.2                                 ! transformation rate of metabolic litter  
                  RLS=.2                                 ! transformation rate of structure litter
            ELSE                                       ! other crops
                  Bio_contrl_fact=Bio_contrl_fact*modelPara(51)
                  RLM=.0405
                  RLS=.0107
            END IF
            ABCO2=.55                                ! Biomass ---> CO2
            A1CO2=.55                                ! non-lignin of structural litter ---> CO2
            RBM=.0164
            HPNC=.1
            XBM=1.
            ! COMPUTE N/C RATIOS
            X1=.1*(NMetabLitt(LD1)+NStructLitt(LD1))/(cropResidu(LD1)+1.E-5)  ! N content
            IF(X1>2.)THEN
                  BMNC=.1                             ! N/C ratio of biomass formed from surface litter  
            ELSE
                  IF(X1>.01)THEN
                        BMNC=1./(20.05-5.0251*X1)       ! APEX-doc Eq. 192
                  ELSE  
                        BMNC=.05
                  END IF    
            END IF    
            HSNC=BMNC/(5.*BMNC+1.)                  ! N/C ratio of slow humus formed from surface microbes
      ELSE                            ! other layers
            ABCO2=.17+.0068*sandFrac(ISL)
            A1CO2=.55 
            RBM=.02
            RLM=.0507
            RLS=.0132
            XBM=.25+.0075*sandFrac(ISL)
            X1=1000.*(NH3_Weight(ISL)+NO3_N_Soil(ISL))/WT(ISL)
            IF(X1>7.15)THEN
                  BMNC=.33
                  HSNC=.083
                  HPNC=.143
            ELSE
                  BMNC=1./(15.-1.678*X1)
                  HSNC=1./(20.-1.119*X1)    ! N/C ratio of slow humus 
                  HPNC=1./(10.-.42*X1)      ! N/C ratio of passive humus
            END IF
      END IF    
      ABP=.003+.00032*clayFrac(ISL)                       ! Biomass ----> passive humus
      SMS(3,ISL)=SMS(3,ISL)+Bio_contrl_fact                ! SMS(3, ISL) : Biological control factor
      ASP=MAX(.001,modelPara(45)-.00009*clayFrac(ISL))   ! Slow humus ----> passive humus
      !     POTENTIAL TRANSFORMATIONS STRUCTURAL LITTER
      X1=RLS*Bio_contrl_fact*EXP(-3.*RLR)                  
      TLSCP=X1*CStructLitt(ISL)                         ! trans of C in structural litter
      TLSLCP=TLSCP*RLR                                     ! trans of C in lignin of structural litter
      TLSLNCP=TLSCP*(1.-RLR)                               ! trans of C in non-lignin compoents of structural litter
      TLSNP=X1*NStructLitt(ISL)                        ! trans of N in non-lignin components of structural litter 
      !     POTENTIAL TRANSFORMATIONS METABOLIC LITTER
      X1=RLM*Bio_contrl_fact
      TLMCP=CMetabLitt(ISL)*X1
      TLMNP=NMetabLitt(ISL)*X1
      !     POTENTIAL TRANSFORMATIONS MICROBIAL BIOMASS
      X1=RBM*Bio_contrl_fact*XBM
      TBMCP=CBiomass(ISL)*X1
      TBMNP=NBiomass(ISL)*X1
      !     POTENTIAL TRANSFORMATIONS SLOW HUMUS
      X1=slowHumusTranRate*Bio_contrl_fact
      THSCP=CSlowHumus(ISL)*X1
      THSNP=NSlowHumus(ISL)*X1
      !     POTENTIAL TRANSFORMATIONS PASSIVE HUMUS
      X1=Bio_contrl_fact*passivHumusTranRate
      THPCP=CPassiveHumus(ISL)*X1
      THPNP=NPassiveHumus(ISL)*X1
      !     ESTIMATE N DEMAND
      A1=1.-A1CO2
      ASX=1.-ASCO2-ASP
      APX=1.-APCO2
      PN1=TLSLNCP*A1*BMNC           ! non-lignin of structural litter ---> Biomass
      PN2=.7*TLSLCP*HSNC            ! lignin of structural litter ---> Slow humus
      PN3=TLMCP*A1*BMNC             ! Metabolic litter ---> Biomass
      PN5=TBMCP*ABP*HPNC
      PN6=TBMCP*(1.-ABP-ABCO2)*HSNC ! Biomass ---> Slow humus
      PN7=THSCP*ASX*BMNC            ! Slow humus ---> Biomass 
      PN8=THSCP*ASP*HPNC            ! Slow humus ---> Passive humus
      PN9=THPCP*APX*BMNC            ! Passive ---> Biomass 
      ! COMPARE SUPPLY AND DEMAND FOR N
      SUM=0.
      CPN1=0.
      CPN2=0.
      CPN3=0.
      CPN4=0.
      CPN5=0.
      X1=PN1+PN2
      IF(TLSNP<X1)THEN
            CPN1=X1-TLSNP
      ELSE
            SUM=SUM+TLSNP-X1
      END IF
      IF(TLMNP<PN3)THEN
            CPN2=PN3-TLMNP
      ELSE
            SUM=SUM+TLMNP-PN3
      END IF
      X1=PN5+PN6
      IF(TBMNP<X1)THEN
            CPN3=X1-TBMNP
      ELSE
            SUM=SUM+TBMNP-X1
      END IF      
      X1=PN7+PN8
      IF(THSNP<X1)THEN
            CPN4=X1-THSNP
      ELSE
            SUM=SUM+THSNP-X1
      END IF
      IF(THPNP<PN9)THEN
            CPN5=PN9-THPNP
      ELSE
            SUM=SUM+THPNP-PN9
      END IF
      !     NH3_Weight(ISL)=NH3_Weight(ISL)+SUM
      WMIN=MAX(1.E-5,NO3_N_Soil(ISL)+SUM)
      DMDN=CPN1+CPN2+CPN3+CPN4+CPN5
      X3=1.
      !     REDUCE DEMAND IF SUPPLY LIMITS
      IF(WMIN<DMDN)X3=WMIN/DMDN
      SMS(5,ISL)=SMS(5,ISL)+X3
      !     ACTUAL TRANSFORMATIONS
      TLSCA=TLSCP*X3
      TLSLCA=TLSLCP*X3
      TLSLNCA=TLSLNCP*X3
      TLSNA=TLSNP*X3
      TLMCA=TLMCP*X3
      TLMNA=TLMNP*X3
      TBMCA=TBMCP*X3
      TBMNA=TBMNP*X3
      THSCA=THSCP*X3
      THSNA=THSNP*X3
      THPCA=THPCP*X3
      THPNA=THPNP*X3
      !     DMDN=DMDN*X3
      SGMN=SGMN+SUM
      XNetNMineralize(ISL)=SUM-DMDN
      !     UPDATE
      IF(XNetNMineralize(ISL)>0.)THEN
            X1=modelPara(96)*XNetNMineralize(ISL)
            NH3_Weight(ISL)=NH3_Weight(ISL)+X1
            NO3_N_Soil(ISL)=NO3_N_Soil(ISL)+XNetNMineralize(ISL)-X1
      ELSE
            X1=NO3_N_Soil(ISL)+XNetNMineralize(ISL)
            IF(X1>0.)THEN
                  NO3_N_Soil(ISL)=X1
            ELSE
                  XNetNMineralize(ISL)=-NO3_N_Soil(ISL)
                  NO3_N_Soil(ISL)=1.E-10    
            END IF   
      END IF    
      DF1=TLSNA
      DF2=TLMNA
      mineralizedN=mineralizedN+XNetNMineralize(ISL)
      SMS(9,ISL)=SMS(9,ISL)+XNetNMineralize(ISL)
      CStructLitt(ISL)=MAX(1.E-10,CStructLitt(ISL)-TLSCA)
      CLgStructLitt(ISL)=MAX(1.E-10,CLgStructLitt(ISL)-TLSLCA)
      NLgStructLitt(ISL)=MAX(1.E-10,CStructLitt(ISL)-CLgStructLitt(ISL))
      CMetabLitt(ISL)=MAX(1.E-10,CMetabLitt(ISL)-TLMCA)
      metabLitt(ISL)=MAX(1.E-10,metabLitt(ISL)-TLMCA/.42)
      lgStructLitt(ISL)=MAX(1.E-10,lgStructLitt(ISL)-TLSLCA/.42)
      structLitt(ISL)=MAX(1.E-10,structLitt(ISL)-TLSCA/.42)
      X3=APX*THPCA+ASX*THSCA+A1*(TLMCA+TLSLNCA)
      CBiomass(ISL)=MAX(1.E-10,CBiomass(ISL)-TBMCA+X3)
      DF3=TBMNA-BMNC*X3
      X1=.7*TLSLCA+TBMCA*(1.-ABP-ABCO2)
      CSlowHumus(ISL)=MAX(1.E-5,CSlowHumus(ISL)-THSCA+X1)
      DF4=THSNA-HSNC*X1
      X1=THSCA*ASP+TBMCA*ABP
      CPassiveHumus(ISL)=MAX(1.E-5,CPassiveHumus(ISL)-THPCA+X1)
      DF5=THPNA-HPNC*X1
      DF6=XZ-NO3_N_Soil(ISL)-NH3_Weight(ISL)
      SMS(10,ISL)=SMS(10,ISL)-DF6
      ADD=DF1+DF2+DF3+DF4+DF5+DF6
      ADF1=ABS(DF1)
      ADF2=ABS(DF2)
      ADF3=ABS(DF3)
      ADF4=ABS(DF4)
      ADF5=ABS(DF5)
      TOT=ADF1+ADF2+ADF3+ADF4+ADF5
      XX=ADD/(TOT+1.E-10)
      NStructLitt(ISL)=MAX(.001,NStructLitt(ISL)-DF1+XX*ADF1)
      NMetabLitt(ISL)=MAX(.001,NMetabLitt(ISL)-DF2+XX*ADF2)
      NBiomass(ISL)=NBiomass(ISL)-DF3+XX*ADF3
      NSlowHumus(ISL)=NSlowHumus(ISL)-DF4+XX*ADF4
      NPassiveHumus(ISL)=NPassiveHumus(ISL)-DF5+XX*ADF5
      !soilResp(ISL)=MAX(0.,.3*TLSLCA+A1CO2*(TLSLNCA+TLMCA)+ABCO2*TBMCA+&
      !ASCO2*THSCA+APCO2*THPCA)
      AD2=CStructLitt(ISL)+CMetabLitt(ISL)+CBiomass(ISL)+CSlowHumus(ISL)+CPassiveHumus(ISL)
      soilResp(ISL)=MAX(0.,AD1-AD2)
      SMM(74,MO)=SMM(74,MO)+soilResp(ISL)
      SMS(8,ISL)=SMS(8,ISL)+soilResp(ISL)
      Tot_SoilRes=Tot_SoilRes+soilResp(ISL)      
      VAR(74)=VAR(74)+soilResp(ISL)
      cropResidu(ISL)=.001*(structLitt(ISL)+metabLitt(ISL))
      HumusMine_RateCons(ISL)=XHN-NSlowHumus(ISL)-NPassiveHumus(ISL)
 
END SUBROUTINE CENTURY 
 
! 4. -------------------------Denitrification and N2O loss ------------------------
SUBROUTINE DenitriN1 
      !     EPIC1102
      !     THIS SUBPROGRAM DEVELOPED BY ARMEN KEMANIAN ESTIMATES DAILY DENITRIFICATION AND N2O LOSSES OF SOIL NO3. 
      !    process: NO3+ ------> N2 and N2O
 
      IMPLICIT NONE
 
      ! local variables
      REAL:: AIRV, X1, X2, X3, ONO3_C, D_F, DNITMX

      weightDenitriN(ISL)=0.
      denitriedN2O=0.
      !     COMPUTE WATER FACTOR
      AIRV=MAX(0.,(Porosity(ISL)-soilWater(ISL))/layerThick)
      X1=0.90+.001*clayFrac(ISL)
      X2=(1.0001-AIRV)/X1
      IF(X2<.8)RETURN
      H2OF(ISL)=1./(1.+X2**(-60))
      !     COMPUTE NITRATE FACTOR
      ONO3_C=MAX(1.E-5,1000.*NO3_N_Soil(ISL)/WT(ISL))
      WNO3F(ISL)=ONO3_C/(ONO3_C+60.)
      !     COMPUTE RESPIRATION FACTOR
      X3=1000.*soilResp(ISL)/WT(ISL)
      CBNF(ISL)=MIN(1.,X3/50.)
      D_F=WNO3F(ISL)*H2OF(ISL)*CBNF(ISL)
      DNITMX=32.
      !weightDenitriN=MIN(modelPara(30),D_F*DNITMX*WT(ISL)/1000.)
      weightDenitriN(ISL)=D_F*DNITMX*WT(ISL)/1000.
      IF(weightDenitriN(ISL)>NO3_N_Soil(ISL))weightDenitriN(ISL)=NO3_N_Soil(ISL)
      !compute N2O as a fraction of weightDenitriN
      denitriedN2O=WNO3F(ISL)*(1.-SQRT(H2OF(ISL)))*(1.-CBNF(ISL)**.25)*weightDenitriN(ISL)
      NO3_N_Soil(ISL)=NO3_N_Soil(ISL)-weightDenitriN(ISL)
END SUBROUTINE DenitriN1
 
! 5 ------------------ NO3 Loss by Denitrification -------------------------------------------------------
SUBROUTINE DenitriN2
      !     EPIC1102
      !     THIS SUBPROGRAM ESTIMATES DAILY LOSS OF NO3 BY DENITRIFICATION.
      !      process: NO3+ ------> N2 and N2O
      IMPLICIT NONE
 
      ! local variables 
      REAL:: AIRV, X1, X2, F_H2O, F 
 
      AIRV=MAX(0.,(Porosity(ISL)-soilWater(ISL))/layerThick)
      X1=0.90+.001*clayFrac(ISL)
      X2=(1.0001-AIRV)/X1
      IF(X2<.8)RETURN
      F_H2O=1./(1.+X2**(-60))
      F=SQRT(Tem_Fact(ISL)*SoilWater_Fact(ISL))*F_H2O
      weightDenitriN(ISL)=NO3_N_Soil(ISL)*F
      IF(weightDenitriN(ISL)>NO3_N_Soil(ISL))weightDenitriN(ISL)=NO3_N_Soil(ISL)
      NO3_N_Soil(ISL)=NO3_N_Soil(ISL)-weightDenitriN(ISL)
      DN2=modelPara(80)*weightDenitriN(ISL)          ! para(80): N2 fraction of total denitrification
      denitriedN2O=weightDenitriN(ISL)-DN2
 
END SUBROUTINE DenitriN2              

!7 ------------------- predict Organic N and Humus loss ------------------
SUBROUTINE OrgN_Loss                       
      !     EPIC1102
      !     THIS SUBPROGRAM PREDICTS DAILY ORGANIC N AND HUMUS LOSS, GIVEN
      !     SOIL LOSS AND ENRICHMENT RATIO.                        
      IMPLICIT NONE
 
      ! local variables
      REAL:: TOT, X1

      TOT=NPassiveHumus(LD1)+NSlowHumus(LD1)+NBiomass(LD1)+NMetabLitt(LD1)+NStructLitt(LD1)
      sedimentLossN=enrichRatioFinal*TOT
      X1=1.-enrichRatioFinal
      NBiomass(LD1)=NBiomass(LD1)*X1
      NSlowHumus(LD1)=NSlowHumus(LD1)*X1
      NPassiveHumus(LD1)=NPassiveHumus(LD1)*X1
      NStructLitt(LD1)=NStructLitt(LD1)*X1
      NMetabLitt(LD1)=NMetabLitt(LD1)*X1

END SUBROUTINE OrgN_Loss  

! 8 ----------------------- predict daily C loss -----------------------
SUBROUTINE C_Loss 
      !     EPIC1102
      !     THIS SUBPROGRAM PREDICTS DAILY C LOSS, GIVEN SOIL LOSS AND
      !     ENRICHMENT RATIO.
      IMPLICIT NONE

      ! local variables
      INTEGER:: I,L
      REAL:: X1, X2, X3, XX, XK, Y1, Y4, QBC, VBC, YBC, YOC, TOT, DK, V, VL, CO, B1, QMAX, &
       C, A, XF, B , Bio_contrl_fact    

      Y1=CBiomass(LD1)
      Y4=percolateFlow(LD1)      ! This is percolation
      QBC=0.
      VBC=0.
      YBC=0.
      YOC=0.
      TOT=CPassiveHumus(LD1)+CSlowHumus(LD1)+CMetabLitt(LD1)+CStructLitt(LD1)
      X1=1.-enrichRatioFinal
      YOC=enrichRatioFinal*TOT
      CSlowHumus(LD1)=CSlowHumus(LD1)*X1
      CPassiveHumus(LD1)=CPassiveHumus(LD1)*X1
      structLitt(LD1)=structLitt(LD1)*X1
      metabLitt(LD1)=metabLitt(LD1)*X1
      lgStructLitt(LD1)=lgStructLitt(LD1)*X1
      CStructLitt(LD1)=CStructLitt(LD1)*X1
      CMetabLitt(LD1)=CMetabLitt(LD1)*X1
      CLgStructLitt(LD1)=CLgStructLitt(LD1)*X1
      NLgStructLitt(LD1)=CStructLitt(LD1)-CLgStructLitt(LD1)    ! ????
      IF(Y1>0.)THEN
            DK=.0001*modelPara(21)*SOC(LD1)                ! DK: m3 t   APEX-doc Eq. 257
            X1=Porosity(LD1)-Wilt_Point(LD1)
            XX=X1+DK
            V=Runoff+Y4             ! total runoff
            X3=0.
            IF(V>0.)THEN
                  X3=Y1*(1.-EXP(-V/XX))
                  CO=X3/(Y4+modelPara(44)*Runoff)     ! para(44)  : ratio of soluble C in runoff to percolate
                  Bio_contrl_fact=modelPara(44)*CO
                  VBC=CO*Y4
                  Y1=Y1-X3
                  QBC=Bio_contrl_fact*Runoff
            END IF    
            ! COMPUTE CBiomass LOSS WITH SEDIMENT
            IF(enrichRatioFinal>0.)THEN
                  Bio_contrl_fact=DK*Y1/XX
                  YBC=enrichRatioFinal*Bio_contrl_fact
            END IF    
      END IF
      CBiomass(LD1)=Y1-YBC
      XX=0.
      DO L=2,Actual_SoilLayers
            ISL=Layer_ID(L)
            Y1=CBiomass(ISL)+VBC
            VBC=0.
            IF(Y1>=.01)THEN
                  V=percolateFlow(ISL)
                  IF(V>0.)THEN
                        IF(soilHorizon(ISL)=='       B'.OR.soilHorizon(ISL)=='       C')THEN
                              DO I=1,3
                                    IF(soilOrder==soilOrders(I))EXIT
                              END DO  
                              IF(I<=3)THEN
                                    SELECT CASE(I)
                                          CASE(1)
                                                B1=2.662
                                                QMAX=10**(B1+.572*LOG10(clayFrac(ISL))-.0602*PH(ISL))
                                          CASE(2)
                                                B1=-206.452
                                                X2=1000.*SOC(ISL)/WT(ISL)
                                                QMAX=B1+.127*X2-46.88*clayFrac(ISL)
                                          CASE(3)
                                                B1=2.141
                                                X1=10.*ironCon(ISL)
                                                QMAX=10**(B1+.403*LOG10(clayFrac(ISL))+.439*LOG10(X1))
                                    END SELECT
                                    XK=10**(-.386-.184*PH(ISL))
                                    VL=soilWater(ISL)
                                    X2=.1*CBiomass(ISL)/(Z(ISL)-XX)
                                    XX=Z(ISL)
                                    C=X2*VL+B1
                                    B=XK*(C-QMAX)-VL
                                    A=-VL*XK
                                    XF=(-B-SQRT(B*B-4.*A*C))/(2.*A)
                                    !RE=XK*QMAX*XF/(1.+XK*XF)-B1
                                    !SBC=.01*RE*V
                                    VBC=.01*XF*V
                                    !Y1=Y1-SBC
                                    !RTO=CSlowHumus(ISL)/(CSlowHumus(ISL)+CPassiveHumus(ISL))
                                    !X1=RTO*SBC
                                    !X2=SBC-X1
                                    !CSlowHumus(ISL)=CSlowHumus(ISL)+X1
                                    !CPassiveHumus(ISL)=CPassiveHumus(ISL)+X2
                              ELSE                                              
                                    VBC=Y1*(1.-EXP(-V/(Porosity(ISL)-Wilt_Point(ISL)+.0001*modelPara(21)*SOC(ISL))))
                              END IF    
                        ELSE                              
                              VBC=Y1*(1.-EXP(-V/(Porosity(ISL)-Wilt_Point(ISL)+.0001*modelPara(21)*SOC(ISL))))
                        END IF    
                  END IF
            END IF
            CBiomass(ISL)=Y1-VBC
      END DO
      SMM(75,MO)=SMM(75,MO)+VBC                ! SMM(75,MO), VAR(75):  Soluble C leached
      VAR(75)=VBC
      SMM(76,MO)=SMM(76,MO)+QBC                ! SMM(76,MO), VAR(76):  Carbon in runoff
      VAR(76)=QBC
      YOC=YOC+YBC
      SMM(77,MO)=SMM(77,MO)+YOC                ! SMM(77,MO), VAR(77):  Carbon Loss with sediment
      VAR(77)=YOC
END SUBROUTINE C_Loss

! --------------------------- 9. NO3 upward caused by soil EP ----------------------------------
SUBROUTINE NO3_EP 
      !     EPIC1102
      !     THIS SUBPROGRAM ESTIMATES UPWARD NO3 MOVEMENT CAUSED BY SOIL EVAPO
      IMPLICIT NONE
 
      ! local variables
      INTEGER:: J
      REAL:: TOT, X1, X2, XX

      IF(NEV==1)RETURN
      TOT=0.
      DO J=NEV,2,-1
            ISL=Layer_ID(J)
            X1=NO3_N_Soil(ISL)-.001*WT(ISL)*modelPara(27)
            IF(X1<=.01)CYCLE
            X2=1.-EXP(-modelPara(62)*soilElement(ISL)/(fracNO3Leach(ISL)*Porosity(ISL)))
            XX=X1*X2
            TOT=TOT+XX
            NO3_N_Soil(ISL)=NO3_N_Soil(ISL)-XX
      END DO
      NO3_N_Soil(LD1)=NO3_N_Soil(LD1)+TOT
END SUBROUTINE NO3_EP

! -------------------------------------10. P upward by soil EP -------------------------------------------
SUBROUTINE PUpward_EP 
      !     EPIC1102
      !     THIS SUBPROGRAM ESTIMATES UPWARD SOLUBLE P MOVEMENT CAUSED BY SOIL
      !     EVAPORATION.
      IMPLICIT NONE
      
      ! local variables
      INTEGER:: J
      REAL:: TOT, XX

      IF(NEV==1)RETURN
      TOT=0.
      DO J=NEV,2,-1
            ISL=Layer_ID(J)
            IF(labileP(ISL)<.001)CYCLE
            XX=labileP(ISL)*MIN(.75, modelPara(43)*soilElement(ISL)/WT(ISL))
            TOT=TOT+XX
            labileP(ISL)=labileP(ISL)-XX
      END DO
      labileP(LD1)=labileP(LD1)+TOT
END SUBROUTINE PUpward_EP

! -------------------------------11. Salt movement by EP ---------------------------------------
SUBROUTINE SaltMove_EP 
      !     EPIC1102
      !     THIS SUBPROGRAM ESTIMATES UPWARD SALT MOVEMENT CAUSED BY SOIL
      !     EVAPORATION
      IMPLICIT NONE
 
      ! local variables
      INTEGER:: J, L
      REAL:: SUM, X1, XX

      IF(NEV==1)RETURN
      J=NEV
      SUM=0.
      DO J=NEV,2,-1
            L=Layer_ID(J)
            X1=saltWeight(L)
            IF(X1<=1.E-5)CYCLE
            XX=MIN(.05*X1,soilElement(L)*X1/(soilWater(L)+soilElement(L)))
            SUM=SUM+XX
            saltWeight(L)=saltWeight(L)-XX
      END DO
      saltWeight(LD1)=saltWeight(LD1)+SUM
END SUBROUTINE SaltMove_EP

!---------------------------- 12. calculate N, P balance --------------------------------
SUBROUTINE Cal_NP_Balance(beginTotN,rainfallN,orgSedimentLostN,surfRunoffNO3,subSurfFlowNO3,vFlowNO3,surfRunoffNO2, &
                          subSurfFlowNO2,vFlowNO2,surfRunoffN2O, subSurfFlowN2O,vFlowN2O, denitriedN,&
                          orgFertN,yieldN, votalizedN, fertNO3,fertNH3,fixedN,burnN,diffusedN2O,finalTotN,KBL)
      !     EPIC1102
      !     THIS SUBPROGRAM COMPUTES THE N & P BALANCES AT THE END OF THE SIMULATION
      !     Using N for example, it can be used for P and K               
      IMPLICIT NONE
      ! global variables
      REAL, INTENT(IN):: beginTotN,rainfallN,orgSedimentLostN,surfRunoffNO3,subSurfFlowNO3,vFlowNO3,surfRunoffNO2, &
                         subSurfFlowNO2,vFlowNO2,surfRunoffN2O, subSurfFlowN2O,vFlowN2O, denitriedN,&
                         orgFertN,yieldN, votalizedN, fertNO3,fertNH3,fixedN,burnN,diffusedN2O,finalTotN
      INTEGER, INTENT(IN):: KBL
      ! local variables
      REAL:: DF, Percent 
	  
      DF=beginTotN + rainfallN + fertNO3 + fertNH3 + orgFertN + fixedN + diffusedN2O &
         -surfRunoffNO3-subSurfFlowNO3-vFlowNO3-surfRunoffNO2-subSurfFlowNO2-vFlowNO2-surfRunoffN2O &
         -subSurfFlowN2O-vFlowN2O-denitriedN-votalizedN-orgSedimentLostN-burnN-yieldN-finalTotN
      
      Percent=100.*DF/finalTotN
      
      SELECT CASE(KBL)
            CASE(1)
                  WRITE(KW(1),'(/T10,A)')'N BALANCE'
                  WRITE(KW(1),1)Percent,DF,beginTotN,rainfallN,orgSedimentLostN,surfRunoffNO3,subSurfFlowNO3,&
                                vFlowNO3,surfRunoffNO2,subSurfFlowNO2,vFlowNO2,surfRunoffN2O,subSurfFlowN2O,vFlowN2O,&
                                denitriedN, yieldN,votalizedN,fertNO3,fertNH3,fixedN,orgFertN,burnN,diffusedN2O,finalTotN
            CASE(2)
                  WRITE(KW(1),'(/T10,A)')'P BALANCE'
                  WRITE(KW(1),2)Percent,DF,beginTotN,orgSedimentLostN,surfRunoffNO3,vFlowNO3,yieldN,fertNO3,orgFertN,finalTotN
            CASE(3)
                  WRITE(KW(1),'(/T10,A)')'K BALANCE'
                  WRITE(KW(1),3)Percent,DF,beginTotN,orgSedimentLostN,surfRunoffNO3,subSurfFlowNO3,vFlowNO3,yieldN,fertNO3,&
                                finalTotN
      END SELECT
 
    1 FORMAT(5X,'PER =',E13.6,2X,'DF  =',E13.6,2X,'BTOT=',E13.6,2X, 'PCP =',E13.6,2X,'Y   =',E13.6,2X,'QNO3=',E13.6&
            /5X,'SNO3=',E13.6,2X,'verticalFlowNO3=',E13.6,2X,'QNO2=',E13.6,2X,'SNO2=',E13.6,2X,'VerticalFlow_NO2=',&
            E13.6,2X,'QN2O=',E13.6/5X,'totN2O=',E13.6,2X,'verticalFlowN2O=',E13.6,2X,'DNIT=',E13.6,2X,'YLD =',E13.6,&
            2X,'VOL =',E13.6,2X,'FNO3=',E13.6/5X,'FNH3=',E13.6,2X,'FIX =',E13.6,2X,'FORG=',E13.6,2X,'BURN=',&
            E13.6,2X,'FN2O=',E13.6,2X,'FTOT=',E13.6)
    2 FORMAT(5X,'PER =',E13.6,2X,'DF  =',E13.6,2X,'BTOT=',E13.6,2X,'Y   =',E13.6,2X,'Q   =',E13.6,2X,'PRK =',E13.6/&
            5X,'YLD =',E13.6,2X,'FPML=',E13.6,2X,'FPO =',E13.6,2X,'ETOT=',E13.6)
    3 FORMAT(5X,'PER =',E13.6,2X,'DF  =',E13.6,2X,'BTOT=',E13.6,2X,'Y   =',E13.6,2X,'Q   =',E13.6,2X,'subsurfLaterFlow =',E13.6&
            /5X,'percolateFlow =',E13.6,2X,'YLD =',E13.6,2X,'FKM =',E13.6,2X,'ETOT=',E13.6)
 
END SUBROUTINE Cal_NP_Balance

! ------------------------------ 13. Calculate C balance -------------------------
SUBROUTINE Cal_C_Balance(beginTotC, orgCLostSediment, solubleCLeach, runoffC, respirC, fertOrgC, nppC, yieldC, burnC, finalTotC)
      !     EPIC1102
      !     THIS SUBPROGRAM COMPUTES THE C BALANCE AT THE END OF THE
      !     SIMULATION
      IMPLICIT NONE
      ! Argument lists
      REAL, INTENT(IN):: orgCLostSediment, solubleCLeach, runoffC, respirC, fertOrgC, nppC, yieldC, burnC, finalTotC 
      REAL, INTENT(INOUT)::beginTotC 
      ! local variables
      REAL:: DF, Percent

      WRITE(KW(1),'(/T10,A)')'C BALANCE'
      
      DF=beginTotC + fertOrgC + nppC -orgCLostSediment-solubleCLeach-runoffC-respirC-burnC-yieldC-finalTotC
      
      Percent=100.*DF/finalTotC
      
      WRITE(KW(1),1)Percent,DF,beginTotC,orgCLostSediment,solubleCLeach,runoffC,respirC,fertOrgC,nppC,yieldC,burnC,finalTotC
      
      beginTotC= finalTotC
 
    1 FORMAT(5X,'PER =',E13.6,2X,'DF  =',E13.6,2X,'BTOT=',E13.6,2X, 'Y   =',E13.6,2X,'PRK =',E13.6,2X,'Q   =',E13.6/&
            5X,'Soil_Respiration=',E13.6,2X,'TFOC=',E13.6,2X,'NPPC=', E13.6,2X,'YLDC=',E13.6,2X,'BURN=',E13.6,&
            2X,'FTOT=',E13.6)
 
END SUBROUTINE Cal_C_Balance

! ------------------------- cal SALT balance ---------------------------
SUBROUTINE Cal_Salt_Balance(SLI,SLF,SLL,SLS,SLQ,SLB,SLE)
      !     EPIC1102
      !     THIS SUBPROGRAM IS THE SALT BALANCE
      IMPLICIT NONE
      ! local variables
      REAL, INTENT(IN):: SLI,SLF,SLL,SLS,SLQ,SLB,SLE
      REAL:: DF, PER
      WRITE(KW(1),2)
      DF=SLB+SLI+SLF-SLL-SLS-SLQ-SLE
      PER=100.*DF/(SLE+.0001)
      WRITE(KW(1),1)PER,DF,SLB,SLI,SLF,SLL,SLS,SLQ,SLE
  
    1 FORMAT(5X,'PER =',E13.6,2X,'DF  =',E13.6,2X,'BTOT=',E13.6,2X,&
     &'IRR =',E13.6,2X,'FERT=',E13.6,2X,'percolateFlow =',E13.6/5X,'subsurfLaterFlow =',E13.6,&
     &2X,'Q   =',E13.6,2X,'FTOT=',E13.6)
    2 FORMAT(/T10,'SALT BALANCE')
END SUBROUTINE Cal_Salt_Balance

END MODULE Nutrient_Module