MODULE Initialize
USE PARM
USE MisceFun_Module
IMPLICIT NONE

CONTAINS
SUBROUTINE Setting_Var
      HEDC=(/"             HUI","             LAI","              RD", "Total_rootWeight","            BIOM",&
             "            PLTC","  AboveGnd_BioMS","            CPHT", "   standCropResi","     actualCropN",&
             "     actualCropP","     actualCropK","     waterStress", "              NS","              PS",&
             "              KS","              TS","              AS", "            SALT","minStress_Factor"/)

      HEDP=(/"            PAPL","            PSRO","      Pest_Leach", "            PSSF","            PSED",&
             "            PDGF","            PDGS","Pest_Intercepted", "    pestSolubity","            PDRN"/)
      HEDS=(/"            TNH3","            TNO3","            PLAB", "            TSOK","            SNOA",&
             "    SoilWater_RZ","   WaterTable_ht","GroundWater_Stor", "StndDedRes_Ognsm","        Crop_Res",&
             "      plawDepSOC","          totSOC","              LS", "              LM","             LSL",&
             "             LSC","             LMC","            LSLC", "            LSNC","             BMC",&
             "             HSC","             HPC","             LSN", "             LMN","             BMN",&
             "             HSN","             HPN","          totSON", "            SALT","            TNO2"/) 
                      
      SID=(/"   1","   2","   3","   4","   5","   6","   7","   8",&
            "   9","  10","  11","  12","  13","  14","  15"," TOT"/)
      soilOrders=(/"ALFISOLS ","MOLLISOLS","ULTISOLS "/)
     
      numEconVar=40;numConVar=4;numPrintVar=100;numDayVar=40;numMonVar=40;numAnnuVar=40
      
      NC=(/0,31,60,91,121,152,182,213,244,274,305,335,366/)    
      
      plantCategoryCode=(/1,2,3,4,5,6,7,8,9,10/)                                        ! Crop category code
	  ! 1 - Warm season annual legume; 
      ! 2 - Cold season annual legume
      ! 3 - Perennial legume
      ! 4 - Warm season annual
      ! 5 - Cold season annual
      ! 6 - Perennial
      ! 7 - Evergreen tree
      ! 8 - Deciduous tree
      ! 9 - Cotton
      ! 10- N-fixing tree
      operationCode=(/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,&  ! Operation code
      23,24,25,26,27/)
      
      KR=(/11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,&
      31,32,33,34,35,36,37,38,39,40/)
 
      WCS=(/25.4,50.8,76.2/)
 
END SUBROUTINE Setting_Var
    
! --------------------------- Initialize Random Values ----------------------------
SUBROUTINE Setting_Random
      
      IX= (/748932582,1985072130,1631331038,67377721,366304404,1094585182, &
            1767585417,1980520317,392682216, 64298628,250756106,1025663860, &
            186056398, 522237216, 213453332, 1651217741, 909094944, 2095891343, &
            203905359, 2001697019, 431442774/)
             
END SUBROUTINE Setting_Random
    
! --------------------------- Allocate variables -----------------------------------
SUBROUTINE Allocate_Parms
      ! 
      ! --------------------------- 1. MFT (20): fertilizer--------------------------
      ALLOCATE(fertName(MFT),KDF(MFT),FCEM(MFT),FCST(MFT),FK(MFT),FN(MFT),FNH3(MFT),FNO(MFT),&              
             FOC(MFT),FP(MFT),FPO(MFT),FSLT(MFT))
      ! --------------------------- 2. MSL (15): soil layers ------------
      ALLOCATE(soilHorizon(MSL),Layer_ID(MSL+1),LORG(MSL),ALS(MSL),labileP(MSL),bulkDensity(MSL),&
             dryBulkDensity(MSL),mineralBulkDensity(MSL),postTillBulkDensity(MSL), flowCoeffP(MSL), &
             CaCO3(MSL),CBNF(MSL),Tem_Fact(MSL),CEC(MSL),CEM(MSL),clayFrac(MSL),initialNO3(MSL),&
             HumusMine_RateCons(MSL),elecCondv(MSL),EQKE(MSL),EQKS(MSL),exchangeK(MSL),&
             ironCon(MSL), fixK(MSL),FreshOrgP_Residu(MSL),H2OF(MSL),lateralFlowCondv(MSL),PH(MSL),&
             stableMineralP(MSL),initialLabileP(MSL), activeMineralP(MSL),Porosity(MSL),sorbRatioP(MSL),&
             XNetNMineralize(MSL),rockFrac(MSL), cropResidu(MSL),sandFrac(MSL), satuateCondv(MSL),&
             soilElement(MSL), siltFrac(MSL),SLTP(MSL),totBase(MSL),solubleK(MSL), &
             subsurfLaterFlow(MSL), percolateFlow(MSL), soilWater(MSL),fracNO3Leach(MSL),soilTem(MSL),&
             SoilWater_Fact(MSL),wiltPoint(MSL),soilSupplyRateK(MSL),soilSupplyRateN(MSL),soilSupplyRateP(MSL),&
             VerticalFlow_NO2(MSL),verticalFlowNO3(MSL),NBiomass(MSL),weightDenitriN(MSL),CPassiveHumus(MSL),&
             NPassiveHumus(MSL),CSlowHumus(MSL),NSlowHumus(MSL),metabLitt(MSL),CMetabLitt(MSL), &
             NMetabLitt(MSL),structLitt(MSL),CStructLitt(MSL),lgStructLitt(MSL),CLgStructLitt(MSL),&
             NLgStructLitt(MSL), NStructLitt(MSL),WNO3F(MSL),SOC(MSL),SON(MSL),&
             SOP(MSL),saltWeight(MSL),WT(MSL),Z(MSL), NPT(6,MSL),pestInSoil(MPS,MSL),SMS(11,MSL),&
             SOIL(20,MSL),RWTX(MSL,MNC),SOL(23,MSL))
      ! ------------------------ 3. MSC (31): crop-related variables -------------------- 
      ALLOCATE(ACO2C(MSC),AFP(MSC),AN2OC(MSC),AO2C(MSC),gasCO2Con(MSC),gasN2OCon(MSC),&
             gasO2Con(MSC),CLCO2(MSC),CLN2O(MSC),CLO2(MSC),DCO2GEN(MSC),DN2G(MSC), &
             DN2OG(MSC),DO2CONS(MSC),DPRC(MSC),DPRN(MSC),DPRO(MSC),DRWX(MSC),&
             EAR(MSC),fieldCapacity(MSC),HKPC(MSC),HKPN(MSC),HKPO(MSC),&
             soilResp(MSC),RWTZ(MSC),Wilt_Point(MSC),SMEA(MSC),SMES(MSC),SOT(MSC),&
             SubSurfFlow_CO2L(MSC),SubSurfFlow_N2O(MSC),SubSurfFlow_O2L(MSC),TPOR(MSC),VerticalFlow_CO2(MSC),&
             VGA(MSC), VGN(MSC),verticalFlowN2O(MSC),VerticalFlow_O2(MSC),VFC(MSC),VWC(MSC),VWP(MSC),&
             CBiomass(MSC), WCO2G(MSC),Weight_CO2L(MSC),WN2O(MSC),WN2OG(MSC),Weight_N2OL(MSC),NO2Weight(MSC),&
             NH3_Weight(MSC),NO3_N_Soil(MSC),WO2G(MSC),Weight_O2L(MSC),ZC(MSC),RWT(MSC,MNC))
      ! --------------------------- 4. NSM(114)  all variables ----------------------------
      ALLOCATE(varName(NSM),SM(NSM),SMY(NSM),VAR(NSM),SMM(NSM,12))
      varName=(/" TMX"," TMN"," RAD","PRCP","SNOF","SNOM","WSPD","RHUM",&                     
      " VPD"," PET","  ET"," PEP","  EP","   Q","  CN"," SSF"," PRK",&               
      "QDRN","IRGA"," QIN","TLGE","TLGW","TLGQ","TLGF","LGIR","LGMI",&               
      "LGMO","  EI"," CVF","USLE","MUSL"," AOF","MUSS","MUST","RUS2",&               
      "RUSL","RUSC"," WK1","RHTT","RRUF","RGRF","  YW"," YON","QNO3",&               
      "SNO3","VNO3"," NMN"," GMN","  DN","NFIX","NITR","AVOL","DRNN",&               
      "  YP"," QAP"," MNP","PRKP","  ER"," FNO","FNO3","FNH3"," FPO",&               
      " FPL"," FSK"," FCO","LIME"," TMP","SW10","SLTI","SLTQ","SLTS",&               
      "SLTF","RSDC","RSPC","CLCH"," CQV"," YOC","YEFK"," QSK"," SSK",&               
      " VSK","SLTV","MUSI","IRDL"," HMN","RNAD","NIMO","FALF"," DN2",&               
      "RLSF"," REK","FULU","DN2O"," FO2","FCO2","CFEM","BURC","BURN",&
      "NPPC"," NEE","FN2O","SNO2","SN2O","VN2O","VNO2","QNO2","QN2O",&
      "UNO3","UNH3","RSSF","DPRK","RSFN","DPKN","DRNP"/)     
      ! --------------------------- 5. MPS (40) Pestiside-related variables----------------- 
      ALLOCATE(pestName(MPS),KDP(MPS),GWPS(MPS),Pest_C_emission(MPS),Pest_Cost(MPS),pestPlantIntercep(MPS),&
              pestHfLifePlant(MPS),pestHfLifeSoil(MPS), Coeff_OrgCarbon(MPS),Pest_Leach(MPS),&
              pestSolubity(MPS),WashOff_Frac(MPS),RSPS(MPS),Pest_Drainged(MPS), PVQ(MPS,90),PVY(MPS,90), &
              SMAP(10,MPS),SMYP(10,MPS),pestVar(12,MPS),SPQ(5,MPS), SPQC(5,MPS),Pest_Sediment(5,MPS), &
              APQ(5,MPS,MYR),APQC(5,MPS,MYR),APY(5,MPS,MYR),AQB(5,MPS,MYR), AYB(5,MPS,MYR),SMMP(12,MPS,12))
              ! MYR : number of simulation years
      ! ---------------------------- 6. MNC (30) number of plants?-----------------------
      ALLOCATE(Crop_Name(MNC),IAP(MNC),IHU(MNC),IYH(MNC),JE(MNC),JP(MNC),JPL(MNC), cropID(MNC), KG(MNC), &
              NCP(MNC),NCR(MNC),NHU(MNC),NHV(MNC),NYLN(MNC), ANA(MNC),AWC(MNC),areationThreshold(MNC),&
              fracYieldK(MNC),fracYieldN(MNC),fracYieldP(MNC),seedCost(MNC),FLF0(MNC),leafWeightFrac(MNC),FNMX(MNC), &
              HUF(MNC),pestDmagFactor(MNC),PSTM(MNC),RDF(MNC),seedRate(MNC),optimalTem(MNC),TFTN(MNC),TFTP(MNC),&
              waterFracYield(MNC),ACET(MNC),AJWA(MNC),AJHI(MNC),cropTolerantAl(MNC),BLYN(MNC),LfWidth(MNC),&
              carbonEmissionSeedWeight(MNC),CPHT(MNC),potentialBioIncrese(MNC),fracGrowSeasonLAIDecline(MNC),&
              totCropBio(MNC),maxLAI(MNC),DMLX(MNC),totCropBioX(MNC),EP(MNC),lintFracCottn(MNC),turnoutFracCottn(MNC),&
              germinateHeatUnit(MNC),maxStomaCond(MNC),HI(MNC),maxCropHeight(MNC), HU(MNC),HUI(MNC),PPL0(MNC),&
              seedYieldPrice(MNC),forageYieldPrice(MNC),PSTF(MNC),factorBioDecline(MNC),RD(MNC),maxRootDep(MNC),&
              Min_Stress_Factor(MNC),factorLAIDecline(MNC),totRootWeight(MNC),RWX(MNC),Current_LAI(MNC),&
              standCropResi(MNC),standDeadResiN(MNC),abvGroundBiom(MNC),AboveGround_Biomass0(MNC),&
              sumActuralEP(MNC),sumPotentialEP(MNC),TCAW(MNC),TCQV(MNC),TCRF(MNC),TCSO(MNC),TCST(MNC),TDM(MNC),&
              TETG(MNC),TFTK(MNC),plantMinT(MNC),THU(MNC),TIRL(MNC),TRA(MNC), TRD(MNC),TVAL(MNC),TVIR(MNC),&
              TYLC(MNC),TYLK(MNC),TYLN(MNC),TYLP(MNC), TYL1(MNC),TYL2(MNC),actualCropK(MNC),UK2(MNC),ULYN(MNC),&
              UNA(MNC),actualCropN(MNC),Optimal_N(MNC),UNMX(MNC),actualCropP(MNC),Optimal_P(MNC),UPMX(MNC),&
              VPDPara2(MNC),VPDThreshold(MNC),RUE(MNC),VPDPara(MNC),WCHT(MNC),WLV(MNC),lowerLimitHI(MNC),&
              WUE(MNC),XDLA0(MNC),XDLAI(MNC),XLAI(MNC),XMTU(MNC),YLD(MNC),potentialHeatUnit(MNC,150),&
              POP(MNC,150),PPLA(MNC,150),totPestControl(MNC,MPS),pointsLAIDevp(2,MNC),pointsFrostDamg(2,MNC),&
              coeffPopuCurve(2,MNC),pointsPlantPopulation(2,MNC),rootPartition(2,MNC),salinityThreshold(2,MNC),&
              CO2EffOnBio(2,MNC),IGMD(10,MNC),IHVD(10,MNC), IPLD(10,MNC),NPSF(10,MNC),STDA(4,MNC),uptakeParaK(4,MNC),&
              uptakeParaN(4,MNC),uptakeParaP(4,MNC),ligninFrac(3,MNC),windEroCrop(3,MNC), &
              cropVar(20,MNC),TotK_Frac(10,MNC), TotN_Frac(10,MNC),TotP_Frac(10,MNC),TPSF(10,MNC),&
              availWaterCrop(10,MNC),CQV(10,MNC),CRF(10,MNC),CSOF(10,MNC), CSTF(10,MNC),DMF(10,MNC), &
              ET_GrowingSeason(10,MNC),HIF(10,MNC),RWF(10,MNC),VIL(10,MNC), VIR(10,MNC),YLCF(10,MNC),YLD1(10,MNC),&
              YLD2(10,MNC),YLKF(10,MNC),YLNF(10,MNC),YLPF(10,MNC), TSFC(7,MNC), growStressFactor(7,MNC), &
              SFMO(7,MNC), SMMC(20,MNC,12), Growing_Stress(7,10,MNC) )
     
      ! ----------------------------- 7. MRO (150) : Crop Rotations?------------------ 
      ALLOCATE(NBC(MRO),numOfPestArr(MRO),numOfTill(MRO),tillDOY(MRO,MNT),LFT(MRO,MNT),LT(MRO,MNT), &
              LYR(MRO,MNT),LPC(MRO,MNT),JH(MRO,MNT),LY(MRO,MNC),NBCX(MRO,MNC),pestKillEffi(MRO,MNT),&
              pestRate(MRO,MNT),fertCFactirArr(MRO,MNT),CND(MRO,MNT),fracHeatUnit(MRO,MNT),HWC(MRO,MNT),&
              QIR(MRO,MNT),baseStockRateArr(MRO,MNT),irrTriggerArr(MRO,MNT),TLMA(MRO,MNT),irrVolumeArr(MRO,MNT),&
              WFA(MRO,MNT))
     
      ! ------------------------------ 8. MNT (100) ----------------- 
      ALLOCATE(equipmentName(MNT), KOMP(MNT),NBE(MNT),NBT(MNT), ICUS(MNT),currentOps(MNT),IHT(MNT), &
              COOP(MNT),COTL(MNT),furrowHeight(MNT),furrowInterval(MNT),EFM(MNT),mixEfficiency(MNT), &              
              PlantRedu_frac(MNT),soilCompactFrac(MNT),Fuel_Use(MNT),harvestEffi(MNT),HMO(MNT), &
              overRideHI(MNT), ridgeHeight(MNT), ridgeInterval(MNT),tillRoughness(MNT),carbonEmiss(MNT), &
              tillDepth(MNT))
    
      ! ------------------------------- 9. NGF (38) ouput file control----------------- 
      ALLOCATE(KFL(NGF))    ! the additional file for testing yield                            
      ALLOCATE(XSP(200,5)) 
      ALLOCATE(NX(200))
END SUBROUTINE Allocate_Parms
 
! ------------------------ initialize other vars --------------------------
SUBROUTINE Initialize_Var(Sulf_con)
      
      REAL, DIMENSION(15):: Sulf_con
      
      ICV=0; ICUS=0
      IDFT(1)=fertTimes
      IDRL=0; IGMD=0; IGO=0; IHRL=0; IHT=0; IHU=0; IHV=0; IHVD=0; GrowingSeason_Days=0; IPY=1
      IPYI=1; IRL=0; ISIX=0; ISX=0; IYR=Year0; JCN=0; JCN0=0; JCN1=0; JE=MNC+1
      JPL=0; varID=0; KDT2=0; KG=0; LC=0; LY=0; LYR=1; LW=1; MXT=0; NBC=0
      NBCX=0; NCP=0; NCR=0; NDT=0; NDF=0; NDP=0; NDT=0; NFA=0; NHU=0; NHV=0 
      Mow_DayLap=0; NPSF=0; NPT=0; NQP=0; NQP0=0; NQP1=0; NWDA=0; NWD0=0; NYLN=0
      ACET=0.; accumRainfall30=100.; AFLG=0.; standLiveDeadPlant=0.; aveWaterFromLagoon=0.; AJHI=0.
      ANA=0.; plawDepLabileP=0.; APQ=0.; APY=0.; AQB=0.; ASW=0.; AYB=0.; BARF=0.; bulkDensity=0.; BLYN=0.
      CaCO3=0.; availWaterCrop=0.; CEC=0.; fertCFactirArr=.05        !CAP=0.; 
      gasO2Con=0.; gasCO2Con=0.; gasN2OCon=0.; initialNO3=0.; CPHT=.01; COST=0.
      CRF=0.; CQV=0.; CSFX=0.; CSO1=0.; CST1=0.; CSTF=0.; groundCover=0.; coverFactorPlantPopu=0.
      abvGroundResidue=0.; CX=1.E-10; CYAV=0.; CYMX=0.;   CYSD=0.; DARF=0.; HumusMine_RateCons=.01
      dikeHeight=0.; DKHL=0.; totCropBio=0.;  totCropBioX=0.; DMF=0.; DN2G=0.; DN2OG=0.; DRWX=0.
      EAR=1.; elecCondv=0.;   EP=0.; ET_GrowingSeason=0.; plawDepFieldCapa=0.; ironCon=0.
      GroundCover_Frac=0.; GroundCover_StandLiveBio=0.; TotK_Frac=0.; TotN_Frac=0.;  TotP_Frac=0.; Fuel_Use=0.
      GS_PlantTranspiration=0.; GSEP=0.; GS_VPD=0.; GSVP=0.; lateralFlowCondv=0.; HIF=0.; HU=0.; HUF=0.; HUI=0.
      plawDepSOC=0.; potentialAvailWater=0.; plawDepSoilWater=0.; PH=0.; potentialHeatUnit=0. 
      PMORF=100.; Porosity=0.; POP=0.; PPL0=0.; PPLA=0.; PRAV=0.; PRB=0.; PRSD=0.
      sorbRatioP=0.; PSTF=1.; PSTM=0.; pestName=' '; pestGrowIndx=0.; PVQ=0.; PVY=0.
      Inflow=0.; QPQB=0.; QPS=0.; RCF=1.; RCM=0.; RD=0.;  RDF=0.; Min_Stress_Factor=0.; Tot_Rainfall=0.
      RidgeHeight2=0.; Ridge_IntervalX=0.; rockFrac=0.; RandRough=.01; cropResidu=0.
      soilResp=0.; baseStockRateArr=0.; RSY=0.; RUSM=0.; totRootWeight=0.; RWF=0.; RWT=0.; RWTX=0.
      RWTZ=0.; RWX=0; soilWaterRZ=0.; sandFrac=0.; SARF=1.E10; satuateCondv=0.
      SET=0.; Growing_Stress=0.; SFMO=0.; SHRL=1.; siltFrac=0.; Current_LAI=0.; SLTP=0.; SM=0.
      SMAP=0.; totBase=0.; Var_GS=0.; SMM=0.; SMMC=0.; SMMP=0.; SMS=0.; SMY=0.
      SMYP=0.; SOIL=0.; SOT=0.; SPQ=0.; Pest_Sediment=0.; SRA=0.; SRD=0.
      SubSurfFlow_NO3=0.; SubSurfFlow_NO2=0.; SubSurfFlow_N2O=0.; SubSurfFlow_CO2L=0.; SubSurfFlow_O2L=0.
      SubSurfFlow_SolubleK=0.; SubSurfFlow_LaibleP=0.; lateralFlow=0.; SSW=0.; soilWater=0.
      standCropResi=deadCropRes0; STDA=0.; standDeadResiN=0.; standDeadResOrganism=0.; StandDeadRes_OrgK=0.
      standDeadResiOrgN=0.; StandDeadRes_OrgP=0.; fracNO3Leach=0.; abvGroundBiom=0.; soilTem=0.
      STV=0.; Sulf_con=0.; totHfHourRain=0.; TAMX=0.; totLabileP=0.; TCAV=0.; TCAW=0.
      TCQV=0.; TCRF=0.; TCMX=0.; TCMN=1.E20; TCSO=0.; TCST=0.; TDM=0.; TEI=0.
      totExchangeK=0.; TET=0.; TETG=0.; totFixK=0.; TFTK=0.; TFTN=0.; TFTP=0.; THK=0.;  THU=0.
      TIRL=0.; TLMF=0.; totActiveP=0.; TMXF=0.; TNOR=0.; ZNO2=0.; totNO3_N=0.; totSOC=0.; TOP=0.
      initialTotSOP=0.; totPestControl=0.; TPOR=0.; TPSF=0.; TQ=0.; TR=0.; TRA=0.; TRD=0.; TRHT=0.
      totCropResi=0.; TSFC=0.; totSolubleK=0.; totSalt=0.; TSN=0.;  snowPackAge=0.; TSR=0.; TSTL=0.;  TSY=0.
      TVAL=0.; TVIR=0.; totSON=0.; TXMN=0.; TXMX=0.; TYC=0.; TYK=0.; TYL1=0.
      TYL2=0.; TYLC=0.; TYLN=0.; TYLP=0.; TYLK=0.; TYN=0.; TYP=0.; TYW=0.
      wiltPoint=0.;  actualCropK=0.; soilSupplyRateN=0.; actualCropN=0.; UNMX=0.; actualCropP=0.; UPMX=0.
      VerticleFlow_LaibleP=0.; cropVar=0.; pestVar=0.; VFC=0.; VIRT=0.; VIL=0.; VIR=0.; verticalFlowNO3=0. 
      verticalFlowN2O=0.; VerticalFlow_CO2=0.; VerticalFlow_O2=0.; VerticleFlow_SolubleK=0.; verticleFlowSalt=0.
      VQ=0.;  VWC=0.; VWP=0.;  VY=0.; VALF1=0.; VerticalFlow_NO2=0.; WCHT=0.; WCO2G=0.; Weight_CO2L=0.; W=0.
      CBiomass=0.; NBiomass=0.; CPassiveHumus=0.; NPassiveHumus=0.; CSlowHumus=0.; NSlowHumus=0.
      CMetabLitt=0.; NMetabLitt=0.; metabLitt=0.; structLitt=0.; NStructLitt=0.
      CStructLitt=0.; lgStructLitt=0.; CLgStructLitt=0.; NLgStructLitt=0.
      WLV=0.; WN2O=0.; WN2OG=0.; Weight_N2OL=0.; NO2Weight=0. ; waterStress=1.; WTN=0.
      NH3_Weight=0.; NO3_N_Soil=0.; WO2G=0.; Weight_O2L=0.; XIM=0.; XMTU=0.; enrichRatioFinal=0.
      YLC=0.; YLCF=0.; YLD=0.; YLD1=0.; YLD2=0.; YLK=0.; YLKF=0.; Yield_N=0.; YLNF=0.; YLP=0.; YLPF=0.; Z=0.
      totCBiomass=0.; totNBiomass=0.; totCPassiveHumus=0.; totNPassiveHumus=0.;  totCSlowHumus=0.; totNSlowHumus=0.
      totMetabLitt=0.; totCMetabLitt=0.; totNMetabLitt=0.; totStructLitt=0.; totCStructLitt=0.
      totLgStructLitt=0.; totCLgStructLitt=0.; totNStructLitt=0.; totNLgStructLitt=0.
      
END SUBROUTINE Initialize_Var
! ------------------------------ Reset ---------------------------------
SUBROUTINE Reset
    INTEGER:: J
    
      IDRL=0; IHRL=0; JE=MNC+1; JPL=0; IRL=0
      IBD= Cal_DOY(Month0,Day0, leapYr)
      JCN=0; JCN0=0; JCN1=0; LW=1; MO=Month0; MO1=MO; NCR=0; NQP=0; NQP0=0; NQP1=0
      NWDA=0; NWD0=0; NYLN=0; APQ=0.; APY=0.; AQB=0.; ASW=0.; AYB=0.; CST1=0. 
      CX=1.E-10; CYAV=0.; CYSD=0.; CYMX=0.; EP=0.; PRAV=0.; PRB=0.; PRSD=0.
      PVQ=0.; PVY=0.; Inflow=0.; QPQB=0.; QPS=0.; RCF=1.; RCM=0.; Min_Stress_Factor=0.; RSY=0.
      SET=0.; SFMO=0.; SM=0.; SMAP=0.; SMM=0.; SMMP=0.; SMY=0.; SMYP=0.; SPQ=0.
      Pest_Sediment=0.;  SRD=0.; STDA=0.;  totHfHourRain=0.; TAMX=0.; TCAV=0.; TCAW=0.
      TCMN=1.E20; TCMX=0.; TCQV=0.; TCRF=0.; TDM=0.; TEI=0.; TET=0.; TETG=0.
      TFTN=0.; TFTP=0.; THU=0.; TQ=0.; TR=0.; TRA=0.; TRD=0.; TRHT=0.;  TSFC=0.
      TSN=0.; TSR=0.; TSTL=0.; TSY=0.; TVIR=0.; TXMX=0.; TXMN=0.; TYC=0.
      TYK=0.; TYL1=0.; TYL2=0.; TYLC=0.; TYLK=0.; TYLN=0.; TYLP=0.; TYN=0.; TYP=0.
      TYW=0.; U10MX=0.; VALF1=0.; pestVar=0.;  W=0.; XIM=0.; YLC=0.; YLK=0.; Yield_N=0.; YLP=0.
      DO J=1,21
          IX(J)=IX0(J)
          randID(J)=randID0(J)
      END DO
      V3=Generate_Random(randID(3))                                                                                                                                               
    
END SUBROUTINE Reset
    
END MODULE Initialize