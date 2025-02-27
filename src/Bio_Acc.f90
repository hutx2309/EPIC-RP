MODULE Bio_Acc
!****************************************************************************
!
!  PURPOSE: A module calculate biomass accumulation considersing the process 
!           of photosynthesis, instead of using the RUE method  
!  (by TXH   06-02-2020)   
! Version 2: excluding r_t    
!****************************************************************************    
USE PARM
implicit none   
! ------------------------ 1. parameters for photosynthesis ------------------
! physical constants
public
real, parameter:: sc = 1367        ! solar constant W/m2 = J/m2/s
real, parameter:: pi = 3.1415926 
real, parameter:: omega = pi/12    ! 
real, parameter:: P = 101325       ! standard air pressure  Pa
real, parameter:: Rgas = 8.314     ! universal gas constant(J/mol/K)  
real, parameter:: Boltz = 5.668E-8 ! Stefan-Boltzmann constant                  J/m2/s/K4
real, parameter:: Latent = 2.4E+6  ! latent heat of water vaporization          J/kg
real, parameter:: V_ca = 1200.0    ! Volumetric heat capacity                   J/m3/degree
real, parameter:: Psy = 0.067      ! psychrometric constant                     kPa/degree

! Parameters for photosynthesis rate 
! See CLM 4.5/5.0 for details
real, parameter:: E_Vcmax = 65330.    !energy of activation for Vcmx(J/mol)
real, parameter:: Ed_Vcmax = 149250
real, parameter:: S_Vcmax = 485

real, parameter:: Rd_Vcmax25 = 0.0089 !ratio of dark respiration to Vcmax at 25oC
real, parameter:: Ea_Rdmax = 46390    !energy of activation for dark respiration(J/mol)
real, parameter:: Ed_Rdmax = 150650
real, parameter:: S_Rdmax = 490

real, parameter:: Ea_Jmax = 43540
real, parameter:: Ed_Jmax = 152040
real, parameter:: S_Jmax = 495

real, parameter:: E_KMC = 79430.      !energy of activation for KMC (J/mol)
real, parameter:: E_KMO = 36380.      !energy of activation for KMO (J/mol)
real, parameter:: KmC25 =  404.9      ! umol / mol    This only applied to C3 in CLM
real, parameter:: KmO25 =  278.4      ! mmol / mol   but different values are used for C4 in CSD  
real, parameter:: Gamma25 = 42.75     ! umol / mol 
real, parameter:: Ea_Gamma = 37830

real, parameter:: s1= 0.3       ! slope of high temperature inhibition function     (0.3)
real, parameter:: s2 = 313.16   ! 1/2 point of high temperature inhibition function (313.16)
real, parameter:: ss3 = 0.2     ! slope of low temperature inhibition function      (0.2)
real, parameter:: s4 = 288.15   ! 1/2 point of low temperature inhibition function  (288.16)
real, parameter:: s5 = 1.3      ! temperature coefficient in gs-a model             (1.3)
real, parameter:: s6 = 328.16   ! temperature coefficient in gs-a model             (328.16)
real, parameter:: B = 0.409     ! biomass conversion efficiency

real, parameter:: Kd = 0.719    ! diffues and scattered diffuse PAR extinction coefficient

real, parameter:: Vcmax25 = 100.7 ! N determined Vcmax at 25 oC (Bonan, 2012; CLM 4.5/5.0)
real, parameter:: Kn = 0.11       ! extinction coefficient of N (Bonan, 2012)
real, parameter:: albg = 0.3      ! soil albedo
! Enviromental Variables
    real:: Latr, E0, SolarDec, T_SR, T_SS, h_SR, h_SS, DayLen, &
           totR, Tot_Sh, Diff_frac, Tss, AVP
    INTEGER:: PreDawn, PostDusk
    real, dimension(24):: AccSolar, AccFrac, Frac, T_air, PAR, NIR
   
! Biomass part 1 --- for radiation calculation --------------
 
    real:: CosTheta, Fdiff, PAR_diff, PAR_dir, NIR_diff, NIR_dir,absPAR_Sun, absPAR_Sha, &
           absRad_tot, absRad_Sha,absRad_Sun
 
    real, dimension(2)::Idir_sun, Idir_sha, Idif_sun, Idif_sha, ref, trans
 
    real:: Kb
    
    real, dimension(3):: Lf2Cp_sun, Lf2Cp_sha  

! Biomass part 2 --- for H20, CO2, heat conductance ---------
    real:: Kkw, f_sun, f_sha,gb_cp, gb_lf, gb_Sun, gb_Sha, rbw_Sha, rbw_Sun, rb_Sun, rb_Sha
           
    
! Biomass part 3 --- for  T_lf, CO2i, Ass_n -----------------
    real:: T_lf, CO2i1, SVPL1, vVPD, SVP,  & 
           Photo1_Sun, Photo1_Sha, &
           SlopeL1, r_swp1_Sun, r_swp1_Sha, Ep_Sun1, NetRd_Sun1, Ep_Sha1, NetRd_Sha1 
           
    real:: del_t_Sun,del_t_Sha, T_lf_Sun, T_lf_Sha, SVPL2_Sun, CO2i2_Sun, SVPL2_Sha, CO2i2_Sha, &
           Photo2_Sun, Photo2_Sha,&
           Slope_lf2_Sun, Slope_lf2_Sha, r_swp2_Sun, r_swp2_Sha, &
           Ep2_Sun, NetRd2_Sun, Ep2_Sha, NetRd2_Sha
! --------------------------- 2. Subroutines and Functions------------------------- 
! *****************************************************************************
! This part contains four main subroutines: 
    ! (1) Two-stream Radiation ;
    ! (2) CO2i calcuation ; 
    ! (3) Instaneous Photosynthesis rate;
    ! (4) Evaporation calculation;
! *******************************************************************************
contains 
    ! ************************ (1) Radiation : Two-stream approach ******************
    ! 
    ! *******************************************************************************
    Subroutine Rad_TwoStream_Dai(lai, chil)
    ! ( Ref Dai et al. 2003 )
    implicit none
    
    real, intent(in)::   lai, &              ! Leaf area index
                         chil                ! Leaf angle distribution factor   
                         !albg, &             ! Soil albedo  
                         !CosTheta            ! Solar zenith angle
    !real, dimension(2), intent(in):: ref, &  ! Leaf reflectance to PAR and NIR
     !                               trans   ! Leaf transmisstance to PAR and NIR 
    
    ! local variables: not accessible to other programs
    real:: phi1, phi2, G_mu, mu, mu_bar, mu_bar2, scat, as, upscat, beta0, be, ce, de, fe, &
           psi, power1, power2, ss1, ss2, p1, p2, p3, p4, f1, f2, h1, h4, sigma 
    
    real:: m1, m2, m3, n1, n2, n3,hh1, hh4, hh2, hh3, hh5, hh6, hh7, hh8, hh9, hh10    
           
    real:: Adir_veg, Tdir_canopy, Adif_veg, Tdif_canopy, Idir_up, Idir_down, Idif_up, Idif_down
 
    integer:: iw
 
     ! equations in CoLM
    
     phi1 = 0.5 - 0.633*chil - 0.33*chil**2
     phi2 = 0.877*(1-2*phi1)
    
     mu = CosTheta
     G_mu = phi1+phi2 * mu

     Kb = G_mu/mu     ! optical depth of the direct beam per unit leaf area
 
    ! intergration of 1/Kb from 0 to 1

if (abs(phi1) > 1e-5 .and. abs(phi2) > 1e-5)  then
   
    mu_bar = 1/phi2*(1-phi1/phi2*log((phi1+phi2)/phi1))
    else if (abs(phi1) < 1e-5) then 
    mu_bar = 1./0.877
    else if (abs(phi2)< 1e-5) then
    mu_bar = 1/(2*phi1)

endif
! squared mu_bar
mu_bar2 = mu_bar**2
 
Do iw = 1, 2      ! PAR and NIR
    
 scat = ref(iw)+ trans(iw)

 as = scat / 2. * G_mu / ( G_mu+ mu * phi2 )*&
      ( 1. - mu * phi1 / ( G_mu + mu * phi2 )*&
      log((G_mu + mu * phi2 + mu * phi1 ) / ( mu * phi1 ) ) )
  
 ! Eq 3.14 in CLM 4.5  Omega^ * Beta^             
 upscat = 0.5 * ( scat + (ref(iw) - trans(iw))*((1.+ chil)*0.5)**2)                                      
 ! Eq. 3.15 
 beta0 = ( 1. + mu_bar * Kb) / ( scat * mu_bar * Kb) * as
 
      be = 1. - scat + upscat                                                   
      ce = upscat                                                              
      de = scat * mu_bar * Kb * beta0                                        
      fe = scat * mu_bar * Kb * ( 1. - beta0 )                               

      psi = sqrt(be*be - ce*ce)/mu_bar                                                   
      power1 = min( psi*lai, 50. )                                            
      power2 = min( Kb*lai, 50. )
      
      ss1 = exp( - power1 )                                                  
      ss2 = exp ( - power2 )   
      
      p1 = be + mu_bar * psi
      p2 = be - mu_bar * psi
      p3 = be + mu_bar * Kb
      p4 = be - mu_bar * Kb

      f1 = 1. - albg*p1/ce
      f2 = 1. - albg*p2/ce
      
      h1 = - ( de * p4 + ce * fe )
      h4 = - ( fe * p3 + ce * de )

      sigma = (mu_bar *Kb)**2 + (ce*ce - be*be)                           

      if (abs(sigma) > 1.e-10 ) then   ! this should be consistant with CLM

         hh1 = h1 / sigma
         hh4 = h4 / sigma
                                                                                
         m1 = f1 * ss1
         m2 = f2 / ss1
         m3 = ( albg - ( hh1 - albg* hh4 ) ) * ss2
 
         n1 = p1 / ce
         n2 = p2 / ce
         n3 = - hh4

         hh2 = (m3*n2 - m2*n3) / (m1*n2 - m2*n1)
         hh3 = (m3*n1 - m1*n3) / (m2*n1 - m1*n2)

         hh5 = hh2 * p1 / ce
         hh6 = hh3 * p2 / ce

         Adir_veg = hh1 + hh2 + hh3                                 
         Tdir_canopy = hh4  * ss2 + hh5 * ss1 + hh6 / ss1                    

         Idir_up = hh1 * (1. - ss2**2) / (2.*Kb)+ hh2 * (1. - ss1*ss2) / (Kb + psi)&
                   + hh3 * (1. - ss2/ss1) / (Kb - psi) 

         Idir_down = hh4 * (1. - ss2**2) / (2.*Kb) + hh5 * (1. - ss1*ss2) / (Kb + psi)&
                     + hh6 * (1. - ss2/ss1) / (Kb - psi)

      else                           

         m1 = f1 * ss1 
         m2 = f2 / ss1 
         m3 = h1 /mu_bar2 * (lai + 1. / (2.*Kb)) * ss2 + albg/ce * (-h1/(2.*Kb)/mu_bar2 * &
              ( p3*lai + p4 / (2.*Kb)) - de) * ss2 + albg * ss2
 
         n1 = p1 / ce
         n2 = p2 / ce
         n3 = ( h1* p4 /(4.*Kb*Kb*mu_bar2) + de) / ce

         hh2 = (m3*n2 - m2*n3)/(m1*n2 - m2*n1)
         hh3 = (m3*n1 - m1*n3)/(m2*n1 - m1*n2)

         hh5 = hh2 * p1/ce
         hh6 = hh3 * p2/ce

         Adir_veg =  - h1/(2.*Kb*mu_bar2) + hh2 + hh3
         Tdir_canopy = 1./ce *( -h1 / (2.*Kb*mu_bar2) * &
                       ( p3*lai + p4 / (2.*Kb) ) - de ) * ss2&
                       + hh5 * ss1 + hh6 / ss1

         Idir_up = (hh2 - h1/(2.*Kb*mu_bar2))*(1.- ss2**2)/(2.*Kb)&
                   + hh3 * lai &
                   + h1/(2.*Kb*mu_bar2) * ( lai*ss2*ss2- (1.- ss2**2)/(2.*Kb))

         Idir_down = (hh5 - (h1*p4/(4.*Kb*Kb*mu_bar)+ de)/ce)*(1. - ss2**2) / (2.*Kb)&
                     + hh6 * lai      &
                     + h1*p3/(ce*4.*Kb*Kb*mu_bar2) * ( lai*ss2*ss2 - (1. - ss2*ss2)/(2.*Kb))     
      endif                           

      Idir_sun(iw) = (1.-scat) * ( 1.-ss2 + 1. / mu_bar * (Idir_up+ Idir_down) ) 
      Idir_sha(iw) = scat * (1.-ss2)+ ( albg*Tdir_canopy + albg*ss2 - Tdir_canopy) - Adir_veg &
               - ( 1. - scat ) /mu_bar * ( Idir_up  + Idir_down ) 
      
     ! diffuse radiation  

      m1 = f1 * ss1
      m2 = f2 / ss1
 
      n1 = p1 / ce
      n2 = p2 / ce
 
      hh7 = -m2 / (m1*n2 - m2*n1)
      hh8 = -m1 / (m2*n1 - m1*n2)

      hh9 = hh7 * p1 / ce
      hh10 = hh8 * p2 / ce

      Adif_veg =  hh7 + hh8                                            
      Tdif_canopy = hh9 * ss1 + hh10 / ss1

      if (abs(sigma)>1.e-6) then 
         Idif_up = hh7 * (1. - ss1*ss2) / (Kb + psi) &
                 + hh8 * (1. - ss2/ss1) / (Kb - psi)
         Idif_down = hh9 * (1. - ss1*ss2) / (Kb + psi) &
                 + hh10 * (1. - ss2/ss1) / (Kb - psi)
      else
         Idif_up   = hh7 * (1. - ss1*ss2) / (Kb + psi) + hh8 * lai
         Idif_down = hh9 * (1. - ss1*ss2) / (Kb + psi) + hh10*lai
      endif

      Idif_sun(iw) = (1.-scat) / mu_bar * (Idif_up + Idif_down)
      Idif_sha(iw) = Tdif_canopy * ( albg -1. ) - ( Adif_veg - 1. ) &
                 - ( 1. - scat ) / mu_bar * ( Idif_up + Idif_down)         
     End Do ! band loop

    End Subroutine Rad_TwoStream_Dai
 
    ! ****************** (2) CO2i : CSD approach ***************************
    ! modified from the CSD model - Yin 2005
    ! **********************************************************************
    Subroutine cal_CO2i(C3C4, CiCa_ratio, CO2a, O2, SVP_lf, CO2i)
           
    implicit none
             real, intent(in):: CiCa_ratio, CO2a, O2
             integer, intent(in):: C3C4
             real, intent(out):: SVP_lf, CO2i
     ! local variables
             real:: VPD_lf, KmC, KmO, Ga_Max, Rd_Vcmax, Ga0, Ga, CiCa
    
    SVP_lf = 0.611 * exp(17.4* T_lf/(T_lf + 239.0))  ! kPa
    VPD_lf = max(0.0, SVP_lf - AVP)

    KmC = KmC25*exp((1./298.-1./(T_lf+273.))*E_KMC/Rgas) ! umol/mol   CO2
    KmO = KmO25*exp((1./298.-1./(T_lf+273.))*E_KMO/Rgas) ! mmol/mol    O2

    ! compensation point in the absence of dark respiration
    Ga_Max = 0.5* exp(-3.3801+5220./(T_lf+273.)/Rgas)*O2*KmC/KmO   ! umol/mol   CO2
    
    ! ---CO2 compensation point (GAMMA)
    ! this can be modifed by CLM
    Rd_Vcmax = Rd_Vcmax25* exp((1./298.-1./(T_lf+273.))*(Ea_Rdmax-E_Vcmax)/Rgas)

    Ga0 = (Ga_Max+Rd_Vcmax*KmC*(1.+O2/KmO))/(1.-Rd_Vcmax) ! for C3 crops

    if (C3C4 == 1) then
    
    Ga = Ga0 
    elseif (C3C4 ==0) then
    
    Ga= Ga0/10.0 
    endif
    !---internal/ambient CO2 ratio, based on data of Morison & Gifford (1983)
    CiCa = 1.-(1.-Ga/CO2a)*(0.14+CiCa_ratio*VPD_lf)
    !---intercellular CO2 concentration
    CO2i = CiCa * CO2a
    
    End Subroutine cal_CO2i
    
    ! ******** (3) Photosynthesis Rate : CLM 4.5/5.0 approach **************
    !
    ! **********************************************************************
    
    Subroutine photo_rate (C3C4, Vcmax25,T_leaf, Lf2Cp, absPAR, O2, CO2i, Photo)
    implicit none
             real, intent(in):: T_leaf, absPAR,O2,CO2i
             real, dimension(3), intent(in):: Lf2Cp
             integer, intent(in):: C3C4
             real, intent(out)::Photo
             real:: Vcmax25
   !local variables
   integer:: C3, C4
   real:: Rt, kp25,Tf, qt, Rd,             & 
          fT_V, fH_V, fh_Vc, fl_Vc, Vcmax, &
          Jmax25, fT_J, fH_J, Jmax,      &
          Rd25, fT_R, fH_R,              &
          Kc, Ko, Gamma, KoKc,           &
          Ac, Aj, Ap, ePAR, Iphi2, Phi2, &
          J, kp, Thetaj, Thetap, detB, Ai, A 
   ! reference : CoLM and CLM 5.0    
  
    C3 = C3C4
    C4 = 1-C3C4
   ! pay attention : Temperature needes to be in K  
    Tf = T_leaf + 273

! *************************** CLM4.5/5.0 ************************
! Ref. Chapter 8 in CLM 4.5
   Rt = 298.15 * Rgas
   Jmax25 = 1.97*Vcmax25
   Rd25 = (0.015 * C3 + 0.025 * C4) * Vcmax25

! only for C4 plants
   kp25 = 20000*Vcmax25  
   qt = (Tf - 298.15)/10  ! Q10   
 
! ---------------- 1. Vcmax ---------------------------
fT_V = exp(E_Vcmax*(1-298.15/Tf)/Rt)
fH_V = (1 + exp((298.15*S_Vcmax - Ed_Vcmax)/Rt))/(1 + exp((Tf*S_Vcmax - Ed_Vcmax)/Rt))

fh_Vc = 1+exp(s1*(Tf-s2))
fl_Vc = 1+exp(ss3*(s4-Tf)) 

Vcmax = Vcmax25*fT_V*fH_V * C3 + C4*Vcmax25 * 2.0**qt/fh_Vc/fl_Vc

! scale to the canopy
Vcmax = Vcmax * Lf2Cp(1)

! Vcmax = Vcmax * WaterStress; % add water stress here

! ----------------2. Jmax ----------------------------
fT_J = exp(Ea_Jmax*(1-298.15/Tf)/Rt)
fH_J = (1 + exp((298.15*S_Jmax - Ed_Jmax)/Rt))/(1 + exp((Tf*S_Jmax - Ed_Jmax)/Rt))

Jmax = Jmax25*fT_J*fH_J
Jmax = Jmax * Lf2Cp(2)
! Jmax = Jmax * WaterStress; % add water stress here

! ----------------3. Rd ------------------------------
fT_R = exp(Ea_Rdmax*(1-298.15/Tf)/Rt)
fH_R = (1 + exp((298.15*S_Rdmax - Ed_Rdmax)/Rt))/(1 + exp((Tf*S_Rdmax - Ed_Rdmax)/Rt))

Rd = Rd25 * fT_R * fH_R * C3 + C4*Rd25*2.0**qt/(1 + exp(s5*(Tf - s6)))
Rd = Rd * Lf2Cp(1)
! Rd = Rd * WaterStress; 
! Rd = Rd * 44*1e-6      !  umol CO2 / m2/s  ---->   g CO2 /m2/s

! =============== 3. cal net assimilation rate ====================

! -------------------------- Ac (umol CO2/m2/s)-------------------- 
! RuBP carboxylase limited rate of carboxylation Ac 
Kc = KmC25 * exp(E_KmC*(1-298.15/Tf)/Rt)               ! umol /mol
Ko = KmO25 * exp(E_KmO*(1-298.15/Tf)/Rt)               ! mmol /mol
Gamma = Gamma25*exp(Ea_Gamma*(1-298.15/Tf)/Rt)         ! umol /mol  

KoKc = Kc *(1+O2/Ko)*C3

Ac = Vcmax * ( CO2i-Gamma ) / ( CO2i + KoKc ) * C3 + Vcmax * C4 

! ----------------------- Aj, umol CO2/m2/s ----------------------------
! The maximum rate of carboxylation allowed by the capacity to regenerate
! RuBP (the light limited rate )
 ePAR = 4.6*absPAR   ! PAR  J/m2/s = W/m2  -----> 4.6 umol photon/J
 
 Iphi2 = 0.5*0.85*ePAR
 Phi2 = 0.7
! the small root of a quadratic function related to Jmax 
 J = (Iphi2 + Jmax -sqrt((Iphi2+ Jmax)**2 - 4*Iphi2*Jmax*Phi2))/2./Phi2
 
 Aj = J * ( CO2i- Gamma ) / ( 4*CO2i+8.*Gamma) * C3 + ePAR*0.05 * C4  
 
! ----------------------- Ap, product-limited rate ---------------------
 kp = kp25*2.0**qt  ! for C4 
 
 ! CO2i is the internal leaf CO2 partial pressure 
 ! P is atmosphere pressure (Pa)
 Ap = 0.5*Vcmax * C3 + kp*CO2i/P/10 * C4
 
! --------------------- co-limitation A --------------------------------
! Ref. Eq . 8.8 in CLM 4.5

    Thetaj = 0.98* C3 + 0.8 * C4
    Thetap = 0.95
    
    detB = max( 0., ( (Ac+Aj)**2 - 4.*Thetaj*Ac*Aj) )
    
    Ai   = (( Ac + Aj ) - sqrt(detB)) / ( 2.*Thetaj)
    
    detB = max( 0., ( (Ap+Ai)**2 - 4.*Thetap*Ap*Ai))
    
    A = ( ( Ai + Ap ) - sqrt(detB) ) / ( 2.* Thetap)  
    
    Photo = max(0.0, 44*1e-6*(A-Rd))  ! net photosynthesis rate  umol CO2/m2/s ---> g CO2/m2/s
    
End Subroutine photo_rate 
    
    ! *************** (4) Conductance : CSD approach ***********************
    !
    ! **********************************************************************
    real function r_swp(Photo, Tf, CO2a, CO2i, rbw)
    implicit none 
        real, intent(in):: Photo, Tf, CO2a, CO2i, rbw
        !local variables
        real:: gc
        
        gc = Photo *(273.+Tf)/0.53717/(CO2a-CO2i)
 
        r_swp = max(1E-10, 1./gc - rbw*1.4)/1.6  
    
    end function r_swp
    
    ! ********************** (5) EP: CSD approach *************************
    !
    ! ********************************************************************** 
    
   Subroutine cal_EP(r_swp, rbw, rb, absRad,Atmo_Tran,f_s, Tf, AVP, Slope, vVPD, Ep, NetRd)
    implicit none
    real,intent(in):: r_swp, rbw, rb, absRad,Atmo_Tran,f_s, Tf, AVP, Slope, vVPD
    real, intent(out):: Ep, NetRd
    ! local variable
    real:: Clear, Blackbody, Longwave_net, Psr, Ptr, Ptd
 
    Clear = max(0., min(1., (Atmo_Tran-0.25)/0.45) )    ! sky clearness
    Blackbody = Boltz*(Tf +273.)**4                     ! black body radiation

    Longwave_net = Blackbody*(0.56-0.079*sqrt(AVP*10.))*(0.1+0.9*Clear)*f_s

    NetRd = absRad - Longwave_net

    !*---intermediate variable related to resistances
    Psr = Psy*(rbw +r_swp)/rb

    !*---radiation-determined term
    Ptr = NetRd*Slope /(Slope+Psr)/Latent;

    !*---vapour pressure-determined term
    
    Ptd = (V_ca*vVPD/rb)/(Slope+Psr)/Latent ! rb : leaf boundary layer resistane to heat (s/m)

    !*---potential evaporation or transpiration
    Ep = max(1.E-10,Ptr+Ptd)
 
 End Subroutine cal_EP
    ! ********************** (6) EP: CSD approach *************************
    !
    ! ********************************************************************** 
    real function del_T(NetRd, Ep, rb) 
    implicit none
       real, intent(in):: NetRd, Ep, rb
 
       del_T = (NetRd -Latent* EP)*rb/V_ca  

      if (del_T < -25) then    
         del_T = -25
      elseif (del_T > 25) then
         del_T = 25
      endif 
 
    end function del_T   
    
! --------------------------- 3. Cal Biomass Accumulation --------------------    
 Subroutine Poten_Growth
 implicit none 
    ! User determined variables 
    ! real:: lat = -27.5    
    ! integer:: DOY = 298 --> Date
    ! real:: lai = 6.0     
    ! real:: Tmax = 30  
    ! real:: Tmin = 20
    ! real:: Wnd = 1.5
    real:: Atmo_Tran = 0.75  ! Atmosphere transparance
    real:: chil = -0.5       ! Leaf angle parameter
    ! real:: Leaf_width = 0.1  ! leaf width  (!! this one should be dynamic ??)  
    real:: O2 = 210          ! mmol / mol   partial pressure of O2
    ! real:: CO2a = 370        ! umol / mol   partial pressure of CO2
    
    real:: Latr, lai, Tmax, Tmin, Wnd, Leaf_width, CO2a, & !qt, Rt, fH_R, fT_R, Rd25, &
           CiCa_ratio  ! ratio of CO2i/CO2a, CO2i interior CO2 within canopy  
    
    integer:: C3C4, C4, C3, t   
 
    ! output variables 
    real, dimension(24):: InsPP, InsEP     ! daily net biomass accumulation and EP
    real:: Bio
    
    ! *************************** Program begin ***************************************
    
    ! *********************************************************************************
    ! interfaces with EPIC
    Latr = YLAT/CLT
    lai = Current_LAI(Crop_Num) 
    Tmax = TMX
    Tmin = TMN
    Wnd = U10
    Leaf_width = LfWidth(Crop_Num)
    CO2a = CO2
    
    C3C4 = 0
    C3 = C3C4
    C4 = 1 - C3C4
    CiCa_ratio = 0.12 * C3 + 0.19 * C4
    
    ref = (/0.11, 0.35/)   ! leaf reflectance of PAR and NIR
    trans = (/0.05, 0.34/) ! leaf transmisstance of PAR and NIR
 
    ! ============================ I. Environmental Factors  ===========================
       
    SolarDec = 23.45*sin(pi*2*(Date+284)/365)*pi/180 ! Solar declination angle in radians
    T_SR = acos(-tan(Latr)*tan(SolarDec))/omega  ! duriation from noon to the hour of sunrise
    T_SS = -acos(-tan(Latr)*tan(SolarDec))/omega ! duriation from noon to the hour of sunset

    h_SR = 12 - T_SR                             ! sunrise hour
    h_SS = 12 - T_SS                             ! sunset hour

    DayLen = 2*acos(-tan(Latr)*tan(SolarDec))/omega
   
    if (SRAD<1.E-5.OR.SRAD>900..OR.KGN(3)==0) then 
    E0 = 1+0.033*cos(2*pi*Date/365)
    ! The maximum Radiation reaches the earth on a clear day (J/m2/day)
    totR = 24*E0*sc*3600*(omega*T_SR*sin(SolarDec)*sin(Latr)+ cos(SolarDec)*cos(Latr)*sin(omega*T_SR))/pi;
    else
    totR = SRAD*1.E6          ! conversion: 1.0 MJ/m2/d = 1000000 J/m2/d 
    endif
    
    PreDawn =floor(12 - T_SR)     
    PostDusk = ceiling(12 - T_SS) 

    Tot_Sh = DayLen*sin(Latr)*sin(SolarDec) + 24*cos(Latr)*cos(SolarDec)*sin(pi*DayLen/24)/pi;
 
    Do t = 12 , PostDusk
    
    if (t < h_SS) then
        
        AccSolar(t) = (t-12)*sin(Latr)*sin(SolarDec) + 12*cos(Latr)*cos(SolarDec)*sin(pi*(t - 12)/12)/pi;
        AccFrac(t)= AccSolar(t)/Tot_Sh;
    
    elseif (t>= h_SS) then
        
        AccSolar(t) = (h_SS-12)*sin(Latr)*sin(SolarDec) + 12*cos(Latr)*cos(SolarDec)*sin(pi*(h_SS - 12)/12)/pi;
        AccFrac(t)= AccSolar(t)/Tot_Sh;
    endif
    
    End Do
    
  ! the fraction of radiation in each hr after predown and before postdusk
 
  Frac(13:PostDusk) = AccFrac(13:PostDusk)- AccFrac(12: PostDusk-1)
  Do t = PreDawn, 12  
   Frac(t)= Frac(24-t+1)
  End Do
  
  ! Diffuse light fraction (FRDF) from atmospheric transmission (Atmo_Tranmiss)
if (Atmo_Tran < 0.22 .OR. Rainfall >0.0 ) then
    
    Diff_frac = 1
elseif (Atmo_Tran > 0.22 .and. Atmo_Tran < 0.35) then
    Diff_frac = 1.-6.4*(Atmo_Tran-0.22)**2 
else
    Diff_frac = 1.47-1.66*Atmo_Tran
endif

 ! Average incoming PAR and NIR within a particular hour (say, 5:00 - 6:00)
 ! Assuming the total solar radiation include 50% PAR and 50% NIR
 ! NIR is not for photosynthesis, but for EP

PAR = 0.5*totR*Frac*Atmo_Tran/3600   !                         (J/m2/s)
NIR = 0.5*totR*Frac*Atmo_Tran/3600

 ! ================== de Wit hourly temperature calculation ================
 
Do t = 1 , 24
    Tss = Tmin + (Tmax - Tmin)*sin(pi*(h_SS - 12 + DayLen/2)/(DayLen + 3))
    if (t >= h_SR .and. t <= h_SS) then
        T_air(t) = Tmin + (Tmax - Tmin)*sin(pi*(t - 12 + DayLen/2)/(DayLen + 3))
        
    elseif (t > h_SS) then
        
        T_air(t) = (Tmin - Tss*exp(-(24-DayLen)/4) + (Tss - Tmin) *exp(-(t-h_SS)/4))/&
            (1- exp(-(24- DayLen)/4))
    elseif (t < h_SR)  then
        T_air(t) = (Tmin - Tss*exp(-(24-DayLen)/4)+ (Tss-Tmin)*exp(-(t+h_SR)/4))/&
            (1- exp(-(24- DayLen)/4))
    endif
    
EndDo

! =========================== vapor pressure ===========================
! SVP = 0.611 * exp(17.4 * T_air ./ (239 + T_air));     % kPa
AVP = 0.611 * exp(17.4 * Tmin / (239 + Tmin))
   
! *************************** Biomass Accumulation Begin **********************
!    by hour
! *****************************************************************************
Do t = 1 , 24
    
    ! ================= (1). Weather condition by hour =============================
    CosTheta = sin(Latr)*sin(SolarDec)+ cos(Latr)*cos(SolarDec)*cos(pi*(t - 12)/12)
    
    if (t >= h_SR .and. t <= h_SS .and. CosTheta > 1.E-4) then
        
      Fdiff = max(Diff_frac,0.15+0.85*(1.- exp(-0.1/CosTheta)))
        
      ! ---incoming diffuse PAR (PARDF) and direct PAR (PARDR)
        PAR_diff = PAR(t) * Fdiff                     ! (J/m2/s)
        PAR_dir = PAR(t) - PAR_diff
      ! *---incoming diffuse NIR (NIRDF) and direct NIR (NIRDR)
        NIR_diff = NIR(t) * Fdiff 
        NIR_dir = NIR(t) - NIR_diff
        
        call Rad_TwoStream_Dai(lai, chil)
         
        absPAR_Sun = PAR_dir * Idir_sun(1) + PAR_diff* Idif_sun(1)
        absPAR_Sha = PAR_diff* Idif_sha(1)
        
        absRad_Sun = absPAR_Sun + NIR_dir * Idir_sun(2) + NIR_diff * Idif_sun(2)
        absRad_Sha = absPAR_Sha + NIR_diff* Idif_sha(2)
        
        absRad_tot = absRad_Sun + absRad_Sha
        
      ! =================== (2). CO2, H2O and Heat resistance =======================
        
        Kkw = 0.5491  ! according to CSD , leaf angle 60,  LAI = 6, scatter = 0.2
        
        !*---fraction of sunlit and shaded components in canopy
        f_sun = 1./Kb/lai*(1.-exp(-Kb*lai))
        f_sha = 1.-f_sun
 
        ! --- scaling-up coefficients from leaf to canopy--2014
        ! cintsha(1) = (1.-exp(-extkn*lai))/extkn - cintsun(1)
        ! Eq. 37a, 37b, 38a, 38b  
        Lf2Cp_sun(1) = (1.-exp(-(Kn + Kb)*lai))/(Kn+ Kb)  ! 0.11 was used by Bonan
        Lf2Cp_sun(2) = (1.-exp(-(Kb+Kd)*lai))/(Kb+Kd)
        Lf2Cp_sun(3) = (1.-exp(-Kb*lai))/Kb
        
        Lf2Cp_sha(1) = (1.-exp(-Kn*lai))/Kn - Lf2Cp_sun(1)
        Lf2Cp_sha(2) = (1.-exp(-Kd*lai))/Kd - Lf2Cp_sun(2)
        Lf2Cp_sha(3) = lai - Lf2Cp_sun(3)
 
        !---boundary layer resistance for canopy, sunlit and shaded leaves
        gb_lf = 0.01*sqrt(WND/Leaf_width)  !Leaf boundary layer conductance for heat transfer  m/s
        !---Canopy boundary layer conductance for heat transfer
        gb_cp = (1.-exp(- 0.5*Kkw *lai))/(0.5*Kkw )*gb_lf  
        ! Kkw: wind extinction coeff in the canopy  ?????????????????
        
        gb_Sun = (1.-exp(-(0.5*Kkw+Kb)*lai))/(0.5*Kkw+Kb)*gb_lf
        gb_Sha = gb_cp - gb_Sun
    
        rb_Sun = 1./gb_Sun     !boundary layer resistance to heat,sunlit part
        rbw_Sun = 0.93*rb_Sun  !boundary layer resistance to H2O, sunlit part
        rb_Sha = 1./gb_Sha     !boundary layer resistance to heat,shaded part
        rbw_Sha = 0.93*rb_Sha  !boundary layer resistance to H2O, shaded part
    
        ! ==================== leaf N contents====================================
        ! N is not considered becasue it was included in Vmax25
        
        ! ================(3) Instaneious Photosynthesis =========================
        ! cal photosynthesis and respiration for sunlit and shaded leaves
        !!! Estimation of leaf temperature is through coupling Ep with A-gs photosynthesis model
        ! However, the CSD method require vapor pressure input, 
        ! the CLM method requires humidity input.
 
        ! --------------------- first round --------------------------------
        ! using Ta to approximate T_lf
        T_lf = T_air(t)
        
        ! **********************  3.1 cal T_lf, CO2i, Ass_n (CSD)*********************
        ! ref: P13 Eq 2, 6, 7, 8  in CSD book
        call cal_CO2i(C3C4, CiCa_ratio, CO2a, O2, SVPL1, CO2i1)
 
        ! ----------------------------3.2 Photosynthesis Rate ------------------------
        !  Ref Dai. 2004 and CLM                       
        
        call photo_rate(C3C4,Vcmax25,T_lf,Lf2Cp_sun,absPAR_Sun,O2,CO2i1,Photo1_Sun)       
        call photo_rate(C3C4,Vcmax25,T_lf,Lf2Cp_sha,absPAR_Sha,O2,CO2i1,Photo1_Sha)
        
         vVPD = max(0., SVPL1 - AVP)
        ! -------------------------------- 3.3 Conductance -----------------------
        ! cal s (Eq.5) in Eq. (2) P11 CSD book
        ! the first round: s is the deviration of saturated VP of T_air
         SlopeL1 = 4158.6 * SVPL1/(T_lf + 239.)**2 
      
        ! --- cal leaf conductance for CO2 (gc) and the stomatal resistance to water (r_swp)
        r_swp1_Sun = r_swp(Photo1_Sun, T_lf, CO2a, CO2i1, rbw_Sun)
        r_swp1_Sha = r_swp(Photo1_Sha, T_lf, CO2a, CO2i1, rbw_Sha)

        !----- cal leaf transpiration using Penman-Monteith equation--------------
        !---net absorbed radiation and Ep
        call cal_EP(r_swp1_Sun, rbw_Sun, rb_Sun, absRad_Sun,Atmo_Tran,f_sun, T_air(t), AVP, SlopeL1, VPD, Ep_Sun1, NetRd_Sun1)
        call cal_EP(r_swp1_Sha, rbw_Sha, rb_Sha, absRad_Sha,Atmo_Tran,f_sha, T_air(t), AVP, SlopeL1, VPD, Ep_Sha1, NetRd_Sha1)

        ! ---- cal the difference between leaf temperature and air temperature
         del_t_Sun = del_T(NetRd_Sun1, Ep_Sun1, rb_Sun)
         del_t_Sha = del_T(NetRd_Sha1, Ep_Sha1, rb_Sha) 

         T_lf_Sun = T_air(t) + del_t_Sun
         T_lf_Sha = T_air(t) + del_t_Sha

         !============ Second round to determine the photosynthesis and T_leaf ====
         ! cal internal CO2 pressure first

         call cal_CO2i(C3C4, CiCa_ratio, CO2a, O2, SVPL2_Sun, CO2i2_Sun)
         call cal_CO2i(C3C4, CiCa_ratio, CO2a, O2, SVPL2_Sha, CO2i2_Sha)

         call photo_rate(C3C4,Vcmax25,T_lf_Sun, Lf2Cp_sun,absPAR_Sun,O2,CO2i2_Sun,Photo2_Sun)
         call photo_rate(C3C4,Vcmax25,T_lf_Sha, Lf2Cp_sha,absPAR_Sha,O2,CO2i2_Sha,Photo2_Sha)

          SVP = 0.611 * exp(17.4* T_air(t)/(T_air(t) + 239.0)) ! kPa    saturated VPD of air

          if (ABS(T_lf_Sun-T_air(t)) > 1.0E-8 .and. ABS(T_lf_Sha - T_air(t)) > 1.0E-8) then
          Slope_lf2_Sun = (SVPL2_Sun - SVP)/(T_lf_Sun - T_air(t))! slope of saturated VPD
          Slope_lf2_Sha = (SVPL2_Sha - SVP)/(T_lf_Sha - T_air(t))
          endif

          !--- cal leaf conductance for CO2 (gc) and the stomatal resistance to water (r_swp)----- 
          r_swp2_Sun = r_swp(Photo2_Sun, T_lf_Sun, CO2a, CO2i2_Sun, rbw_Sun) 
          r_swp2_Sha = r_swp(Photo2_Sha, T_lf_Sha, CO2a, CO2i2_Sha, rbw_Sha) 

          !--- cal leaf transpiration using Penman-Monteith equation--------------
          !---net absorbed radiation
          call cal_EP(r_swp2_Sun, rbw_Sun, rb_Sun, absRad_Sun,Atmo_Tran, f_sun, T_lf_Sun, AVP, Slope_lf2_Sun, &
                      vVPD,Ep2_Sun, NetRd2_Sun)
          call cal_EP(r_swp2_Sha, rbw_Sha, rb_Sha, absRad_Sha,Atmo_Tran, f_sha, T_lf_Sha, AVP, Slope_lf2_Sha, &
                      vVPD,Ep2_Sha, NetRd2_Sha)

          ! ***************** End of T_lf, An, CO2i *********************************

           InsPP(t) = (Photo2_Sun + Photo2_Sha)*3600   ! g CO2 /m2/hr    daily net accumulation

           InsEP(t) = (Ep2_Sun + Ep2_Sha)*3600
    !else    ! at night
    !       qt = (t_air(t)+273 - 298.15)/10  ! q10   
    !       rt = 298.15 * rgas
    !       rd25 = (0.015 * c3 + 0.025 * c4) * vcmax25
    !       
    !       ft_r = exp(ea_rdmax*(1-298.15/(t_air(t)+273))/rt)
    !       fh_r = (1 + exp((298.15*s_rdmax - ed_rdmax)/rt))/(1 + exp((t_air(t)+273)*s_rdmax - ed_rdmax)/rt)
    !
    !       respir(t) = rd25 * ft_r * fh_r * c3 + rd25*2.0**qt/(1 + exp(s5*((t_air(t)+273) - s6)))*c4
    !
    !       respir(t) = respir(t) * lai* 3600*44*1e-6   ! g co2 /m2/hr              nightly consumption
    !    
    ENDIF
 
    End Do
    
   Bio = sum(InsPP)*0.409* 0.01  ! g CO2 /m2 --> g Biomass/m2 --> t Biomass/ha   0.409 might be like HI
   potentialBioIncrese(Crop_Num)=max(1.E-5,Bio)
   totCropBio(Crop_Num)=totCropBio(Crop_Num)+potentialBioIncrese(Crop_Num)
   ! WRITE(KW(51),"(1X, I4, 1X, I4,  3F8.3 )") IYR, Date, absPAR_Sun + absPAR_Sha, Bio, totCropBio(Crop_Num)
 !  print *, "Daily biomass accumulation:"
 !  print *, Bio, ' g /m2/day'
 !  pause
 
 End Subroutine Poten_Growth
    
    
END MODULE Bio_Acc 