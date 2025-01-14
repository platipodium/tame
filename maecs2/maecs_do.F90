!#include "fabm_driver.h"

!> @brief This is the main routine where right-hand-sides are calculated
!> @details
!> NOTE: Although this subroutine looks as if it's not a part of any module,
!! it is temporarily included in the fabm_hzg_maecs module (inside maecs.F90)
!! when compiling the documentation, such that the subroutine is documented
!! under the 'Data Type Documentation' chapter, where the in-body docs are also listed
!>
!> **Phytoplankton Equations**
!> \n We distinguish between mass state variables
!! (in units of carbon, nitrogen, & phosphorus) and property state variables.
!! \lref{For a textual narration and equations, see sec.,sec:ModStr,.}\n
!>  Current 'traits' are:
!> - nitrogen allocated to rubisco [-] (frac_Rub)
!> - Chla content of chloroplasts [chl-a/chl-C] (theta)
!>
!> **General code structure:**
!> 1. Calculation of quotas, internal states, potential rates
!> 2. Calculation of fluxes, mass exchange rates &  rates of change of traits variables
!> 3. Assign mass exchange rates ('rhs(j,i)')
!> 4. Assign rates of change of 'traits' property variables
!>
!> \n **Detailed Descriptions:**
! @todo: althougth HIDE_IN_BODY_DOCS=NO, the body-documentation is not included! HAS TO BE FIXED
! @todo: add equations
! @todo: why UNIT instead of secs_pr_day?
! @todo: 'sensitivities' does not seem to be a proper name choice. Something more intuitive?
subroutine maecs_do(self,_ARGUMENTS_DO_)

use fabm_types
use maecs_types
use maecs_functions
use maecs_primprod
use maecs_grazing

! !INPUT PARAMETERS:
 class (type_hzg_maecs),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
type (type_maecs_rhs)    :: rhsv
type (type_maecs_phy)    :: phy   ! phytoplankton type containing state and trait information
type (type_maecs_zoo)    :: zoo   ! zooplankton type
type (type_maecs_om)     :: dom, det, nut, uptake, exud, lossZ, floppZ, nquot
type (type_maecs_env)    :: env
type (type_maecs_switch) :: mswitch
type (type_maecs_traitdyn)::acclim
type (type_maecs_sensitivities) :: sens

! --- LOCAL MODEL VARIABLES:
integer  :: i, j, iz, ihour, iloop
real(rk) :: reminT, degradT       ! Temp dependent remineralisation and hydrolysis rates
! --- QUOTA and FRACTIONS
real(rk) :: phys_status, dQN_dt, dRchl_phyC_dt=0.0_rk ! []

! --- ZOOPLANKTON GRAZING, INGESTION, MORTALITY, RESPIRATION...
real(rk) :: graz_rate   ! carbon-specific grazing rate                          [d^{-1}]
real(rk) :: zoo_respC   ! temperature dependent carbon respiration rate  [mmolC m^{-3} d^{-1}]
real(rk) :: zoo_mort
real(rk) :: decay       ! pigment-specific decay rate                          [d^{-1}]
real(rk) :: denitrate   ! pelagic N-loss by denitrification, emulating benthic pool and suboxic micro-environments
real(rk) :: deporate    ! pelagic, "volumetric" deposition, slowly refueling N-losses
real(rk) :: qualPOM, ddegN, ddegP     !  POM quality -> degradation
real(rk) :: qualDOM, dremN, dremP     !  DOM quality -> degradation
real(rk) :: secs_pr_day = 86400.0_rk
! --- AGGREGATION
real(rk) :: aggreg_rate ! aggregation among phytoplankton and between phytoplankton & detritus [d^{-1}]
real(rk) :: viral_rate ! loss rate due to viral/parasite infections  [d^{-1}]
real(rk) :: vir_lysis = 0.1_rk  !assumes that virally infected cells equally fuel POM and DOM pools
logical  :: out = .true.
!   if(36000.eq.secondsofday .and. mod(julianday,1).eq.0 .and. outn) out=.true.
real(rk) :: pdet, no3
real(rk) :: att, fa, fts, relmort, dq_dt, dQP_dt
real(rk) :: QP_phy_max, rqn
real(rk) :: det_prod, dom_dep, nh3f
real(rk) :: radsP,Oxicminlim,Denitrilim,Anoxiclim,Rescale,rP
real(rk) :: flO2, flODU, aPO4
real(rk),parameter :: relaxO2=0.04_rk
real(rk),parameter :: T0 = 288.15_rk ! reference Temperature fixed to 15 degC
real(rk),parameter :: Q10b = 1.5_rk
real(rk) :: Cprod, Nprod, Pprod
real(rk) :: poc, doy, zmax, sal, add_aggreg_rate
real(rk) :: AnoxicMin,Denitrific,OxicMin,Nitri,OduDepo,OduOx,pDepo, Anammox
real(rk) :: prodO2, rhochl, uptNH4, uptNO3, uptchl, uptN, respphyto,faeces, min_Cmass
real(rk) :: a_lit, ksat_graz, food
real(rk) :: vir_max = 3.0_rk
real(rk) :: vird, dvir_dt, infect, vdilg, vrepl, vadap, vmort, virf, vire
logical  :: IsCritical = .false. ! phyC and phyN below reasonable range ?
#define _KAI_ 2
#define _MARKUS_ 1
! #define debugMK
#define _DEBUG_ 0
! #define UNIT / 86400
#define UNIT *1.1574074074E-5_rk
#define HALFQ 0.5
 _LOOP_BEGIN_

#if _DEBUG_
write(*,'(A)') 'begin DO'
#endif

! First retrieve current (local) state  variable values
!#S_GET
!---------- GET for each state variable ----------
  _GET_(self%id_nutN, nut%N)  ! Dissolved Inorganic Nitrogen DIN in mmol-N/m**3
  _GET_(self%id_phyC, phy%C)  ! Phytplankton Carbon in mmol-C/m**3
  _GET_(self%id_phyN, phy%N)  ! Phytplankton Nitrogen in mmol-N/m**3
  _GET_(self%id_detC, det%C)  ! Detritus Carbon in mmol-C/m**3
  _GET_(self%id_detN, det%N)  ! Detritus Nitrogen in mmol-N/m**3
  _GET_(self%id_domC, dom%C)  ! Dissolved Organic Carbon in mmol-C/m**3
  _GET_(self%id_domN, dom%N)  ! Dissolved Organic Nitrogen in mmol-N/m**3
if (self%RubiscoOn) then
      _GET_(self%id_Rub, phy%Rub)  ! fraction of Rubisco in -
end if
if (self%PhotoacclimOn) then
      _GET_(self%id_chl, phy%chl)  ! Chl in mg-Chl/m**3
end if
if (self%PhosphorusOn) then
      _GET_(self%id_nutP, nut%P)  ! Dissolved Inorganic Phosphorus DIP in mmol-P/m**3
      _GET_(self%id_phyP, phy%P)  ! Phytplankton Phosphorus in mmol-P/m**3
      _GET_(self%id_detP, det%P)  ! Detritus Phosphorus in mmol-P/m**3
      _GET_(self%id_domP, dom%P)  ! Dissolved Organic Phosphorus in mmol-P/m**3
end if
if (self%SiliconOn) then
      _GET_(self%id_nutS, nut%Si)  ! Dissolved Inorganic Silicon Si in mmol-Si/m**3
      _GET_(self%id_phyS, phy%Si)  ! Phytplankton Silicon in mmol-Si/m**3
      _GET_(self%id_detS, det%Si)  ! Detritus Silicon in mmol-Si/m**3
end if
if (self%GrazingOn) then
      _GET_(self%id_zooC, zoo%C)  ! Zooplankton Carbon in mmol-C/m**3
end if
if (self%BioOxyOn) then
      _GET_(self%id_nh3, env%nh3)  ! dissolved ammonium in mmolN/m**3
      _GET_(self%id_oxy, env%oxy)  ! dissolved oxygen in mmolO2/m**3
      _GET_(self%id_odu, env%odu)  ! dissolved reduced substances in mmolO2/m**3
end if
if (self%VirusOn) then
      _GET_(self%id_vir, phy%vir)  ! Virus C density in cells in -
end if
if (self%NResOn) then
      _GET_(self%id_RNit, env%RNit)  ! N-reservoir in mmol-N/m**3
end if
!#E_GET
if (.not. self%GrazingOn) then
      zoo%C = self%zooC_initial
      zoo%N = self%zooC_initial * self%const_NC_zoo
      zoo%P = self%zooC_initial * self%const_PC_zoo
end if

!S_GED
_GET_(self%id_temp, env%temp)  ! water temperature
_GET_(self%id_par, env%par)    ! light photosynthetically active radiation

_SET_DIAGNOSTIC_(self%id_dPAR,env%par)         !average Photosynthetically_Active_Radiation_

if (self%ChemostatOn) then
  if (_AVAILABLE_(self%id_CO2)) then
    _GET_(self%id_CO2, env%CO2)  ! CO2
  else
    env%CO2 = 0.0_rk ! todo: throw an error, if necessary dependency cannot be found
  end if
end if

#ifdef debugMK
if ( nut%N /= nut%N .or. nut%N<0.0_rk) then; write(0,*) 'ERROR: maecs_do#161 nut%N = ',nut%N; stop; endif
if ( phy%C /= phy%C .or. phy%C<0.0_rk) then; write(0,*) 'ERROR: maecs_do#162 phy%C = ',phy%C; stop; endif
if ( phy%N /= phy%N .or. phy%N<0.0_rk) then; write(0,*) 'ERROR: maecs_do#163 phy%N = ',phy%N; stop; endif
if ( det%C /= det%C .or. det%C<0.0_rk) then; write(0,*) 'ERROR: maecs_do#164 det%C = ',det%C; stop; endif
if ( det%N /= det%N .or. det%N<0.0_rk) then; write(0,*) 'ERROR: maecs_do#165 det%N = ',det%N; stop; endif
if ( dom%C /= dom%C .or. dom%C<0.0_rk) then; write(0,*) 'ERROR: maecs_do#166 dom%C = ',dom%C; stop; endif
if ( dom%N /= dom%N .or. dom%N<0.0_rk) then; write(0,*) 'ERROR: maecs_do#167 dom%N = ',dom%N; stop; endif
if (self%RubiscoOn) then
  if ( phy%Rub /= phy%Rub .or. phy%Rub<0.0_rk ) then; write(0,*) 'ERROR: maecs_do#169 phy%Rub = ',phy%Rub; stop; endif
end if
if (self%PhotoacclimOn) then
  if ( phy%chl /= phy%chl .or. phy%chl<0.0_rk ) then; write(0,*) 'ERROR: maecs_do#172 phy%chl = ',phy%chl; stop; endif
end if
if (self%PhosphorusOn) then
  if ( nut%P /= nut%P .or. nut%P<0.0_rk ) then; write(0,*) 'ERROR: maecs_do#175 nut%P = ',nut%P; stop; endif
  if ( phy%P /= phy%P .or. phy%P<0.0_rk ) then; write(0,*) 'ERROR: maecs_do#176 phy%P = ',phy%P; stop; endif
  if ( det%P /= det%P .or. det%P<0.0_rk ) then; write(0,*) 'ERROR: maecs_do#177 det%P = ',det%P; stop; endif
  if ( dom%P /= dom%P .or. dom%P<0.0_rk ) then; write(0,*) 'ERROR: maecs_do#178 dom%P = ',dom%P; stop; endif
end if
if (self%SiliconOn) then
  if ( nut%Si /= nut%Si .or. nut%Si<0.0_rk ) then; write(0,*) 'ERROR: maecs_do#181 nut%Si = ',nut%Si; stop; endif
  if ( phy%Si /= phy%Si .or. phy%Si<0.0_rk ) then; write(0,*) 'ERROR: maecs_do#182 phy%Si = ',phy%Si; stop; endif
  if ( det%Si /= det%Si .or. det%Si<0.0_rk ) then; write(0,*) 'ERROR: maecs_do#183 det%Si = ',det%Si; stop; endif
end if
if (self%GrazingOn) then
  if ( zoo%C /= zoo%C .or. zoo%C<0.0_rk ) then; write(0,*) 'ERROR: maecs_do#186 zoo%C = ',zoo%C; stop; endif
end if
if (self%BioOxyOn) then
  if ( env%nh3 /= env%nh3 ) then; write(0,*) 'ERROR: maecs_do#189 env%nh3 = ',env%nh3; stop; endif
  if ( env%oxy /= env%oxy ) then; write(0,*) 'ERROR: maecs_do#190 env%oxy = ',env%oxy; stop; endif
  if ( env%odu /= env%odu ) then; write(0,*) 'ERROR: maecs_do#191 env%odu = ',env%odu; stop; endif
end if
if (self%VirusOn) then
  if ( phy%vir /= phy%vir .or. phy%vir<0.0_rk ) then; write(0,*) 'ERROR: maecs_do#194 phy%vir = ',phy%vir; stop; endif
end if
if (self%NResOn) then
  if ( env%RNit /= env%RNit ) then; write(0,*) 'ERROR: maecs_do#197 env%RNit = ',env%RNit; stop; endif
end if
if ( env%temp /= env%temp ) then; write(0,*) 'ERROR: maecs_do#199 env%temp = ',env%temp; stop; endif
if ( env%par /= env%par   ) then; write(0,*) 'ERROR: maecs_do#200 env%par = ',env%par; stop; endif
#endif


!E_GED  ! list outcommented due to different usage of zmax and doy (see light extinction)

! @ingroup main
!> @fn fabm_hzg_maecs::maecs_do ()
!> 1. Calculation of quotas, internal states, potential rates
!>   - call min_mass with method=_KAI_, store phy\%C and \%N in phy\%reg
!>   - call calc_internal_states: retrieve phy\%Q\%X, phy\%theta, phy\%frac\%X
!>   - if PhotoacclimOn=.false., calculate:
!>     - phy\%chl=phy\%C * self\%frac_chl_ini
!>     - phy\%frac\%theta = self\%frac_chl_ini * self\%itheta_max
!>     - phy\%theta= self\%frac_chl_ini / (self\%frac_Rub_ini * phy\%relQ\%N**self\%sigma)
!>   - call calc_sensitivities: retrieve potential rates: @f$f_T@f$, sens\%upt\_pot\%C (=LH), sens\%upt\_pot\%X (= @f$ V_X @f$), sens\%P\_max
!> @todo: min_mass correction of phy%\C and phy\%N at this stage requires specification of threshold values. What about back-calculating phy\%reg\%N from the smooth_small corrected phy\%Q\%N?

! --- checking and correcting extremely low state values  ------------
call min_mass(self,phy, min_Cmass, IsCritical, method=_KAI_) ! minimal reasonable Phy-C and -Nitrogen

if(self%maxVal .gt. 0.0d0) then
  if(phy%chl .gt. self%maxVal .or. phy%Rub .gt. self%maxVal) IsCritical=.true.
endif
if(IsCritical .and. env%par*self%frac_PAR * self%alpha .gt. 0.5d0) IsCritical=.false.

if(IsCritical .and. .not. self%ChemostatOn ) then
!if(IsCritical .or. (self%ChemostatOn .and. (nut%P .lt. self%small_finite .or. nut%N .lt. self%small_finite) )) then
  rhsv%nutN=0.0d0
  rhsv%nutP=0.0d0
  rhsv%nutS=0.0d0
  rhsv%phyC=0.0d0
  rhsv%phyN=0.0d0
  rhsv%phyP=0.0d0
  rhsv%phyS=0.0d0
  rhsv%zooC=0.0d0
  rhsv%detC=0.0d0
  rhsv%detN=0.0d0
  rhsv%detP=0.0d0
  rhsv%detS=0.0d0
  rhsv%domC=0.0d0
  rhsv%domN=0.0d0
  rhsv%domP=0.0d0
  rhsv%RNit=0.0d0
  rhsv%Rub=0.0d0
  rhsv%chl=0.0d0
  rhsv%nh3=0.0d0
  rhsv%oxy=0.0d0
  rhsv%odu=0.0d0
  rhsv%vir=0.0d0
else
! --- stoichiometry of autotrophs (calculating QN_phy, frac_R, theta, and QP_phy)
call calc_internal_states(self,phy,det,dom,zoo)

! trait-hack: QmxaP ~ QN
if(self%mort_ODU .gt. 0.99) then
  rqn  = phy%relQ%N
  if(rqn .gt. 1.0) rqn = 1.0
  if(rqn .lt. 0.05) rqn = 0.05

  QP_phy_max = self%QP_phy_0+ rqn*(self%QP_phy_max*self%QP_phy_0)
  phy%relQ%P = ( phy%Q%P - self%QP_phy_0 )/( rqn * self%QP_phy_max)
endif


if (.not. self%PhotoacclimOn) then
   phy%chl         = phy%C * self%frac_chl_ini   ! total Chl mg-CHL/m3
   phy%frac%theta  = self%frac_chl_ini * self%itheta_max
   phy%theta       = self%frac_chl_ini /(self%frac_Rub_ini*phy%relQ%N**self%sigma)
! g-CHL/mol-C*m3
! write (*,'(A,2(F10.3))') 'theta:',phy%relQ%N**self%sigma,phy%theta
end if

call calc_sensitivities(self,sens,phy,env,nut,acclim)

!if (IsCritical .and. .false. ) then
!  write (*,'(A,4(F10.3))') 'fR=',phy%reg%C,phy%frac%Rub,self%small_finite + self%rel_chloropl_min,phy%frac%NutUpt
!end if
!if (phy%chl .lt. 0.01d0) then
!  phy%theta     = phy%chl / (phy%rel_chloropl * phy%reg%C)   ! trait variable
!  phy%frac%theta= phy%theta * phy%rel_chloropl * self%itheta_max ! []     no
!  write (*,'(A,3(F10.3))') 'fT=',phy%chl,phy%reg%C,phy%Rub end if
!write (*,'(A,4(F10.3))') 'PAR, T, th, P =',env%par,env%temp,phy%theta, sens%upt_pot%C

!> @fn fabm_hzg_maecs::maecs_do ()
!> 2. Calculation of fluxes, mass exchange rates &  rates of change of traits variables
!! & Specify rates of change of traits variables
!>   - call maecs_primprod::photosynthesis(): this is where everything happens!
!>   - if GrazingOn:
!>     - graz\_rate=rate retrieved from call maecs_grazing::grazing()
!>     - lossZ\%X=lossZNut\%X, floppZ\%X=lossZDet\%X retrieved from call maecs_grazing::grazing_losses()
!>     - calculate graz_rate retr= graz_rate * zoo\%C and zoo_mort
!>   - calc. aggreg_rate @f$ = \mathrm{phi\_agg}*(1-e^{-0.02*\mathrm{dom\%C}}) * phy\%N *det\%N @f$
!>     - future work: aggreg_rate=f(size), \latexonly (see section \ref{sec:partagg}) \endlatexonly \n
!>   - calc. degradT=self\%hydrol * @f$ f_T @f$ and reminT=self\%remin * @f$ f_T @f$
!> @todo: aggregation equation: where does 0.02 come from? are the results sensitive to this par?
!> @todo: specific graz_rate becomes pop. grazing rate. Do this at the rhs calculations
!> @todo: graz_rate: no temperature modification: forgotten?

! --- ALGAL GROWTH and EXUDATION RATES, physiological trait dynamics ----------------
#ifdef debugMK
   if ( exud%C/= exud%C .or. abs(exud%C)>huge(0.0_rk)) then
     write(0,*) 'ERROR: maecs_do#302 INF/NaN detected exud%C = ',exud%C; stop;
   endif
   if ( exud%N/= exud%N .or. abs(exud%N)>huge(0.0_rk)) then
     write(0,*) 'ERROR: maecs_do#305 INF/NaN detected exud%N = ',exud%N; stop;
   endif
   if ( self%PhosphorusOn ) then
     if ( exud%P/= exud%P .or. abs(exud%P)>huge(0.0_rk) ) then
       write(0,*) 'ERROR: maecs_do#309 INF/NaN detected exud%P = ',exud%P; stop;
     endif
   endif
#endif
call photosynthesis(self,sens,phy,nut,uptake,exud,acclim)

! ----------------       grazing        -------------------------------------------
if (self%GrazingOn) then

  ksat_graz = self%k_grazC

!  --- quadratic closure term
  relmort = 1.0d0
  fts = 1.0d0 ! relevance of microoo grazing
  _GET_(self%id_attpar, att)
  if (self%GrazTurbOn .ge. 0) then
!    _GET_(self%id_attf_dep, att_f)
!  _GET_GLOBAL_ (self%id_doy,doy) !day of year
   _GET_(self%id_attpar, att)
   _SET_DIAGNOSTIC_(self%id_datt,att)
   _SET_DIAGNOSTIC_(self%id_vphys, -1)
   select case (self%GrazTurbOn)
     case (0)
      _GET_GLOBAL_ (self%id_doy,doy) !day of year
      _GET_HORIZONTAL_(self%id_zmax, zmax)  ! max depth
       fa = self%zm_fa_inf +(1.0_rk-self%zm_fa_inf)/(1+exp(0.16*(zmax-24.0)))
       relmort = fa + fa*self%zm_fa_delmax*sens%f_T2*0.25*(1-sin(2*(doy+45)*Pi/365.0))**2
       ksat_graz = fa * self%k_grazC
     case (1)
       fa = 1.0_rk/(1.0_rk+exp(2*(att-self%zm_fa_inf)))
       relmort = 1.0_rk + self%zm_fa_delmax* fa
     case (2)
       relmort = 1.0d0 + self%zm_fa_delmax/(att+self%zm_fa_inf)
     case (3,4)
       relmort = 1.0d0 + sens%f_T2*self%zm_fa_delmax/(att+self%zm_fa_inf) ! assumes greater fish/larvae abundance in summer
     case (5,6)
      _GET_GLOBAL_ (self%id_doy,doy) !day of year
       relmort=1.0d0 + self%zm_fa_delmax*sens%f_T2*0.5*(1-sin(2*(doy+75)*Pi/365.0))
     case (7)
      _GET_GLOBAL_ (self%id_doy,doy) !day of year
      _GET_HORIZONTAL_(self%id_zmax, zmax)  ! max depth
    !f(z)=sigmoidal function of depth with an upper plateau (100%) at 0-10 m and a lower (10%) for 30+
       fz=0.5/(1+exp(-zmax*0.5_rk+self%a_fz))
       relmort=1.0d0 + self%zm_fa_delmax*sens%f_T2*0.5*(1-fz*sin(2*(doy+75)*Pi/365.0))
     case (8)
      _GET_GLOBAL_ (self%id_doy,doy) !day of year
      _GET_(self%id_sal,sal) ! salinity
      _GET_HORIZONTAL_(self%id_zmax, zmax)  ! max depth
       fts = 1._rk/(1+exp(self%zm_fa_inf*(sal+self%mort_ODU)))
       fa = 0.01_rk+1.0_rk*max(fts, 1.0_rk/(1+exp(0.4*self%zm_fa_inf*(zmax-25.0))))
 !      if(sal .lt. 0.0_rk) sal = 0.0_rk
 !      if(sal .gt. 40.0_rk) sal = 40.0_rk
       if (self%zm_fa_inf .gt. 1E-4) fa =  fa + 1.0_rk/(1+exp(self%zm_fa_inf*(sal-12)))
 ! fa =  fa + 0.4_rk/(1+exp(self%zm_fa_inf*(sal-12)))
!   fts = self%zm_fa_delmax*sens%f_T2*0.25*(1-sin(2*(doy+25)*Pi/365.0))**2 ! seasonal increase in top-predation
       fts = self%zm_fa_delmax*sens%f_T2*0.5*(1-sin(2*(doy+0)*Pi/365.0)) ! seasonal increase in top-predation
!       relmort = fa *(1+fts) + 0.5*fts
       relmort = fa*(0.1*zoo%C+fts)
!       ksat_graz = (0.5*(1+fa)*(1+fts)+0.0_rk) * self%k_grazC
       ksat_graz = (fa+0.5_rk) * self%k_grazC
       _SET_DIAGNOSTIC_(self%id_vphys, fa)       !average Temporary_diagnostic_
! relevance of microzoo/mixotrophy grazing increases at low DIP
!      fts = 1.0d0 + 0.5_rk/(1+exp(10*(nut%P - 0.3_rk)))
    end select
  end if !self%GrazTurbOn .gt. 0
  zoo_mort   = self%mort_zoo * relmort* sens%f_T**self%fT_exp_mort ! * zoo%C
  if (self%GrazTurbOn .eq. 4 .or. self%GrazTurbOn .gt. 5) zoo_mort   = zoo_mort + self%mort_zoo

!  if (self%GrazTurbOn .eq. 0)  zoo_mort   = zoo_mort*
!!  write (*,'(A,4(F11.3))') 'Zm=',att,relmort,zoo%C,zoo_mort

! feeding threshold and grazer disturbance
  if (phy%C .lt. 20*det%C/(20+det%C) ) then
    food=0.0_rk
  else
    food=phy%C
  end if
! calls grazing function accounting for homeostasis
  call grazing(self%g_max * sens%f_T2,ksat_graz,food,graz_rate) ! * fts
  zoo%feeding = graz_rate
  zoo_respC   = self%basal_resp_zoo * sens%f_T  !  basal respiration of grazers
  nquot       = type_maecs_om(1.0_rk, phy%Q%N, phy%Q%P, phy%Q%Si )
  mswitch     = type_maecs_switch(self%PhosphorusOn,self%SiliconOn,.true. )

! --- calculates zooplankton loss rates (excretion->Nut, floppy+egestion->Det), specific to C
  call grazing_losses(zoo,zoo_respC,nquot,lossZ,floppZ, mswitch)
!  --- transform from specific to bulk grazing rate
  graz_rate   = graz_rate * zoo%C

                                  !isP, isSi, isTotIng

else
  graz_rate   = 0.0_rk
!  if (self%ChemostatOn .and. .not. IsCritical) graz_rate = 0.2*phy%C
  zoo_mort    = 0.0_rk
  lossZ       = type_maecs_om(0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk)
  floppZ      = type_maecs_om(0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk)
end if


! --- phytoplankton aggregation -------------------------------------------------
! If biovolume is primarily determined by the nitrogen content, also for detritus
!aggreg_rate = self%phi_agg * dom%C * (phy%N + det%N)                    ! [d^{-1}]

!_GET_(self%id_fracR,phys_status )

dom_dep     = self%agg_doc*dom%C/(1.0_rk+self%agg_doc*dom%C)
aggreg_rate = min(self%phi_agg * dom_dep * (phy%N + det%N),0.1_rk)
!         vS * exp(-4*phys_status )                ! [d^{-1}]
!aggreg_rate = aggreg_rate * exp(-4*phy%rel_phys ) TODO: DOM quality as proxy for TEP

! additional mortality due to H2S stress (EC_50 :55 mmol-H2S/m3; Kuester2005  or for Skel cost 3 mmol-H2S/m3, Breteler1991)
if (self%BioOxyOn) then
   aggreg_rate = aggreg_rate + self%mort_ODU* env%odu
endif
#ifdef debugMK
   if ( aggreg_rate/= aggreg_rate .or. abs(aggreg_rate)>huge(0.0_rk) .or. aggreg_rate < 0.0_rk ) then
     write(0,*) 'ERROR: dom%C        = ',dom%C
     write(0,*) 'ERROR: self%agg_doc = ',self%agg_doc
     write(0,*) 'ERROR: dom_dep      = ',dom_dep
     write(0,*) 'ERROR: self%phi_agg = ',self%phi_agg
     write(0,*) 'ERROR: phy%N        = ',phy%N
     write(0,*) 'ERROR: det%N        = ',det%N
     if (self%BioOxyOn) then
       write(0,*) 'ERROR: self%mort_ODU = ',self%mort_ODU
       write(0,*) 'ERROR: env%odu       = ',env%odu
     endif
     !-----------------------------------------------------------
     write(0,*) 'ERROR: maecs_do#446 INF/NaN/below_zero detected aggreg_rate = ',aggreg_rate; stop;
!      write(0,'(''ERROR: maecs_do#'',i4.0,'' INF/NaN/below_zero detected aggreg_rate = '',e20.10)') &
!                    __LINE__, aggreg_rate; stop;
   endif
#endif

! --- phytoplankton viral losses --------------------------------
if (self%vir_loss .gt. self%small_finite .or. self%VirusOn ) then
  if (self%VirusOn ) then
    ! ---  dynamic approach --------------------------------
    vird = phy%vir/phy%reg%C   ! density relative to host biomass
!    viral_rate = self%vir_loss * vird
! non-linear impact on host mortality with threshold
   ! self%vir_infect: steepness factor
    vire = exp(self%vir_infect*(1.0_rk-vird))
    virf = 1.0_rk/(1.0_rk + vire)     ! smooth step function
    viral_rate = self%vir_loss * virf                  ! host mortality
  else
    ! ---  static approach --------------------------------
    if (uptake%C .gt. self%small_finite .and. self%vir_mu*sens%f_T .gt. self%vir_loss) then  ! infection only at daytime
      viral_rate  = 0.1 - (uptake%C - infect)/(self%vir_mu*sens%f_T-self%vir_loss)
    ! TODO: update with current dyn version
    ! TODO: check for netagive rates!
    else
      viral_rate  = 0.0_rk
    endif
    viral_rate  = viral_rate * self%vir_loss
  endif
  ! if (viral_rate .lt. 0.0_rk) viral_rate = 0.0_rk
  ! a faction of virally infected cells goes into dead material
  aggreg_rate = aggreg_rate + (1.0_rk-vir_lysis)* viral_rate

  ! the remainder is lysis (for all elements)
  exud%C      = exud%C + vir_lysis* viral_rate
  exud%N      = exud%N + vir_lysis* viral_rate * phy%Q%N ! the remainder is lysis.N
  if (self%PhosphorusOn) then
    exud%P     = exud%P + vir_lysis* viral_rate * phy%Q%P ! the remainder is lysis.P
  endif
endif
#ifdef debugMK
   if ( exud%C/= exud%C .or. abs(exud%C)>huge(0.0_rk)) then
     write(0,*) 'ERROR: vir_lysis   = ',vir_lysis
     write(0,*) 'ERROR: viral_rate  = ',viral_rate
     !-----------------------------------------------------------
     write(0,*) 'ERROR: maecs_do#488 INF/NaN detected exud%C = ',exud%C; stop;
   endif
#endif
!_____________________________________________________________________________
!
!      turnover of long-term N-reservoir (denitrification + wet N-deposition)
if (self%NResOn) then
!pelagic N-loss by denitrification, emulating benthic pool and suboxic micro-environments
 denitrate = self%denit * 4 * sens%f_T * (1.0d0 - exp(-det%N/self%PON_denit)) * det%N

!pelagic, "volumetric" deposition, slowly refueling N-losses
 deporate  = self%denit * exp(-4*sens%f_T) * env%RNit

 rhsv%RNit = denitrate - deporate
else
  deporate  = 0.0d0
  denitrate = 0.0d0
endif


!> @fn fabm_hzg_maecs::maecs_do ()
!> 3. Assign mass exchange rates ('rhs(j,i)')
!>   - phyC= uptake - dil - exud - aggreg\_rate - graz\_rate
!> @todo: add the rhs equations


! right hand side of ODE (rhs)
!__________________________________________________________________________
!
! PHYTOPLANKTON C
rhsv%phyC = uptake%C              * phy%C &
           - self%dil             * phy%C &
           - exud%C               * phy%C &  !TODO: move loss rates to mu, also checking for
           - aggreg_rate          * phy%C &  !      trait dependencies
           - graz_rate
! write (*,'(A,3(F9.4))') 'flxc=',uptake%C, phy%C,nut%P
! write (*,'(A,3(F9.4))') 'c=',phy%chl,phy%Rub, phy%C
#ifdef debugMK
   if ( rhsv%phyC/= rhsv%phyC .or. abs(rhsv%phyC)>huge(0.0_rk)) then
     write(0,*) 'ERROR: uptake%C    = ',uptake%C
     write(0,*) 'ERROR: exud%C      = ',exud%C
     write(0,*) 'ERROR: aggreg_rate = ',aggreg_rate
     write(0,*) 'ERROR: graz_rate   = ',graz_rate
     !-----------------------------------------------------------
     write(0,*) 'ERROR: maecs_do#519 INF/NaN detected rhsv%phyC = ',rhsv%phyC; stop;
   endif
#endif

!_____________________________________________________________________________
!
! PHYTOPLANKTON N
rhsv%phyN =  uptake%N             * phy%C &
           - exud%N               * phy%C &
           - aggreg_rate          * phy%N &
           - self%dil             * phy%N &
           - graz_rate * phy%Q%N       ! Q_N is bounded for safety reasons: TODO change grazing units to N
!rhsv%phyN = 0.0_rk
#ifdef debugMK
   if ( rhsv%phyN/= rhsv%phyN .or. abs(rhsv%phyN)>huge(0.0_rk)) then
     write(0,*) 'ERROR: uptake%N    = ',uptake%N
     write(0,*) 'ERROR: exud%N      = ',exud%N
     write(0,*) 'ERROR: aggreg_rate = ',aggreg_rate
     write(0,*) 'ERROR: graz_rate   = ',graz_rate
     write(0,*) 'ERROR: phy%Q%N     = ',phy%Q%N
     !-----------------------------------------------------------
     write(0,*) 'ERROR: maecs_do#540 INF/NaN detected rhsv%phyN = ',rhsv%phyN; stop;
   endif
#endif

!_____________________________________________________________________________

!> @fn fabm_hzg_maecs::maecs_do ()
!> 4. Assign rates of change of 'traits' property variables
!>    - if PhotoacclimOn: rhsv%chl=A + B
!>      + A = rhsv\%phyC*phy\%theta*rel_chlorpl [= gC/m3/d * gchl/gchlorpl-C * gchloropl-C/gC]
!>      + B = dRchl/dtheta* dtheta/dt + dRchl/dfracR*dfracR/dt+dRchl/dQN * dQN/dt
!>      + all terms in B except dQN/dt are calculated in maecs_primprod::photosynthesis()
!>      + dQN/dt = (rhsv\%phyN* phyC - rhsv\*phyC*phyN )/(phyC^2)
!>    - if RubiscoOn: rhsv%Rub=A+B
!>      + A = rhsv\%phyC * phyRub/phyC
!>      + B = dfracR_dt is calculated in maecs_primprod::photosynthesis()
! viral density is treated like a trait (because specific to PhyC)
 if (self%VirusOn ) then
   poc   = zoo%C + 3*det%C + phy%reg%C
!!  poc   = dom%P + det%P + phy%P + self%small_finite
!!  poc   = dom%C + det%C + phy%C + self%small_finite
           ! encounter prob. free virus conc around infected cell
           ! average distance ~ C^-1/3 + Gaussian/diffusive spots
!  infect  = infect * exp(-1.0_rk/((phy%reg%C/10)**0.667_rk) )
!  vdilg = 0*(uptake%C - 0*exud%C)      ! dilution of viral concentration by new biomass production
!  if (vdilg .lt. 0.0_rk) vdilg = 0.0_rk

! cross section by non-algae particles (bacteria+colloids)

! a_lit = att - (self%a_spm*(det%C+zoo%C)+ self%a_doc*dom%C + self%a_phyc*phy%C + self%a_chl*phy%chl)
!  a_lit = abs(a_lit/(self%a_phyc + self%small_finite))
 ! viral replication
 if (self%vir_mu .gt. 0.0_rk ) then
  vrepl = self%vir_mu *phy%relQ%N**2/(HALFQ**2+phy%relQ%N**2)! non-linear dependence on stoichiometry
  if (self%PhosphorusOn) vrepl = vrepl * phy%relQ%P**2/(HALFQ**2+phy%relQ%P**2)
      ! depends on host P-stoichiometry (Wilson et al 1996, Clasen&Elser 2007)
 else
!  vrepl = -self%vir_mu *  phy%relQ%N ! linear dependence on stoichiometry
  !if (self%PhosphorusOn) vrepl = vrepl * phy%relQ%P
  vrepl = -self%vir_mu/(1.0_rk+ exp(-4.5*(phy%relQ%N-1.0_rk)))   ! linear dependence on stoichiometry
  if (self%PhosphorusOn) vrepl = vrepl/(1.0_rk+ exp(-4.5*(phy%relQ%P-1.0_rk)))
  if (self%GrazTurbOn .eq. 8) vrepl = vrepl*exp(-2*fa)  ! further reduce viral growth in coastal waters

 endif   !phy%C* (1.0_rk+phy%reg%C/self%vir_phyC)
  vrepl = vrepl * sens%f_T *phy%C**2/poc   ! cross section interception
 ! vrepl = vrepl * sens%f_T *phy%C
 ! vrepl = vrepl * phy%C
  vrepl = vrepl/(1.0_rk+ exp(-self%vir_infect*(vir_max-vird)))   !capacity reached

! viral removal by preferential decline of more infected hosts
!  vadap = 0.0_rk
  vadap = self%vir_loss * virf**2 * self%vir_infect *vire  ! marginal host loss due to infection
!  vadap = vadap * self%vir_phyC/(phy%reg%C+self%vir_phyC) *smooth_small(vir_max-vird,1.0_rk) !self%small
  vadap = vadap * exp(-phy%C/self%vir_phyC)*smooth_small(vir_max-vird,1.0_rk)
  vadap = vadap * vird/(vird+self%vir_spor_C)
 ! pathogenic diversity
! death and spore formation of viral cells
  vmort = self%vir_spor_r * sens%f_T**2 * vird/(vird+self%vir_spor_C)

  dvir_dt =  (vrepl - vadap - vmort) *phy%vir
  rhsv%vir = rhsv%phyC * phy%vir/phy%reg%C + dvir_dt

!  acclim%fac1= infect
!  acclim%fac2= -vdilg
  acclim%fac3= vrepl

 endif

 if (self%PhosphorusOn) then
  ! ---  PHYTOPLANKTON P
   rhsv%phyP = uptake%P              * phy%C    &
              - exud%P               * phy%C    &
              - self%dil             * phy%P    &
              - aggreg_rate          * phy%P    &
              - graz_rate            * phy%Q%P
   dQP_dt    = (rhsv%phyP * phy%C - rhsv%phyC * phy%P) / (phy%reg%C*phy%reg%C)
 endif

!if (abs(phy%C) .gt. 1d-4) then
 if (self%PhotoacclimOn ) then ! check for too small biomasses %chl

! PHYTOPLANKTON CHLa
     ! note that theta*rel_chloropl in units [mg Chla (mmol C)^{-1}]
!   dQN_dt        = (rhsv%phyN * phy%reg%C - rhsv%phyC * phy%reg%N) / (phy%reg%C*phy%reg%C)
   dQN_dt        = (rhsv%phyN * phy%C - rhsv%phyC * phy%N) / (phy%reg%C*phy%reg%C)
! TODO: dangerous to work with RHS instead of net uptake rates (mortality has no physiological effect)
#ifdef debugMK
   if ( dQN_dt/= dQN_dt .or. abs(dQN_dt)>huge(0.0_rk)) then
     write(0,*) 'ERROR: Term1 = ',(rhsv%phyN * phy%C - rhsv%phyC * phy%N)
     write(0,*) 'ERROR: rhsv%phyN = ', rhsv%phyN
     write(0,*) 'ERROR: phy%C     = ', phy%C
     write(0,*) 'ERROR: rhsv%phyC = ', rhsv%phyC
     write(0,*) 'ERROR: phy%N     = ', phy%N
     write(0,*) 'ERROR: Term2 = ',(phy%reg%C*phy%reg%C)
     !-----------------------------------------------------------
     write(0,*) 'ERROR: maecs_do#635 INF/NaN detected dQN_dt = ',dQN_dt; stop;
   endif
#endif
   dRchl_phyC_dt =  acclim%dRchl_dtheta * acclim%dtheta_dt   &
                  + acclim%dRchl_dfracR * acclim%dfracR_dt   &
                  + acclim%dRchl_dQN    * dQN_dt
#ifdef debugMK
  if ( dRchl_phyC_dt/= dRchl_phyC_dt .or. abs(dRchl_phyC_dt)>huge(0.0_rk)) then
    write(0,*) 'ERROR: Term1 = ',acclim%dRchl_dtheta * acclim%dtheta_dt
    write(0,*) 'ERROR: Term2 = ',acclim%dRchl_dfracR * acclim%dfracR_dt
    write(0,*) 'ERROR: Term3 = ',acclim%dRchl_dQN    * dQN_dt
    !-----------------------------------------------------------------------------
    write(0,*) 'ERROR: maecs_do#647 INF/NaN detected dRchl_phyC_dt = ',dRchl_phyC_dt; stop;
  endif
#endif

! pigment decay to relieve from artificially high pigm:C ratios at very low phyC
   decay = self%decay_pigm * (phy%theta*phy%rel_chloropl*self%itheta_max-0.5d0)**5
!   decay = self%decay_pigm * (phy%frac%theta-0.5d0)**5

  ! surge release at unrealistic partitioning (at sigma=1 and Q->Q0)
   if(phy%frac%theta+phy%frac%Rub .gt. 0.98d0) decay = decay + 1.0d0/(1.0d0+exp(-20*(phy%frac%theta+phy%frac%Rub-1.0d0)))

   rhsv%chl = phy%theta * phy%frac%Rub * phy%relQ%N**self%sigma * rhsv%phyC  &
               + dRchl_phyC_dt * phy%C  - decay * phy%chl

#ifdef debugMK
if ( rhsv%chl /= rhsv%chl ) then;
!   if ( dQN_dt/= dQN_dt .or. abs(dQN_dt)>huge(0.0_rk)) then
!      write(0,*) 'ERROR: Term1 = ',(rhsv%phyN * phy%C - rhsv%phyC * phy%N)
!      write(0,*) 'ERROR: Term2 = ',(phy%reg%C*phy%reg%C)
!      !-----------------------------------------------------------------------------
!      write(0,*) 'ERROR: INF/NaN detected dQN_dt = ',dQN_dt
!    endif
!    if ( dRchl_phyC_dt/= dRchl_phyC_dt .or. abs(dRchl_phyC_dt)>huge(0.0_rk)) then
!     write(0,*) 'ERROR: Term1 = ',acclim%dRchl_dtheta * acclim%dtheta_dt
!     write(0,*) 'ERROR: Term2 = ',acclim%dRchl_dfracR * acclim%dfracR_dt
!     write(0,*) 'ERROR: Term3 = ',acclim%dRchl_dQN    * dQN_dt
!     !-----------------------------------------------------------------------------
!     write(0,*) 'ERROR: INF/NaN detected dRchl_phyC_dt = ',dRchl_phyC_dt
!   endif

  if ( phy%theta    /= phy%theta     ) write(0,*) 'ERROR: NaN detected phy%theta     = ',phy%theta
  if ( phy%frac%Rub /= phy%frac%Rub  ) write(0,*) 'ERROR: NaN detected phy%frac%Rub  = ',phy%frac%Rub
  if ( phy%relQ%N   /= phy%relQ%N    ) write(0,*) 'ERROR: NaN detected phy%relQ%N    = ',phy%relQ%N
  if ( self%sigma   /= self%sigma    ) write(0,*) 'ERROR: NaN detected self%sigma    = ',self%sigma
  if ( dRchl_phyC_dt/= dRchl_phyC_dt ) write(0,*) 'ERROR: NaN detected dRchl_phyC_dt = ',dRchl_phyC_dt
  if ( phy%C        /= phy%C         ) write(0,*) 'ERROR: NaN detected phy%C         = ',phy%C
  if ( decay        /= decay         ) write(0,*) 'ERROR: NaN detected decay         = ',decay
  if ( phy%chl      /= phy%chl       ) write(0,*) 'ERROR: NaN detected phy%chl       = ',phy%chl

  if ( abs(phy%theta    ) > huge(0.0_rk) ) write(0,*) 'ERROR: INF detected phy%theta     = ',phy%theta
  if ( abs(phy%frac%Rub ) > huge(0.0_rk) ) write(0,*) 'ERROR: INF detected phy%frac%Rub  = ',phy%frac%Rub
  if ( abs(phy%relQ%N   ) > huge(0.0_rk) ) write(0,*) 'ERROR: INF detected phy%relQ%N    = ',phy%relQ%N
  if ( abs(self%sigma   ) > huge(0.0_rk) ) write(0,*) 'ERROR: INF detected self%sigma    = ',self%sigma
  if ( abs(dRchl_phyC_dt) > huge(0.0_rk) ) write(0,*) 'ERROR: INF detected dRchl_phyC_dt = ',dRchl_phyC_dt
  if ( abs(phy%C        ) > huge(0.0_rk) ) write(0,*) 'ERROR: INF detected phy%C         = ',phy%C
  if ( abs(decay        ) > huge(0.0_rk) ) write(0,*) 'ERROR: INF detected decay         = ',decay
  if ( abs(phy%chl      ) > huge(0.0_rk) ) write(0,*) 'ERROR: INF detected phy%chl       = ',phy%chl

  write(0,*) 'ERROR: Term1                  = ',phy%theta * phy%frac%Rub * phy%relQ%N**self%sigma * rhsv%phyC
  write(0,*) 'ERROR: phy%relQ%N**self%sigma = ',phy%relQ%N**self%sigma
  write(0,*) 'ERROR: dRchl_phyC_dt * phy%C  = ',dRchl_phyC_dt * phy%C
  write(0,*) 'ERROR: decay * phy%chl        = ',decay * phy%chl

  !-----------------------------------------------------------
  write(0,*) 'ERROR: maecs_do#701 rhsv%chl = ',rhsv%chl; stop;
endif
#endif

!if (rhsv%chl .lt. -200.d0) then
!  phy%theta     = phy%chl / (phy%rel_chloropl * phy%reg%C)   ! trait variable
!  phy%frac%theta= phy%theta * phy%rel_chloropl * self%itheta_max ! []     no
!  write (*,'(A,4(F14.3))') 'dChl=',rhsv%chl,phy%chl,dRchl_phyC_dt,rhsv%phyC,dQN_dt
!  write (*,'(A,6(F14.3))') 'aa=',acclim%dRchl_dtheta,acclim%dtheta_dt,acclim%dRchl_dfracR,acclim%dfracR_dt,acclim%dRchl_dQN, dQN_dt
!  write (*,'(A,2(F14.3))') 'dmdx=',acclim%fac1,acclim%fac2
!end if

!write (*,'(A,4(F10.3))') 'rhs chl=', phy%theta * phy%frac%Rub * phy%relQ%N**self%sigma * rhsv%phyC,dRchl_phyC_dt * phy%reg%C*1E1,phy%relQ%N**self%sigma,phy%theta

!_____________________________________________ _________________________________
 end if ! PhotoacclimOn

 if (self%RubiscoOn) then
!   decay = self%decay_pigm * (exp(phy%frac%Rub)-1.1d0)
   decay = self%decay_pigm * (phy%frac%Rub-0.5d0)**5
   rhsv%Rub  = phy%Rub/phy%reg%C * rhsv%phyC + acclim%dfracR_dt * phy%C - decay*phy%Rub
 end if
!else  rhsv%Rub  = 0.0d0  rhsv%chl  = 0.0d0
!endif !if (abs(phy%C) .gt. 1d-4)

!________________________________________________________________________________
!
!  additional aggregation mediated by TEP production following C-overconsumption and release
  dq_dt = 0.0_rk
! dq_dt=min(dQN_dt/self%QN_phy_max, dQP_dt/self%QP_phy_max)*(1.0_rk-exp(-0.1*env%par)
! dq_dt=dQP_dt/(1E-4+phy%Q%P) + dQN_dt/(1E-4+phy%Q%N)
!  fag = 1.0_rk/(1+exp(10*(dq_dt)))
!  add_aggreg_rate = fag * min(self%phi_agg * (phy%N + det%N),0.25_rk)
  add_aggreg_rate = 0.0_rk
  aggreg_rate = aggreg_rate + add_aggreg_rate
  rhsv%phyN =  rhsv%phyN  - add_aggreg_rate * phy%N
  rhsv%phyC =  rhsv%phyC  - add_aggreg_rate * phy%C
  _SET_DIAGNOSTIC_(self%id_pPads,dq_dt )       !average Temporary_diagnostic_
! _SET_DIAGNOSTIC_(self%id_datt, fag)

!________________________________________________________________________________
!
! ZOOPLANKTON zoo%feeding
if (self%GrazingOn) then
   rhsv%zooC   =  zoo%yield * graz_rate       &
                - zoo_mort          * zoo%C   &
                - self%dil          * zoo%C   &
                - lossZ%C           * zoo%C
else
   rhsv%zooC      = 0.0_rk
end if

! ------------------------------------------------------------------
!  ---  POM&DOM quality, relative to max N-quota of phytoplankton
if (self%remNP .gt. -1E-5) then ! base variant: quality correlates with detritus N:C
   qualPOM     = (1.0d0-self%Nqual) + self%Nqual * det%N /(det%C + self%small_finite)  * self%iK_QN
   qualDOM     = (1.0d0-self%Nqual) + self%Nqual * dom%N /(dom%C + self%small_finite)  * self%iK_QN
else                      ! alternative: quality correlates with detritus P:C
   qualPOM     = (1.0d0-self%Nqual) + self%Nqual * det%P /(det%C + self%small_finite)  * self%iK_QP
   qualDOM     = (1.0d0-self%Nqual) + self%Nqual * dom%P /(dom%C + self%small_finite)  * self%iK_QP
end if

!  ---  hydrolysis & remineralisation rate (temp dependent)
degradT     = self%hydrol * sens%f_T * qualPOM
reminT      = self%remin  * sens%f_T * qualDOM

!  ---  hydrolysis & remineralisation depend on quality, here propto N/C quota of OM
!  acceleration: rate difference for N-pool
if (self%remNP .gt. -1E-5) then ! base variant: quality correlates with detritus N:C
   ddegN       = self%hydrol * sens%f_T * max(1.0d0 - qualPOM, 0.0d0)
   ddegP       = self%remNP * ddegN
   dremN       = self%remin * sens%f_T * max(1.0d0 - qualDOM, 0.0d0)
   dremP       = self%remNP * dremN
else                      ! alternative: quality correlates with detritus P:C
   ddegP       = self%hydrol * sens%f_T * max(1.0d0 - qualPOM, 0.0d0)
   ddegN       = -self%remNP * ddegP
   dremP       = self%remin * sens%f_T * max(1.0d0 - qualDOM, 0.0d0)
   dremN       = -self%remNP * dremP
end if
!________________________________________________________________________________
!
!  --- DETRITUS C
det_prod    = floppZ%C              * zoo%C   &
             + aggreg_rate          * phy%C   &
             + zoo_mort             * zoo%C

rhsv%detC   = det_prod                        &
             - self%dil             * det%C   &
             - degradT              * det%C

!________________________________________________________________________________
!
!  --- DETRITUS N
rhsv%detN   = floppZ%N              * zoo%C   &
             + aggreg_rate          * phy%N   &
             - self%dil             * det%N   &
             + zoo_mort             * zoo%N   &
             - (degradT + ddegN)    * det%N   &
             - denitrate

!________________________________________________________________________________
!
!  --- DOC
 Cprod      = reminT                * dom%C
rhsv%domC   = exud%C                * phy%C   &
             + degradT              * det%C   &
             - self%dil             * dom%C   &
             - Cprod
!________________________________________________________________________________
!
!  --- DON
Nprod       = (reminT + dremN)    * dom%N
rhsv%domN   = exud%N                * phy%C   &
             + (degradT + ddegN)    * det%N   &
             - self%dil             * dom%N   &
             - Nprod
!________________________________________________________________________________
!
! DIC
!if (self%BioCarbochemOn) then
!  rhsv%dic     = -uptake%grossC      * phy%C   &
!                + uptake%lossC       * phy%C   &
!                + reminT             * dom%C   &
!                + lossZ%C            * zoo%C
!
!_SET_ODE_(self%id_dic,rhsv%dic UNIT)
!end if
!________________________________________________________________________________
!
!  --- DIN
rhsv%nutN   = -uptake%N            * phy%C    &
             + Nprod                          &
             + lossZ%N             * zoo%C    &
             + self%dil * (self%nutN_initial - nut%N) &
             + deporate

if(self%ChemostatOn .and. self%remin .lt. 0.0001d0) rhsv%nutN = 0.0d0

!________________________________________________________________________________
!
if (self%PhosphorusOn) then
  ! ---  PHYTOPLANKTON P
   rhsv%phyP = rhsv%phyP- add_aggreg_rate          * phy%P

  !  --- DETRITUS P
   rhsv%detP = floppZ%P              * zoo%C    &
              + aggreg_rate          * phy%P    &
              - self%dil             * det%P    &
              + zoo_mort             * zoo%P    &
              - (degradT + ddegP)    * det%P    ! quality enhances P remin

  !  --- DOP
   Pprod     = (reminT + dremP)      * dom%P
   rhsv%domP = exud%P                * phy%C    &
              + (degradT + ddegP)    * det%P    &
              - self%dil             * dom%P    &
              - Pprod

  !  --- DIP
   rhsv%nutP = - uptake%P            * phy%C    &
              + Pprod                           &
              + lossZ%P              * zoo%C    &
              + self%dil * (self%nutP_initial - nut%P)

   if(self%ChemostatOn .and. self%remin .lt. 0.0001d0) rhsv%nutP = 0.0d0

end if

!________________________________________________________________________________
!
if (self%SiliconOn) then
  ! ---  PHYTOPLANKTON Si
   rhsv%phyS = uptake%Si              * phy%C    &
              - self%dil             * phy%Si    &
              - aggreg_rate          * phy%Si    &
              - graz_rate            * phy%Q%Si
  !  --- DETRITUS Si
   rhsv%detS = floppZ%Si              * zoo%C    &
              + aggreg_rate          * phy%Si    &
              - self%dil             * det%Si    &
              - degradT              * det%Si
  !  --- Dissolved Si
   rhsv%nutS = - uptake%Si            * phy%C    &
              + degradT              * det%Si    &
              + lossZ%Si              * zoo%C    &
              + self%dil * (self%nutS_initial - nut%Si)

end if

!---------- RHS for BGC/diagensis model part ----------
! Fortran 2003 version of OMEXDIA+P biogeochemical model
! The OMEXDIA+P+MPB model is based on the OMEXDIA model (see Soetard et al. 1996a)
! P-cycle is added by kai wirtz

if (self%BioOxyOn) then

! ---------- temperature    TODO: retrieve from existing temp variables
   f_T    = sens%f_T

! ---------- manages overlapping state variables
   no3    = max(nut%N - env%nh3,  0.0d0)
! ---------- remineralisation limitations
   Oxicminlim = env%oxy/(env%oxy+self%ksO2oxic+relaxO2*(env%nh3+env%odu))
   Denitrilim = (1.0_rk-env%oxy/(env%oxy+self%kinO2denit)) * no3/(no3+self%ksNO3denit)
   Anoxiclim  = (1.0_rk-env%oxy/(env%oxy+self%kinO2anox)) * (1.0_rk-no3/(no3+self%kinNO3anox))
   Rescale    = 1.0_rk/(Oxicminlim+Denitrilim+Anoxiclim)

! extra-omexdia P -dynamics
  if (self%PhosphorusOn) then
!   PO4-adsorption ceases when critical capacity is reached
!   [FeS] approximated by ODU
!   po4    = nut%P
    radsP      = self%rPAds * degradT * nut%P * max(env%odu,self%PAdsODU)
    rhsv%nutP  = rhsv%nutP - radsP
    rhsv%detP  = rhsv%detP + radsP
!   rP     = self%rFast * (1.0_rk - Oxicminlim)
!   Pprod  = rP * pdet
  endif

! Oxic mineralisation, denitrification, anoxic mineralisation
! then the mineralisation rates
   OxicMin    = Cprod*Oxicminlim*Rescale        ! oxic mineralisation
   Denitrific = Cprod*Denitrilim*Rescale        ! Denitrification
   AnoxicMin  = Cprod*Anoxiclim *Rescale        ! anoxic mineralisation

! reoxidation and ODU deposition
   Nitri      = f_T * self%rnit   * env%nh3 * env%oxy/(env%oxy + self%ksO2nitri &
                  + relaxO2*(dom%C + env%odu))
   OduOx      = f_T * self%rODUox * env%odu * env%oxy/(env%oxy + self%ksO2oduox &
                  + relaxO2*(env%nh3 + dom%C))

!  pDepo      = min(1.0_rk,0.233_rk*(wDepo)**0.336_rk )
   pDepo      = 0.0_rk
   OduDepo    = AnoxicMin * pDepo

!  dynamics of env%oxy ~ dissolved oxygen
  rhsv%oxy    = - OxicMin - 2.0_rk* Nitri - OduOx &
                - lossZ%C * zoo%C + uptake%C * phy%C  &
                + self%dil * (self%oxy_initial - env%oxy)

!  dynamics of odu ~ dissolved reduced substances
  rhsv%odu    = (AnoxicMin - OduOx - OduDepo)
!  dynamics of no3 ~ dissolved nitrate
!!  rhsv%no3 = (-0.8_rk*Denitrific + Nitri - uptNO3)

! Anammox: NH3 oxidation by nitrite, here related to NO3
  Anammox    = self%rAnammox * Anoxiclim *Rescale * env%nh3 * no3/(no3+self%ksNO3denit)

! preference for NH3 in DIN-uptake of autotrophs
!  nh3f        = 1.0d0 -exp(-5*env%nh3/smooth_small(nut%N,self%small))
!  nh3f        = 1.0d0/(1.0d0+exp(-30*(env%nh3/(nut%N + env%nh3+0.01d0 )-0.5d0)))
  nh3f        = env%nh3/(nut%N+ env%nh3+0.1d0)
!  dynamics of nh3 ~ dissolved ammonium
  rhsv%nh3    = Nprod - Nitri + lossZ%N * zoo%C  & !/ (1.0_rk + self%NH3Ads)
               - nh3f*uptake%N * phy%C &!env%nh3/(nut%N+self%small) *
               + self%dil * (self%nh3_initial - env%nh3) &
               - Anammox - max(env%nh3 - nut%N,  0.0d0)

!  if(env%nh3 .lt. 0.1d0) write (*,'(A,11(F10.3))') 'nh3=', env%nh3,nut%N,rhsv%nh3,nh3f,- Nitri,uptake%N,- nh3f*uptake%N*phy%C,self%dil*(self%nh3_initial-env%nh3),(0.8d0 * Denitrific+Anammox)*1E5

  rhsv%nutN   = rhsv%nutN - 0.8d0 * Denitrific - Anammox

!  dynamics of pdet ~ detritus-P
!    rhsv%pdet = (radsP - f_T * Pprod)
!  dynamics of po4 ~ dissolved phosphate
!  rhsv%po4 = (f_T * Pprod - radsP)
end if !BioOxyOn

end if !IsCritical

!#S_ODE
!---------- ODE for each state variable ----------
  _SET_ODE_(self%id_nutN, rhsv%nutN UNIT)
  _SET_ODE_(self%id_phyC, rhsv%phyC UNIT)
  _SET_ODE_(self%id_phyN, rhsv%phyN UNIT)
  _SET_ODE_(self%id_detC, rhsv%detC UNIT)
  _SET_ODE_(self%id_detN, rhsv%detN UNIT)
  _SET_ODE_(self%id_domC, rhsv%domC UNIT)
  _SET_ODE_(self%id_domN, rhsv%domN UNIT)
if (self%RubiscoOn) then
      _SET_ODE_(self%id_Rub, rhsv%Rub UNIT)
end if
if (self%PhotoacclimOn) then
      _SET_ODE_(self%id_chl, rhsv%chl UNIT)
end if
if (self%PhosphorusOn) then
      _SET_ODE_(self%id_nutP, rhsv%nutP UNIT)
      _SET_ODE_(self%id_phyP, rhsv%phyP UNIT)
      _SET_ODE_(self%id_detP, rhsv%detP UNIT)
      _SET_ODE_(self%id_domP, rhsv%domP UNIT)
end if
if (self%SiliconOn) then
      _SET_ODE_(self%id_nutS, rhsv%nutS UNIT)
      _SET_ODE_(self%id_phyS, rhsv%phyS UNIT)
      _SET_ODE_(self%id_detS, rhsv%detS UNIT)
end if
if (self%GrazingOn) then
      _SET_ODE_(self%id_zooC, rhsv%zooC UNIT)
end if
if (self%BioOxyOn) then
      _SET_ODE_(self%id_nh3, rhsv%nh3 UNIT)
      _SET_ODE_(self%id_oxy, rhsv%oxy UNIT)
      _SET_ODE_(self%id_odu, rhsv%odu UNIT)
end if
if (self%VirusOn) then
      _SET_ODE_(self%id_vir, rhsv%vir UNIT)
end if
if (self%NResOn) then
      _SET_ODE_(self%id_RNit, rhsv%RNit UNIT)
end if
!#E_ODE

! artifical, serial nutrient input to illustrate co-limitation dynamics in 0D
!if (self%ChemostatOn) then
!  _GET_GLOBAL_ (self%id_doy,doy) !day of year
!  select case (doy)
!           case (89:92)
!            _SET_ODE_(self%id_nutP, 5*(1.0-nut%P) UNIT)
!           case (29:32)
!            _SET_ODE_(self%id_nutN, 5*(16.0-nut%N) UNIT)
!           case (59:62)
!            _SET_ODE_(self%id_nutS, 5*(16.0-nut%Si) UNIT)
!  end select
!endif

!_SET_DIAGNOSTIC_(self%id_vphys, exp(-self%sink_phys*phy%relQ%N * phy%relQ%P))       !average

! experimental formulation for emulating P-adsorption at particles in the water column and at the bottom interface

!_GET_HORIZONTAL_(self%id_o2flux, flO2)
!_GET_HORIZONTAL_(self%id_oduflux, flODU)
!_GET_HORIZONTAL_(self%id_zmax, zmax)  ! max depth

!aPO4 = (flODU-flO2)/(zmax+self%small)

!_SET_DIAGNOSTIC_(self%id_vphys, aPO4)       !average Temporary_diagnostic_

!________________________________________________________________________________
! set diag variables, mostly from PrimProd module

!define _REPLNAN_(X) X !changes back to original code
#define _REPLNAN_(X) nan_num(X)

if (self%kwFzmaxMeth > 0) then
  !> Back-calculate background attenuation (kw) consistent with get_light subroutine
  a_lit = att - (self%a_spm*(det%C+zoo%C)+ self%a_doc*dom%C  &
    + self%a_phyc*phy%C + self%a_chl*phy%chl)
  _SET_DIAGNOSTIC_(self%id_bgatt, _REPLNAN_(a_lit))
endif

!#S_DIA
if (self%DebugDiagOn) then
  _SET_DIAGNOSTIC_(self%id_tmp, _REPLNAN_(acclim%tmp))       !average Temporary_diagnostic_
!  _SET_DIAGNOSTIC_(self%id_fac4, _REPLNAN_(dRchl_phyC_dt))   !average Auxiliary_diagnostic_
!  _SET_DIAGNOSTIC_(self%id_fac4, _REPLNAN_(rhsv%phyC * phy%vir/phy%reg%C))   !average Auxiliary_diagnostic_
  _SET_DIAGNOSTIC_(self%id_fac4, _REPLNAN_(- vadap ))  !average Auxiliary_diagnostic
  _SET_DIAGNOSTIC_(self%id_fac5, _REPLNAN_( - vmort)) !average Auxiliary_diagnostic_
!  _SET_DIAGNOSTIC_(self%id_fac5, _REPLNAN_(sens%f_T2)) !average Auxiliary_diagnostic_
!  _SET_DIAGNOSTIC_(self%id_fac5, _REPLNAN_(acclim%dRchl_dfracR*acclim%dfracR_dt)) !average Auxiliary_diagnostic_
  _SET_DIAGNOSTIC_(self%id_fac3, _REPLNAN_(acclim%fac3)) !average Auxiliary_diagnostic_dRchl_dtheta*acclim%dtheta_dt
  _SET_DIAGNOSTIC_(self%id_fac1, _REPLNAN_(acclim%fac1))     !average dtheta_dt_due_to_flex_theta_
  _SET_DIAGNOSTIC_(self%id_fac2, _REPLNAN_(acclim%fac2))     !average dtheta_dt_due_to_grad_theta_
end if
if (self%BGC0DDiagOn) then
  _SET_DIAGNOSTIC_(self%id_GPPR, _REPLNAN_(phy%gpp*phy%C))   !average gross_primary_production_
!  _SET_DIAGNOSTIC_(self%id_Denitr, _REPLNAN_(0.8*Denitrific)) !average denitrification_rate_
  _SET_DIAGNOSTIC_(self%id_DNP, _REPLNAN_(nut%N/(nut%P+self%small))) !average DIN:DIP_ratio_
  _SET_DIAGNOSTIC_(self%id_QNP, _REPLNAN_(phy%Q%N/phy%Q%P))  !average N:P_ratio_
  _SET_DIAGNOSTIC_(self%id_qualPOM, _REPLNAN_(qualPOM))      !average Quality_of_POM_
  _SET_DIAGNOSTIC_(self%id_qualDOM, _REPLNAN_(qualDOM))      !average Quality_of_DOM_
!  _SET_DIAGNOSTIC_(self%id_no3, _REPLNAN_(no3))              !average Nitrate_
end if
if (self%PhysiolDiagOn) then
  _SET_DIAGNOSTIC_(self%id_chl2C, _REPLNAN_(phy%theta*phy%rel_chloropl/12)) !average chlorophyll:carbon_ratio_=_chl-a/chloroplast-C_*_chloroplast-C/phy-molC_*_1molC/12gC_
  _SET_DIAGNOSTIC_(self%id_Theta, _REPLNAN_(phy%theta))      !average Theta_
  _SET_DIAGNOSTIC_(self%id_fracR, _REPLNAN_(phy%frac%Rub))   !average Rubisco_fract._allocation_
  _SET_DIAGNOSTIC_(self%id_fracT, _REPLNAN_(phy%frac%theta)) !average LHC_fract._allocation_
  _SET_DIAGNOSTIC_(self%id_fracNU, _REPLNAN_(phy%frac%NutUpt)) !average Nut._Uptake_fract._allocation_
  _SET_DIAGNOSTIC_(self%id_QN, _REPLNAN_(phy%Q%N))           !average N:C_ratio_
  _SET_DIAGNOSTIC_(self%id_QP, _REPLNAN_(phy%Q%P))           !average P:C_ratio_
  _SET_DIAGNOSTIC_(self%id_QSi, _REPLNAN_(phy%Q%Si))         !average Si:C_ratio_
  _SET_DIAGNOSTIC_(self%id_aVN, _REPLNAN_(acclim%aV%N))      !average N-uptake_activity_
  _SET_DIAGNOSTIC_(self%id_aVP, _REPLNAN_(acclim%aV%P))      !average P-uptake_activity_
  _SET_DIAGNOSTIC_(self%id_aVSi, _REPLNAN_(acclim%aV%Si))    !average Si-uptake_activity_
  _SET_DIAGNOSTIC_(self%id_faN, _REPLNAN_(acclim%fA%N))      !average N-uptake_affinity_allocation_
  _SET_DIAGNOSTIC_(self%id_faP, _REPLNAN_(acclim%fA%P))      !average P-uptake_affinity_allocation_
  _SET_DIAGNOSTIC_(self%id_faSi, _REPLNAN_(acclim%fA%Si))    !average Si-uptake_affinity_allocation_
  _SET_DIAGNOSTIC_(self%id_rQN, _REPLNAN_(phy%relQ%N))       !average Relative_N-Quota_
  _SET_DIAGNOSTIC_(self%id_rQP, _REPLNAN_(phy%relQ%P))       !average Relative_P-Quota_
  _SET_DIAGNOSTIC_(self%id_rQSi, _REPLNAN_(phy%relQ%Si))     !average Relative_Si-Quota_
end if
if (self%RateDiagOn) then
  _SET_DIAGNOSTIC_(self%id_phyUR, _REPLNAN_(uptake%C))       !average Phytoplankton_C_Uptake_Rate_
  _SET_DIAGNOSTIC_(self%id_phyRER, _REPLNAN_(-phy%resp))     !average Phytoplankton_Respiration_Rate_
  _SET_DIAGNOSTIC_(self%id_phyELR, _REPLNAN_(-exud%C))       !average Phytoplankton_Exudation_Loss_Rate_
  _SET_DIAGNOSTIC_(self%id_phyALR, _REPLNAN_(-aggreg_rate))  !average Phytoplankton_Aggregation_Loss_Rate_
  _SET_DIAGNOSTIC_(self%id_phyVLR, _REPLNAN_(-viral_rate))   !average Phytoplankton_Viral_Loss_Rate_
  _SET_DIAGNOSTIC_(self%id_phyGLR, _REPLNAN_(-graz_rate/phy%reg%C)) !average Phytoplankton_Grazing_Loss_Rate_
!  _SET_DIAGNOSTIC_(self%id_vsinkr, _REPLNAN_(exp(-self%sink_phys*phy%relQ%N*phy%relQ%P))) !average Relative_Sinking_Rate_
  _SET_DIAGNOSTIC_(self%id_zoomort, _REPLNAN_(zoo_mort))     !average Zooplankton_Mortality_Rate_
end if
!#E_DIA


#if _DEBUG_
write(*,'(A)') 'end DO'
#endif

  _LOOP_END_

end subroutine maecs_do

!> @brief handles vertical movement for depth-varying movement rates
!> @details phyto sinking rate depends on the nutritional state, so for each node:
!! \n \f$ phy\%relQ \f$ obtained by calling calc_internal_states(self,phy,det,dom,zoo)
!! \n then \f$ phyQstat=phy\%relQ\%N * phy\%relQ\%P \f$
!! \n finally, vs_phy = maecs_functions::sinking(self\%vS_phy, phyQstat, vs_phy)
subroutine maecs_get_vertical_movement(self,_ARGUMENTS_GET_VERTICAL_MOVEMENT_)

use maecs_functions
use maecs_types

implicit none
!
! !INPUT PARAMETERS:
 class (type_hzg_maecs),intent(in) :: self
_DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_
 !   REALTYPE, intent(in)              ::vstokes
type (type_maecs_phy):: phy !< maecs phytoplankton type
type (type_maecs_zoo) :: zoo
type (type_maecs_om):: det
type (type_maecs_om):: dom
logical  :: IsCritical = .false. ! phyC and phyN below reasonable range ?

!
! !LOCAL VARIABLES:
REALTYPE    :: phyQstat, ef, vs_phy,vs_det, phyEner, minPigm,min_Cmass, minc, zmax, par, vs_zoo
!REALTYPE    :: aggf, agge=16.d0
REALTYPE, parameter :: secs_pr_day = 86400.d0
!EOP
!-----------------------------------------------------------------------
!BOC

!define _REPLNAN_(X) X !changes back to original code
#define _REPLNAN_(X) nan_num(X)

_FABM_LOOP_BEGIN_

#if _DEBUG_
write(*,'(A)') 'begin vert_move'
#endif

   ! Retrieve phtoplankton state

   !Retrieve the 'phyQstat' directly as a diagnostic variable: does not work yet.
   !fabm_get_bulk_diagnostic_data(self%id_phyqstat,phyQstatD) !where, phyQstat=relQ%N*relQ%P
   !_GET_(self%id_phyqstat,phyQstatD)

   !Calculate manually
   _GET_(self%id_phyC, phy%C)  ! Phytplankton Carbon in mmol-C/m**3
   _GET_(self%id_phyN, phy%N)  ! Phytplankton Nitrogen in mmol-N/m**3
   _GET_(self%id_detC, det%C)  ! Detritus Nitrogen in mmol-N/m**3
!   _GET_(self%id_detN, det%N)  ! Detritus Nitrogen in mmol-N/m**3
   _GET_(self%id_domC, dom%C)  ! DONitrogen in mmol-N/m**3
    if (self%PhosphorusOn) then
!      _GET_(self%id_detP, det%P)  ! Detritus Phosphorus in mmol-P/m**3
      _GET_(self%id_phyP, phy%P)  ! Phytplankton Phosphorus in mmol-P/m**3
    endif
 !
 ! ****** SINKING AS A FUNCTION OF INTERNAL STATES *************
 !
    call min_mass(self,phy, min_Cmass, IsCritical, method=_KAI_)
    call calc_internal_states(self,phy,det,dom,zoo)
    !write (*,'(A,2(F10.3))') 'phy%relQ%N, phy%relQ%P=', phy%relQ%N, phy%relQ%P

    phyQstat = 1.0_rk
! infected cells sink faster to save population
    if (self%vir_loss .gt. self%small_finite .or. self%VirusOn ) then
      _GET_(self%id_vir, phy%vir)  ! Virus C density in cells in -
      phyQstat = phyQstat * 1.0_rk/(1.0_rk + exp(-self%vir_infect*(1.0_rk-phy%vir/(phy%C+self%small_finite))))
! threshold virus with multi-stage replication
    endif
! non-linear response of sinking speed to nutrient limitation
    if(self%sink_phys .gt. 0) then ! weak non-linearity
      phyQstat = phyQstat * phy%relQ%N**2/(HALFQ**2+phy%relQ%N**2)
      if (self%PhosphorusOn) then
        phyQstat =phyQstat * phy%relQ%P**2/(HALFQ**2+phy%relQ%P**2)
      end if
     ! Calculate sinking speed
      vs_phy = -self%vS_phy * exp( -self%sink_phys * phyQstat)
    else   ! strong non-linearity
      phyQstat = phyQstat * 1.0_rk/(1.0_rk + exp(self%sink_phys*(phy%relQ%N-HALFQ)))
      if (self%PhosphorusOn) then
        phyQstat =phyQstat * 1.0_rk/(1.0_rk + exp(self%sink_phys*(phy%relQ%P-HALFQ)))
      end if
      vs_phy = -self%vS_phy * (1.0_rk-phyQstat)
   endif
   ! nutrient limitation ; TODO check product rule and add other elements such as Si
   ! phyQstat = phy%relQ%N * phy%relQ%P

   ! energy limitation ; TODO check function and quantity
!   phyEner  = phy%gpp / (self%V_NC_max*self%zeta_CN)
! smoothed minimum of energy and nutrient limitation;
!   phyQstat = phyQstat - smooth_small(phyQstat - phyEner, self%small)
!  call sinking(self%vS_phy, phyQstat, vs_phy)

! ascending of top-conditioned cells
   if(self%genMeth .gt. 0) then
     vs_phy = vs_phy + self%vS_phy * exp(-3.0d0+self%genMeth*0.2d0)
   elseif(self%genMeth .lt. 0) then
     vs_phy = vs_phy - self%vS_phy * exp(-3.0d0-self%genMeth*0.2d0)
   endif
   !if (self%RateDiagOn) then
   ! write (*,'(A)',advance='no') '' ! Silly Fix to 'NETCDF: Numeric conversion not representable' problem ??
   !  _SET_DIAGNOSTIC_(self%id_vsinkr, _REPLNAN_(vs_phy)) !average Relative Sinking Velocity
   !end if
   !CONSTANT SINKING
   !vs_phy = self%vS_phy

   _SET_DIAGNOSTIC_(self%id_vsinkr, _REPLNAN_(vs_phy)) !average Relative_Sinking_Rate_

  _GET_HORIZONTAL_(self%id_zmax, zmax)  ! max depth
! increase vDet in shallow water:co-agulation with lithogenic particles
!   vs_det = vs_det * (1.0_rk+ 4*exp(-(zmax/20.0_rk)**2))
   vs_phy = vs_phy / secs_pr_day
   !write (*,'(A,2(F10.3))') 'phyQstat, vs_phy=', phyQstat, vs_phy
!   vs_det = -self%vS_det*aggf/secs_pr_day
   vs_det = -self%vS_det / secs_pr_day

   vs_det = vs_det * (0.01*det%C)**0.38_rk

! slowing down of vertical velocities at high and very low concentration to smooth numerical problems in shallow, pesitional boxes
   ef     = 20_rk/(1+zmax)
   vs_det = vs_det * 1.0_rk/(1.0_rk+((0.003*det%C +0.02*det%N + 100*self%small_finite/(det%C+self%small_finite))*ef )**4 )
! additional slowdown in very shallow waters
!   if (vs_det*secs_pr_day .gt. 1.0_rk .and. zmax .lt. 8.0_rk) then
!     ef = exp(2*(zmax-5.0_rk))
!     vs_det = vs_det * (1.0_rk+ef)/(3.0_rk+ef)
!   endif
   vs_phy = vs_phy * 1.0_rk/(1.0_rk+((0.002*phy%C + 100*self%small_finite/(phy%C+self%small_finite))*ef )**4 )

  !set the rates
   _SET_VERTICAL_MOVEMENT_(self%id_detC,vs_det)
   _SET_VERTICAL_MOVEMENT_(self%id_detN,vs_det)
   _SET_VERTICAL_MOVEMENT_(self%id_phyN,vs_phy)
   _SET_VERTICAL_MOVEMENT_(self%id_phyC,vs_phy)
!   if (ef .lt. 1.0_rk) then
!    _GET_(self%id_par, par)  ! light photosynthetically active radiatio
!    vs_zoo = 22.0_rk * (1.0_rk-2.0_rk/(1.0_rk + exp(1.0d0-par)))
!    _SET_VERTICAL_MOVEMENT_(self%id_zooC,vs_zoo / secs_pr_day)
!   else
!    _SET_VERTICAL_MOVEMENT_(self%id_zooC,0.0_rk)
!   endif
   if (self%PhosphorusOn) then
      _SET_VERTICAL_MOVEMENT_(self%id_phyP,vs_phy)
      _SET_VERTICAL_MOVEMENT_(self%id_detP,vs_det)
   end if
   if (self%SiliconOn) then
      _SET_VERTICAL_MOVEMENT_(self%id_phyS,vs_phy)
      _SET_VERTICAL_MOVEMENT_(self%id_detS,vs_det)
   end if
   if (self%PhotoacclimOn) then
      _SET_VERTICAL_MOVEMENT_(self%id_chl, vs_phy)
      _SET_VERTICAL_MOVEMENT_(self%id_Rub, vs_phy)
   end if
   if (self%VirusOn) then
      _SET_VERTICAL_MOVEMENT_(self%id_vir, vs_phy)
   end if

#if _DEBUG_
write(*,'(A)') 'end vert_move'
#endif

_FABM_LOOP_END_

end subroutine maecs_get_vertical_movement

subroutine maecs_do_surface(self,_ARGUMENTS_DO_SURFACE_)
   use maecs_functions

   class (type_hzg_maecs), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_

   real(rk) :: tot_vi_C,tot_vi_N,tot_vi_P,tot_vi_S, O2flux,O2airbl,oxy,tot_vi_GPPR,tot_vi_Denitr

!define _REPLNAN_(X) X !changes back to original code
#define _REPLNAN_(X) nan_num(X)

   _HORIZONTAL_LOOP_BEGIN_

#if _DEBUG_
write(*,'(A)') 'begin surface_DO'
#endif

      !if (self%BGC2DDiagOn) then
        !_GET_HORIZONTAL_(self%id_GPPR_vertint,tot_vi_GPPR)
        !_SET_HORIZONTAL_DIAGNOSTIC_(self%id_GPPR_vertint_diag,_REPLNAN_(tot_vi_GPPR))
        !if (self%BioOxyOn) then
        !  _GET_HORIZONTAL_(self%id_Denitr_vertint,tot_vi_Denitr)
        !  _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Denitr_vertint_diag, _REPLNAN_(tot_vi_Denitr))
        !end if
      !end if

      if (self%Budget2DDiagOn) then
      _GET_HORIZONTAL_(self%id_totN_vertint,tot_vi_N)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_totN_vertint_diag,_REPLNAN_(tot_vi_N))
      _GET_HORIZONTAL_(self%id_totC_vertint,tot_vi_C)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_totC_vertint_diag,_REPLNAN_(tot_vi_C))
      if (self%PhosphorusOn) then
         _GET_HORIZONTAL_(self%id_totP_vertint,tot_vi_P)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_totP_vertint_diag,_REPLNAN_(tot_vi_P))
      end if
      if (self%SiliconOn) then
         _GET_HORIZONTAL_(self%id_totS_vertint,tot_vi_S)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_totS_vertint_diag,_REPLNAN_(tot_vi_S))
      end if
      end if

! --- wet and dry deposition of NO3
      _SET_SURFACE_EXCHANGE_(self%id_nutN, self%N_depo UNIT)
! --- atmospheric deposition of PO4
      if (self%PhosphorusOn) then
         _SET_SURFACE_EXCHANGE_(self%id_nutP, self%P_depo UNIT)
      end if

! --- oxygen flux between sea water and air -----
      if (self%BioOxyOn) then
! O2 flux across the boundary layer
! O2airbl is the saturation concentration of O2
! airsea_ex is the average diffusivity coefficient (m2/sec) divided by the thickness of the boundary layer.
! for O2 in mmol m-3, the rate of exchange in mmol m-2 s-1).
! Positive values imply a flux into the water, negative: out of the water.
         O2airbl = self%O2_sat
!        _GET_HORIZONTAL_(self%id_O2airbl, O2airbl)! boundary layer dissolved oxygen in mmolO2/m**3
        _GET_(self%id_oxy, oxy)   ! sea water dissolved oxygen in mmolO2/m**3

        O2flux  = self%ex_airsea * (O2airbl - oxy)!
        _SET_SURFACE_EXCHANGE_(self%id_oxy, O2flux )
        if (self%BGC2DDiagOn) then
          _SET_HORIZONTAL_DIAGNOSTIC_(self%id_O2flux_diag, _REPLNAN_(O2flux)) ! converts mmol/m2.s to mmol/m2.d
        end if
      endif

#if _DEBUG_
write(*,'(A)') 'end surface_DO'
#endif

   _HORIZONTAL_LOOP_END_

end subroutine maecs_do_surface
! potential entries in maecs_deps.lst;but might work only using GOTM-input scheme
! O2airbl	mmol-C/m**3 horizontal_dependency O2airbl surface_molecular_oxygen_partial_pressure_difference_between_sea_water_and_air #BioOxyOn
!N2air	mmol-C/m**3 horizontal_dependency N2air mole_concentration_of_atomic_nitrogen_in_air
!N2flux	mmol-N/m**2/d horizontal_diagnostic_variable  'N2flux','mmol-N/m**2/d','nitrogen_flux_between_sea_water_and_air'
