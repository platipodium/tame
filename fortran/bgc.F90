#include "fabm_driver.h"
!define _REPLNAN_(X) X !changes back to original code
#define _REPLNAN_(X) nan_num(X)
! converts biological unit d-1 into physical FABM/driver unit s-1 for RHS
#define UNIT *1.1574074074E-5_rk  

!----------------------------------------
!	tame/bgc
!
!> @brief ** BioGeoChemical Equations **
!> This is the main BGC routine where right-hand-sides are calculated
!! for organic material (Det, DOM) with arbitrary element units such carbon, nitrogen, & phosphorus
!!  and arbitrary dissolved chemical species such as NO3, NH4, or CO2 as state variables.
!----------------------------------------
! !INTERFACE:
module tame_bgc

use fabm_types
use tame_types
use tame_functions

implicit none

 private
 ! !PUBLIC_DERIVED_TYPES:
 type,extends(type_base_model),public :: type_tame_bgc
    type (type_state_variable_id)      :: id_var(NUM_ELEM*2+NUM_CHEM) ! TODO : flexible num of DOM & POM 
    type (type_dependency_id) :: id_par,id_temp
  !	type (type_horizontal_dependency_id) :: id_taub
  	type (type_diagnostic_variable_id) :: id_din,id_rate !id_chla,id_GPP,id_NPP
    real(rk) :: remineral,hydrolysis,alloc_N,Nqual,CNref,DenitKNO3,denit,T_ref,rq10,dil
    integer :: tlim
  contains
    procedure :: initialize
    procedure :: do
  ! procedure :: do_bottom
  ! procedure :: get_sinking_rate
end type

integer :: det_index(NUM_ELEM), dom_index(NUM_ELEM), dix_index(NUM_CHEM)
type (type_tame_chemical) :: dix
integer :: num_chemicals= NUM_CHEM
integer :: num_elements = NUM_ELEM !,parameter
integer :: TransIndex_DOMDIX(NUM_ELEM), TransIndex2_DOMDIX(1,2)

!----------------------------------------
contains
!----------------------------------------!BOP
!
! !IROUTINE: Initialise the tame/bgc model
! Reading namelist and registration of variables with FABM
!
subroutine initialize(self,configunit)
 class (type_tame_bgc), intent(inout), target :: self
 integer,		intent(in)		:: configunit
 integer :: i, i0, n
 character :: elem
 call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
 call self%register_dependency(self%id_temp, standard_variables%temperature)

 !call self%register_state_dependency(self%id_phy, 'phy','','' ) !, scale_factor=days_per_sec
 call self%get_parameter(self%remineral, 'remineral','d-1','DOM remineralisation rate', default=0.1_rk )
 call self%get_parameter(self%hydrolysis, 'hydrolysis','d-1','detritus hydrolysis rate', default=0.05_rk )
 call self%get_parameter(self%alloc_N, 'alloc_N','-','nh4 - NO3 product ratio remineralisation', default=0.5_rk )
 call self%get_parameter(self%Nqual, 'Nqual','-','OM fraction w quality prop to N:Cratio ', default=1.0_rk )
 call self%get_parameter(self%CNref, 'CNref','Redfield','POM quality relative to carbon : nitrogen ratio (mol C/mol N)', default=6.625_rk )
 call self%get_parameter(self%DenitKNO3, 'DenitKNO3','mmol N/m3','half-saturation NO3 denitrification', default=1.0_rk )
 call self%get_parameter(self%denit, 'denit','1/d','pelagic denitrification rate', default=0.01_rk )
 call self%get_parameter(self%T_ref, 'T_ref','Kelvin','reference temperature', default=293.0_rk )
 call self%get_parameter(self%rq10, 'rq10','-','temperature dependence Q10', default=0.175_rk )
 !call self%get_parameter(self%dil, 'dil','-','dilution rate', default=0.0_rk)
 call self%get_parameter(self%tlim, 'tlim','0: none, 1: flagellate-style, 2: cyanobacteria-style','temperature limitation of growth', default=0 )
 !call self%register_dependency(self%id_taub, standard_variables%bottom_stress)
 call self%register_diagnostic_variable(self%id_din, 'DIN', 'mmol-N m-3', 'dissolved inorganic nitrogen')
 call self%register_diagnostic_variable(self%id_rate, 'rate', 'd-1', 'rate')
 !call self%register_diagnostic_variable(self%id_NPP, 'NPP',  'mmol/m3/d',   'net primary production')

do i = 1,num_chemicals !
    call self%register_state_variable(self%id_var(i), chemicals(i),'mmol m-3',chemicals(i))
! *,chemicals(i)
end do
i0 = num_chemicals

! set indices of element vectors and pointers
do i = 1,num_elements !
  det_index(i) = i0+2*i-1
  dom_index(i) = i0+2*i
  elem = ElementList(i:i)
  call self%register_state_variable(self%id_var(det_index(i)), 'det_' // elem,'mmol-' // elem // ' m-3','Detritus ' // trim(ElementName(i)))
  call self%register_state_variable(self%id_var(dom_index(i)), 'dom_' // elem,'mmol-' // elem // ' m-3','Dissolved Organic ' // trim(ElementName(i)))
  ! print *,det_index(i), ElementList(i:i),dom_index(i)
end do

end subroutine initialize
!----------------------------------------
!

! !IROUTINE: Right hand sides of tame/bgc model
!
! !INTERFACE:
 subroutine do(self,_ARGUMENTS_DO_)
 _DECLARE_ARGUMENTS_DO_
! !LOCAL VARIABLES:
 class (type_tame_bgc),intent(in) :: self

  real(rk) :: remin_chemical(NUM_ELEM),qualDetv(NUM_ELEM),qualDOMv(NUM_ELEM) !, target
  real(rk) :: dom_element(NUM_ELEM), det_element(NUM_ELEM), dix_chemical(NUM_CHEM)
  type (type_tame_elem)     :: det,dom
  real(rk) :: rhs(NUM_CHEM+2*NUM_ELEM)
	real(rk) :: par, temp, ddix, denitrate, nitrate, qualDet, qualDOM
!  real(rk) :: phy,zoo
  type (type_tame_env)      :: env
 !type (type_tame_switch) :: mswitch
  type (type_tame_sensitivities) :: sens
!type (stoich_pointer), dimension(5)::elem ! struct-pointer addressing elements wthin loops
! --- LOCAL MODEL VARIABLES:
  integer  :: i, j, i0, Index_Det_No_NorC, Index_DOX_No_NorC
  real(rk) :: hydrol_rate, remin_rate, remineral , hydrolysis  ! Temp dependent remineralisation and hydrolysis rates
  real(rk) :: aggreg_rate ! particle aggregation
  logical  :: out = .true.
!   if(36000.eq.secondsofday .and. mod(julianday,1).eq.0 .and. outn) out=.true.
! The following is the inverse of seconds_per_day 1/86400
! 
! transfer matrix of indices for DOX to produced DIX/chemical  (e.g. IndexOf_DOP->IndexOf_PO4) TODO: move to tame_types?
TransIndex_DOMDIX(1) = 0    ! C: no chemical if CO" is not resolved, see "chemicals" above
TransIndex_DOMDIX(2) = -1   ! N: -1 partitioned between NO3-chemical 1 and NH4-chemical 2
TransIndex_DOMDIX(3) = 3    ! P: 3rd chemcal PO4
! partitioning from DON break-down to N-species (NO3, NH4,..) TODO: move to tame_types?
TransIndex2_DOMDIX(1,1) = 1 ! partitioning of DON to 1st and 2nd chemical (NO3, NH4)
TransIndex2_DOMDIX(1,2) = 2 !

do i = 1,NUM_CHEM !
   call set_chem_pointer(dix,dix_chemical,chemicals(i), i)
!   print *,i,chemicals(i),dix%index%NO3,dix%index%NH4,dix_chemical(i)
   dix_index(i) = i
end do

i0 = NUM_CHEM
do i = 1,num_elements !
  ! internally link the element resolving vectors of OM  
  ! print *,ElementList(i:i),i
  call set_pointer(det,det_element,ElementList(i:i), i)
  call set_pointer(dom,dom_element,ElementList(i:i), i)
  ! print *,det_index(i),det%index%C
  det_index(i) = i0+2*i-1
  dom_index(i) = i0+2*i
end do

! TODO: add and solve pointer issue
!!if (associated(det%index%P)) Index_Det_No_NorC = det%index%P
!!if (associated(dom%index%P)) Index_DOX_No_NorC = dom%index%P
 !allocate(remin_chemical(num_chemicals), stat=rc) allocate(qualDetv(num_elements), stat=rc)

 _LOOP_BEGIN_

! First retrieve current (local) state  variable values
!---------- GET for each state variable ----------
do i = 1,NUM_CHEM ! e.g., CO2, NO3, NH4 (PO4)
!  if (_AVAILABLE_(self%id_dix(i))) then ddix
   _GET_(self%id_var(i), dix_chemical(i))  ! Dissolved Inorganic Nutrient DIX in mmol-X/m**3
!  end if
end do
i0=i
!  retrieve OM variables for each element
do i = 1,num_elements ! e.g., C, N, P (Si, Fe)
  _GET_(self%id_var(det_index(i)), det_element(i))  ! Detritus Organics in mmol-C/m**3
  _GET_(self%id_var(dom_index(i)), dom_element(i))  ! Dissolved Organics in mmol-C/m**3
!   print *,'det_',ElementList(i:i),det_element(i)
end do

!---------- get ambient conditions ----------
_GET_(self%id_temp, env%temp)  ! water temperature
!_GET_(self%id_temp, ddix)  ! water temperaturedummy

call calc_sensitivities(sens,env,self%rq10,self%T_ref)
! call photosynthesis(self,sens,phy,nut,uptake,exud,acclim)

!___________________________________________________________________
!  ---  POM&DOM quality, relative to Refield ?
! TODO: merge POM and DOM !
!  Nqual = 1 full N:C dependency   0: only fresh material
qualDet   = (1.0_rk-self%Nqual) + self%Nqual * det%N /(det%C + small) * self%CNref
qualDOM   = (1.0_rk-self%Nqual) + self%Nqual * dom%N /(dom%C + small) * self%CNref
! distribute  preferential degradation rate to elements (N:fast; C:quality dep; others: intermediate)

! TODO: compress
if (index(ElementList,'N') .gt. 0)  then 
  qualDetv(det%index%N) = 1.0_rk
  qualDOMv(dom%index%N)  = 1.0_rk
endif
if (index(ElementList,'C') .gt. 0)  then 
  qualDetv(det%index%C) = qualDet
  qualDOMv(dom%index%C)  = qualDOM
endif
if (index(ElementList,'P') .gt. 0)  then 
  qualDetv(det%index%P) = (1.0_rk + qualDet)/2
  qualDOMv(dom%index%P)  = (1.0_rk + qualDOM)/2
endif
!_____________________________________________________________________________
!     denitrification
! pelagic N-loss by denitrification, emulating benthic pool and suboxic micro-environments
! check for the carbon index
! calculate substrate (all OM)
if (self%denit .gt. 0.0_rk) then
!! if (associated(det%C)) then ! POC Glud LO 2015 (suboxic spots in particles)
  ! calculate oxidant (NO3)
  if (associated(dix%NO3)) then ! TODO: add nitrite NO2
    nitrate = dix%NO3
  elseif (associated(dix%DIN)) then
    nitrate = dix%DIN    ! lower denitrication self%denit
  else
    nitrate = 0.1_rk ! TODO: replace by SMALL
  endif
  denitrate = self%denit * sens%f_T * det%C * qualDet * (1.0_rk-exp(-nitrate/self%DenitKNO3))
!else
!  denitrate = 0.0_rk
endif

!  ---  hydrolysis & remineralisation rate (temp dependent)
hydrol_rate  = self%hydrolysis * sens%f_T
remin_rate   = self%remineral  * sens%f_T
!print *,'remin_rate=',remin_rate,sens%f_T,self%rq10,env%temp,self%T_ref

!________________________________________________________________________________
!
!  --- DETRITUS C
! TODO: link to other modules
! get det_prod(i)  =   !(floppZ%C + zoo_mort) * zoo%C + aggreg_rate * phy%C
! get dom_prod(i)  =   ! exud%C * phy%C

do i = 1,num_elements ! e.g., N  ( C, Si, Fe, P)
  hydrolysis = hydrol_rate * qualDetv(i) * det_element(i)
  remineral  = remin_rate  * qualDOMv(i) * dom_element(i)
  rhs(det_index(i)) =  - hydrolysis
  rhs(dom_index(i)) =  + hydrolysis - remineral

  ! transfer matrix of remineralised DOX to DIX
  j = TransIndex_DOMDIX(i)
  !print *,i,j,':',remineral,':',remin_rate, qualDOMv(i) ,dom_element(i)

  if (j .gt. 0) then
     remin_chemical(j) = remineral
  elseif (j .lt. 0) then ! partitioning between NO3 and NH4
     remin_chemical(TransIndex2_DOMDIX(-j,1)) = remineral * self%alloc_N
     remin_chemical(TransIndex2_DOMDIX(-j,2)) = remineral * (1.0_rk - self%alloc_N)
 !    print *,-j,TransIndex2_DOMDIX(-j,1),TransIndex2_DOMDIX(-j,2),remineral,self%alloc_N

  endif
end do
! add denitrification of POC(!) Glud et al LO 2015 (suboxic spots in particles)
if(self%denit .gt. 0.0_rk) then ! TODO: associated(det%C) .AND. 
    j = det_index(det%index%C)
    rhs(j) = rhs(j) - denitrate
endif
!Index_DetN Index_DON Index_DetC Index_DOC Index_Det_No_NorC Index_DOX_No_NorC
! TODO: link to other modules
! GET nut_prod(i)  =   !  -uptake%N * phy%C + lossZ%N * zoo%C

! here, nutrients are only remineralised (e.g., uptake in tame_phy)
do i = 1,num_chemicals
  rhs(dix_index(i)) = remin_chemical(i) !+ nut_prod(i)
!  print *,dix_index(i),' rhs=',remin_chemical(i)
end do

!  chemostat mode
#if 0
if (self%dil .gt. 0.0_rk) then
  do i = 1,num_elements ! e.g., N  ( C, Si, Fe, P)
    rhs(det_index(i)) = rhs(det_index(i)) - self%dil * det_element(i)
    rhs(dom_index(i)) = rhs(dom_index(i)) - self%dil * dom_element(i)
  end do
  do i = 1,num_chemicals
    j = dix_index(i)
    !Error: ‘chemical0’ at (1) is not a member of the ‘type_tame_chemical’ structure; did you mean ‘chemical’?
    !/Users/Lemmen/devel/fabm/generalized-aquatic-ecosystem-model/fortran/bgc.F90:266:122:
    ! todo KAI
    rhs(j) = rhs(j) + self%dil * (dix_chemical0(i) - dix_chemical(i))
  end do
endif
#endif

! tell FABM about right hand sides ....
do i = 1,dom_index(NUM_ELEM)
  _ADD_SOURCE_(self%id_var(i), rhs(i) UNIT)
!  print *,i,dom_index(NUM_ELEM)
end do

_SET_DIAGNOSTIC_(self%id_din, dix%NO3+dix%NH4)    !average
_SET_DIAGNOSTIC_(self%id_rate, remin_rate)       

!_SET_DIAGNOSTIC_(self%id_vphys, exp(-self%sink_phys*phy%relQ%N * phy%relQ%P))       !average
! experimental formulation for emulating P-adsorption at particles in the water column and at the bottom interface
!_GET_HORIZONTAL_(self%id_o2flux, flO2)!_GET_HORIZONTAL_(self%id_oduflux, flODU)!_GET_HORIZONTAL_(self%id_zmax, zmax)  ! max depth
!aPO4 = (flODU-flO2)/(zmax+self%small)
!_SET_DIAGNOSTIC_(self%id_vphys, aPO4)       !average Temporary_diagnostic_
!________________________________________________________________________________

!if (self%BGC0DDiagOn) then
!  _SET_DIAGNOSTIC_(self%id_qualDOM, _REPLNAN_(qualDOM))      !average Quality_of_DOM_
!end if

_LOOP_END_
end subroutine do

end module tame_bgc
