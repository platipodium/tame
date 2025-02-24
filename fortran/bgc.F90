#include "fabm_driver.h"

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
    type (type_state_variable_id)      :: id_var(NUM_ELEM*2+NUM_CHEM)
    type (type_dependency_id) :: id_par,id_temp
  !	type (type_horizontal_dependency_id) :: id_taub
  !	type (type_diagnostic_variable_id) :: id_chla,id_GPP,id_NPP
    real(rk) :: remineral,hydrolysis,alloc_N,Nqual,CNref,DenitKno3,denit,T_ref,rq10,dil
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
integer :: TransIndex_DOMDIX(NUM_ELEM), TransIndex_DON

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
 character(len=3) :: chemicals(NUM_CHEM) = (/'no3','nh4','po4'/)

 call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
 call self%register_dependency(self%id_temp, standard_variables%temperature)

 !call self%register_state_dependency(self%id_phy, 'phy','','' )
 call self%get_parameter(self%remineral, 'remineral','1/d','DOM remineralisation rate', default=0.1_rk , scale_factor=days_per_sec)
 call self%get_parameter(self%hydrolysis, 'hydrolysis','1/d','detritus hydrolysis rate', default=0.05_rk , scale_factor=days_per_sec)
 call self%get_parameter(self%alloc_N, 'alloc_N','-','nh4 - no3 product ratio remineralisation', default=0.5_rk )
 call self%get_parameter(self%Nqual, 'Nqual','-','OM fraction w quality prop to N:Cratio ', default=1.0_rk )
 call self%get_parameter(self%CNref, 'CNref','Redfield','POM quality relative to carbon : nitrogen ratio (mol C/mol N)', default=6.625_rk )
 call self%get_parameter(self%DenitKno3, 'DenitKno3','mmol N/m3','half-saturation no3 denitrification', default=1.0_rk )
 call self%get_parameter(self%denit, 'denit','1/d','pelagic denitrification rate', default=0.01_rk , scale_factor=days_per_sec)
 call self%get_parameter(self%T_ref, 'T_ref','Kelvin','reference temperature', default=293.0_rk )
 call self%get_parameter(self%rq10, 'rq10','-','temperature dependence Q10', default=0.175_rk )
 !call self%get_parameter(self%dil, 'dil','-','dilution rate', default=0.0_rk)
 call self%get_parameter(self%tlim, 'tlim','0: none, 1: flagellate-style, 2: cyanobacteria-style','temperature limitation of growth', default=0 )
  !call self%register_dependency(self%id_taub, standard_variables%bottom_stress)
 call self%register_diagnostic_variable(self%id_din, 'DIN', 'mmol-N m-3', 'dissolved inorganic nitrogen')
! call self%register_diagnostic_variable(self%id_chla, 'chla', 'mg chl a/m3', 'chlorophyll concentration')
 !call self%register_diagnostic_variable(self%id_NPP, 'NPP',  'mmol/m3/d',   'net primary production')

do i = 1,num_chemicals !
    call self%register_state_variable(self%id_var(i), chemicals(i),'dummy unit','dummy long name')
    print *,chemicals(i)
end do
i0 = i
! partitioning of DIN production from DON between N-species (NO3, NH4,..)
TransIndex_DOMDIX(1) = -1 ! N:1 partitioned between NO3-chemical 1 and NH4-chemical 2
TransIndex_DON = 2
TransIndex_DOMDIX(2) = 3 !

Index_DIN
do i = 1,num_chemicals !
   
chemicals
  if (associated(dix%no3)) then ! TODO: add nitrite NO2
    nitrate = dix%no3
  elseif (associated(dix%din)) then

! set indices of element vectors and pointers
do i = 1,num_elements !
!  call set_pointer(det,ElementList(i:i), i)
!  call set_pointer(dom,ElementList(i:i), i)
!  det_index(i) = i0+2*i-1
!  dom_index(i) = i0+2*i
  call self%register_state_variable(self%id_var(det_index(i)), 'det_' // ElementList(i:i),'dummy unit','dummy long name')
  call self%register_state_variable(self%id_var(dom_index(i)), 'dom_' // ElementList(i:i),'dummy unit','dummy long name')
  ! ERROR: Program received signal SIGSEGV, Segmentation fault.
  ! #3  0x00005555555ffd9e in __tame_bgc_MOD_initialize ()

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

!  real(kind=rk), allocatable :: remin_chemical(:),qualDetv(:),qualDOMv(:) !, target
!  real(kind=rk), allocatable :: det_prod(:),dom_prod(:),nut_prod(:)
  real(rk) :: remin_chemical(NUM_ELEM),qualDetv(NUM_ELEM),qualDOMv(NUM_ELEM) !, target
  !real(rk) :: det_prod(NUM_ELEM),dom_prod(NUM_ELEM),nut_prod(NUM_ELEM)
  real(rk) :: dom_element(NUM_ELEM), det_element(NUM_ELEM)
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
  real(rk) :: remineral , hydrolysis  ! Temp dependent remineralisation and hydrolysis rates
  real(rk) :: aggreg_rate ! particle aggregation
  logical  :: out = .true.
!   if(36000.eq.secondsofday .and. mod(julianday,1).eq.0 .and. outn) out=.true.

! The following is the inverse of seconds_per_day 1/86400
#define UNIT *1.1574074074E-5_rk

! 
do i = 1,num_chemicals !
    dix_index(i) = i
end do

i0 = num_chemicals
do i = 1,num_elements !
  ! internally link the element resolving vectors of OM  
  call set_pointer(det,det_element,ElementList(i:i), i)
  call set_pointer(dom,dom_element,ElementList(i:i), i)
  det_index(i) = i0+2*i-1
  dom_index(i) = i0+2*i
end do

! TODO: add and solve pointer issue
!!if (associated(det%index%P)) Index_Det_No_NorC = det%index%P
!!if (associated(dom%index%P)) Index_DOX_No_NorC = dom%index%P

 !allocate(remin_chemical(num_chemicals), stat=rc) allocate(qualDetv(num_elements), stat=rc)
 !allocate(qualDOMv(num_elements), stat=rc) allocate(det_prod(num_elements), stat=rc, source=0.0_rk)
 !allocate(dom_prod(num_elements), stat=rc, source=0.0_rk) allocate(nut_prod(num_elements), stat=rc, source=0.0_rk)
 !allocate(det_element(num_elements), stat=rc)

 _LOOP_BEGIN_

! First retrieve current (local) state  variable values
!---------- GET for each state variable ----------
do i = 1,num_chemicals ! e.g., CO2, NO3, NH4 (PO4)
!  if (_AVAILABLE_(self%id_dix(i))) then ddix
   _GET_(self%id_var(i), dix%chemical(i))  ! Dissolved Inorganic Nutrient DIX in mmol-X/m**3
!  end if
end do
i0=i
!  retrieve OM variables for each element
do i = 1,num_elements ! e.g., N  ( C, Si, Fe, P)
  _GET_(self%id_var(det_index(i)), det_element(i))  ! Detritus Organics in mmol-C/m**3
  _GET_(self%id_var(dom_index(i)), dom_element(i))  ! Dissolved Organics in mmol-C/m**3
end do
!do i = 1,num_phyclass ! e.g., DIA, FLA, CYA  (or 1,5,30,100 µm)
!  _GET_(self%id_phy(i), phy%class(i))  ! Detritus Organics in mmol-C/m**3
!end do
!---------- get ambient conditions ----------
_GET_(self%id_temp, env%temp)  ! water temperature
!_GET_(self%id_temp, ddix)  ! water temperaturedummy
! _GET_(self%id_par, env%par)    ! light photosynthetically active radiation

!_SET_DIAGNOSTIC_(self%id_PAR_diag,env%par)         !average Photosynthetically_Active_Radiation_

call calc_sensitivities(sens,env,self%rq10,self%T_ref)
! call photosynthesis(self,sens,phy,nut,uptake,exud,acclim)

!___________________________________________________________________
!  ---  POM&DOM quality, relative to Refield ?
! TODO: merge POM and DOM !
!  Nqual = 1 full N:C dependency   0: only fresh material
qualDet   = (1.0_rk-self%Nqual) + self%Nqual * det%N /(det%C + small) * self%CNref
qualDOM   = (1.0_rk-self%Nqual) + self%Nqual * dom%N /(dom%C + small) * self%CNref
! distribute  preferential degradation rate to elements (N:fast; C:quality dep; others: intermediate)

! TODO: add and solve pointer issue
!if (associated(det%index%N))  qualDetv(det%index%N) = 1.0_rk
!if (associated(dom%index%N))  qualDOMv(dom%index%N)  = 1.0_rk
!if (associated(det%index%C))  qualDetv(det%index%C) = qualDet
!if (associated(dom%index%C))  qualDOMv(dom%index%C)  = qualDOM
!if (associated(Index_Det_No_NorC)) qualDetv(Index_Det_No_NorC) = (1.0_rk + qualDet)/2
!if (associated(Index_DOX_No_NorC)) qualDOMv(Index_DOX_No_NorC) = (1.0_rk + qualDOM)/2
!_____________________________________________________________________________
!     denitrification
! pelagic N-loss by denitrification, emulating benthic pool and suboxic micro-environments
! check for the carbon index
! calculate substrate (all OM)
!if (associated(dom%C)) sum_OM = sum_OM + dom%C* qualDOM !

!! if (associated(det%C)) then ! POC Glud LO 2015 (suboxic spots in particles)
  ! calculate oxidant (NO3)
  if (associated(dix%no3)) then ! TODO: add nitrite NO2
    nitrate = dix%no3
  elseif (associated(dix%din)) then
    nitrate = dix%din    ! lower denitrication self%denit
  else
    nitrate = 0.1_rk ! TODO: replace by SMALL
  endif
  denitrate = self%denit * sens%f_T * det%C * qualDet * (1.0_rk-exp(-nitrate/self%DenitKNO3))
!else
!  denitrate = 0.0_rk
!endif

!  ---  hydrolysis & remineralisation rate (temp dependent)
hydrolysis  = self%hydrolysis * sens%f_T
remineral   = self%remineral  * sens%f_T

!________________________________________________________________________________
!
!  --- DETRITUS C
! TODO: link to other modules
! get det_prod(i)  =   !(floppZ%C + zoo_mort) * zoo%C + aggreg_rate * phy%C
! get dom_prod(i)  =   ! exud%C * phy%C

do i = 1,num_elements ! e.g., N  ( C, Si, Fe, P)
  hydrolysis = self%hydrolysis * qualDetv(i) * det_element(i)
  remineral  = self%remineral  * qualDOMv(i) * dom_element(i)
  rhs(det_index(i)) =  - hydrolysis
  rhs(dom_index(i)) =  + hydrolysis - remineral

  ! transfer matrix of remineralised DOX to DIX
  j = TransIndex_DOMDIX(i)
  if (j .gt. 0) then
     remin_chemical(j) = remineral
  elseif (j .lt. 0) then ! partitioning between NO3 and NH4
     remin_chemical(-j) = remineral * self%alloc_N
     remin_chemical(TransIndex_DON) = remineral * (1.0_rk - self%alloc_N)
  endif
end do
! add denitrification of POC(!) Glud et al LO 2015 (suboxic spots in particles)
if(associated(det%C)) then
    j = det_index(det%index%C)
    rhs(j) = rhs(j) - denitrate
endif
!Index_DetN Index_DON Index_DetC Index_DOC Index_Det_No_NorC Index_DOX_No_NorC
! TODO: link to other modules
! GET nut_prod(i)  =   !  -uptake%N * phy%C + lossZ%N * zoo%C

! here, nutrients are only remineralised (e.g., uptake in tame_phy)
do i = 1,num_chemicals
  rhs(dix_index(i)) = remin_chemical(i) !+ nut_prod(i)
end do

!  chemostat mode
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
    !rhs(j) = rhs(j) + self%dil * (dix%chemical0(i) - dix%chemical(i))
  end do
endif

! tell FABM about right hand sides ....
do i = 1,dom_index(NUM_ELEM)
  _ADD_SOURCE_(self%id_var(i), rhs(i) UNIT)
end do


if Index_DIN  
do i = 1,dom_index(NUM_ELEM)
  _ADD_SOURCE_(self%id_var(i), rhs(i) UNIT)
end do

_SET_DIAGNOSTIC_(self%id_din, )       !average

!_SET_DIAGNOSTIC_(self%id_vphys, exp(-self%sink_phys*phy%relQ%N * phy%relQ%P))       !average
! experimental formulation for emulating P-adsorption at particles in the water column and at the bottom interface
!_GET_HORIZONTAL_(self%id_o2flux, flO2)!_GET_HORIZONTAL_(self%id_oduflux, flODU)!_GET_HORIZONTAL_(self%id_zmax, zmax)  ! max depth
!aPO4 = (flODU-flO2)/(zmax+self%small)
!_SET_DIAGNOSTIC_(self%id_vphys, aPO4)       !average Temporary_diagnostic_
!________________________________________________________________________________

!define _REPLNAN_(X) X !changes back to original code
#define _REPLNAN_(X) nan_num(X)

!if (self%BGC0DDiagOn) then
!  _SET_DIAGNOSTIC_(self%id_qualDet, _REPLNAN_(qualDet))      !average Quality_of_POM_
!  _SET_DIAGNOSTIC_(self%id_qualDOM, _REPLNAN_(qualDOM))      !average Quality_of_DOM_
!end if

_LOOP_END_
end subroutine do

end module tame_bgc
