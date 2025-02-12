#include "fabm_driver.h"
# define NUM_ELEM 3
# define NUM_CHEM 3
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
  real(rk) :: remineral,hydrolysis,alloc_N,Nqual,CNref,DenitKno3,denit,T_ref,rq10
  integer :: tlim

  type (type_tame_elem)     :: det,dom
  type (type_tame_chemical) :: dix
  integer :: num_chemicals= NUM_CHEM
  integer :: num_elements = NUM_ELEM !,parameter
  integer :: det_vec_index(NUM_ELEM), dom_vec_index(NUM_ELEM), dix_vec_index(NUM_CHEM)
  character(len = 6) ::  element= 'CNPSF'
  integer :: TransIndex_DOMDIX(NUM_ELEM), TransIndex_DON

contains
	procedure :: initialize
	procedure :: do
 	procedure :: set_element
 	procedure :: set_pointer
! procedure :: do_bottom
! procedure :: get_sinking_rate
end type
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

 real(rk),parameter :: secs_per_day = 86400._rk
 real(rk),parameter :: small = 1.E-4_rk
 integer :: i, i0, n
 character(len=3) :: chemicals(NUM_CHEM)

 chemicals(1)='no3'
 chemicals(2)='nh4'
 chemicals(3)='po4'
 !real(rk), pointer :: carbon_ptr(number_of_carbon_elements) => null

!call self%register_state_variable(id_dix(1), 'carbon')
! TODO rewrite as loop 
!call self%register_state_variable(id_dix(3), 'po4')
do i = 1,self%num_chemicals ! 
    call self%register_state_variable(self%id_var(i), chemicals(i))
    self%dix_vec_index(i) = i
    print *,chemicals(i)
end do
i0 = i
!dix%no3 => dix%chemical(1)

! partitioning of DIN production from DON between N-species (NO3, NH4,..)
self%TransIndex_DOMDIX(1) = -1 ! N:1 partitioned between NO3-chemical 1 and NH4-chemical 2
self%TransIndex_DON = 2
self%TransIndex_DOMDIX(2) = 3 ! 

! set indices of element vectors and pointers
do i = 1,self%num_elements ! 
  call set_pointer(self%det,self%element(i:i), i)
  call set_pointer(self%dom,self%element(i:i), i)
  self%det_vec_index(i) = i0+2*i-1
  self%dom_vec_index(i) = i0+2*i
  call self%register_state_variable(self%id_var(self%det_vec_index(i)), 'det%' // self%element(i:i))
  call self%register_state_variable(self%id_var(self%dom_vec_index(i)), 'dom%' // self%element(i:i))
end do
 !call self%register_state_dependency(self%id_phy, 'phy','','' )
 call self%get_parameter(self%remineral, 'remineral','1/d','DOM remineralisation rate', default=0.1 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%hydrolysis, 'hydrolysis','1/d','detritus hydrolysis rate', default=0.05 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%alloc_N, 'alloc_N','-','nh4 - no3 product ratio remineralisation', default=0.5 )
 call self%get_parameter(self%Nqual, 'Nqual','-','OM fraction w quality prop to N:Cratio ', default=1. )
 call self%get_parameter(self%CNref, 'CNref','Redfield','POM quality relative to carbon : nitrogen ratio (mol C/mol N)', default=6.625 )
 call self%get_parameter(self%DenitKno3, 'DenitKno3','mmol N/m3','half-saturation no3 denitrification', default=1. )
 call self%get_parameter(self%denit, 'denit','1/d','pelagic denitrification rate', default=0.01 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%T_ref, 'T_ref','Kelvin','reference temperature', default=293. )
 call self%get_parameter(self%rq10, 'rq10','-','temperature dependence Q10', default=0.175 )
 call self%get_parameter(self%tlim, 'tlim','0: none, 1: flagellate-style, 2: cyanobacteria-style','temperature limitation of growth', default=0 )
 call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
 call self%register_dependency(self%id_temp, standard_variables%temperature)
 !call self%register_dependency(self%id_taub, standard_variables%bottom_stress)
 !call self%register_diagnostic_variable(self%id_chla, 'chla', 'mg chl a/m3', 'chlorophyll concentration')
 !call self%register_diagnostic_variable(self%id_NPP, 'NPP',  'mmol/m3/d',   'net primary production')

end subroutine initialize
!----------------------------------------
!

! !IROUTINE: Right hand sides of tame/bgc model
!
! !INTERFACE:
 subroutine do(self,_ARGUMENTS_DO_)
! !LOCAL VARIABLES:
 class (type_tame_bgc),intent(in) :: self

!  real(kind=rk), allocatable :: remin_chemical(:),qualDetv(:),qualDOMv(:) !, target
!  real(kind=rk), allocatable :: det_prod(:),dom_prod(:),nut_prod(:)
!  type (type_tame_elem), allocatable :: rhs_det,rhs_dom
  real(rk) :: remin_chemical(NUM_ELEM),qualDetv(NUM_ELEM),qualDOMv(NUM_ELEM) !, target
  real(rk) :: det_prod(NUM_ELEM),dom_prod(NUM_ELEM),nut_prod(NUM_ELEM)
  real(rk) :: rhsv(NUM_CHEM+2*NUM_ELEM)
	real(rk) :: par, temp
!  real(rk) :: phy,zoo
  type (type_tame_env)      :: env
 !type (type_tame_switch) :: mswitch
  type (type_tame_sensitivities) :: sens
!type (stoich_pointer), dimension(5)::elem ! struct-pointer addressing elements wthin loops
! --- LOCAL MODEL VARIABLES:
  integer  :: i, j, Index_NO3, Index_Det_No_NorC, Index_DOX_No_NorC
  real(rk) :: remineral_rate , hydrolysis_rate  ! Temp dependent remineralisation and hydrolysis rates
  real(rk) :: aggreg_rate ! particle aggregation 
  logical  :: out = .true.
!   if(36000.eq.secondsofday .and. mod(julianday,1).eq.0 .and. outn) out=.true.
#define UNIT *1.1574074074E-5_rk ! 1/86400

if (associated(self%det%index%P)) Index_Det_No_NorC = self%det%index%P 
if (associated(self%dom%index%P)) Index_DOX_No_NorC = self%dom%index%P 

 !allocate(remin_chemical(self%num_chemicals), stat=rc) allocate(qualDetv(self%num_elements), stat=rc)
 !allocate(qualDOMv(self%num_elements), stat=rc) allocate(det_prod(self%num_elements), stat=rc, source=0.0_rk)
 !allocate(dom_prod(self%num_elements), stat=rc, source=0.0_rk) allocate(nut_prod(self%num_elements), stat=rc, source=0.0_rk)
 !allocate(det%element(self%num_elements), stat=rc)
 
 _LOOP_BEGIN_

! First retrieve current (local) state  variable values
!---------- GET for each state variable ----------
do i = 1,self%num_chemicals ! e.g., CO2, NO3, NH4 (PO4)
!  if (_AVAILABLE_(self%id_dix(i))) then
   _GET_(self%id_var(i), self%dix%chemical(i))  ! Dissolved Inorganic Nutrient DIX in mmol-X/m**3
!  end if
end do
i0=i
!  retrieve OM variables for each element 
do i = 1,self%num_elements ! e.g., N  ( C, Si, Fe, P)
  _GET_(self%id_var(self%det_vec_index(i)), self%det%element(i))  ! Detritus Organics in mmol-C/m**3
  _GET_(self%id_var(self%dom_vec_index(i)), self%dom%element(i))  ! Dissolved Organics in mmol-C/m**3
end do
!do i = 1,num_phyclass ! e.g., DIA, FLA, CYA  (or 1,5,30,100 µm)
!  _GET_(self%id_phy(i), phy%class(i))  ! Detritus Organics in mmol-C/m**3
!end do
!---------- get ambient conditions ----------
_GET_(self%id_temp, env%temp)  ! water temperature
! _GET_(self%id_par, env%par)    ! light photosynthetically active radiation

!_SET_DIAGNOSTIC_(self%id_PAR_diag,env%par)         !average Photosynthetically_Active_Radiation_

call calc_sensitivities(sens,env,self%rq10,self%T_ref)
! call photosynthesis(self,sens,phy,nut,uptake,exud,acclim)

!___________________________________________________________________
!  ---  POM&DOM quality, relative to Refield ?
! TODO: merge POM and DOM !
!  Nqual = 1 full N:C dependency   0: only fresh material
qualDet   = (1.0d0-self%Nqual) + self%Nqual * self%det%N /(self%det%C + self%small) * self%CNref
qualDOM   = (1.0d0-self%Nqual) + self%Nqual * self%dom%N /(self%dom%C + self%small) * self%CNref
! distribute  preferential degradation rate to elements (N:fast; C:quality dep; others: intermediate)
if (associated(self%det%index%N))  qualDetv(self%det%index%N) = 1.0d0
if (associated(self%dom%index%N))  qualDOMv(self%dom%index%N)  = 1.0d0
if (associated(self%det%index%C))  qualDetv(self%det%index%C) = qualDet
if (associated(self%dom%index%C))  qualDOMv(self%dom%index%C)  = qualDOM
if (associated(Index_Det_No_NorC)) qualDetv(Index_Det_No_NorC) = (1.0d0 + qualDet)/2
if (associated(Index_DOX_No_NorC)) qualDOMv(Index_DOX_No_NorC) = (1.0d0 + qualDOM)/2
!_____________________________________________________________________________
!     denitrification 
! pelagic N-loss by denitrification, emulating benthic pool and suboxic micro-environments
! check for the carbon index 
! calculate substrate (all OM)
!if (associated(self%dom%C)) sum_OM = sum_OM + dom%C* qualDOM ! 

if (associated(self%det%C)) then ! POC Glud LO 2015 (suboxic spots in particles) 
  ! calculate oxidant (NO3)
  if (associated(self%dix%no3)) then ! TODO: add nitrite NO2
    nitrate = self%dix%no3
    Index_NO3 = 1
  elseif (associated(self%dix%din)) then
    nitrate = self%dix%din    ! lower denitrication self%denit
    Index_NO3 = 1
  else
    nitrate = 0.1_rk ! TODO: replace by SMALL
  endif
  denitrate = self%denit * sens%f_T * self%det%C * qualDet * (1.0d0-exp(-nitrate/self%DenitKNO3)) 
elseif
  denitrate = 0.0_rk
endif

!  ---  hydrolysis & remineralisation rate (temp dependent)
hydrolysis  = self%hydrolysis * sens%f_T 
remineral   = self%remineral  * sens%f_T 

!________________________________________________________________________________
!
!  --- DETRITUS C
! TODO: link to other modules
! get det_prod(i)  =   !(floppZ%C + zoo_mort) * zoo%C + aggreg_rate * phy%C  
! get dom_prod(i)  =   ! exud%C * phy%C 

do i = 1,self%num_elements ! e.g., N  ( C, Si, Fe, P)
  hydrolysis_rate = self%hydrolysis * qualDetv(i) * self%det%element(i:i)
  remineral_rate  = self%remineral  * qualDOMv(i) * self%dom%element(i:i)
  rhsv(self%det_vec_index(i)) = det_prod(i) - hydrolysis_rate
  rhsv(self%dom_vec_index(i)) = dom_prod(i) + hydrolysis_rate - remineral_rate
  
  ! transfer matrix of remineralised DOX to DIX
  j = self%TransIndex_DOMDIX(i)
  if (j .gt. 0) then
     remin_chemical(j) = remineral_rate
  elseif (j .lt. 0) then ! partitioning between NO3 and NH4
     remin_chemical(-j) = remineral_rate * self%alloc_N
     remin_chemical(self%TransIndex_DON) = remineral_rate * (1.0_rk - self%alloc_N)
  endif
end do
! add denitrification of POC(!) Glud et al LO 2015 (suboxic spots in particles)
if(associated(self%det%index%C)) rhsv(self%det_vec_index(i)) = rhsv(self%det_vec_index(i)) + denitrate 

!Index_DetN Index_DON Index_DetC Index_DOC Index_Det_No_NorC Index_DOX_No_NorC
! TODO: link to other modules
! GET nut_prod(i)  =   !  -uptake%N * phy%C + lossZ%N * zoo%C 

do i = 1,self%num_chemicals 
  rhsv(self%dix_vec_index(i)) = remin_chemical(i) + nut_prod(i) 
end do
!if(associated(self%Index_NO3)) rhsv%nut%chemical(self%Index_NO3) = rhsv%nut%chemical(self%Index_NO3) - 0.8_rk*denitrate 
!if(associated(self%Index_NO3)) rhsv%nut%chemical(self%Index_NO3) = rhsv%nut%chemical(self%Index_NO3) - 0.8_rk*denitrate 

!  chemostat mode 
if (_AVAILABLE_(self%dil)) then
  if (self%dil .gt. 0.0d0) then
    do i = 1,self%num_elements ! e.g., N  ( C, Si, Fe, P)
      rhsv(self%det_vec_index(i)) = rhsv(self%det_vec_index(i)) - self%dil * self%det%element(i)
      rhsv(self%dom_vec_index(i)) = rhsv(self%dom_vec_index(i)) - self%dil * self%dom%element(i)
    end do
    do i = 1,self%num_chemicals 
      rhsv(self%dix_vec_index(i)) = rhsv(self%dix_vec_index(i)) + self%dil * (self%dix%chemical0(i) - self%dix%chemical(i))
    end do
  endif
endif

! tell FABM about right hand sides ....
do i = 1,dom_vec_index(NUM_ELEM) 
  _ADD_SOURCE_(self%id_var(i), rhsv(i) UNIT)
end do

!_SET_DIAGNOSTIC_(self%id_vphys, exp(-self%sink_phys*phy%relQ%N * phy%relQ%P))       !average
! experimental formulation for emulating P-adsorption at particles in the water column and at the bottom interface
!_GET_HORIZONTAL_(self%id_o2flux, flO2)!_GET_HORIZONTAL_(self%id_oduflux, flODU)!_GET_HORIZONTAL_(self%id_zmax, zmax)  ! max depth
!aPO4 = (flODU-flO2)/(zmax+self%small)
!_SET_DIAGNOSTIC_(self%id_vphys, aPO4)       !average Temporary_diagnostic_
!________________________________________________________________________________

!define _REPLNAN_(X) X !changes back to original code
#define _REPLNAN_(X) nan_num(X)

if (self%BGC0DDiagOn) then
  _SET_DIAGNOSTIC_(self%id_qualDet, _REPLNAN_(qualDet))      !average Quality_of_POM_
  _SET_DIAGNOSTIC_(self%id_qualDOM, _REPLNAN_(qualDOM))      !average Quality_of_DOM_
end if

_LOOP_END_

end subroutine do

! set indices of elements vectors 

subroutine set_element(object, name, value)
  type(type_tame_elem), intent(inout) :: object
  character(*), intent(in) :: name
  integer, intent(in) :: value
  select case (name)
  case ('C')  ; object%C = value
  case ('N') ; object%N = value
  case ('P')  ; object%P = value
  case ('S')  ; object%Si = value
  case ('F')  ; object%Fe = value
  end select
end subroutine set_element

! set pointer to indexed elements vector 
subroutine set_pointer(object, name, value)
  type(type_tame_elem), intent(inout) :: object
  character(*), intent(in) :: name
  integer, intent(in) :: value
  select case (name)
  case ('C') 
    object%C => object%element(value)
    object%index%C = value
  case ('N')
    object%N => object%element(value) 
    object%index%N = value
  case ('P')
    object%P => object%element(value) 
    object%index%P = value
  case ('Si')
    object%Si => object%element(value) 
    object%index%Si = value
  case ('Fe')
    object%Fe => object%element(value) 
    object%index%Fe = value
  end select
end subroutine set_pointer

