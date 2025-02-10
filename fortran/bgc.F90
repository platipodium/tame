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


 real(kind=rk), allocatable, target ::  id_dix(:)
!! type (type_tame_elem_id)  :: id_det,id_dom
!! type (type_tame_elem_index)  :: Index_Det,Index_DOM
!	type (type_state_variable_id) :: id_no3,id_nh4,id_o2,id_po4
!	type (type_state_variable_id) :: id_phy,id_zoo
!!	type (type_dependency_id) :: id_par,id_temp
!	type (type_horizontal_dependency_id) :: id_taub
!	type (type_diagnostic_variable_id) :: id_chla,id_GPP,id_NPP
!! real(rk) :: remineral,hydrolysis,alloc_N,Nqual,CNref,DenitKno3,denit,T_ref,rq10
!! integer :: tlim

 contains
	procedure :: initialize
	procedure :: do
 	procedure :: set_element
 	procedure :: set_pointer
! procedure :: do_bottom
! procedure :: get_sinking_rate
 end type
!EOP
!----------------------------------------
 CONTAINS
!----------------------------------------!BOP
!
! !IROUTINE: Initialise the tame/bgc model
!
! !INTERFACE:
 subroutine initialize(self,configunit)
!
! !DESCRIPTION:
! Reading namelist and registration of variables with FABM
!
! !INPUT PARAMETERS: class(type_tame_bgc),intent(inout),target :: self
 integer,		intent(in)		:: configunit
 #define NUM_ELEM 3

 real(rk),parameter :: secs_per_day = 86400._rk
 real(rk),parameter :: small = 1.E-4_rk
 integer :: TransIndex_DOMDIX(NUM_ELEM), TransIndex_DON
 integer :: i, Index_Det_No_NorC, Index_DOX_No_NorC
 
 !real(rk), pointer :: carbon_ptr(number_of_carbon_elements) => null
 integer,parameter :: num_chemicals= 3
 integer,parameter :: num_elements = NUM_ELEM
 character(len = 6) ::  element
 element = 'CNPSF'
 allocate(id_dix(num_chemicals), stat=rc)
 !if rc /= 0 goto 99
 allocate(dix%chemical(num_chemicals), stat=rc)
 !if rc /= 0 goto 99

!call self%register_state_variable(id_dix(1), 'carbon')
! TODO rewrite as loop 
call self%register_state_variable(id_dix(1), 'no3')
call self%register_state_variable(id_dix(2), 'nh4')
call self%register_state_variable(id_dix(3), 'po4')

! TODO rewrite as loop 
dix%no3 => dix%chemical(1)
dix%nh4 => dix%chemical(2)
dix%po4 => dix%chemical(3)

!carbon_id_ptr(1) => id_dix(1)
! partitioning of DIN production from DON between N-species (NO3, NH4,..)
TransIndex_DOMDIX(1) = -1 ! N:1 partitioned between NO3-chemical 1 and NH4-chemical 2
TransIndex_DON = 2
TransIndex_DOMDIX(2) = 3 ! 

! set indices of elements vectors and pointers
do i = 1,num_elements ! 
  set_element(Index_Det,element(i), i)
  set_element(Index_DOM,element(i), i)
  set_pointer(id_det,element(i), i)
  set_pointer(id_dom,element(i), i)
end do

! Index_Det_No_NorC, Index_DOX_No_NorC

!allocate(id_det%element(num_elements), stat=rc) allocate(id_dom%element(num_elements), stat=rc) !if rc /= 0 goto 99
!id_det%C => id_det%element(1) id_det%N => id_det%element(2)
!id_dom%C => id_dom%element(1) id_dom%N => id_dom%element(2)

 !call self%register_state_variable(self%, 'no3','mmol N/m3','nitrate', default=30.0 )
 !call self%register_state_variable(self%id_nh4, 'nh4','mmol N/m3','ammonium', default=2.0 )
 ! TODO rewrite as loop  ?? ::
 call self%register_state_variable(self%id_det%C, 'det%C','mmol C/m3','detritus carbon', default=0.0 )
 call self%register_state_variable(self%id_det%N, 'det%N','mmol N/m3','detritus nitrogen', default=0.0 )
 call self%register_state_variable(self%id_det%P, 'det%P','mmol P/m3','detritus phosphorus', default=0.0 )
 call self%register_state_variable(self%id_dom%C, 'dom%C','mmol C/m3','Dissolved Organic Carbon', default=0.0 )
 call self%register_state_variable(self%id_dom%N, 'dom%N','mmol N/m3','Dissolved Organic Nitrogen', default=0.0 )
 call self%register_state_variable(self%id_dom%P, 'dom%P','mmol P/m3','Dissolved Organic Phosphorus', default=0.0 )
 !call self%register_state_variable(self%id_o2, 'o2','mmol O2/m3','oxygen', default=280.0 )
 !call self%register_state_variable(self%id_po4, 'po4','mmol P/m3','phosphate', default=2.0 )
 !call self%register_state_dependency(self%id_phy, 'phy','','' )
 !call self%register_state_dependency(self%id_zoo, 'zoo','','' )
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
 class (type_hzg_tame),intent(in) :: self

!  real(kind=rk), allocatable :: remin_chemical(:),qualDetv(:),qualDOMv(:) !, target
!  real(kind=rk), allocatable :: det_prod(:),dom_prod(:),nut_prod(:)
!  type (type_tame_elem), allocatable :: rhs_det,rhs_dom
  real(rk) :: remin_chemical(NUM_ELEM),qualDetv(NUM_ELEM),qualDOMv(NUM_ELEM) !, target
  real(rk) :: det_prod(NUM_ELEM),dom_prod(NUM_ELEM),nut_prod(NUM_ELEM)  
	real(rk) :: par, temp
!  real(rk) :: phy,zoo
  type (type_tame_chemicals):: dix, rhs_nut(NUM_ELEM)
  type (type_tame_elem)     :: dom, det, rhs_det(NUM_ELEM), rhs_dom(NUM_ELEM)
  type (type_tame_env)      :: env
 !type (type_tame_switch) :: mswitch
  type (type_tame_sensitivities) :: sens
!type (stoich_pointer), dimension(5)::elem ! struct-pointer addressing elements wthin loops
! --- LOCAL MODEL VARIABLES:
  integer  :: i, j, Index_NO3
  real(rk) :: remineral_rate , hydrolysis_rate  ! Temp dependent remineralisation and hydrolysis rates
  real(rk) :: aggreg_rate ! particle aggregation 
  logical  :: out = .true.
!   if(36000.eq.secondsofday .and. mod(julianday,1).eq.0 .and. outn) out=.true.
#define UNIT *1.1574074074E-5_rk ! 1/86400

 !allocate(remin_chemical(self%num_chemicals), stat=rc) allocate(qualDetv(self%num_elements), stat=rc)
 !allocate(qualDOMv(self%num_elements), stat=rc) allocate(det_prod(self%num_elements), stat=rc, source=0.0_rk)
 !allocate(dom_prod(self%num_elements), stat=rc, source=0.0_rk) allocate(nut_prod(self%num_elements), stat=rc, source=0.0_rk)

 allocate(det%element(self%num_elements), stat=rc)
 !if rc /= 0 goto 99
 allocate(dom%element(self%num_elements), stat=rc)
 !if rc /= 0 goto 99
 
 ! set indices of elements vectors and pointers
do i = 1,self%num_elements ! 
  set_element(Index_Det,element(i), i)
  set_element(Index_DOM,element(i), i)
  set_pointer(det,element(i), i)
  set_pointer(dom,element(i), i)
end do
! det%C => det%element(1) det%N => det%element(2)
! dom%C => dom%element(1) dom%N => dom%element(2)

 _LOOP_BEGIN_

! First retrieve current (local) state  variable values
!---------- GET for each state variable ----------
do i = 1,self%num_chemicals ! e.g., CO2, NO3, NH4 (PO4)
!  if (_AVAILABLE_(self%id_dix(i))) then
   _GET_(self%id_dix(i), dix%chemical(i))  ! Dissolved Inorganic Nutrient DIX in mmol-X/m**3
!  end if
end do
! TODO: create array of OM classes
do i = 1,self%num_elements ! e.g., N  ( C, Si, Fe, P)
  _GET_(self%id_det%element(i), det%element(i))  ! Detritus Organics in mmol-C/m**3
  _GET_(self%id_dom%element(i), dom%element(i))  ! Dissolved Organics in mmol-C/m**3
end do
!do i = 1,num_phyclass ! e.g., DIA, FLA, CYA  (or 1,5,30,100 µm)
!  _GET_(self%id_phy(i), phy%class(i))  ! Detritus Organics in mmol-C/m**3
!end do
!---------- GET ambient physico-stoich conditions ----------
_GET_(self%id_temp, env%temp)  ! water temperature
! _GET_(self%id_par, env%par)    ! light photosynthetically active radiation

!_SET_DIAGNOSTIC_(self%id_PAR_diag,env%par)         !average Photosynthetically_Active_Radiation_

call calc_sensitivities(self,sens,env)
! call photosynthesis(self,sens,phy,nut,uptake,exud,acclim)

!___________________________________________________________________
!  ---  POM&DOM quality, relative to Refield ?
! TODO: merge POM and DOM !
!  Nqual = 1 full N:C dependency   0: only fresh material
qualDet   = (1.0d0-self%Nqual) + self%Nqual * det%N /(det%C + self%small) * self%CNref
qualDOM   = (1.0d0-self%Nqual) + self%Nqual * dom%N /(dom%C + self%small) * self%CNref
! distribute  preferential degradation rate to elements (N:fast; C:quality dep; others: intermediate)
if (associated(self%Index_Det%N))  qualDetv(self%Index_Det%N) = 1.0d0
if (associated(self%Index_DOM%N))   qualDOMv(self%Index_DOM%N)  = 1.0d0
if (associated(self%Index_Det%C))  qualDetv(self%Index_DetC) = qualDet
if (associated(self%Index_DOM%C))   qualDOMv(self%Index_DOC)  = qualDOM
if (associated(self%Index_Det_No_NorC)) qualDetv(self%Index_Det_No_NorC) = (1.0d0 + qualDet)/2
if (associated(self%Index_DOX_No_NorC)) qualDOMv(self%Index_DOX_No_NorC) = (1.0d0 + qualDOM)/2
!_____________________________________________________________________________
!     denitrification 
! pelagic N-loss by denitrification, emulating benthic pool and suboxic micro-environments
! check for the carbon index 
! calculate substrate (all OM)
!if (associated(self%dom%C)) sum_OM = sum_OM + dom%C* qualDOM ! 

if (associated(self%det%C)) then ! POC Glud LO 2015 (suboxic spots in particles) 
  ! calculate oxidant (NO3)
  if (associated(self%dix%no3)) then ! TODO: add nitrite NO2
    nitrate = dix%no3
    Index_NO3 = 1
  elseif (associated(self%dix%din)) then
    nitrate = dix%din    ! lower denitrication self%denit
    Index_NO3 = 1
  else
    nitrate = 0.1_rk ! TODO: replace by SMALL
  endif
  denitrate = self%denit * sens%f_T * det%C * qualDet * (1.0d0-exp(-nitrate/self%DenitKNO3)) 
elseif
  denitrate = 0.0_rk
endif

!  ---  hydrolysis & remineralisation rate (temp dependent)
hydrolysis  = self%hydrolysis * sens%f_T 
remineral   = self%remineral  * sens%f_T 

!________________________________________________________________________________
!
!  --- DETRITUS C
! _GET_(self%id_, phy%mort, zoo%mort)
! TODO: link to other modules
! GET det_prod(i)  =   !(floppZ%C + zoo_mort) * zoo%C + aggreg_rate * phy%C  
! GET dom_prod(i)  =   ! exud%C * phy%C 

do i = 1,num_elements ! e.g., N  ( C, Si, Fe, P)
  hydrolysis_rate = self%hydrolysis * qualDetv(i) * det%element(i)
  remineral_rate  = self%remineral  * qualDOMv(i) * dom%element(i)
  rhs_det%element(i) = det_prod(i) - hydrolysis_rate
  rhs_dom%element(i) = dom_prod(i) + hydrolysis_rate - remineral_rate
  
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
if(associated(self%Index_Det%C)) rhsv%det%element(self%Index_Det%C) = rhsv%det%element(self%Index_Det%C) + denitrate 

!Index_DetN Index_DON Index_DetC Index_DOC Index_Det_No_NorC Index_DOX_No_NorC
! TODO: link to other modules
! GET nut_prod(i)  =   !  -uptake%N * phy%C + lossZ%N * zoo%C 

do i = 1,self%num_chemicals 
  rhs_nut%chemical(i) = remin_chemical(i) + nut_prod(i) 
end do
!if(associated(self%Index_NO3)) rhsv%nut%chemical(self%Index_NO3) = rhsv%nut%chemical(self%Index_NO3) - 0.8_rk*denitrate 
!if(associated(self%Index_NO3)) rhsv%nut%chemical(self%Index_NO3) = rhsv%nut%chemical(self%Index_NO3) - 0.8_rk*denitrate 

!  chemostat mode 
if (_AVAILABLE_(self%dil)) then
  if (self%dil .gt. 0.0d0) then
    do i = 1,num_elements ! e.g., N  ( C, Si, Fe, P)
      rhsv%det%element(i) = rhsv%det%element(i) - self%dil * det%element(i)
      rhsv%dom%elem(i) = rhsv%dom%elem(i) - self%dil * dom%elem(i)
    end do
    do i = 1,num_nutrients 
      rhsv%nut%chemical(i) = rhsv%nut%chemical(i) - self%dil * (self%dix%chemical0(i) - self%dix%chemical(i))
    end do
  endif
endif

! tell FABM about right hand sides ....
do i = 1,num_elements 
  _SET_ODE_(self%id_det(i), rhsv%det%element(i) UNIT)
  _SET_ODE_(self%id_dom(i), rhsv%dom%elem(i) UNIT)
end do
do i = 1,num_nutrients 
  _SET_ODE_(self%id_dix(i), rhsv%nut%chemical(i) UNIT)
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

end subroutine tame_do

! set indices of elements vectors 

subroutine set_element(object, name, value)
  type(t), intent(inout) :: object
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
  type(t), intent(inout) :: object
  character(*), intent(in) :: name
  integer, intent(in) :: value
  select case (name)
  case ('C') ; object%C => object%element(value)
  case ('N') ; object%N => object%element(value) 
  case ('P') ; object%P => object%element(value) 
  case ('S') ; object%Si => object%element(value) 
  case ('F') ; object%Fe => object%element(value) 
  end select
end subroutine set_pointer

