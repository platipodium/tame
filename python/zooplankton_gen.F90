#include "fabm_driver.h"
!----------------------------------------
!
!	selma/zooplankton
!
! here some comments
!----------------------------------------
!
!----------------------------------------
! !INTERFACE:
 MODULE selma_zooplankton
!
! !USES:
 use fabm_types

 implicit none

 private
!
 ! !PUBLIC_DERIVED_TYPES:
 type,extends(type_base_model),public :: type_selma_zooplankton
	type (type_state_variable_id) :: id_c
	type (type_state_variable_id) :: id_prey1,id_prey2,id_prey3,id_aa,id_po,id_dd,id_o2
	type (type_dependency_id) :: id_par,id_temp
	type (type_horizontal_dependency_id) :: id_taub
	type (type_diagnostic_variable_id) :: id_chla,id_GPP,id_NPP
	real(rk) :: id_c0,id_rfr,id_rfc,id_pref3,id_nue,id_sigma_b,id_iv,id_graz,id_toptz,id_zcl1
	integer :: id_nprey

 contains
	procedure :: initialize
	procedure :: do
 end type
!EOP
!----------------------------------------
 CONTAINS

!----------------------------------------!BOP
!
! !IROUTINE: Initialise the selma/zooplankton model
!
! !INTERFACE:
 subroutine initialize(self,configunit)
!
! !DESCRIPTION:
! Reading namelist and registration of variables with FABM
!
! !INPUT PARAMETERS: class(type_selma_zooplankton),intent(inout),target :: self
 integer,		intent(in)		:: configunit

 real(rk),parameter :: secs_per_day = 86400._rk

 call self%register_state_variable(self%id_c, 'c','mmol N/m3','concentration', default=0.001 , vertical_movement=wz)
 call self%register_state_dependency(self%id_prey1, 'prey1','mmol N/m3','prey 1' )
 call self%register_state_dependency(self%id_prey2, 'prey2','mmol N/m3','prey 2' )
 call self%register_state_dependency(self%id_prey3, 'prey3','mmol N/m3','prey 3' )
 call self%register_state_dependency(self%id_aa, 'aa','mmol N/m3','ammonium' )
 call self%register_state_dependency(self%id_po, 'po','mmol P/m3','phosphate' )
 call self%register_state_dependency(self%id_dd, 'dd','mmol N/m3','detritus' )
 call self%register_state_dependency(self%id_o2, 'o2','mmol O2/m3','oxygen' )
 call self%get_parameter(self%id_c0, 'c0','mmol N/m3','background concentration, default = 0.0', default=0.001 )
 call self%get_parameter(self%id_rfr, 'rfr','mol P/mol N','phosphorus : nitrogen ratio, default = 0.0625', default=0.0625 )
 call self%get_parameter(self%id_rfc, 'rfc','mol C/mol N','carbon : nitrogen ratio, default = 6.625', default=6.625 )
 call self%get_parameter(self%id_nprey, 'nprey','','number of prey, default = 1', default=3 )
 call self%get_parameter(self%id_pref3, 'pref3','-','preference for prey 3, default = 1.0', default=0.5 )
 call self%get_parameter(self%id_nue, 'nue','m3/d/mmol N','respiration rate, default = 0.01', default=0.01 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%id_sigma_b, 'sigma_b','m3/d/mmol N','mortality rate, default = 0.03', default=0.03 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%id_iv, 'iv','1/(mmol N/m3','Ivlev constant, quadratic3), default = 1.2', default=1.2 )
 call self%get_parameter(self%id_graz, 'graz','1/d','grazing rate, default = 0.5', default=0.5 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%id_toptz, 'toptz','deg C','optimal temperature for grazing, default = 20.0', default=20.0 )
 call self%get_parameter(self%id_zcl1, 'zcl1','-','closure parameter, default = 50.0', default=50.0 )
 call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
 call self%register_dependency(self%id_temp, standard_variables%temperature)
 call self%register_dependency(self%id_taub, standard_variables%bottom_stress)
 call self%register_diagnostic_variable(self%id_chla, 'chla', 'mg chl a/m3', 'chlorophyll concentration')
 call self%register_diagnostic_variable(self%id_GPP, 'GPP',  'mmol/m3/d',   'gross primary production')
 call self%register_diagnostic_variable(self%id_NPP, 'NPP',  'mmol/m3/d',   'net primary production')

 end subroutine initialize
!EOC

!----------------------------------------!BOP
!

! !IROUTINE: Right hand sides of selma/zooplankton model
!
! !INTERFACE:
 subroutine do(self,_ARGUMENTS_DO_)
!
! !LOCAL VARIABLES:
	real(rk) :: par, temp
  real(rk) :: c
  real(rk) :: prey1,prey2,prey3,aa,po,dd,o2
! real(rk) :: 
  real(rk),parameter :: secs_per_day = 86400._rk

! Enter spatial_loops (if any)
 _LOOP_BEGIN_

  _GET_(self%id_prey1, prey1)		! prey 1
  _GET_(self%id_prey2, prey2)		! prey 2
  _GET_(self%id_prey3, prey3)		! prey 3
  _GET_(self%id_aa, aa)		! ammonium
  _GET_(self%id_po, po)		! phosphate
  _GET_(self%id_dd, dd)		! detritus
  _GET_(self%id_par, par)		! downwelling_photosynthetic_radiative_flux
  _GET_(self%id_temp, temp)		! temperature

!----------------------------------------


!----------------------------------------
! _SET_ODE_(self%id_c,  )
! _SET_ODE_(self%id_,  )
! _SET_ODE_(self%id_prey1,  )
!  _SET_DIAGNOSTIC_(self%id_chla, )		! 'mg chl a/m3', 'chlorophyll concentration'
!  _SET_DIAGNOSTIC_(self%id_GPP, )		!  'mmol/m3/d',   'gross primary production'
!  _SET_DIAGNOSTIC_(self%id_NPP, )		!  'mmol/m3/d',   'net primary production'

! Leave spatial loops (if any)
 _LOOP_END_

 END subroutine do
!EOC

!----------------------------------------
!----------------------------------------
END  MODULE selma_zooplankton

!----------------------------------------
