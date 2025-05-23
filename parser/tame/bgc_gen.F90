! SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
!
! SPDX-License-Identifier: Apache-2.0

#include "fabm_driver.h"
!----------------------------------------
!
!	tame/bgc
!
! here some comments
!----------------------------------------
!
!----------------------------------------
! !INTERFACE:
 MODULE tame_bgc
!
! !USES:
 use fabm_types

 implicit none

 private
!
 ! !PUBLIC_DERIVED_TYPES:
 type,extends(type_base_model),public :: type_tame_bgc
	type (type_state_variable_id) :: id_no3,id_nh4,id_det%C,id_det%N,id_dom%C,id_dom%N,id_o2,id_po4
	type (type_state_variable_id) :: id_phy,id_zoo
	type (type_dependency_id) :: id_par,id_temp
	type (type_horizontal_dependency_id) :: id_taub
	type (type_diagnostic_variable_id) :: id_chla,id_GPP,id_NPP
	real(rk) :: id_remineral,id_hydrolysis,id_alloc_N,id_Nqual,id_CNref,id_DenitKno3,id_denit,id_T_ref,id_rq10
	integer :: id_tlim

 contains
	procedure :: initialize
	procedure :: do
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

 real(rk),parameter :: secs_per_day = 86400._rk

 call self%register_state_variable(self%id_no3, 'no3','mmol N/m3','nitrate', default=30.0 )
 call self%register_state_variable(self%id_nh4, 'nh4','mmol N/m3','ammonium', default=2.0 )
 call self%register_state_variable(self%id_det%C, 'det%C','mmol C/m3','detritus carbon', default=0.0 )
 call self%register_state_variable(self%id_det%N, 'det%N','mmol N/m3','detritus nitrogen', default=0.0 )
 call self%register_state_variable(self%id_dom%C, 'dom%C','mmol C/m3','Dissolved Organic Carbon', default=0.0 )
 call self%register_state_variable(self%id_dom%N, 'dom%N','mmol N/m3','Dissolved Organic Nitrogen', default=0.0 )
 call self%register_state_variable(self%id_o2, 'o2','mmol O2/m3','oxygen', default=280.0 )
 call self%register_state_variable(self%id_po4, 'po4','mmol P/m3','phosphate', default=2.0 )
 call self%register_state_dependency(self%id_phy, 'phy','','' )
 call self%register_state_dependency(self%id_zoo, 'zoo','','' )
 call self%get_parameter(self%id_remineral, 'remineral','1/d','DOM remineralisation rate', default=0.1 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%id_hydrolysis, 'hydrolysis','1/d','detritus hydrolysis rate, default = 0.003', default=0.05 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%id_alloc_N, 'alloc_N','','nh4 - no3 product ratio remineralisation', default=0.5 )
 call self%get_parameter(self%id_Nqual, 'Nqual','','OM fraction w quality prop to N:Cratio ', default=1. )
 call self%get_parameter(self%id_CNref, 'CNref','Redfield','POM quality relative to carbon : nitrogen ratio (mol C/mol N), default = 6.625', default=6.625 )
 call self%get_parameter(self%id_DenitKno3, 'DenitKno3','mmol N/m3','half-saturation no3 denitrification, default = 0.25', default=1. )
 call self%get_parameter(self%id_denit, 'denit','1/d','pelagic denitrification rate', default=0.01 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%id_T_ref, 'T_ref','Kelvin','reference temperature', default=293. )
 call self%get_parameter(self%id_rq10, 'rq10','','temperature dependence Q10', default=0.175 )
 call self%get_parameter(self%id_tlim, 'tlim','0: none, 1: flagellate-style, 2: cyanobacteria-style','temperature limitation of growth, default = 0', default=0 )
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

! !IROUTINE: Right hand sides of tame/bgc model
!
! !INTERFACE:
 subroutine do(self,_ARGUMENTS_DO_)
!
! !LOCAL VARIABLES:
	real(rk) :: par, temp
  real(rk) :: no3,nh4,det%C,det%N,dom%C,dom%N,o2,po4
  real(rk) :: phy,zoo
! real(rk) ::
  real(rk),parameter :: secs_per_day = 86400._rk

! Enter spatial_loops (if any)
 _LOOP_BEGIN_

  _GET_(self%id_no3, no3)		! nitrate
  _GET_(self%id_nh4, nh4)		! ammonium
  _GET_(self%id_det%C, det%C)		! detritus carbon
  _GET_(self%id_det%N, det%N)		! detritus nitrogen
  _GET_(self%id_dom%C, dom%C)		! Dissolved Organic Carbon
  _GET_(self%id_dom%N, dom%N)		! Dissolved Organic Nitrogen
  _GET_(self%id_o2, o2)		! oxygen
  _GET_(self%id_phy, phy)		!
  _GET_(self%id_par, par)		! downwelling_photosynthetic_radiative_flux
  _GET_(self%id_temp, temp)		! temperature

!----------------------------------------


!----------------------------------------
! _SET_ODE_(self%id_no3,  )
! _SET_ODE_(self%id_nh4,  )
! _SET_ODE_(self%id_phy,  )
!  _SET_DIAGNOSTIC_(self%id_chla, )		! 'mg chl a/m3', 'chlorophyll concentration'
!  _SET_DIAGNOSTIC_(self%id_GPP, )		!  'mmol/m3/d',   'gross primary production'
!  _SET_DIAGNOSTIC_(self%id_NPP, )		!  'mmol/m3/d',   'net primary production'

! Leave spatial loops (if any)
 _LOOP_END_

 END subroutine do
!EOC

!----------------------------------------
!----------------------------------------
END  MODULE tame_bgc

!----------------------------------------
