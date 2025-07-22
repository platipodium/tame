! SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
!
! SPDX-License-Identifier: Apache-2.0

#include "fabm_driver.h"
!----------------------------------------
!
!	tame/tame
!
! here some comments
!----------------------------------------
!
!----------------------------------------
! !INTERFACE:
 MODULE tame_tame
!
! !USES:
 use fabm_types

 implicit none

 private
!
 ! !PUBLIC_DERIVED_TYPES:
 type,extends(type_base_model),public :: type_tame_tame
	type (type_state_variable_id) :: id_c
	type (type_dependency_id) :: id_par,id_temp
	type (type_horizontal_dependency_id) :: id_taub
	type (type_diagnostic_variable_id) :: id_chla,id_GPP,id_NPP
	real(rk) :: id_wdz,id_wz,id_sedrate,id_tau_crit,id_kc
	integer :: id_env_type

 contains
	procedure :: initialize
	procedure :: do
 end type
!EOP
!----------------------------------------
 CONTAINS

!----------------------------------------!BOP
!
! !IROUTINE: Initialise the tame/tame model
!
! !INTERFACE:
 subroutine initialize(self,configunit)
!
! !DESCRIPTION:
! Reading namelist and registration of variables with FABM
!
! !INPUT PARAMETERS: class(type_tame_tame),intent(inout),target :: self
 integer,		intent(in)		:: configunit

 real(rk),parameter :: secs_per_day = 86400._rk

 call self%register_state_variable(self%id_c, 'c','mmol N/m3','concentration', default=0.001 , vertical_movement=wz)
 call self%get_parameter(self%id_env_type, 'env_type','Define environment type, either fresh or marine','(Define environment type, either fresh or marine), default = marine, default = marine', default=marine )
 call self%get_parameter(self%id_wdz, 'wdz','positive: upwards/floating, negative: downwards/sm/d','vertical velocity of detritus (m/d), default = -4.5', default=-4.5 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%id_wz, 'wz','m/d','vertical velocity [positive: upwards/floating, negative: downwards/sinking], default = 0.0', default=-0.5 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%id_sedrate, 'sedrate','m/d','sedimentation rate, default = 0.0 2.25', default=2.25 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%id_tau_crit, 'tau_crit','N/m2','critical shear stress, default = 0.07', default=0.07 )
 call self%get_parameter(self%id_kc, 'kc','m2/mmol N','specific light attenuation of detritus', default=0.5 )
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

! !IROUTINE: Right hand sides of tame/tame model
!
! !INTERFACE:
 subroutine do(self,_ARGUMENTS_DO_)
!
! !LOCAL VARIABLES:
	real(rk) :: par, temp
  real(rk) :: c
! real(rk) ::
  real(rk),parameter :: secs_per_day = 86400._rk

! Enter spatial_loops (if any)
 _LOOP_BEGIN_

  _GET_(self%id_par, par)		! downwelling_photosynthetic_radiative_flux
  _GET_(self%id_temp, temp)		! temperature

!----------------------------------------


!----------------------------------------
! _SET_ODE_(self%id_c,  )
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
END  MODULE tame_tame

!----------------------------------------
