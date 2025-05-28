# SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
#
# SPDX-License-Identifier: Apache-2.0

#include "fabm_driver.h"
!----------------------------------------
!
!	selma/phytoplankton
!
! here some comments
!----------------------------------------
!
!----------------------------------------
! !INTERFACE:
 MODULE selma_phytoplankton
!
! !USES:
 use fabm_types

 implicit none

 private
!
 ! !PUBLIC_DERIVED_TYPES:
 type,extends(type_base_model),public :: type_selma_phytoplankton
	type (type_state_variable_id) :: id_c
	type (type_state_variable_id) :: id_aa,id_nn,id_o2,id_po,id_dd
	type (type_bottom_state_variable_id) :: id_fl
	type (type_dependency_id) :: id_par,id_temp
	type (type_horizontal_dependency_id) :: id_taub
	type (type_diagnostic_variable_id) :: id_chla,id_GPP,id_NPP
	real(rk) :: id_c0,id_rfr,id_rfc,id_imin,id_alpha,id_r0,id_nb,id_deltao,id_Yc,id_wz,id_kc,id_sedrate,id_tau_crit,id_tll
	integer :: id_tlim

 contains
	procedure :: initialize
	procedure :: do
	procedure :: do_bottom
 end type
!EOP
!----------------------------------------
 CONTAINS

!----------------------------------------!BOP
!
! !IROUTINE: Initialise the selma/phytoplankton model
!
! !INTERFACE:
 subroutine initialize(self,configunit)
!
! !DESCRIPTION:
! Reading namelist and registration of variables with FABM
!
! !INPUT PARAMETERS: class(type_selma_phytoplankton),intent(inout),target :: self
 integer,		intent(in)		:: configunit

 real(rk),parameter :: secs_per_day = 86400._rk

 call self%register_state_variable(self%id_c, 'c','mmol N/m3','concentration', default=0.001 , vertical_movement=wz)
 call self%register_state_dependency(self%id_aa, 'aa','mmol N/m3','ammonium' )
 call self%register_state_dependency(self%id_nn, 'nn','mmol N/m3','nitrate' )
 call self%register_state_dependency(self%id_o2, 'o2','mmol O2/m3','oxygen' )
 call self%register_state_dependency(self%id_po, 'po','mmol P/m3','phosphate' )
 call self%register_state_dependency(self%id_dd, 'dd','mmol N/m3','detritus' )
 call self%register_state_dependency(self%id_fl, 'fl','mmol N/m2','fluff' )
 call self%get_parameter(self%id_c0, 'c0','mmol N/m3','background concentration, default = 0.0', default=0.001 )
 call self%get_parameter(self%id_rfr, 'rfr','mol P/mol N','phosphorus : nitrogen ratio, default = 0.0625', default=0.0625 )
 call self%get_parameter(self%id_rfc, 'rfc','mol C/mol N','carbon : nitrogen ratio, default = 6.625', default=6.625 )
 call self%get_parameter(self%id_imin, 'imin','W/m2','minimal optimal light radiation, default = 50.0', default=35.0 )
 call self%get_parameter(self%id_alpha, 'alpha','mmol N/m3','half-saturation for nutrient uptake, default = 0.25', default=0.25 )
 call self%get_parameter(self%id_r0, 'r0','1/d','maximum growth rate, default = 1.3', default=1.3 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%id_nb, 'nb','1/d','excretion rate, default = 0.01', default=0.01 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%id_deltao, 'deltao','1/d','mortality rate, default = 0.02', default=0.02 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%id_Yc, 'Yc','mmol C/mg Chl a','carbon : chlorophyll a ratio, default = 6.25', default=6.25 )
 call self%get_parameter(self%id_wz, 'wz','m/d','vertical velocity [positive: upwards/floating, negative: downwards/sinking], default = 0.0', default=-0.5 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%id_kc, 'kc','m2/mmol N','specific light attenuation', default=0.5 )
 call self%get_parameter(self%id_sedrate, 'sedrate','m/d','sedimentation rate, default = 0.0 2.25', default=2.25 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%id_tau_crit, 'tau_crit','N/m2','critical shear stress, default = 0.07', default=0.07 )
 call self%get_parameter(self%id_tlim, 'tlim','0: none, 1: flagellate-style, 2: cyanobacteria-stdegrees C^2','temperature limitation of growth, default = 0', default=0 )
 call self%get_parameter(self%id_tll, 'tll','degrees C^2','half-saturation temperature, squared, default = 100.0', default=100.0 )
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

! !IROUTINE: Right hand sides of selma/phytoplankton model
!
! !INTERFACE:
 subroutine do(self,_ARGUMENTS_DO_)
!
! !LOCAL VARIABLES:
	real(rk) :: par, temp
  real(rk) :: c
  real(rk) :: aa,nn,o2,po,dd,fl
! real(rk) ::
  real(rk),parameter :: secs_per_day = 86400._rk

! Enter spatial_loops (if any)
 _LOOP_BEGIN_

  _GET_(self%id_aa, aa)		! ammonium
  _GET_(self%id_nn, nn)		! nitrate
  _GET_(self%id_o2, o2)		! oxygen
  _GET_(self%id_po, po)		! phosphate
  _GET_(self%id_dd, dd)		! detritus
  _GET_(self%id_par, par)		! downwelling_photosynthetic_radiative_flux
  _GET_(self%id_temp, temp)		! temperature

!----------------------------------------


!----------------------------------------
! _SET_ODE_(self%id_c,  )
! _SET_ODE_(self%id_,  )
! _SET_ODE_(self%id_aa,  )
!  _SET_DIAGNOSTIC_(self%id_chla, )		! 'mg chl a/m3', 'chlorophyll concentration'
!  _SET_DIAGNOSTIC_(self%id_GPP, )		!  'mmol/m3/d',   'gross primary production'
!  _SET_DIAGNOSTIC_(self%id_NPP, )		!  'mmol/m3/d',   'net primary production'

! Leave spatial loops (if any)
 _LOOP_END_

 END subroutine do
!EOC

!----------------------------------------
!BOP
!

! !IROUTINE: Right hand sides of benthic selma/phytoplankton model
!
! !INTERFACE:
 subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
 class(type_selma_phytoplankton),intent(in) :: self
 _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !LOCAL VARIABLES:
 real(rk) :: r0
 real(rk) :: taub

!EOP
!----------------------------------------
!BOC
 CONTAINS

!----------------------------------------!BOP
!
! Enter spatial loops over the horizontal domain (if any).
 _HORIZONTAL_LOOP_BEGIN_
 _GET_HORIZONTAL_(self%id_taub, taub)		! bottom_stress
! _SET_BOTTOM_ODE_(self%id_fl,  )

! Leave spatial loops over the horizontal domain (if any).
_HORIZONTAL_LOOP_END_

end subroutine do_bottom
!EOC
!----------------------------------------

!----------------------------------------
END  MODULE selma_phytoplankton

!----------------------------------------
