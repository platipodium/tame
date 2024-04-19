#include "fabm_driver.h"
!----------------------------------------
!
!	selma/selma
!
! here some comments
!----------------------------------------
!
!----------------------------------------
! !INTERFACE:
 MODULE selma_selma
!
! !USES:
 use fabm_types

 implicit none

 private
!
 ! !PUBLIC_DERIVED_TYPES:
 type,extends(type_base_model),public :: type_selma_selma
	type (type_state_variable_id) :: id_nn,id_aa,id_dd,id_o2,id_po
	type (type_dependency_id) :: id_par,id_temp
	type (type_horizontal_dependency_id) :: id_taub
	type (type_diagnostic_variable_id) :: id_chla,id_GPP,id_NPP
	real(rk) :: id_wdz,id_wpo4,id_dn,id_dn_sed,id_kc,id_q10_rec,id_ade_r0,id_alphaade,id_q10_recs,id_tau_crit,id_sedrate,id_erorate,id_sedratepo4,id_eroratepo4,id_po4ret,id_pburialrate,id_fl_burialrate,id_pliberationrate,id_ipo4th,id_maxsed,id_br0,id_fds,id_pvel
	integer :: id_env_type,id_newflux

 contains
	procedure :: initialize
	procedure :: do
 end type
!EOP
!----------------------------------------
 CONTAINS

!----------------------------------------!BOP
!
! !IROUTINE: Initialise the selma/selma model
!
! !INTERFACE:
 subroutine initialize(self,configunit)
!
! !DESCRIPTION:
! Reading namelist and registration of variables with FABM
!
! !INPUT PARAMETERS: class(type_selma_selma),intent(inout),target :: self
 integer,		intent(in)		:: configunit

 real(rk),parameter :: secs_per_day = 86400._rk

 call self%register_state_variable(self%id_nn, 'nn','mmol N/m3','nitrate', default=20.0 )
 call self%register_state_variable(self%id_aa, 'aa','mmol N/m3','ammonium', default=1.0 )
 call self%register_state_variable(self%id_dd, 'dd','mmol N/m3','detritus', default=0.0 )
 call self%register_state_variable(self%id_o2, 'o2','mmol O2/m3','oxygen', default=280.0 )
 call self%register_state_variable(self%id_po, 'po','mmol P/m3','phosphate', default=4.0 )
 call self%get_parameter(self%id_env_type, 'env_type','Define environment type, either fresh or marine','(Define environment type, either fresh or marine), default = marine, default = marine', default=marine )
 call self%get_parameter(self%id_wdz, 'wdz','positive: upwards/floating, negative: downwards/spositive: upwards/floating, negative: downwards/s1/d','vertical velocity of detritus (m/d), default = -4.5', default=-4.5 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%id_wpo4, 'wpo4','positive: upwards/floating, negative: downwards/s1/d','vertical velocity of suspended P-Fe (m/d), default = -1.0', default=-1.0 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%id_dn, 'dn','1/d','detritus mineralization rate, default = 0.003', default=0.003 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%id_dn_sed, 'dn_sed','1/d','sediment mineralization rate, default = 0.002', default=0.002 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%id_kc, 'kc','m2/mmol N','specific light attenuation of detritus', default=0.5 )
 call self%get_parameter(self%id_q10_rec, 'q10_rec','1/K','temperature dependence of detritus remineralization, default = 0.15', default=0.15 )
 call self%get_parameter(self%id_ade_r0, 'ade_r0','1/d','maximum chemoautolithotrophic denitrification rate, default = 0.1', default=0.1 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%id_alphaade, 'alphaade','mmol N/m3','half-saturation constant for chemoautolithotrophic denitrification, default = 1.0', default=1.0 )
 call self%get_parameter(self%id_q10_recs, 'q10_recs','1/K','temperature dependence of sediment remineralization, default = 0.175', default=0.175 )
 call self%get_parameter(self%id_tau_crit, 'tau_crit','N/m2','critical shear stress, default = 0.07', default=0.07 )
 call self%get_parameter(self%id_sedrate, 'sedrate','m/d','detritus sedimentation rate, default = 2.25', default=2.25 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%id_erorate, 'erorate','1/d','sediment erosion rate, default = 6.0', default=6.0 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%id_sedratepo4, 'sedratepo4','m/d','P-Fe sedimentation rate, default = 0.5', default=0.5 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%id_eroratepo4, 'eroratepo4','1/d','P-Fe erosion rate, default = 6.0', default=6.0 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%id_po4ret, 'po4ret','-','phosphate retention rate, oxic sediments, default = 0.18', default=0.18 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%id_pburialrate, 'pburialrate','1/d','phosphate burial rate, default = 0.007', default=0.007 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%id_fl_burialrate, 'fl_burialrate','1/d','sediment burial rate, default = 0.001', default=0.001 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%id_pliberationrate, 'pliberationrate','1/d','phosphate liberation rate, anoxic sediments, default = 0.1', default=0.1 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%id_ipo4th, 'ipo4th','mmol P/m2','maximum phosphorus density available for burial, default = 100.0', default=100.0 )
 call self%get_parameter(self%id_maxsed, 'maxsed','mmol N/m2','maximum active sediment density, default = 1000.0', default=1000.0 )
 call self%get_parameter(self%id_br0, 'br0','1/d','bioresuspension rate, default = 0.03', default=0.03 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%id_fds, 'fds','-','fraction of sediment remineralization fueled by denitrification, default = 0.7', default=0.7 )
 call self%get_parameter(self%id_pvel, 'pvel','m/d','piston velocity, default = 5.0', default=5.0 , scale_factor=1.0_rk/secs_per_day)
 call self%get_parameter(self%id_newflux, 'newflux','-','oxygen flux type, default = 2', default=2 )
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

! !IROUTINE: Right hand sides of selma/selma model
!
! !INTERFACE:
 subroutine do(self,_ARGUMENTS_DO_)
!
! !LOCAL VARIABLES:
	real(rk) :: par, temp
  real(rk) :: nn,aa,dd,o2,po
! real(rk) :: 
  real(rk),parameter :: secs_per_day = 86400._rk

! Enter spatial_loops (if any)
 _LOOP_BEGIN_

  _GET_(self%id_nn, nn)		! nitrate
  _GET_(self%id_aa, aa)		! ammonium
  _GET_(self%id_dd, dd)		! detritus
  _GET_(self%id_o2, o2)		! oxygen
  _GET_(self%id_par, par)		! downwelling_photosynthetic_radiative_flux
  _GET_(self%id_temp, temp)		! temperature

!----------------------------------------


!----------------------------------------
! _SET_ODE_(self%id_nn,  )
! _SET_ODE_(self%id_aa,  )
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
END  MODULE selma_selma

!----------------------------------------
