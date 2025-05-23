! SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
!
! SPDX-License-Identifier: Apache-2.0

#include "fabm_driver.h"

module tame_npzdep

   use fabm_types
   use fabm_particle

   implicit none

   private

   !type, extends(type_base_model), public :: type_tame_npzdep
   type, extends(type_particle_model), public :: type_tame_npzdep
      ! Variable identifiers
      type (type_state_variable_id)      :: id_n, id_p, id_z, id_d
      type (type_state_variable_id)      :: id_dic
      type (type_dependency_id)          :: id_par
      type (type_surface_dependency_id)  :: id_I_0
      type (type_diagnostic_variable_id) :: id_PPR, id_NPR, id_dPAR
      type (type_horizontal_dependency_id)          :: id_bedstress,id_wdepth     !new
      type (type_horizontal_diagnostic_variable_id) :: id_w_bot    !new
      ! type (type_horizontal_diagnostic_variable_id),allocatable,dimension(:) :: id_ndep    !new
      type (type_bottom_state_variable_id),allocatable :: id_ndep(:)

      ! Model parameters
      real(rk) :: p0, z0, kc, i_min, rmax, gmax, iv, alpha, rpn, rzn, rdn, rpdu, rpdl, rzd
      real(rk) :: dic_per_n
   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_ppdd
   end type

contains

   subroutine initialize(self, configunit)
      class (type_tame_npzdep), intent(inout), target :: self
      integer,                intent(in)            :: configunit

      real(rk), parameter :: d_per_s = 1.0_rk / 86400.0_rk
      real(rk) :: w_p, w_d, vel_crit    !new

      ! Store parameter values in our own derived type
      ! NB: all rates must be provided in values per day in te configuration file,
      ! and are converted here to values per second by specifying scale_factor=d_per_s.
      call self%get_parameter(self%p0,'p0','mmol m-3','background phytoplankton concentration ',default=0.0225_rk)
      call self%get_parameter(self%z0,'z0','mmol m-3','background zooplankton concentration',default=0.0225_rk)
      call self%get_parameter(self%kc,'kc','m2 mmol-1','specific light extinction of phytoplankton and detritus',default=0.03_rk)
      call self%get_parameter(self%i_min,'i_min','W m-2','minimum light intensity in euphotic zone',default=25.0_rk)
      call self%get_parameter(self%rmax,'rmax','d-1','maximum specific growth rate of phytoplankton',default=1.0_rk,scale_factor=d_per_s)
      call self%get_parameter(self%gmax,'gmax','d-1','maximum specific grazing rate of zooplankton',default=0.5_rk,scale_factor=d_per_s)
      call self%get_parameter(self%iv,'iv','m3 mmol-1','Ivlev grazing constant',default=1.1_rk)
      call self%get_parameter(self%alpha,'alpha','mmol m-3','half-saturation nutrient concentration for phytoplankton',default=0.3_rk)
      call self%get_parameter(self%rpn,'rpn','d-1','loss rate of phytoplankton to nutrients',default=0.01_rk,scale_factor=d_per_s)
      call self%get_parameter(self%rzn,'rzn','d-1','loss rate of zooplankton to nutrients',default=0.01_rk,scale_factor=d_per_s)
      call self%get_parameter(self%rdn,'rdn','d-1','detritus remineralization rate',default=0.003_rk,scale_factor=d_per_s)
      call self%get_parameter(self%rpdu,'rpdu','d-1','phytoplankton mortality in euphotic zone',default=0.02_rk,scale_factor=d_per_s)
      call self%get_parameter(self%rpdl,'rpdl','d-1','phytoplankton mortality below euphotic zone',default=0.1_rk,scale_factor=d_per_s)
      call self%get_parameter(self%rzd,'rzd','d-1','zooplankton mortality',default=0.02_rk,scale_factor=d_per_s)
      call self%get_parameter(self%dic_per_n,'dic_per_n','-','molar C:N ratio of biomass',default=106.0_rk/16.0_rk)
      call self%get_parameter(w_p,'w_p','m d-1','vertical velocity of phytoplankton (<0 for sinking)',default=-1.0_rk, scale_factor=d_per_s)
      call self%get_parameter(w_d,'w_d','m d-1','vertical velocity of detritus  (<0 for sinking)',default=-5.0_rk,scale_factor=d_per_s)
      call self%get_parameter(vel_crit,'vel_crit','m/s','critical bed shear velocity for deposition',default=0.01_rk)    !new


      ! Register state variables
      call self%register_state_variable(self%id_n, 'nut', 'mmol m-3', 'nutrients',    4.5_rk, minimum=0.0_rk, no_river_dilution=.true.)
      call self%register_state_variable(self%id_p, 'phy', 'mmol m-3', 'phytoplankton',0.0_rk, minimum=0.0_rk, vertical_movement=w_p)
      call self%register_state_variable(self%id_z, 'zoo', 'mmol m-3', 'zooplankton',  0.0_rk, minimum=0.0_rk)
      call self%register_state_variable(self%id_d, 'det', 'mmol m-3', 'detritus',     4.5_rk, minimum=0.0_rk, vertical_movement=w_d)

      ! Register the contribution of all state variables to total nitrogen
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_n)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_p)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_z)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_d)


      ! Link to pools in which to deposit matter.    !new
      allocate(self%id_ndep(1))
      call self%register_state_dependency(self%id_ndep(1),'deposition_target','N/m^2','target pool for det deposition')

      !call self%request_coupling_to_model(self%id_ndep(1),'deposition_target',self%id_d)

      ! Register optional link to external DIC pool
      call self%register_state_dependency(self%id_dic, 'dic', 'mmol m-3', 'total dissolved inorganic carbon', required=.false.)

      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_PPR,  'PPR', 'mmol m-3 d-1', 'gross primary production rate')
      call self%register_diagnostic_variable(self%id_NPR,  'NPR', 'mmol m-3 d-1', 'net community production rate')
      call self%register_diagnostic_variable(self%id_dPAR, 'PAR', 'W m-2',        'photosynthetically active radiation')

      ! Register environmental dependencies
      call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
      call self%register_dependency(self%id_I_0, standard_variables%surface_downwelling_photosynthetic_radiative_flux)

      ! Let phytoplankton (including background concentration) and detritus contribute to light attentuation
      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, self%id_p, scale_factor=self%kc)
      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, self%id_d, scale_factor=self%kc)
      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, self%p0 * self%kc)
   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_tame_npzdep), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk)            :: n, p, z, d, par, I_0
      real(rk)            :: iopt, rpd, primprod, g, dn
      real(rk), parameter :: secs_pr_day = 86400.0_rk

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         ! Retrieve current (local) state variable values.
         _GET_(self%id_n,n) ! nutrient
         _GET_(self%id_p,p) ! phytoplankton
         _GET_(self%id_z,z) ! zooplankton
         _GET_(self%id_d,d) ! detritus

         ! Retrieve current environmental conditions.
         _GET_(self%id_par,par)          ! local photosynthetically active radiation
         _GET_SURFACE_(self%id_I_0,I_0)  ! surface photosynthetically active radiation

         ! Light acclimation formulation based on surface light intensity.
         iopt = max(0.25*I_0,self%I_min)

         ! Loss rate of phytoplankton to detritus depends on local light intensity.
         if (par >= self%I_min) then
            rpd = self%rpdu
         else
            rpd = self%rpdl
         end if

         ! Define some intermediate quantities that will be reused multiple times.
         primprod = fnp(self%rmax, self%alpha, n, p + self%p0, par, iopt)
         g = fpz(self%gmax, self%iv, p, z + self%z0)
         dn = - primprod + self%rpn*p + self%rzn*z + self%rdn*d

         ! Set temporal derivatives
         _ADD_SOURCE_(self%id_n,dn)
         _ADD_SOURCE_(self%id_p,primprod - g - self%rpn*p - rpd*p)
         _ADD_SOURCE_(self%id_z,g - self%rzn*z - self%rzd*z)
         _ADD_SOURCE_(self%id_d,rpd*p + self%rzd*z - self%rdn*d)

         ! If an externally maintained DIC pool is present, change the DIC pool according to the
         ! the change in nutrients (assuming constant C:N ratio)
         if (_AVAILABLE_(self%id_dic)) _ADD_SOURCE_(self%id_dic,self%dic_per_n*dn)

         ! Export diagnostic variables
         _SET_DIAGNOSTIC_(self%id_dPAR,par)
         _SET_DIAGNOSTIC_(self%id_PPR,primprod*secs_pr_day)
         _SET_DIAGNOSTIC_(self%id_NPR,(primprod - self%rpn*p)*secs_pr_day)

      ! Leave spatial loops (if any)
      _LOOP_END_
   end subroutine do

   subroutine do_ppdd(self, _ARGUMENTS_DO_PPDD_)
      class (type_tame_npzdep), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_PPDD_

      real(rk)            :: n, p, z, d, par, I_0
      real(rk)            :: iopt, rpd, dn, primprod
      real(rk), parameter :: secs_pr_day = 86400.0_rk

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         ! Retrieve current (local) state variable values.
         _GET_(self%id_n,n) ! nutrient
         _GET_(self%id_p,p) ! phytoplankton
         _GET_(self%id_z,z) ! zooplankton
         _GET_(self%id_d,d) ! detritus

         ! Retrieve current environmental conditions.
         _GET_(self%id_par,par)          ! local photosynthetically active radiation
         _GET_SURFACE_(self%id_I_0,I_0)  ! surface photosynthetically active radiation

         ! Light acclimation formulation based on surface light intensity.
         iopt = max(0.25*I_0,self%I_min)

         ! Loss rate of phytoplankton to detritus depends on local light intensity.
         if (par >= self%I_min) then
            rpd = self%rpdu
         else
            rpd = self%rpdl
         end if

         ! Rate of primary production will be reused multiple times - calculate it once.
         primprod = fnp(self%rmax, self%alpha, n, p + self%p0, par, iopt)

         ! Assign destruction rates to different elements of the destruction matrix.
         ! By assigning with _SET_DD_SYM_(i,j,val) as opposed to _SET_DD_(i,j,val),
         ! assignments to dd(i,j) are automatically assigned to pp(j,i) as well.
         _SET_DD_SYM_(self%id_n,self%id_p,primprod)                           ! snp
         _SET_DD_SYM_(self%id_p,self%id_z,fpz(self%gmax,self%iv,p,z+self%z0)) ! spz
         _SET_DD_SYM_(self%id_p,self%id_n,self%rpn*p)                         ! spn
         _SET_DD_SYM_(self%id_z,self%id_n,self%rzn*z)                         ! szn
         _SET_DD_SYM_(self%id_d,self%id_n,self%rdn*d)                         ! sdn
         _SET_DD_SYM_(self%id_p,self%id_d,rpd*p)                              ! spd
         _SET_DD_SYM_(self%id_z,self%id_d,self%rzd*z)                         ! szd

         ! If an externally maintained DIC pool is present, change the DIC pool according to the
         ! the change in nutrients (assuming constant C:N ratio)
         dn = - primprod + self%rpn*p + self%rzn*z + self%rdn*d
         if (_AVAILABLE_(self%id_dic)) _SET_PP_(self%id_dic,self%id_dic,self%dic_per_n*dn)

         ! Export diagnostic variables
         _SET_DIAGNOSTIC_(self%id_dPAR,par)
         _SET_DIAGNOSTIC_(self%id_PPR,primprod*secs_pr_day)
         _SET_DIAGNOSTIC_(self%id_NPR,(primprod-self%rpn*p)*secs_pr_day)

      ! Leave spatial loops (if any)
      _LOOP_END_
   end subroutine do_ppdd

   ! Phytoplankton growth limited by light and nutrient availability
   elemental real(rk) function fnp(rmax, alpha, n, p, par, iopt)
      real(rk), intent(in) :: rmax, alpha, n, p, par, iopt

      fnp = rmax * par / iopt * exp(1.0_rk - par / iopt) * n / (alpha + n) * p
   end function fnp

   ! Ivlev formulation for zooplankton grazing on phytoplankton
   elemental real(rk) function fpz(gmax, iv, p, z)
      real(rk), intent(in) :: gmax, iv, p, z

      fpz = gmax * (1.0_rk - exp(-iv * iv * p * p)) * z
   end function fpz

end module tame_npzdep
