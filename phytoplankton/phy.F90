#include "fabm_driver.h"

!module examples_npzd_phy ! TAME phytoplankton module
module tame_phytoplankton

use fabm_types
use tame_types
use tame_functions

   implicit none

   private

   !type, extends(type_base_model), public :: type_examples_npzd_phy
   type, extends(type_base_model), public :: type_tame_phytoplankton

      ! Variable identifiers
      type (type_state_variable_id)      :: id_phytoplankton     ! Phytoplankton biomass
      type (type_state_variable_id)      :: id_n     ! N


      type (type_dependency_id)          :: id_par   ! PAR light

      !type (type_dependency_id)          :: id_grazing

      !type (type_dependency_id)          :: id_n     ! Nutrient
      !type (type_dependency_id)          :: id_d     ! Day of the year

      !type (type_surface_dependency_id)  :: id_I_0   ! Surface irradiance
      !type (type_dependency_id)          :: id_z     ! Zooplankton

      ! Model parameters
      real(rk) :: p0
      real(rk) :: rmax                ! Growth metabolism parameters
      real(rk) :: gamma               ! Light exploitation
      real(rk) :: s0                  ! Sinking
      real(rk) :: Vpotn           ! Uptake parameters
      real(rk) :: Kn           ! Uptake parameters
      real(rk) :: resp                ! Respiration parameters
      !real(rk) :: amax, pred_d        ! Grazing parameters (attack rate, escape with size and half saturation)
      !real(rk) :: z ! Zooplankton
      real(rk) :: n0 ! Background N
      real(rk) :: nremin ! N remineralization

   contains
      procedure :: initialize
      procedure :: do
      !procedure :: do_ppdd
   end type

contains

   subroutine initialize(self, configunit)
      class (type_examples_npzd_phy), intent(inout), target :: self
      integer,                        intent(in)            :: configunit

      real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk

      ! Store parameter values in our own derived type
      ! NB: all rates must be provided in values per day and are converted here to values per second.
      call self%get_parameter(self%rmax,    'rmax',    'd-1',          'maximum production rate',            default=8.0_rk, scale_factor=d_per_s)
      call self%get_parameter(self%gamma,   'gamma',   'microE-1 m-2', 'light absorption rate scaling',      default=1.0_rk)
      call self%get_parameter(self%p0,    'p0',    'mmol m-3',  'background concentration ',                default=0.0225_rk)

      call self%get_parameter(self%s0,      's0',      'd-1',          'default sinking rate',               default=0.1_rk, scale_factor=d_per_s)

      call self%get_parameter(self%Vpotn,   'Vpotn',   'mmol.d-1',      'N potential uptake rate',            default=8.0_rk, scale_factor=d_per_s)
      call self%get_parameter(self%Kn,   'Kn',   'mmol m-3',      'N demi saturation',            default=8.0_rk)

      call self%get_parameter(self%n0,       'n0',       'mmol m-3',     'Background Nitrogen',                default=1.0_rk)
      call self%get_parameter(self%nremin,       'nremin',       'd-1',     'Nitrogen remineralization rate',     default=0.003_rk, scale_factor=d_per_s)

      call self%get_parameter(self%resp,     'resp',     'mmol',         'carbon cost per nitrogen uptake',    default=0.2_rk)

      !call self%get_parameter(self%amax,    'amax',    'd-1',          'maximum attack rate',                default=0.2_rk, scale_factor=d_per_s)
      !call self%get_parameter(self%pred_d,  'pred_d',  'na',           'predator half saturation',           default=3.0_rk)
      !call self%get_parameter(self%z,       'z',       'mmol m-3',     'predators',                          default=0.5_rk)

      ! Register state variables
      call self%register_state_variable(self%id_phytoplankton,    'phytoplankton', 'mmol m-3', 'concentration', 1.0_rk, minimum=0.0_rk)
      call self%register_state_variable(self%id_n,     'n', 'mmol m-3', 'concentration', 5.0_rk, minimum=0.0_rk)

      ! Register environmental dependencies
      !call self%register_dependency(self%id_grazing, "grazing", 'd-1', 'grazing pressure', required = .false.)!, scale_factor = d_per_s) ! Zooplankton activity
      call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
      !call self%register_dependency(self%id_n,   standard_variables%total_nitrogen)

   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_examples_npzd_phy), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk)            :: phytoplankton, par, d, n
      real(rk)            :: production, respiration, sinking!, grazing
      real(rk)            :: ninput
      real(rk), parameter :: secs_pr_day = 86400.0_rk

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         ! Retrieve current (local) state variable values.
         _GET_(self%id_phytoplankton, phytoplankton)         ! phytoplankton

         !if ( _AVAILABLE_(self%id_grazing) ) then
         !   _GET_(self%id_grazing, grazing)
         !else
         !   grazing = 0
         !end if

         ! Retrieve current environmental conditions.
         _GET_(self%id_par,par)          ! local photosynthetically active radiation
         !_GET_(self%id_n, n)

         ! Nitrogen concentration
         _GET_(self%id_n, n)

         ! Production
         production = light_absorb(self%rmax, self%gamma, par) * self%rmax * uptake(self%Kn, n)

         ! Losses
         respiration = self%resp * uptake(self%Kn, n) * self%Vpotn
         sinking = self%s0
         !grazing = holling2(self%amax, self%pred_d, self%z, p)

         ! Nitrogen dynamics
         ninput = sinking * (phytoplankton + self%p0) * self%nremin - self%Vpotn * uptake(self%Kn, n) * (phytoplankton + self%p0) + (self%n0 - n) * self%nremin

         ! Set temporal derivatives
         _ADD_SOURCE_(self%id_phytoplankton, (production  - sinking - respiration) * (phytoplankton + self%p0) )
         _ADD_SOURCE_(self%id_n, ninput )

      ! Leave spatial loops (if any)
      _LOOP_END_
   end subroutine do

! ---------- Model Functions ------------ !

   elemental real(rk) function light_absorb(rmax, gamma, Ik)
      real(rk), intent(in) :: rmax, gamma, Ik

      light_absorb = 1 - exp( -gamma * Ik / rmax )
   end function light_absorb

   elemental real(rk) function uptake(Kn, n)
      real(rk), intent(in) :: Kn, n

      uptake = n / ( Kn + n )
   end function uptake

!end module examples_npzd_phy
end module tame_phytplankton

!-----------------------------------------------------------------------
! Copyright Bolding & Bruggeman ApS - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
