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
      type (type_state_variable_id)      :: id_din, id_po4     ! Nutrients

      
      type (type_dependency_id)          :: id_par   ! PAR light
      
      !type (type_dependency_id)          :: id_grazing 

      !type (type_dependency_id)          :: id_n     ! Nutrient

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
      real(rk) :: n0 ! Background N
      real(rk) :: nremin ! N remineralization 

      real(rk) :: uptake_rate(NUM_NUTRIENT) ! Vector of uptake rates
      real(rk) :: affinity(NUM_NUTRIENT) ! Vector of uptake rates
      !real(rk) :: uptake_chemicals(NUM_GROWTH_CHEM) ! Vector of uptake chemicals

   contains
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class (type_tame_phytoplankton), intent(inout), target :: self
      integer,                        intent(in)            :: configunit
      integer :: i ! Indice dummy

      !real(rk), parameter :: days_per_sec = 1.0_rk/86400.0_rk
      integer, parameter :: counter = 0

      ! Store parameter values in our own derived type
      ! NB: all rates must be provided in values per day and are converted here to values per second.
      call self%get_parameter(self%rmax,    'rmax',    'd-1',          'maximum production rate',            default=8.0_rk, scale_factor=days_per_sec)
      call self%get_parameter(self%gamma,   'gamma',   'microE-1 m-2', 'light absorption rate scaling',      default=1.0_rk)
      call self%get_parameter(self%p0,    'p0',    'mmol m-3',  'background concentration ',                default=0.0225_rk)

      call self%get_parameter(self%s0,      's0',      'd-1',          'default sinking rate',               default=0.1_rk, scale_factor=days_per_sec)

      call self%get_parameter(self%Vpotn,   'Vpotn',   'mmol.d-1',      'N potential uptake rate',            default=8.0_rk, scale_factor=days_per_sec)
      call self%get_parameter(self%Kn,   'Kn',   'mmol m-3',      'N demi saturation',            default=8.0_rk)

      call self%get_parameter(self%n0,       'n0',       'mmol m-3',     'Background Nitrogen',                default=1.0_rk)

      call self%get_parameter(self%resp,     'resp',     'mmol',         'carbon cost per nitrogen uptake',    default=0.2_rk)

      do i = 1, NUM_NUTRIENT
         call self%get_parameter(self%uptake_rate(i), 'dummy', 'mmol', 'dummy uptake',    default=0.2_rk)
         call self%get_parameter(self%affinity(i), 'dummy', 'mmol', 'dummy affinity',    default=0.2_rk)
!call self%get_parameter(self%uptake_chemicals(i), 'dummy', 'mmol', 'dummy name',    default=0.2_rk)
      end do

      ! Register state variables
      call self%register_state_variable(self%id_phytoplankton,    'phytoplankton', 'mmol m-3', 'concentration', 1.0_rk, minimum=0.0_rk)

      ! Register environmental dependencies
      !call self%register_dependency(self%id_grazing, "grazing", 'd-1', 'grazing pressure', required = .false.)!, scale_factor = days_per_sec) ! Zooplankton activity
      call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)

      !do i = 1, NUM_GROWTH_CHEM
      !   call self%get_parameter(self%uptake_rates(i), 'dummy', 'mmol', 'dummy uptake',    default=0.2_rk)
      !   !call self%get_parameter(self%uptake_chemicals(i), 'dummy', 'mmol', 'dummy name',    default=0.2_rk)
      !end do
      !call self%register_dependency(self%id_no3,     'no3', 'mmol m-3', 'concentration', required =.false.)
      !call self%register_dependency(self%id_nh4,     'nh4', 'mmol m-3', 'concentration', required =.false.)
      call self%register_state_dependency(self%id_din,  'din', 'mmol m-3', 'DIN concentration', required =.true.) ! Needs to be turned into a loop
      call self%register_state_dependency(self%id_po4,  'po4', 'mmol m-3', 'po4 concentration', required =.true.)

   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_tame_phytoplankton), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk)            :: phytoplankton, par
      real(rk)            :: production, respiration, sinking
      real(rk)            :: uptake(NUM_NUTRIENT)
      real(rk)            :: nutrient(NUM_NUTRIENT)
      real(rk)            :: exudation(NUM_NUTRIENT)
      integer             :: i ! Indice dummy

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         ! Retrieve current (local) state variable values.
         _GET_(self%id_phytoplankton, phytoplankton)         ! phytoplankton

         !do i = 1, NUM_NUTRIENT
         !_GET_( self%id_nutrient(i), nutrient(i) ) ! Nutrient target for later
         !end do

         _GET_( self%id_din, nutrient(1) )
         _GET_( self%id_po4, nutrient(2) ) 


         ! Retrieve current environmental conditions.
         _GET_(self%id_par,par)          ! local photosynthetically active radiation

         ! Nutrient uptake
         !do i = 1, NUM_NUTRIENT
         !uptake(i) = nutrient_uptake( nutient(i), self%uptake_rate(i), self%affinity(i) ) ! Nutrient target
         !end do
         do i = 1,NUM_NUTRIENT
            uptake(i) = nutrient_uptake( nutrient(i), self%uptake_rate(i), self%affinity(i) ) / stoichiometry(i)  ! Add the nutrient limitation law for phytoplankton
         end do

         ! Production
         production = light_absorb(self%rmax, self%gamma, par) * self%rmax * minval( uptake )

         do i = 1,NUM_NUTRIENT
            exudation(i) =  ( uptake(i) - minval( uptake ) ) * stoichiometry(i) ! Add DOX as a dependency
         end do

         ! Losses 
         respiration = self%resp * uptake(1) * self%Vpotn ! C loss for DIN-uptake
         sinking = self%s0

         ! Nutrient dynamics
         !_ADD_SOURCE_(self%id_nut(i), uptake(i) * stoichiometry(i) ) ! Nutrients sink 
         !_ADD_SOURCE_(self%id_dox(i), exudation(i) ) 
         ! TO DO : Need to take away N, but DIN is a diagnostic of bgc.F90. How to know whether to tak NO3 or NH4 away?

         ! Set temporal derivatives
         _ADD_SOURCE_(self%id_phytoplankton, (production  - sinking - respiration) * (phytoplankton + self%p0) )
         
         ! Exhudation to DOM

      ! Leave spatial loops (if any)
      _LOOP_END_
   end subroutine do

! ---------- Model Functions ------------ !

   elemental real(rk) function light_absorb(rmax, gamma, Ik)       
      real(rk), intent(in) :: rmax, gamma, Ik

      light_absorb = 1 - exp( -gamma * Ik / rmax )
   end function light_absorb

   elemental real(rk) function nutrient_uptake(Kn, n, p)       
      real(rk), intent(in) :: Kn, n, p

      nutrient_uptake = n / ( Kn + n )
   end function nutrient_uptake
   
!end module examples_npzd_phy
end module tame_phytoplankton