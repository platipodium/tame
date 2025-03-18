#include "fabm_driver.h"
! converts biological unit d-1 into physical FABM/driver unit s-1 for RHS!
! physical drivers usually have s-1 as a natural time unit.  So we
! need to convert the biological time unit d-1 into s-1.
! We agree to make all rate constants and parameters in fabm.yaml in d-1.
! Conversion occurs at fabm host interaciton (_ADD_SOURCE_)
! We use `secs_per_day` and `days_per_sec` for this conversion.

!module examples_npzd_phy ! TAME phytoplankton module
module tame_phytoplankton

use fabm_types
! for demonstration only
use tame_types, only : secs_per_day => tame_secs_per_day , days_per_sec, small
use tame_types
use tame_functions
   implicit none
   private
   !type, extends(type_base_model), public :: type_examples_npzd_phy
   type, extends(type_base_model), public :: type_tame_phytoplankton
      ! Variable identifiers
      type (type_state_variable_id)      :: id_phytoplankton     ! Phytoplankton biomass
!      type (type_state_variable_id)      :: id_no3, id_nh4, id_po4     ! id_din Nutrients
      type (type_state_variable_id)      :: id_var(NUM_CHEM+2*NUM_ELEM) ! TODO : flexible num of DOM & POM

      type (type_dependency_id)          :: id_par   ! PAR light
      type (type_diagnostic_variable_id) :: id_nut,id_nut2,id_rate
      !type (type_dependency_id)          :: id_grazing
      !type (type_dependency_id)          :: id_n     ! Nutrient
      !type (type_surface_dependency_id)  :: id_I_0   ! Surface irradiance
      !type (type_dependency_id)          :: id_z     ! Zooplankton
      ! Model parameters
      real(rk) :: p0
      real(rk) :: rmax                ! Growth metabolism parameters
      real(rk) :: gamma               ! Light exploitation
      real(rk) :: s0                  ! Sinking
      real(rk) :: resp                ! Respiration parameters
      real(rk) :: K_P, K_N            !
      real(rk) :: nut_limitation(NUM_NUTRIENT) ! Vector of limitation degree
      real(rk) :: HalfSatNut(NUM_NUTRIENT) ! Vector of half-saturations
      !real(rk) :: uptake_chemicals(NUM_GROWTH_CHEM) ! Vector of uptake chemicals

   contains
      procedure :: initialize
      procedure :: do
   end type

   integer :: det_index(NUM_ELEM), dom_index(NUM_ELEM) ! TODO: global setting!

contains
   subroutine initialize(self, configunit)
      class (type_tame_phytoplankton), intent(inout), target :: self
      integer,                        intent(in)            :: configunit
      integer :: i ! Indice dummy
      !real(rk), parameter :: days_per_sec = 1.0_rk/86400.0_rk
      integer, parameter :: counter = 0
      character :: elem

      ! Store parameter values in our own derived type
      ! NB: all rates must be provided in values per day and are converted here to values per second.
      call self%get_parameter(self%rmax, 'rmax', 'd-1',         'maximum production rate',  default=2.5_rk) !, scale_factor=days_per_sec
      call self%get_parameter(self%gamma,'gamma','microE-1 m-2','light absorption scaling', default=1.0_rk)
      call self%get_parameter(self%p0,   'p0',   'mmol m-3',    'background concentration ',default=0.0225_rk)
      call self%get_parameter(self%s0,   's0',   'd-1',         'default sinking rate',     default=0._rk) !, scale_factor=days_per_sec
      call self%get_parameter(self%K_P,  'K_P',  'mmol m-3',    'P half-saturation',        default=0.4_rk)
      call self%get_parameter(self%K_N,  'K_N',  'mmol m-3',    'N half-saturation',        default=4.0_rk)
      call self%get_parameter(self%resp, 'resp', 'mmol',        'carbon cost per nitrogen uptake',    default=0.2_rk)
      ! TODO redesign with transparent indices
      ! Also redesign get_parameter call from auto-generated parameter name like
      ! 'K_'//(self%halfsat(i)%name)
      self%HalfSatNut(1) = self%K_N
      self%HalfSatNut(2) = self%K_P

      ! Register state variables
      call self%register_state_variable(self%id_phytoplankton,'phytoplankton', 'mmol m-3', 'concentration', 1.0_rk, minimum=0.0_rk)

      ! Register environmental dependencies
      !call self%register_dependency(self%id_grazing, "grazing", 'd-1', 'grazing pressure', required = .false.)!, scale_factor = days_per_sec) ! Zooplankton activity
      call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)

      !do i = 1, NUM_GROWTH_CHEM
      !   call self%get_parameter(self%nut_limitations(i), 'dummy', 'mmol', 'dummy uptake',    default=0.2_rk)
      !   !call self%get_parameter(self%uptake_chemicals(i), 'dummy', 'mmol', 'dummy name',    default=0.2_rk)
      !end do
      !call self%register_dependency(self%id_no3,     'no3', 'mmol m-3', 'concentration', required =.false.)
!      call self%register_state_dependency(self%id_no3,  'NO3', 'mmol m-3', 'NO3 concentration', required =.true.) ! Needs to be turned into a loop
!      call self%register_state_dependency(self%id_nh4,  'NH4', 'mmol m-3', 'NH4 concentration', required =.true.) ! Needs to be turned into a loop
!      call self%register_state_dependency(self%id_po4,  'PO4', 'mmol m-3', 'PO4 concentration', required =.true.)
      do i = 1,NUM_CHEM !
         call self%register_state_dependency(self%id_var(i), chemicals(i),'mmol m-3',chemicals(i))
      end do
      !  retrieve OM variables for each element
      do i = 1,NUM_ELEM ! e.g., C, N, P (Si, Fe)
         det_index(i) = NUM_CHEM+2*i-1
         dom_index(i) = NUM_CHEM+2*i
         elem = ElementList(i:i)
         call self%register_state_dependency(self%id_var(det_index(i)), 'det_' // elem,'mmol-' // elem // ' m-3','Detritus ' // trim(ElementName(i)))
         call self%register_state_dependency(self%id_var(dom_index(i)), 'dom_' // elem,'mmol-' // elem // ' m-3','Dissolved Organic ' // trim(ElementName(i)))
  ! print *,det_index(i), ElementList(i:i),dom_index(i)
      end do
      call self%register_diagnostic_variable(self%id_nut, 'nut1','mmol-? m-3', 'nutrient related')
      call self%register_diagnostic_variable(self%id_nut2,'nut2','mmol-? m-3', 'nutrient related')
      call self%register_diagnostic_variable(self%id_rate,'rate','mmol-? m-3', 'nutrient related')

   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_tame_phytoplankton), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk)            :: phytoplankton, par, din, func,rhs_phy,nh4_part,no3_part
      real(rk)            :: production, respiration, sinking, new, nut_lim_tot
      real(rk)            :: nutrient_lim(NUM_NUTRIENT), exudation(NUM_NUTRIENT)
      real(rk)            :: nutrient(NUM_NUTRIENT), din_no3, din_nh4
      real(rk)            :: dix_chemical(NUM_CHEM), part(NUM_CHEM)
      logical             :: ncrit(NUM_CHEM)
      integer             :: i ! Index
      real(rk)            :: chem_stoichiometry(NUM_CHEM)=(/0.0625, 0.0625, 0.0094/) ! Redfield TODO

      ! Enter spatial loops (if any)
      !print *,stoichiometry
      _LOOP_BEGIN_
         ! Retrieve current (local) state variable values.
         _GET_(self%id_phytoplankton, phytoplankton)     ! phytoplankton carbon

         !do i = 1, NUM_NUTRIENT
         !_GET_( self%id_nutrient(i), nutrient(i) ) ! Nutrient target for later
         !end do

         do i = 1,NUM_CHEM ! e.g., CO2, NO3, NH4 (PO4)
            _GET_(self%id_var(i), dix_chemical(i))  ! Dissolved Inorganic Nutrient DIX in mmol-X/m**3
         end do
         ! TODO replace by TransIndex_DOMDIX, TransIndex2_DOMDIX, which should be set globally
         nutrient(1) = dix_chemical(1) + dix_chemical(2)
         nutrient(2) = dix_chemical(3)
         ! should be a vector * matrix multiplication, fortran MATMUL

         ! Retrieve current environmental conditions.
         _GET_(self%id_par,par)          ! local photosynthetically active radiation
         ! Nutrient uptake
         !do i = 1, NUM_NUTRIENT
         !uptake(i) = limitation( nutient(i), self%nut_limitation(i), self%affinity(i) ) ! Nutrient target
         !end do
         nut_lim_tot = 0._rk
         do i = 1,NUM_NUTRIENT
            nutrient_lim(i) = nutrient(i)/(self%HalfSatNut(i)+nutrient(i))
            nut_lim_tot = nut_lim_tot + 1.0_rk/(nutrient_lim(i) + small)
!            nutrient_lim(i) = limitation( self%affinity(i)*nutrient(i)) !/ chem_stoichiometry(i)  ! Add the nutrient limitation law for phytoplankton
         end do
         nut_lim_tot = 1.0_rk/nut_lim_tot
      !   _SET_DIAGNOSTIC_(self%id_nut, nutrient_lim(1) )
         ! Production
         func = 1.0_rk - exp( -self%gamma * par / self%rmax )
         production = self%rmax * func * nut_lim_tot
!         production = self%rmax * light_absorb(self%rmax, self%gamma, par) !* minval( nutrient_lim )
      !   _SET_DIAGNOSTIC_(self%id_rate, production )

         !do i = 1,NUM_NUTRIENT
         !   exudation(i) =  ( uptake(i) - minval( uptake ) ) * chem_stoichiometry(i) ! Add DOX as a dependency
         !end do

         ! Losses
         respiration = self%resp  ! C loss for DIN-uptake
         sinking = self%s0

         ! Nutrient dynamics
         ! 3 cases for NO3+NH4 partitioning in DIN usage
         ! both low: 50%, both high: relational, one low: 0+100%
         part = 1._rk
         do i = 1,2
           ncrit(i) = (dix_chemical(i) .LT. self%HalfSatNut(1)/10 )
           if (ncrit(i)) part(i) = 0._rk
         end do
         if (ncrit(1) .AND. ncrit(2))  part(1) = 0.5_rk
         if (.NOT.(ncrit(1)) .AND. .NOT.(ncrit(2)))  part(1) = dix_chemical(1)/nutrient(1)
         part(2) = 1._rk - part(1)

!         _SET_DIAGNOSTIC_(self%id_nut, part(1)*new * chem_stoichiometry(1) )

         ! nutrient uptake = new production times Redfield chem_stoichiometry -> passed to BGC DIX variables
         new = production * phytoplankton
         _SET_DIAGNOSTIC_(self%id_rate,  new)

         do i = 1,NUM_CHEM
            _ADD_SOURCE_(self%id_var(i), -part(i)*new* chem_stoichiometry(i) UNIT) ! Nutrients sink
         end do
         ! temporal derivative for phytoplankton C
         rhs_phy = (production  - sinking - respiration) * (phytoplankton + self%p0)
         _ADD_SOURCE_(self%id_phytoplankton, rhs_phy  UNIT)

         ! Exudation to DOM (proportional to C-respiration)
         new = respiration * phytoplankton
         do i = 1,NUM_ELEM
            if (ElementList(i:i) .NE. 'C') _ADD_SOURCE_(self%id_var(dom_index(i)), new*stoichiometry(i) UNIT) ! Nutrients sink
         end do
         _SET_DIAGNOSTIC_(self%id_nut, new*stoichiometry(2) )
         _SET_DIAGNOSTIC_(self%id_nut2,new * stoichiometry(3) )

         ! sinking to POM
         new = sinking * phytoplankton

         do i = 1,NUM_ELEM  ! C, N, P (Si, Fe)
            _ADD_SOURCE_(self%id_var(det_index(i)), new*stoichiometry(i) UNIT) !
         end do

      ! Leave spatial loops (if any)
      _LOOP_END_
   end subroutine do

! ---------- Model Functions ------------ !

   elemental real(rk) function light_absorb(rmax, gamma, Ik)
      real(rk), intent(in) :: rmax, gamma, Ik
      light_absorb = 1 - exp( -gamma * Ik / rmax )
   end function light_absorb

   elemental real(rk) function limitation( Nut, Aff)
      real(rk), intent(in) :: Aff, Nut
      limitation = Aff*Nut/(Aff*Nut + 1)
   end function limitation

!end module examples_npzd_phy
end module tame_phytoplankton
