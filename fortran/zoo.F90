#include "fabm_driver.h"
!! converts biological unit d-1 into physical FABM/driver unit s-1 for RHS
!@todo: avoid preprocessor calls and use parameters
#define UNIT *1.1574074074E-5_rk 

!! TAME zooplankton--any generic predator--(fixed stoichiometry) module
module tame_zooplankton
   use fabm_types
   use tame_types
   use tame_functions
   implicit none

   private
   type, extends(type_base_model), public :: type_tame_zooplankton
      !! State variables ids
      type (type_state_variable_id) :: id_biomass ! biomass
      type (type_state_variable_id) :: id_var(NUM_CHEM+2*NUM_ELEM) ! nutrients for exudation

      !! Dependency ids
      type (type_dependency_id) :: id_par ! light
      type (type_dependency_id) :: id_prey ! prey

      !! Diagnostic variable ids
      type (type_diagnostic_variable_id) :: id_nut, id_nut2, id_rate

      !! Model parameters
      real(rk) :: p0 ! background concentration
      real(rk) :: rmax ! maximum ingestion rate
      real(rk) :: gamma ! ivlev coefficient
      real(rk) :: s0 ! Sinking rate
      real(rk) :: resp ! Respiration rate

      contains
      procedure :: initialize
      procedure :: do
   end type

   integer :: det_index(NUM_ELEM), dom_index(NUM_ELEM) ! TODO: global setting!

   contains

   subroutine initialize(self, configunit)
      class (type_tame_zooplankton), intent(inout), target :: self
      integer, intent(in) :: configunit
      integer :: i ! dummy index
      !real(rk), parameter :: days_per_sec = 1.0_rk/86400.0_rk
      integer, parameter :: counter = 0 !@what
      character :: elem

      !! Store parameter values in our own derived type
      !! all rates must be provided in values per day and are converted here to values per second.
      !! @todo: add min and max values
      call self%get_parameter(self%rmax, 'rmax', 'd-1', 'maximum ingestion rate', default=2.5_rk) 
      call self%get_parameter(self%gamma,'gamma','mmol-C m-3', 'Ivlev coeficient', default=1.0_rk)
      call self%get_parameter(self%p0, 'p0', 'mmol-C m-3', 'background concentration ',default=0.0225_rk)
      call self%get_parameter(self%s0, 's0', 'm d-1', 'default sinking velocity', default=0._rk) 
      call self%get_parameter(self%resp, 'resp', 'd-1', 'respiration rate', default=0._rk) 

      !! Register state variables
      call self%register_state_variable(self%id_biomass,'biomass', 'mmol-C m-3', 'concentration', 1.0_rk, minimum=0.0_rk)

      !! Register environmental dependencies
      call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
      
      !! Register external dependencieS
      call self%register_state_dependency(self%id_prey, prey,'mmol-C m-3')

      do i = 1,NUM_CHEM !
         call self%register_state_dependency(self%id_var(i), chemicals(i),'mmol m-3',chemicals(i))
      end do

      !! retrieve OM variables for each element
      do i = 1,NUM_ELEM ! e.g., C, N, P (Si, Fe)
         det_index(i) = NUM_CHEM+2*i-1
         dom_index(i) = NUM_CHEM+2*i
         elem = ElementList(i:i)
         call self%register_state_dependency(self%id_var(det_index(i)), 'det_' // elem,'mmol-' // elem // ' m-3','Detritus ' // trim(ElementName(i)))
         call self%register_state_dependency(self%id_var(dom_index(i)), 'dom_' // elem,'mmol-' // elem // ' m-3','Dissolved Organic ' // trim(ElementName(i)))
         ! print *,det_index(i), ElementList(i:i),dom_index(i)
      end do

      ! call self%register_diagnostic_variable(self%id_nut, 'nut1','mmol-? m-3', 'nutrient related')
      ! call self%register_diagnostic_variable(self%id_nut2,'nut2','mmol-? m-3', 'nutrient related')
      ! call self%register_diagnostic_variable(self%id_rate,'rate','mmol-? m-3', 'nutrient related')
   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_tame_zooplankton), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_
      integer :: i ! dummy index
      real(rk) :: biomass, prey
      real(rk) :: production, respiration, sinking, new, nut_lim_tot
      real(rk) :: exudation(NUM_NUTRIENT)
      real(rk) :: nutrient(NUM_NUTRIENT), din_no3, din_nh4
      real(rk) :: dix_chemical(NUM_CHEM), part(NUM_CHEM)
      logical :: ncrit(NUM_CHEM)

      real(rk) :: chem_stoichiometry(NUM_CHEM)=(/0.0625, 0.0625, 0.0094/) ! Redfield TODO !@rodo: change to predator stoichiometry

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_
      _GET_(self%id_biomass, biomass) ! biomass carbon

      !@what: is zooplankton exudating to nutrient pool or to DOM pool?
      do i = 1, NUM_NUTRIENT
         _GET_( self%id_nutrient(i), nutrient(i) ) ! Nutrient target for later
      end do

      do i = 1,NUM_CHEM ! e.g., CO2, NO3, NH4 (PO4)
         _GET_(self%id_var(i), dix_chemical(i)) 
      end do

      ! TODO replace by TransIndex_DOMDIX, TransIndex2_DOMDIX, which should be set globally
      nutrient(1) = dix_chemical(1) + dix_chemical(2)
      nutrient(2) = dix_chemical(3)

      !! Retrieve current environmental conditions.
      _GET_(self%id_par,par) ! local photosynthetically active radiation
      _GET_(self%id_prey,prey) ! local photosynthetically active radiation

      !! Predation
      func = 1.0_rk - exp(-self%gamma*prey/ self%rmax )
      production = self%rmax * func 
      ! _SET_DIAGNOSTIC_(self%id_rate, production )

      !! Exhudation
      !@what: is zooplankton exudating to nutrient pool or to DOM pool
      do i = 1,NUM_NUTRIENT
         !!exudation(i) = ( uptake(i) - minval( uptake ) ) * chem_stoichiometry(i) ! Add DOX as a dependency
      end do

      !! Losses
      respiration = self%resp !@todo: include other physiological loss terms here

      !! adding predator and prey sources
      _ADD_SOURCE_(self%id_prey, -production*biomass UNIT) ! Prey is grazed on
      _ADD_SOURCE_(self%id_biomass, (production - sinking - respiration) * (biomass + self%p0) UNIT)

      !! Exudation to DOM (proportional to C-respiration)
      new = respiration * biomass
      do i = 1,NUM_ELEM
         !if (ElementList(i:i) .NE. 'C') _ADD_SOURCE_(self%id_var(dom_index(i)), new*stoichiometry(i) UNIT) ! Nutrients sink
      end do
      _SET_DIAGNOSTIC_(self%id_nut, new*stoichiometry(2) )
      _SET_DIAGNOSTIC_(self%id_nut2,new * stoichiometry(3) )

      ! sinking to POM 
      new = sinking * biomass

      !! @what: should not be sinking be managed by physical driver? I see this useful probably only on 0d setups. A flag is required!
      do i = 1,NUM_ELEM ! C, N, P (Si, Fe)
         _ADD_SOURCE_(self%id_var(det_index(i)), new*stoichiometry(i) UNIT) ! 
      end do 

      ! Leave spatial loops (if any)
      _LOOP_END_
   end subroutine do

end module tame_zooplankton
