! SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
!
! SPDX-License-Identifier: Apache-2.0

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

      !! Dependency ids
      !type (type_dependency_id) :: id_Q(NUM_ELEM) ! Prey stoichiometry
      type (type_state_variable_id) :: id_prey ! prey

      !! Diagnostic variable ids
      !type (type_diagnostic_variable_id) :: id_nut, id_nut2, id_rate
      type (type_diagnostic_variable_id) :: id_dummy

      !! Model parameters
      real(rk) :: max_ingestion, saturation ! maximum ingestion rate, food saturation
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
      integer, parameter :: counter = 0 !@what
      character :: elem

      !! Store parameter values in our own derived type
      !! all rates must be provided in values per day and are converted here to values per second.
      !! @todo: add min and max values
      call self%get_parameter(self%saturation, 'saturation', 'mmol-C m-3', 'grazing saturation', default=2.5_rk)
      call self%get_parameter(self%max_ingestion, 'max_ingestion', 'd-1', 'maximum ingestion rate', default=2.5_rk)
      call self%get_parameter(self%resp, 'resp', 'd-1', 'respiration rate', default=0._rk)

      !! Register state variables
      call self%register_state_variable(self%id_biomass,'biomass', 'mmol-C m-3', 'concentration', 1.0_rk, minimum=0.0_rk)

      !! Register environmental dependencies

      !! Register external dependencies
      call self%register_state_dependency(self%id_prey, 'prey','mmol-C m-3', 'prey source')

      !! Retrieve prey stoichiometric composition (if FlexStoich)
      !do i = 1,NUM_ELEM
      !   elem = ElementList(i:i)
      !   if (elem .NE. 'C') then
      !      call self%register_dependency(self%id_Q(i), 'PreyQ_' // elem,'mol-' // elem // ' mol-C-1', elem // ':C-quota', require = .false.)
      !   endif
      !end do ! Prey is supposed to be in carbon, so only extracting other elements, if available

      call self%register_diagnostic_variable(self%id_dummy, 'dummy','', '')

      ! call self%register_diagnostic_variable(self%id_nut, 'nut1','mmol-? m-3', 'nutrient related')
      ! call self%register_diagnostic_variable(self%id_nut2,'nut2','mmol-? m-3', 'nutrient related')
      ! call self%register_diagnostic_variable(self%id_rate,'rate','mmol-? m-3', 'nutrient related')
   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_tame_zooplankton), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_
      integer :: i ! dummy index
      real(rk) :: biomass, prey, func
      real(rk) :: production, respiration, new
      real(rk) :: exudation(NUM_NUTRIENT)

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_
      _GET_(self%id_biomass, biomass) ! biomass carbon

      !@what: is zooplankton exudating to nutrient pool or to DOM pool?
      !do i = 1, NUM_NUTRIENT
      !   _GET_( self%id_nutrient(i), nutrient(i) ) ! Nutrient target for later
      !end do
      ! Excretion to ammonia NH4+?

      !! Retrieve current environmental conditions.
      _GET_(self%id_prey,prey)

      !! Predation
      func = 1.0_rk - exp(- prey / self%saturation )
      production = self%max_ingestion * func
      ! _SET_DIAGNOSTIC_(self%id_rate, production )

      !! Exhudation
      !@what: is zooplankton exudating to nutrient pool or to DOM pool
      !do i = 1,NUM_NUTRIENT
         !!exudation(i) = ( uptake(i) - minval( uptake ) ) * chem_stoichiometry(i) ! Add DOX as a dependency
      !end do

      !! Losses
      respiration = self%resp !@todo: include other physiological loss terms here

      _ADD_SOURCE_(self%id_prey, -production*biomass * days_per_sec ) ! UNIT Prey is grazed on
      _ADD_SOURCE_(self%id_biomass, (production - respiration) * biomass * days_per_sec ) ! UNIT
      _SET_DIAGNOSTIC_(self%id_dummy, production)

      !! Exudation to DOM (proportional to C-respiration)
      !new = respiration * biomass
      !do i = 1,NUM_ELEM
         !if (ElementList(i:i) .NE. 'C') _ADD_SOURCE_(self%id_var(dom_index(i)), new*stoichiometry(i) UNIT) ! Nutrients sink
      !end do
      !_SET_DIAGNOSTIC_(self%id_nut, new*stoichiometry(2) )
      !_SET_DIAGNOSTIC_(self%id_nut2, new*stoichiometry(3) )

      ! sinking to POM
      !new = sinking * biomass

      !! @what: should not be sinking be managed by physical driver? I see this useful probably only on 0d setups. A flag is required!
      !do i = 1,NUM_ELEM ! C, N, P (Si, Fe)
      !   _ADD_SOURCE_(self%id_var(det_index(i)), new*stoichiometry(i) UNIT) !
      !end do

      ! Leave spatial loops (if any)
      _LOOP_END_
   end subroutine do

end module tame_zooplankton
