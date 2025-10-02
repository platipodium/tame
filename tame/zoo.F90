! SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
! SPDX-FileContributor Thomas Imbert <thomas.imbert@hereon.de>
! SPDX-FileContributor Ovidio Garcia <ovidio.garcia@hereon.de>
! SPDX-FileContributor Carsten Lemmen <carsten.lemmen@hereon.de>
! SPDX-License-Identifier: Apache-2.0

#include "fabm_driver.h"

!! TAME zooplankton--any generic predator--(fixed stoichiometry) module
module tame_zooplankton
   use fabm_types
   use tame_types
   use tame_functions
   implicit none

   private
   type, extends(type_base_model), public :: type_tame_zooplankton

      !! State variables ids
      type (type_state_variable_id) :: id_zooplankton_C ! zooplankton_C
      !! Dependency ids
      type (type_dependency_id)     :: id_prey_Q(NUM_ELEM) ! Prey and stoichiometry
      type(type_dependency_id)      :: id_temp ! for temperature dependence
      type (type_state_variable_id) :: id_prey, id_dom(NUM_ELEM), id_det(NUM_ELEM) ! prey
      !! Diagnostic variable ids
      type (type_diagnostic_variable_id) :: id_dummye, id_dummya, id_dummy_N, id_dummy_P, id_rhs_phyC
      type (type_diagnostic_variable_id) :: id_zoo_elem(NUM_ELEM) ! Zooplankton stoichiometry
      !! Model parameters
      real(rk) :: max_ingestion, saturation, sloppy ! maximum ingestion rate, food saturation
      real(rk) :: resp, Q10 ! respiration rate

      contains
      procedure :: initialize
      procedure :: do
   end type

   integer :: det_index(NUM_ELEM), dom_index(NUM_ELEM) ! TODO: global setting!

   contains
!! ADD ELEM TO AGGREGATE VARIABLE (like in total_X in bgc)
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
      call self%get_parameter(self%sloppy, 'sloppy', '', 'sloppy feeding', default=0._rk)
      call self%get_parameter(self%Q10, 'Q10', '--', 'temperature sensitivity (ref temperature 20 degC)', default=2.2_rk)

      !! Register state variables
      call self%register_state_variable(self%id_zooplankton_C,'zooplankton_C', 'mmol-C m-3', 'concentration', 1.0_rk, minimum=0.0_rk)

      !! Register environmental dependencies
      call self%register_dependency(self%id_temp, standard_variables%temperature)

      !! Register external dependencies
      call self%register_state_dependency(self%id_prey, 'prey','mmol-C m-3', 'prey source')
      do i = 1,NUM_ELEM
         elem = ElementList(i:i)

         call self%register_state_dependency(self%id_det(i), 'det_' // elem,'mol-' // elem // ' mol-C-1', 'Detritus ' // elem)
         call self%register_state_dependency(self%id_dom(i), 'dom_' // elem,'mol-' // elem // ' mol-C-1', 'Dissolved Organic ' // elem)
         if (elem .NE. 'C') then
            ! Element stoichiometry in Zooplankton
            call self%register_diagnostic_variable(self%id_zoo_elem(i), 'zoo_' // elem,'mol-' // elem // ' m-3', 'zooplankton ' // elem)
            call self%add_to_aggregate_variable(type_universal_standard_variable(name='total_' // trim(ElementName(i)), units='mmol-' // elem // ' m-3', &
              aggregate_variable=.true., conserved=.true.), self%id_zoo_elem(i))
      !! Retrieve prey stoichiometric composition (if FlexStoich)
            call self%register_dependency(self%id_prey_Q(i), 'PreyQ_' // elem,'mol-' // elem // ' mol-C-1', elem // ':C-quota prey')
         endif
      end do ! Prey is supposed to be in carbon, so only extracting other elements, if available

      call self%register_diagnostic_variable(self%id_rhs_phyC, 'rhs_phyC','mol-C m-3 d-1', 'phy-C grazing rate')
      call self%register_diagnostic_variable(self%id_dummye, 'dummye','', '')
      call self%register_diagnostic_variable(self%id_dummya, 'dummya','', '')
      call self%register_diagnostic_variable(self%id_dummy_N, 'dummy_N','', '')
      call self%register_diagnostic_variable(self%id_dummy_P, 'dummy_P','', '')

   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_tame_zooplankton), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_
      integer :: i ! dummy index
      real(rk) :: zooplankton_C, prey, ivlev
      real(rk) :: ingestion, feeding, respiration, new, growth
      real(rk) :: elem_Q(NUM_ELEM), Excess_C_upt(NUM_ELEM) !
      real(rk) :: excretion(NUM_ELEM)            ! exudation of C,N,P
      real(rk) :: temp, temp_factor
      character :: elem

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_
      _GET_(self%id_zooplankton_C, zooplankton_C) ! zooplankton_C carbon

      !! Retrieve nutrient quotas of prey
      do i = 2, NUM_ELEM
         _GET_( self%id_prey_Q(i), elem_Q(i) ) ! elem_Q target for later
      end do

      !! Retrieve current environmental conditions.
      _GET_(self%id_prey,prey)
      _GET_(self%id_temp, temp) ! temperature
      temp_factor = self%Q10**(0.1_rk*(temp-20.0_rk))

      !! Predation
      ivlev     =  1.0_rk - exp(- prey / self%saturation )
      feeding   =  temp_factor * self%max_ingestion * ivlev
      ingestion =  (1.0_rk - self%sloppy) * feeding 

      !! respiratory losses : TODO: link to feeding activity
      respiration = temp_factor * self%resp !@todo: include other physiological loss terms here

      Excess_C_upt(1) = - respiration
      do i = 2,NUM_ELEM
         Excess_C_upt(i) = (1._rk - elem_Q(i)/zoo_stoichiometry(i)) * ingestion - respiration 
      end do
      excretion(1) = max(0.0, maxval(Excess_C_upt)) 
      Excess_C_upt = Excess_C_upt - excretion(1)
      excretion(2:NUM_ELEM)  = zoo_stoichiometry(2:NUM_ELEM) * (-Excess_C_upt(2:NUM_ELEM))

      ! C-based growth rate of grazer 
      growth = ingestion - respiration - excretion(1) ! 
      _ADD_SOURCE_(self%id_zooplankton_C, growth * zooplankton_C * days_per_sec ) ! Zooplankton growth
      !! TODO: Zooplankton can store C in fat, should not excrete it directly, maybe as a state variable?

      ! C-based loss rate of grazed prey 
       _ADD_SOURCE_(self%id_prey, -feeding * zooplankton_C * days_per_sec ) ! Prey is ingested
      ! additional storage as diag since the rate is needed for mass balance calculation in the prey (phy) module

      _SET_DIAGNOSTIC_(self%id_rhs_phyC, -feeding * zooplankton_C )


      ! Resolving the DOM and POM targets with biogeochemistry
      ! Excretion to ammonia NH4+?
      do i = 1, NUM_ELEM
         _ADD_SOURCE_( self%id_dom(i), excretion(i) * zooplankton_C * days_per_sec ) ! DOM target for later
         _ADD_SOURCE_( self%id_det(i), self%sloppy * feeding* zoo_stoichiometry(i)* zooplankton_C * days_per_sec ) ! POM target for later
         if (i .gt. 1) _SET_DIAGNOSTIC_(self%id_zoo_elem(i), zoo_stoichiometry(i) * zooplankton_C )
      end do

      _SET_DIAGNOSTIC_(self%id_dummya, Excess_C_upt(1) )! maxval(excretion) )
      _SET_DIAGNOSTIC_(self%id_dummye, excretion(1))
      _SET_DIAGNOSTIC_(self%id_dummy_N,  excretion(2))
      _SET_DIAGNOSTIC_(self%id_dummy_P,  excretion(3) )

      ! sinking to POM
      !new = sinking * zooplankton_C

      ! Leave spatial loops (if any)
      _LOOP_END_
   end subroutine do

end module tame_zooplankton
