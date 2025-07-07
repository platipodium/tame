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
      type (type_state_variable_id) :: id_biomass ! biomass

      !! Dependency ids
      type (type_dependency_id) :: id_prey_Q(NUM_ELEM) ! Prey stoichiometry
      type(type_dependency_id) :: id_temp ! for temperature dependence
      type (type_state_variable_id) :: id_prey, id_dom_(NUM_ELEM), id_det_(NUM_ELEM) ! prey

      !! Diagnostic variable ids
      !type (type_diagnostic_variable_id) :: id_nut, id_nut2, id_rate
      type (type_diagnostic_variable_id) :: id_dummy, id_dummy_sloppy

      !! Model parameters
      real(rk) :: max_ingestion, saturation, sloppy ! maximum ingestion_rate rate, food saturation
      real(rk) :: resp, Q10 ! respiration_rate rate

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
      call self%register_state_variable(self%id_biomass,'biomass', 'mmol-C m-3', 'concentration', 1.0_rk, minimum=0.0_rk)

      !! Register environmental dependencies
      call self%register_dependency(self%id_temp, standard_variables%temperature)

      !! Register external dependencies
      call self%register_state_dependency(self%id_prey, 'prey','mmol-C m-3', 'prey source')
      do i = 1,NUM_ELEM
         elem = ElementList(i:i)
         if (elem .NE. 'C') then
            call self%register_state_dependency(self%id_dom_(i), 'dom_' // elem,'mol-' // elem // ' mol-C-1', 'Dissolved Organic' // elem)

            call self%register_state_dependency(self%id_det_(i), 'det_' // elem,'mol-' // elem // ' mol-C-1', 'Particulate Organic' // elem)

      !! Retrieve prey stoichiometric composition (if FlexStoich)
            call self%register_dependency(self%id_prey_Q(i), 'PreyQ_' // elem,'mol-' // elem // ' mol-C-1', elem // ':C-quota')

         endif
      end do ! Prey is supposed to be in carbon, so only extracting other elements, if available

      call self%register_diagnostic_variable(self%id_dummy, 'dummy','', '')
      call self%register_diagnostic_variable(self%id_dummy_sloppy, 'dummy_sloppy','', '')

      ! call self%register_diagnostic_variable(self%id_nut, 'nut1','mmol-? m-3', 'nutrient related')
      ! call self%register_diagnostic_variable(self%id_nut2,'nut2','mmol-? m-3', 'nutrient related')
      ! call self%register_diagnostic_variable(self%id_rate,'rate','mmol-? m-3', 'nutrient related')
   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_tame_zooplankton), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_
      integer :: i ! dummy index
      real(rk) :: biomass, prey, ivlev
      real(rk) :: ingestion_rate, respiration_rate, new, growth_rate, assimilation_rate
      real(rk) :: exudation_rate(NUM_ELEM), elem_Q(NUM_ELEM), Q_diff(NUM_ELEM) ! elem_Q is the prey quotas
      real(rk) :: C_excess(NUM_ELEM), elem_assimilation(NUM_ELEM)              ! assimilation or excess of nutrients
      real(rk) :: sloppy_feeding(NUM_ELEM), excretion_rate(NUM_ELEM)
      real(rk) :: temp, temp_factor
      character :: elem

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_
      _GET_(self%id_biomass, biomass) ! biomass carbon

      !@what: is zooplankton exudating to elem_Q pool or to DOM pool?
      do i = 1, NUM_ELEM
         elem = ElementList(i:i)
         if (elem .NE. 'C') then
            _GET_( self%id_prey_Q(i), elem_Q(i) ) ! elem_Q target for later
         endif
      end do
      ! Excretion to ammonia NH4+?

      !! Retrieve current environmental conditions.
      _GET_(self%id_prey,prey)
      _GET_(self%id_temp, temp) ! temperature

      temp_factor = self%Q10**(0.1+rk*(temp-20.0_rk))

      !! Predation
      ivlev = 1.0_rk - exp(- prey / self%saturation )
      ingestion_rate =  temp_factor * self%max_ingestion * ivlev
      ! _SET_DIAGNOSTIC_(self%id_rate, ingestion_rate )

      !! Losses
      respiration_rate = temp_factor * self%resp !@todo: include other physiological loss terms here

      !! Exhudation
      !@what: is zooplankton exudating to elem_Q pool or to DOM pool
      !do i = 1,NUM_ELEM
      !   elem = ElementList(i:i)
      !   if (elem .NE. 'C') then
            !exudation_rate(i) = max(0.0_rk, elem_Q(i) -  zoo_fixed_stoichiometry(i) ) * ingestion_rate*biomass * self%sloppy  !
            ! Should exudate C if elem is limiting?

            ! Difference in stoichiometry of the prey, related to the copepods'
      !      Q_diff(i) = (elem_Q(i) -  zoo_fixed_stoichiometry(i) + respiration_rate) / zoo_fixed_stoichiometry(i)

      !      if (Q_diff(i) < 0 ) then ! Element is limiting
      !         C_excess(i) = -Q_diff(i) * prey ! C excess in relation to Quota_X
      !      endif
            ! @ Will need to excrete excess nutrient, if one is limiting

      !      elem_diff_to_limiting = Q_diff - min( Q_diff )
      !      elem_excess = elem_diff_to_limiting * prey * elem_Q ! Proportion of elements in excess, in relation to the most limiting

            !sloppy_feeding(i) = ingestion_rate*biomass * (1 - self%sloppy) * elem_Q(i)
      !   endif
      !end do

      ! Resolving the stoechiometry in zooplankton
      do i = 1,NUM_ELEM
         elem = ElementList(i:i)
         if (elem .NE. 'C') then

            ! Relative difference in the stoichiometry of the prey
            Q_diff(i) = (elem_Q(i) -  zoo_fixed_stoichiometry(i) + respiration_rate) / zoo_fixed_stoichiometry(i)

            if (Q_diff(i) < 0 ) then ! Element is limiting
               C_excess(i) = -Q_diff(i) ! C excess in the prey, in %
               elem_assimilation(i) = 1._rk ! Element assimilation, in %

            else ! The element is in excess
               C_excess(i) = 0._rk                             ! Only a part of the element is excreeted to be in balamce with C
               elem_assimilation(i) = 1._rk / ( Q_diff(i) + 1._rk ) ! % of element assimilated, the rest is excreeted

            endif

         endif
      end do

      if ( maxval(C_excess) > 0 ) then ! One of the elements is limiting, balancing with the other nutrients
         elem_assimilation = elem_assimilation * (1._rk - maxval(C_excess) + C_excess) ! % Difference in stoichiometry to the limiting element
      endif

      ! Growth rate of the copepod
      growth_rate = ingestion_rate * ( 1._rk - maxval(C_excess) ) ! Carbon that is not excreted is used

      _ADD_SOURCE_(self%id_prey, -ingestion_rate*biomass * days_per_sec ) ! Prey is ingested
      _ADD_SOURCE_(self%id_biomass, (growth_rate - respiration_rate) * biomass * days_per_sec ) ! Zooplankton growth

      ! Excretion and ingestion of elements
      excretion_rate = 1._rk - assimilation_rate  ! What is not assimilated is excreted
      excretion_rate(1) = maxval(C_excess)       ! Carbon excretion

      !! Zooplankton can stor C in fat, should not excrete it directly, maybe as a state variable?
      ! As an state variable, controlled by a boolean (IF_STORAGE T/F)

      ! Resolving the DOM and POM targets with biogeochemistry
      _SET_DIAGNOSTIC_(self%id_dummy, exudation_rate(2))
      _SET_DIAGNOSTIC_(self%id_dummy_sloppy, sloppy_feeding(2))

      do i = 1, NUM_ELEM
         elem = ElementList(i:i)
         if (elem .NE. 'C') then
            _ADD_SOURCE_( self%id_dom_(i), excretion_rate(i) * ingestion_rate * prey * elem_Q(i) ) !* days_per_sec DOM target for later
            _ADD_SOURCE_( self%id_det_(i), sloppy_feeding(i) ) ! POM target for later

         !else
         !   _ADD_SOURCE_( self%id_dom_(i), ( excretion_rate(i) * ingestion_rate + respiration_rate ) * prey  ) ! DOM target for later
         endif                                                                   ! DIC
      end do

      !! exudation_rate to DOM (proportional to C-respiration_rate)
      !new = respiration_rate * biomass
      !do i = 1,NUM_ELEM
         !if (ElementList(i:i) .NE. 'C') _ADD_SOURCE_(self%id_var(dom_index(i)), new*stoichiometry(i) UNIT) ! elem_Qs sink
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
