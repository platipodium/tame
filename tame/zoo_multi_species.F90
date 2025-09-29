! SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
! SPDX-FileContributor Thomas Imbert <thomas.imbert@hereon.de>
! SPDX-FileContributor Ovidio Garcia <ovidio.garcia@hereon.de>
! SPDX-FileContributor Carsten Lemmen <carsten.lemmen@hereon.de>
! SPDX-License-Identifier: Apache-2.0

#include "fabm_driver.h"

!! TAME zooplankton--any generic predator--(fixed stoichiometry) module
module tame_zooplanktonmulti
   use fabm_types
   use tame_types
   use tame_functions
   implicit none

   public NUM_ZOO
   integer, parameter :: NUM_ZOO = 1 !! Number of Zoo PFTs modelled
   character(len=3) :: register_name

   private
   type, extends(type_base_model), public :: type_tame_zooplanktonmulti

      !! State variables ids
      type (type_state_variable_id) :: id_zooplankton_C(NUM_ZOO) ! zooplankton_C

      !! Dependency ids
      type (type_dependency_id) :: id_prey_Q(NUM_ELEM) ! Prey and stoichiometry
      type(type_dependency_id) :: id_temp ! for temperature dependence
      type (type_state_variable_id) :: id_prey, id_dom_(NUM_ELEM), id_det_(NUM_ELEM) ! prey

      !! Diagnostic variable ids
      type (type_diagnostic_variable_id) :: id_dummye, id_dummya, id_dummy_N, id_dummy_P
      type (type_diagnostic_variable_id) :: id_zoo_elem(NUM_ELEM) ! Zooplankton stoichiometry

      !! Model parameters
      !real(rk) :: max_ingestion, saturation, sloppy ! maximum ingestion_rate rate, food saturation
      !real(rk) :: resp, Q10 ! respiration_rate rate

      real(rk) :: max_ingestion(NUM_ZOO), saturation(NUM_ZOO)!, sloppy(NUM_ZOO) 
      real(rk) :: resp(NUM_ZOO), Q10(NUM_ZOO) 

      contains
      procedure :: initialize
      procedure :: do
   end type

   integer :: det_index(NUM_ELEM), dom_index(NUM_ELEM) ! TODO: global setting!

   contains
!! ADD ELEM TO AGGREGATE VARIABLE (like in total_X in bgc)
   subroutine initialize(self, configunit)
      class (type_tame_zooplanktonmulti), intent(inout), target :: self
      integer, intent(in) :: configunit
      integer :: i ! dummy index
      integer, parameter :: counter = 0 !@what
      character :: elem

      !! Store parameter values in our own derived type
      !! all rates must be provided in values per day and are converted here to values per second.
      !! @todo: add min and max values

      !! Register parameters for the different PFTs
      do i = 1,NUM_ZOO
         write(register_name,'(i3)') i
         register_name = adjustl(register_name)

         call self%get_parameter(self%saturation(i), 'saturation_'//trim(register_name), 'mmol-C m-3', 'grazing saturation', default=2.5_rk)
         call self%get_parameter(self%max_ingestion(i), 'max_ingestion_'//trim(register_name), 'd-1', 'maximum ingestion rate', default=2.5_rk)
         call self%get_parameter(self%resp(i), 'resp_'//trim(register_name), 'd-1', 'respiration rate', default=0._rk)
         !call self%get_parameter(self%sloppy, 'sloppy', '', 'sloppy feeding', default=0._rk)
         call self%get_parameter(self%Q10(i), 'Q10_'//trim(register_name), '--', 'temperature sensitivity (ref temperature 20 degC)', default=2.2_rk)

         !! Register state variables
         call self%register_state_variable(self%id_zooplankton_C(i),'zooplankton_C_'//trim(register_name), 'mmol-C m-3', 'concentration', 1.0_rk, minimum=0.0_rk)
      end do

      !! Register environmental dependencies
      call self%register_dependency(self%id_temp, standard_variables%temperature)

      !! Register external dependencies
      call self%register_state_dependency(self%id_prey, 'prey','mmol-C m-3', 'prey source')
      do i = 1,NUM_ELEM
         elem = ElementList(i:i)
                     
         call self%register_state_dependency(self%id_det_(i), 'det_' // elem,'mol-' // elem // ' mol-C-1', 'Particulate Organic' // elem)

         if (elem .NE. 'C') then
            call self%register_state_dependency(self%id_dom_(i), 'dom_' // elem,'mol-' // elem // ' mol-C-1', 'Dissolved Organic' // elem)

            ! Element stoichiometry in Zooplankton
            call self%register_diagnostic_variable(self%id_zoo_elem(i), 'zoo_' // elem,'mol-' // elem // ' m-3', 'zooplankton ' // elem)
            call self%add_to_aggregate_variable(type_universal_standard_variable(name='total_' // trim(ElementName(i)), units='mmol-' // elem // ' m-3', &
              aggregate_variable=.true., conserved=.true.), self%id_zoo_elem(i))

      !! Retrieve prey stoichiometric composition (if FlexStoich)
            call self%register_dependency(self%id_prey_Q(i), 'PreyQ_' // elem,'mol-' // elem // ' mol-C-1', elem // ':C-quota')

         endif
      end do ! Prey is supposed to be in carbon, so only extracting other elements, if available

      call self%register_diagnostic_variable(self%id_dummye, 'dummye','', '')
      call self%register_diagnostic_variable(self%id_dummya, 'dummya','', '')
      
      call self%register_diagnostic_variable(self%id_dummy_N, 'dummy_N','', '')
      call self%register_diagnostic_variable(self%id_dummy_P, 'dummy_P','', '')

   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_) 
      class (type_tame_zooplanktonmulti), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_
      integer :: i,j ! dummy index
      real(rk) :: zooplankton_C, prey, ivlev, C_consumption
      real(rk) :: ingestion_rate, respiration_rate, growth_rate

      !! Variables to store intermediary results of the model for the stoichiometry
      real(rk) :: exudation_rate(NUM_ELEM), elem_Q(NUM_ELEM), Q_diff(NUM_ELEM) ! elem_Q is the prey quotas
      real(rk) :: C_excess(NUM_ELEM), elem_assimilation(NUM_ELEM)              ! assimilation or excess of nutrients
      real(rk) :: zoo_stoich_scale!(NUM_ELEM)                                  ! adjust the zooplankton stoichiometry when respiration > 0
      real(rk) :: sloppy_feeding(NUM_ELEM), excretion_rate(NUM_ELEM)
      real(rk) :: temp, temp_factor
      character :: elem

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

      !! Retrieve current environmental conditions, independant of PFTs
      _GET_(self%id_prey,prey)
      _GET_(self%id_temp, temp) ! temperature

      !! Stoichiometry composition of prey
      !@what: is zooplankton exudating to elem_Q pool or to DOM pool?
      do i = 1, NUM_ELEM
         elem = ElementList(i:i)
         if (elem .NE. 'C') then
            _GET_( self%id_prey_Q(i), elem_Q(i) ) ! elem_Q target for later
         endif
      end do
      ! Excretion to ammonia NH4+?

      do i = 1,NUM_ZOO

        _GET_(self%id_zooplankton_C(i), zooplankton_C) ! zooplankton_C carbon

        temp_factor = self%Q10(i)**(0.1_rk*(temp-20.0_rk))

        !! Predation
        ivlev = 1.0_rk - exp(- prey / self%saturation(i) )
        ingestion_rate =  temp_factor * self%max_ingestion(i) * ivlev

        !! Losses
        respiration_rate = temp_factor * self%resp(i) !@todo: include other physiological loss terms here

        !! Resolving the stoechiometry in zooplankton
        C_consumption = respiration_rate / ingestion_rate

        !! Resolve the stoichiometric balance between the zooplankton and the prey
        do j = 1,NUM_ELEM
            elem = ElementList(i:i)
            if (elem .NE. 'C') then

                ! Relative difference in the stoichiometry of the prey
                !Q_diff(i) = (elem_Q(i) -  zoo_fixed_stoichiometry(i) / ( 1._rk + respiration_rate ) ) / zoo_fixed_stoichiometry(i)
                
                !Q_diff(i) = (elem_Q(i) -  zoo_fixed_stoichiometry(i)) / zoo_fixed_stoichiometry(i) + respiration_rate / ingestion_rate

                !zoo_stoich_scale = zoo_fixed_stoichiometry(i) / (1._rk + respiration_rate / ingestion_rate)
                !Q_diff(i) = (elem_Q(i) -  zoo_stoich_scale) / zoo_stoich_scale !+ respiration_rate / ingestion_rate
                
                if ( C_consumption .GE. 1._rk ) then ! Simplification if respiration exceeds grazing
                C_excess(j) = 0._rk
                elem_assimilation(j) = 0._rk ! 0% assimilation as all the C ingested is used for maintenance

                else ! Simplification if respiration is not exceeding

                zoo_stoich_scale = elem_Q(j) / (1._rk - respiration_rate / ingestion_rate) ! R / I of the C will be consumed for maintenance costs
                Q_diff(j) = (zoo_stoich_scale - zoo_fixed_stoichiometry(j)) / zoo_fixed_stoichiometry(j) !+ respiration_rate / ingestion_rate

                if (Q_diff(j) < 0._rk ) then        ! Element is limiting
                    C_excess(j) = -Q_diff(j)     ! C excess in the prey, in %
                    elem_assimilation(j) = 1._rk ! Element assimilation, in %

                else ! The element is in excess
                    C_excess(j) = 0._rk                                  ! Only a part of the element is excreted to be in balance with C
                    elem_assimilation(j) = 1._rk / ( Q_diff(j) + 1._rk ) ! % of element assimilated, the rest is excreted

                endif

                endif 

            endif
        end do

        if ( maxval(C_excess) > 0._rk ) then ! One of the elements is limiting, balancing with the other nutrients
            elem_assimilation = elem_assimilation * (1._rk - maxval(C_excess) + C_excess) ! % Difference in stoichiometry to the limiting element
        endif

        ! Growth rate of the copepod
        growth_rate = ingestion_rate * ( 1._rk - maxval(C_excess) ) ! Carbon that is not excreted is used

        _ADD_SOURCE_(self%id_prey, -ingestion_rate * zooplankton_C * days_per_sec ) ! Prey is ingested
        _ADD_SOURCE_(self%id_zooplankton_C(i), (growth_rate - respiration_rate) * zooplankton_C * days_per_sec ) ! Zooplankton growth

        ! Excretion and ingestion of elements
        excretion_rate = 1._rk - elem_assimilation  ! What is not assimilated is excreted
        excretion_rate(1) = maxval(C_excess)        ! Carbon excretion

        !! Zooplankton can store C in fat, should not excrete it directly, maybe as a state variable?
        ! As a state variable, controlled by a boolean (IF_STORAGE T/F)

        ! Resolving the DOM and POM targets with biogeochemistry
      
        do j = 1, NUM_ELEM
            elem = ElementList(i:i)

            if (growth_rate - respiration_rate < 0._rk) then ! Zooplankton starves and dies
                _ADD_SOURCE_( self%id_det_(j), -(growth_rate - respiration_rate) * zooplankton_C * zoo_fixed_stoichiometry(j) * days_per_sec ) ! Fecal pellets: POC 
            endif

            if (elem .NE. 'C') then
                _ADD_SOURCE_( self%id_dom_(j), excretion_rate(j) * zooplankton_C * ingestion_rate * elem_Q(j) * days_per_sec ) ! DOM target for later
                
                !_ADD_SOURCE_( self%id_det_(i), sloppy_feeding(i) ) ! POM target for later
                _SET_DIAGNOSTIC_(self%id_zoo_elem(j), zoo_fixed_stoichiometry(j) * zooplankton_C ) 

            else
                _ADD_SOURCE_( self%id_det_(j), ingestion_rate * maxval(C_excess) * zooplankton_C * days_per_sec ) ! Fecal pellets: POC 
            endif                                                             ! DIC ?

        end do


      !! Retrieve diagnostics 
      _SET_DIAGNOSTIC_(self%id_dummya, 1._rk - maxval(C_excess) )! maxval(C_excess) )  
      _SET_DIAGNOSTIC_(self%id_dummye, maxval(elem_assimilation) )
      !_SET_DIAGNOSTIC_(self%id_dummy_N, ingestion_rate * elem_Q(2) * elem_assimilation(2) / ( (growth_rate - respiration_rate) * zoo_fixed_stoichiometry(2) ) )
      _SET_DIAGNOSTIC_(self%id_dummy_N, respiration_rate )
      _SET_DIAGNOSTIC_(self%id_dummy_P, ingestion_rate * elem_Q(3) * elem_assimilation(3) / ( (growth_rate - respiration_rate) * zoo_fixed_stoichiometry(3) ) )

      ! sinking to POM
      !new = sinking * zooplankton_C

      !! @what: should not be sinking be managed by physical driver? I see this useful probably only on 0d setups. A flag is required!
      !do i = 1,NUM_ELEM ! C, N, P (Si, Fe)
      !   _ADD_SOURCE_(self%id_var(det_index(i)), new*stoichiometry(i) UNIT) !
      !end do
      
      end do

      ! Leave spatial loops (if any)
      _LOOP_END_
   end subroutine do

end module tame_zooplanktonmulti
