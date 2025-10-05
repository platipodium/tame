! SPDX-License-Identifier: CC0-1.0
! SPDX-License-Identifier: Apache-2.0
! SPDX-FileCopyrightText: 2024-2025 Helmholtz-Zentrum hereon GmbH
! SPDX-FileContributor Kai Wirtz <kai.wirtz@hereon.de>
! SPDX-FileContributor Thomas Imbert <thomas.imbert@hereon.de>
! SPDX-FileContributor Ovidio Garcia <ovidio.garcia@hereon.de>
! SPDX-FileContributor Carsten Lemmen <carsten.lemmen@hereon.de>
! TAME life module: phytoplankton with flexible stoichiometry and zooplankton

#include "fabm_driver.h"

module tame_life

use fabm_types
use fabm_expressions
use tame_types
use tame_functions
use tame_stoich_functions
   implicit none
   private
   
   type, extends(type_base_model), public :: type_tame_life
      ! ========== Phytoplankton Variable identifiers ==========
      type (type_state_variable_id)      :: id_phytoplankton_C     ! Phytoplankton biomass
      type (type_state_variable_id)      :: id_var(NUM_CHEM+2*NUM_ELEM) ! DOM & POM
      type (type_dependency_id)          :: id_par, id_temp, id_Phy_X_old(NUM_ELEM)
      type (type_diagnostic_variable_id) :: id_nut, id_nut2, id_din, id_rate, id_day_of_year
      type (type_diagnostic_variable_id) :: id_Q(NUM_ELEM), id_dQ_dt(NUM_ELEM), id_phy_elem(NUM_ELEM)
      type (type_global_dependency_id)   :: id_doy
      
      ! ========== Zooplankton Variable identifiers ==========
      type (type_state_variable_id)      :: id_zooplankton_C ! zooplankton_C
      type (type_diagnostic_variable_id) :: id_dummye, id_dummya, id_dummy_N, id_dummy_P
      type (type_diagnostic_variable_id) :: id_zoo_elem(NUM_ELEM) ! Zooplankton stoichiometry
      
      ! ========== Shared diagnostic for information exchange ==========
      type (type_diagnostic_variable_id) :: id_rhs_phyC ! grazing rate (calculated in zoo, used in phy)
      
      ! ========== Phytoplankton Model parameters ==========
      real(rk) :: num_dum        ! generic dummy parameter to control numerical settings
      real(rk) :: qmort          ! density dep mortality
      real(rk) :: rmax           ! Growth metabolism parameters
      real(rk) :: gamma          ! Light exploitation
      real(rk) :: s0             ! Sinking
      real(rk) :: resp           ! Respiration parameters
      real(rk) :: K_P, K_N       ! Half-saturation constants
      real(rk) :: nut_limitation(NUM_NUTRIENT) ! Vector of limitation degree
      real(rk) :: HalfSatNut(NUM_NUTRIENT) ! Vector of half-saturations
      logical  :: FlexStoich
      
      ! ========== Zooplankton Model parameters ==========
      real(rk) :: max_ingestion, saturation, sloppy ! maximum ingestion rate, food saturation
      real(rk) :: resp_zoo, Q10 ! zooplankton respiration rate and Q10

   contains
      procedure :: initialize
      procedure :: do
   end type

   integer :: det_index(NUM_ELEM), dom_index(NUM_ELEM) ! TODO: global setting!
   integer :: num_chem_of_nut(NUM_NUTRIENT), share_nut_chemindex(NUM_NUTRIENT,NUM_CHEM)

contains

   subroutine initialize(self, configunit)
      class (type_tame_life), intent(inout), target :: self
      integer, intent(in) :: configunit
      integer :: i, j ! Indice dummy
      integer, parameter :: counter = 0
      character :: elem

      ! ========== Phytoplankton Parameters ==========
      call self%get_parameter(self%rmax, 'rmax', 'd-1', 'maximum production rate', default=2.5_rk)
      call self%get_parameter(self%gamma, 'gamma', 'microE-1 m-2', 'light absorption scaling', default=1.0_rk)
      call self%get_parameter(self%num_dum, 'num_dum', '-', 'dummy parameter', default=1._rk)
      call self%get_parameter(self%qmort, 'qmort', 'mmol-C-1 m3 d-1', 'density dep mortality', default=0.0_rk)
      call self%get_parameter(self%s0, 's0', 'd-1', 'default sinking rate', default=0._rk)
      call self%get_parameter(self%K_P, 'K_P', 'mmol m-3', 'P half-saturation', default=0.4_rk)
      call self%get_parameter(self%K_N, 'K_N', 'mmol m-3', 'N half-saturation', default=4.0_rk)
      call self%get_parameter(self%resp, 'resp', 'mmol', 'carbon cost per nitrogen uptake', default=0.2_rk)
      call self%get_parameter(self%FlexStoich, 'FlexStoich', '', 'Is FlexStoich?', default=.true.)
      
      ! ========== Zooplankton Parameters ==========
      call self%get_parameter(self%saturation, 'saturation', 'mmol-C m-3', 'grazing saturation', default=2.5_rk)
      call self%get_parameter(self%max_ingestion, 'max_ingestion', 'd-1', 'maximum ingestion rate', default=2.5_rk)
      call self%get_parameter(self%resp_zoo, 'resp_zoo', 'd-1', 'zooplankton respiration rate', default=0._rk)
      call self%get_parameter(self%sloppy, 'sloppy', '', 'sloppy feeding', default=0._rk)
      call self%get_parameter(self%Q10, 'Q10', '--', 'temperature sensitivity (ref temperature 20 degC)', default=2.2_rk)
      
      dphyXdt_crit = self%num_dum * fixed_stoichiometry
      
      self%HalfSatNut(1) = self%K_N
      self%HalfSatNut(2) = self%K_P

      ! ========== Register state variables ==========
      call self%register_state_variable(self%id_phytoplankton_C, 'phytoplankton_C', 'mmol m-3', 'concentration', 1.0_rk, minimum=0.0_rk)
      call self%register_state_variable(self%id_zooplankton_C, 'zooplankton_C', 'mmol-C m-3', 'concentration', 1.0_rk, minimum=0.0_rk)

      ! ========== Register environmental dependencies (shared) ==========
      call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
      call self%register_dependency(self%id_temp, standard_variables%temperature)
      call self%register_dependency(self%id_doy, standard_variables%number_of_days_since_start_of_the_year)

      ! ========== Register chemical dependencies ==========
      do i = 1, NUM_CHEM
         call self%register_state_dependency(self%id_var(i), chemicals(i), 'mmol m-3', chemicals(i))
      end do

      ! ========== Register OM variables for each element ==========
      do i = 1, NUM_ELEM ! e.g., C, N, P (Si, Fe)
         det_index(i) = NUM_CHEM + 2*i - 1
         dom_index(i) = NUM_CHEM + 2*i
         elem = ElementList(i:i)
         
         call self%register_state_dependency(self%id_var(det_index(i)), 'det_' // elem, 'mmol-' // elem // ' m-3', 'Detritus ' // trim(ElementName(i)))
         call self%register_state_dependency(self%id_var(dom_index(i)), 'dom_' // elem, 'mmol-' // elem // ' m-3', 'Dissolved Organic ' // trim(ElementName(i)))
         call self%register_dependency(self%id_Phy_X_old(i), 'old_Phy_' // elem, 'mol-' // elem // ' mol-C-1', 'previous phy' // elem) 

         if (elem .NE. 'C') then
            ! Phytoplankton element diagnostics
            call self%register_diagnostic_variable(self%id_Q(i), 'Q_' // elem, 'mol-' // elem // ' mol-C-1', elem // ':C-quota')
            call self%register_diagnostic_variable(self%id_dQ_dt(i), 'dQ_dt_' // elem, 'mol-' // elem // ' mol-C-1 d-1', 'change in ' // elem // ':C-quota')
            call self%register_diagnostic_variable(self%id_phy_elem(i), 'phy_' // elem, 'mol-' // elem // ' m-3', 'phytoplankton ' // elem)
            call self%add_to_aggregate_variable(type_universal_standard_variable(name='total_' // trim(ElementName(i)), units='mmol-' // elem // ' m-3', &
              aggregate_variable=.true., conserved=.true.), self%id_phy_elem(i))
            
            ! Zooplankton element diagnostics
            call self%register_diagnostic_variable(self%id_zoo_elem(i), 'zoo_' // elem, 'mol-' // elem // ' m-3', 'zooplankton ' // elem)
            call self%add_to_aggregate_variable(type_universal_standard_variable(name='total_' // trim(ElementName(i)), units='mmol-' // elem // ' m-3', &
              aggregate_variable=.true., conserved=.true.), self%id_zoo_elem(i))

         endif
      end do
      ! ==== partitioning of chemical-nutrient uptake by phy  ======
      num_chem_of_nut = 0 
      do i = 1, NUM_CHEM
         j = chem2nut(i)
         num_chem_of_nut(j) = num_chem_of_nut(j) + 1 
         share_nut_chemindex(j,num_chem_of_nut(j)) = i 
      end do
    
      ! ========== Shared diagnostic for phy-zoo coupling ==========
      call self%register_diagnostic_variable(self%id_rhs_phyC, 'rhs_phyC', 'mmol-C m-3 d-1', 'phy-C grazing rate')
      
      ! ========== Other diagnostic variables ==========
      call self%register_diagnostic_variable(self%id_din, 'din', 'mmol-N m-3', 'DIN', output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_nut, 'nut1', 'mmol-? m-3', 'nutrient related', output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_nut2, 'nut2', 'mmol-? m-3', 'nutrient related', output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_rate, 'rate', 'mmol-? m-3', 'dummy')
      call self%register_diagnostic_variable(self%id_day_of_year, 'day_of_year', 'd', 'day_of_year')
      call self%register_diagnostic_variable(self%id_dummye, 'dummye', '', '')
      call self%register_diagnostic_variable(self%id_dummya, 'dummya', '', '')
      call self%register_diagnostic_variable(self%id_dummy_N, 'dummy_N', '', '')
      call self%register_diagnostic_variable(self%id_dummy_P, 'dummy_P', '', '')

   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_tame_life), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) ::  par, temp,  dphy, dphyt, rdphy, dphyXdt, dtime=120., rhs_phy, doy, doy0
      real(rk) :: production, sinking, mort, chem_change, loss, exud, nut_lim, nut_lim_tot
      real(rk) :: exudation(NUM_NUTRIENT), nutrient(NUM_NUTRIENT), nut_change(NUM_NUTRIENT)
      real(rk) :: quota(NUM_ELEM), phy_X_old(NUM_ELEM), phy_X_change(NUM_ELEM), dix_chemical(NUM_CHEM), part(NUM_CHEM)
      logical  :: ncrit(NUM_CHEM)
      integer  :: i, j, ni, ie ! Indices
      real(rk) :: phytoplankton_C, zooplankton_C, prey, ext_rhs_phyC, ivlev, temp_factor
      real(rk) :: ingestion, feeding, respiration, rhs_zoo, din, phy_X, func, resp_hetero,sloppy_feed, sum_part
      real(rk) :: prey_Q(NUM_ELEM), Excess_C_upt(NUM_ELEM),excretion(NUM_ELEM)

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_
   ! ===============================================================
   ! get environmental and own states
   ! ==============================================================
     
      _GET_GLOBAL_(self%id_doy, doy)       ! day of year
      _SET_DIAGNOSTIC_(self%id_day_of_year, doy)

      ! Retrieve current (local) state variable values.
      _GET_(self%id_phytoplankton_C, phytoplankton_C)     ! phytoplankton carbon
      _GET_(self%id_zooplankton_C, zooplankton_C) ! zooplankton_C carbon
      !! Retrieve current chemical conditions
      do i = 1, NUM_CHEM ! e.g., CO2, NO3, NH4 (PO4)
         _GET_(self%id_var(i), dix_chemical(i))  ! Dissolved Inorganic Nutrient DIX in mmol-X/m**3
      end do
      ! Retrieve current environmental conditions.
      _GET_(self%id_par, par)          ! local photosynthetically active radiation
      _GET_(self%id_temp, temp)        ! water temperature
      _GET_(self%id_Phy_X_old(1), doy0)

      temp_factor = self%Q10**(0.1_rk*(temp - 20.0_rk))

   ! ===============================================================
   ! nutrient limitation of autotrophs
   ! ==============================================================

      ! TODO replace by TransIndex_DOMDIX, TransIndex2_DOMDIX, which should be set globally
      ! calculate nutrients from chemicals (e.g., DIC=CO2+HCO3, DIN=NO3+NH4)
      nutrient = 0.0_rk
      do i = 1, NUM_CHEM ! e.g., CO2, NO3, NH4 (PO4)
         nutrient(chem2nut(i)) = nutrient(chem2nut(i)) + dix_chemical(i)
      end do
      _SET_DIAGNOSTIC_(self%id_din, nutrient(1))

      ! effect of nutrient co-limitation
      !   here: depends on ambient concentration, not on internal quota! 
      ! TODO: check and re-formulate sum-rule!
      nut_lim_tot = 0._rk
      do i = 1, NUM_NUTRIENT
         nut_lim = nutrient(i)/(self%HalfSatNut(i) + nutrient(i))
         nut_lim_tot = nut_lim_tot + 1.0_rk/(nut_lim + small)  ! sum rule for co-limitation
      end do
      nut_lim_tot = 1.0_rk/nut_lim_tot
      _SET_DIAGNOSTIC_(self%id_nut, nut_lim_tot*100.0_rk)

   ! ===============================================================
   ! net growth rate (photosynthesis) of autotrophs - without grazing
   ! ==============================================================
      ! Production
      func = 1.0_rk - exp(-self%gamma * par / self%rmax)
      production = self%rmax * func * nut_lim_tot
      !_SET_DIAGNOSTIC_(self%id_rate, production)

      ! All losses: respiration + sinking
      respiration = self%resp                     ! C loss for DIN-uptake
      ! sinking     = self%s0                       ! TODO - move to FABM's set_settling
      mort        = self%qmort * phytoplankton_C  ! density-dependent mortality (e.g., virus)

      ! combine to get the temporal derivative for phytoplankton C
      rhs_phy = (production  - respiration - mort) * phytoplankton_C

   ! ===============================================================
   !  internal element ratios
   !              as response function of T, PAR, ambient nutrients 
   ! ==============================================================
      ! Set quota either as flexible or constant (Redfield)
      if (self%FlexStoich) then
         ! Flexible regulation of non-Redfield stoichiometry (C:N:P)
         do i = 1, NUM_NUTRIENT           
            j  = nut2othernut(i) ! index of complementary, co-limiting nutrient
            ie = nut2elem(i)
            quota(ie) = calc_quota(nutrient(i), nutrient(j), par, temp, i, j)
         end do

         if (doy0 .ge. 0 .AND. doy0 .lt. 367) then 
            do i = 1, NUM_ELEM
               if (ElementList(i:i) .NE. 'C') then
                  _GET_(self%id_Phy_X_old(i), phy_X_old(i))
                  phy_X = phytoplankton_C*quota(i)
                  dtime = (doy - doy0 + 0.001_rk*days_per_sec)
                  dphyXdt = (phy_X - phy_X_old(i))/dtime
                  rdphy = dphyXdt/dphyXdt_crit(i)
                  if (abs(rdphy) .gt. 0.0001_rk) then
                     phy_X_change(i) = dphyXdt_crit(i) * (-1._rk + 2._rk/(1._rk + exp(-2*rdphy)))
                     quota(i) = (phy_X_old(i) + phy_X_change(i)*dtime) / (phytoplankton_C + small)
                  else
                     phy_X_change(i) = dphyXdt
                  endif
               else
                  quota(i) = fixed_stoichiometry(i)
               endif
            end do
         else ! initial period with memory
            phy_X_change = 0.0_rk
         endif
      else
         do i = 1, NUM_ELEM
            quota(i) = fixed_stoichiometry(i)
            phy_X_change(i) = 0.0_rk
         end do
      endif

   ! ===============================================================
   !  grazing by herbivores
   ! ==============================================================
      prey   = phytoplankton_C ! Using phytoplankton as prey
      prey_Q = quota

   ! functional response
      ivlev     = 1.0_rk - exp(-prey / self%saturation)
      feeding   = temp_factor * self%max_ingestion * ivlev
      ingestion = (1.0_rk - self%sloppy) * feeding 

      !! Respiratory losses of heterotrophs
      resp_hetero = temp_factor * self%resp_zoo
      
      ! Store grazing rate for phytoplankton stoich. regulation
      ext_rhs_phyC = -feeding * zooplankton_C
      ! Also store as diagnostic for output
      _SET_DIAGNOSTIC_(self%id_rhs_phyC, ext_rhs_phyC)

      ! virtual C-uptake to balance stoichiometry changes by predation
      Excess_C_upt(1) = -resp_hetero
      do i = 2, NUM_ELEM
         Excess_C_upt(i) = (1._rk - prey_Q(i)/zoo_stoichiometry(i)) * ingestion - resp_hetero 
      end do
      ! the maximal imbalance determines C-exudation, if any
      excretion(1) = max(0.0, maxval(Excess_C_upt)) 
      ! all other exudation rates of chemicals are rescaled accordingly
      !  .. and should be negative  
      Excess_C_upt = Excess_C_upt - excretion(1)
      excretion(2:NUM_ELEM) = zoo_stoichiometry(2:NUM_ELEM) * (-Excess_C_upt(2:NUM_ELEM))
      ! C-based growth rate of grazer 
      rhs_zoo = (ingestion - resp_hetero - excretion(1)) * zooplankton_C 

   ! =============================================================
   !  feedback of phytoplankton nutrient uptake
   ! =============================================================
      ! avoiding too strong draw-down in a rare element, if another element  
      !   of the same nutrient is relatively high
      ! e.g., NO3+NH4 partitioning in DIN uptake 
      part = 1.0_rk
      do j = 1, NUM_NUTRIENT
         ! more than one chemical per nutrient -> partitioning
         if (num_chem_of_nut(j) > 1) then
            sum_part = 0._rk
            do ni = 1, num_chem_of_nut(j)
              i = share_nut_chemindex(j,ni)
              ! reduce partitioning at low concentration
              part(i) = 1.0_rk - exp(-dix_chemical(i)/nut_minval(j))
              sum_part = sum_part + part(i)
            end do
            ! re-normalize partitioning coefficients
            do ni = 1, num_chem_of_nut(j)
              i = share_nut_chemindex(j,ni)
              part(i) = part(i)/(sum_part+ small)
            end do
         endif
      end do
      !  estimate nutrient demand by quota changes based on old phy_X
      do i = 1, NUM_CHEM
         j = chem2elem(i)
         ! detected change in phyC minus known changes
         dphyt = phy_X_change(j) - (rhs_phy + ext_rhs_phyC) * quota(j)
         ! total change: uptake by growth assimilation and quota change
         chem_change = -part(i)*(production * quota(j)* phytoplankton_C + dphyt) 
         _ADD_SOURCE_(self%id_var(i), chem_change * days_per_sec) 
      end do

      _SET_DIAGNOSTIC_(self%id_phy_elem(1), doy)
      do i = 2, NUM_ELEM
         _SET_DIAGNOSTIC_(self%id_Q(i), quota(i))
         _SET_DIAGNOSTIC_(self%id_dQ_dt(i), phy_X_change(i))
         _SET_DIAGNOSTIC_(self%id_phy_elem(i), quota(i) * phytoplankton_C)
      end do

      ! C-based growth rate of grazer 
      _ADD_SOURCE_(self%id_zooplankton_C, rhs_zoo * days_per_sec)
     ! Final change in phy C 
      _ADD_SOURCE_(self%id_phytoplankton_C, (rhs_phy + ext_rhs_phyC) * days_per_sec)
   
      ! Resolving the DOM and POM targets with biogeochemistry
      ! C-respiration (to DIC) TODO: add DIC and 
      !     corresponding element flux to DOM
      exud = respiration * phytoplankton_C
      do i = 1, NUM_ELEM
         _ADD_SOURCE_(self%id_var(dom_index(i)), excretion(i) * zooplankton_C * days_per_sec)
         ! Mortality and sloppy feeding part to detrital POM
         sloppy_feed = self%sloppy * feeding  * zooplankton_C
         loss = mort *  phytoplankton_C 
         _ADD_SOURCE_(self%id_var(det_index(i)), (loss + sloppy_feed) * quota(i) *days_per_sec) !

         if (i .gt. 1) then
            ! Exudation to DOM (proportional to C-respiration)
            _ADD_SOURCE_(self%id_var(dom_index(i)), exud*quota(i) * days_per_sec)

            _SET_DIAGNOSTIC_(self%id_zoo_elem(i), zoo_stoichiometry(i) * zooplankton_C)
         endif 
      end do
      _SET_DIAGNOSTIC_(self%id_nut2, Excess_C_upt(3))
      _SET_DIAGNOSTIC_(self%id_rate, Excess_C_upt(1))
      _SET_DIAGNOSTIC_(self%id_dummya, Excess_C_upt(2))
      _SET_DIAGNOSTIC_(self%id_dummye, excretion(1))
!      _SET_DIAGNOSTIC_(self%id_dummy_N, excretion(2))
!      _SET_DIAGNOSTIC_(self%id_dummy_P, excretion(3))
      _SET_DIAGNOSTIC_(self%id_dummy_N, resp_hetero)
      _SET_DIAGNOSTIC_(self%id_dummy_P, rhs_zoo)
      
   !rhs_zoo = (ingestion - resp_hetero - excretion(1)
  ! Leave spatial loops (if any)
      _LOOP_END_
   end subroutine do
   ! ============================================================================
   ! MODEL FUNCTIONS
   ! ============================================================================
   elemental real(rk) function light_absorb(rmax, gamma, Ik)
      real(rk), intent(in) :: rmax, gamma, Ik
      light_absorb = 1 - exp(-gamma * Ik / rmax)
   end function light_absorb

   elemental real(rk) function limitation(Nut, Aff)
      real(rk), intent(in) :: Aff, Nut
      limitation = Aff*Nut/(Aff*Nut + 1)
   end function limitation

end module tame_life