! SPDX-License-Identifier: CC0-1.0
! SPDX-FileCopyrightText: 2024-2025 Helmholtz-Zentrum hereon GmbH
! SPDX-FileContributor Kai Wirtz <kai.wirtz@hereon.de>
! SPDX-FileContributor Thomas Imbert <thomas.imbert@hereon.de>
! TAME phytoplankton module with flexible stoichiometry partly based on MAECS (Wirtz et al)
!

#include "fabm_driver.h"
module tame_phyto

use fabm_types
use fabm_expressions
use tame_types
use tame_functions
use tame_stoich_functions
   implicit none
   private
   !type, extends(type_base_model), public :: type_examples_npzd_phy
   type, extends(type_base_model), public :: type_tame_phyto
      ! Variable identifiers
      type (type_state_variable_id)      :: id_phytoplankton_C     ! Phytoplankton biomass
!      type (type_state_variable_id)      :: id_no3, id_nh4, id_po4     ! id_din Nutrients
      type (type_state_variable_id)      :: id_var(NUM_CHEM+2*NUM_ELEM) ! TODO : flexible num of DOM & POM
     ! type_interior_standard_variable(name='DIN',units='mmol-N m-3')
      type (type_dependency_id)   :: id_par, id_temp,id_nut_change(NUM_CHEM),id_Phy_X_old(NUM_ELEM),Phy_X(NUM_ELEM),id_rhs_phy_old
      type (type_diagnostic_variable_id) :: id_nut,id_nut2, id_din,id_rate,id_day_of_year, id_Q(NUM_ELEM),id_dQ_dt(NUM_ELEM),id_phy_elem(NUM_ELEM), id_rhs_phyC
      type (type_global_dependency_id)     :: id_doy
      !type (type_surface_dependency_id)  :: id_I_0   ! Surface irradiance
      ! Model parameters
      real(rk) :: num_dum        ! generic dummy parameter to control numerical settings
      real(rk) :: qmort          ! density dep mortality
      real(rk) :: rmax           ! Growth metabolism parameters
      real(rk) :: gamma          ! Light exploitation
      real(rk) :: s0             ! Sinking
      real(rk) :: resp           ! Respiration parameters
      real(rk) :: K_P, K_N       !
      real(rk) :: nut_limitation(NUM_NUTRIENT) ! Vector of limitation degree
      real(rk) :: HalfSatNut(NUM_NUTRIENT) ! Vector of half-saturations
      logical :: FlexStoich
   contains
      procedure :: initialize
      procedure :: do
   end type

   integer :: det_index(NUM_ELEM), dom_index(NUM_ELEM) ! TODO: global setting!

contains
   subroutine initialize(self, configunit)
      class (type_tame_phyto), intent(inout), target :: self
      integer,                        intent(in)     :: configunit
      integer :: i ! Indice dummy
      integer, parameter :: counter = 0
      character :: elem

      ! Store parameter values in our own derived type
      ! NB: all rates must be provided in values per day and are converted here to values per second.
      call self%get_parameter(self%rmax, 'rmax', 'd-1',         'maximum production rate',  default=2.5_rk) !, scale_factor=days_per_sec
      call self%get_parameter(self%gamma,'gamma','microE-1 m-2','light absorption scaling', default=1.0_rk)
      call self%get_parameter(self%num_dum,   'num_dum',   '-', 'dummy parameter',default=1._rk)
      call self%get_parameter(self%qmort,   'qmort',   'mmol-C-1 m3 d-1',    'density dep mortality',default=0.0_rk)
      call self%get_parameter(self%s0,   's0',   'd-1',         'default sinking rate',     default=0._rk) !, scale_factor=days_per_sec
      call self%get_parameter(self%K_P,  'K_P',  'mmol m-3',    'P half-saturation',        default=0.4_rk)
      call self%get_parameter(self%K_N,  'K_N',  'mmol m-3',    'N half-saturation',        default=4.0_rk)
      call self%get_parameter(self%resp, 'resp', 'mmol',        'carbon cost per nitrogen uptake',    default=0.2_rk)
      call self%get_parameter(self%FlexStoich,  'FlexStoich',  '',    'Is FlexStoich?',        default=.true.)
      
      dphyXdt_crit = self%num_dum * fixed_stoichiometry

      ! TODO redesign with transparent indices
      ! Also redesign get_parameter call from auto-generated parameter name like
      ! 'K_'//(self%halfsat(i)%name)
      self%HalfSatNut(1) = self%K_N
      self%HalfSatNut(2) = self%K_P

      !call self%register_dependency(self%id_din,'din', 'mmol-N m-3','dummy')
      ! Register state variables
      call self%register_state_variable(self%id_phytoplankton_C,'phytoplankton_C', 'mmol m-3', 'concentration', 1.0_rk, minimum=0.0_rk)

      ! Register environmental dependencies
      !call self%register_dependency(self%id_prey, "grazing", 'd-1', 'grazing pressure', required = .false.)!, scale_factor = days_per_sec) ! Zooplankton activity
      call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
      call self%register_dependency(self%id_temp, standard_variables%temperature)
      call self%register_dependency(self%id_doy,standard_variables%number_of_days_since_start_of_the_year)

      do i = 1,NUM_CHEM !
         call self%register_state_dependency(self%id_var(i), chemicals(i),'mmol m-3',chemicals(i))
      end do

      call self%register_dependency(self%id_rhs_phy_old, 'old_rhs_phy', 'mmol-C m-3 d-1', 'previous dphyC_dt') 
      call self%register_diagnostic_variable(self%id_rhs_phyC, 'rhs_phyC','mmol-C m-3 d-1', 'current dphyC_dt')

      !  retrieve OM variables for each element
      do i = 1,NUM_ELEM ! e.g., C, N, P (Si, Fe)
         det_index(i) = NUM_CHEM+2*i-1  ! TODO: move to tame_types
         dom_index(i) = NUM_CHEM+2*i
         elem = ElementList(i:i)
         call self%register_state_dependency(self%id_var(det_index(i)), 'det_' // elem,'mmol-' // elem // ' m-3','Detritus ' // trim(ElementName(i)))
         call self%register_state_dependency(self%id_var(dom_index(i)), 'dom_' // elem,'mmol-' // elem // ' m-3','Dissolved Organic ' // trim(ElementName(i)))
  ! print *,det_index(i), ElementList(i:i),dom_index(i)
         !call self%register_dependency(self%Phy_X(i), 'old_Phy_' // elem,'mol-' // elem // ' mol-C-1', elem // ':C-quota')
         !call self%register_dependency(self%id_Phy_X_old(i), temporal_mean(self%Phy_X(i), period=900.0_rk, resolution=60.0_rk)) ! period=900.0_rk, resolution=60.0_rk
        ! call self%register_dependency(self%id_Phy_X_old(i), last_state(self%Phy_X(i))) 
         call self%register_dependency(self%id_Phy_X_old(i), 'old_Phy_' // elem,'mol-' // elem // ' mol-C-1', 'previous phy' // elem) 

         if (elem .NE. 'C') then  ! here only non-carbon elements as Q_C=1 and phytoplankton biomass assumed to be in carbon units 
            call self%register_diagnostic_variable(self%id_Q(i), 'Q_' // elem,'mol-' // elem // ' mol-C-1', elem // ':C-quota')
            call self%register_diagnostic_variable(self%id_dQ_dt(i), 'dQ_dt_' // elem,'mol-' // elem // ' mol-C-1 d-1', 'change in ' // elem // ':C-quota')
            call self%register_diagnostic_variable(self%id_phy_elem(i), 'phy_' // elem,'mol-' // elem // ' m-3', 'phytoplankton ' // elem)
 !!           call self%register_dependency(self%id_Phy_X_old(i), temporal_mean(self%id_Q(i), period=2._rk*60._rk, resolution=60._rk))
!            call self%register_dependency(self%id_Phy_X_old(i), 'old_Q_' // elem,'mol-' // elem // ' mol-C-1', elem // ':C-quota')
!         if (elem == 'P') call self%register_dependency(self%id_Phy_X_old(3), int_change_in_phosphorus)

            call self%add_to_aggregate_variable(type_universal_standard_variable(name='total_' // trim(ElementName(i)), units='mmol-' // elem // ' m-3', &
              aggregate_variable=.true., conserved=.true.), self%id_phy_elem(i))
!            call self%add_to_aggregate_variable(type_universal_standard_variable(name='total_' // trim(ElementName(i)), units='mmol-' // elem // ' m-3', &
!              aggregate_variable=.true., conserved=.true.), self%id_phytoplankton_C, scale_factor = fixed_stoichiometry(i))
         endif
      end do
      call self%register_diagnostic_variable(self%id_din, 'din','mmol-N m-3', 'DIN',output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_nut, 'nut1','mmol-? m-3', 'nutrient related',output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_nut2,'nut2','mmol-? m-3', 'nutrient related',output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_rate,'rate','mmol-? m-3', 'dummy')
      call self%register_diagnostic_variable(self%id_day_of_year,'day_of_year', 'd','day_of_year')

   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_tame_phyto), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_
      real(rk)     :: phytoplankton_C, par, temp, din, nut, phy_X, func, dphy,dphyt,rdphy, dphyXdt, dtime=120., rhs_phy, doy, doy0
      real(rk)     :: production, respiration, sinking, new, mort, quota_2, loss,sgn, nut_lim_tot, old_rhs_phyC, avg_rhs_phyC
      real(rk)     :: chem_change, delta_q,ANutC(NUM_NUTRIENT,NUM_NUTRIENT),iANutC(NUM_NUTRIENT,NUM_NUTRIENT)
      real(rk)     :: nutrient_lim(NUM_NUTRIENT), exudation(NUM_NUTRIENT)
      real(rk)     :: nutrient(NUM_NUTRIENT),rhs_nut(NUM_NUTRIENT),nut_change(NUM_NUTRIENT)
      real(rk)     :: quota(NUM_ELEM), phy_X_old(NUM_ELEM), phy_X_change(NUM_ELEM)
      real(rk)     :: dNut_dt0(NUM_NUTRIENT), quota2, dNut, dQ_dNut(NUM_ELEM,NUM_NUTRIENT)
      real(rk)     :: dix_chemical(NUM_CHEM),total_elem(NUM_ELEM), part(NUM_CHEM), part_safe(NUM_CHEM),rhs_chem(NUM_CHEM)
      logical      :: ncrit(NUM_CHEM), IsPhosporus
      integer      :: i,j,ni,ie ! Indices

      ! Enter spatial loops (if any)
      !print *,stoichiometry
      _LOOP_BEGIN_

         _GET_GLOBAL_(self%id_doy,doy)       ! day of year
         _SET_DIAGNOSTIC_(self%id_day_of_year,doy)

         ! Retrieve current (local) state variable values.
         _GET_(self%id_phytoplankton_C, phytoplankton_C)     ! phytoplankton carbon

         do i = 1,NUM_CHEM ! e.g., CO2, NO3, NH4 (PO4)
            _GET_(self%id_var(i), dix_chemical(i))  ! Dissolved Inorganic Nutrient DIX in mmol-X/m**3
         end do

         ! Retrieve current environmental conditions.
         _GET_(self%id_par,par)          ! local photosynthetically active radiation
         _GET_(self%id_par,temp)         ! water temperature
         _GET_(self%id_Phy_X_old(1), doy0)
         _GET_(self%id_rhs_phy_old, old_rhs_phyC)

         ! TODO replace by TransIndex_DOMDIX, TransIndex2_DOMDIX, which should be set globally
         nutrient = 0.0_rk
         do i = 1,NUM_CHEM ! e.g., CO2, NO3, NH4 (PO4)
            nutrient(chem2nut(i)) = nutrient(chem2nut(i))+ dix_chemical(i)
         end do
         _SET_DIAGNOSTIC_(self%id_din,nutrient(1))

         ! TODO: check and re-formulate sum-rule!
         nut_lim_tot = 0._rk
         do i = 1,NUM_NUTRIENT
            nutrient_lim(i) = nutrient(i)/(self%HalfSatNut(i)+nutrient(i))
            nut_lim_tot = nut_lim_tot + 1.0_rk/(nutrient_lim(i) + small)  ! sum rule for co-limitation, see Wirtz & Kerimoglu 2016
!            nutrient_lim(i) = limitation( self%affinity(i)*nutrient(i)) !/ chem_stoichiometry(i)  ! Add the nutrient limitation law for phytoplankton
         end do
         nut_lim_tot = 1.0_rk/nut_lim_tot
         _SET_DIAGNOSTIC_(self%id_rate, nut_lim_tot*100.0_rk  ) ! 

 !        _SET_DIAGNOSTIC_(self%id_nut2, nut_lim_tot*1.E2 )
         ! Production
         func = 1.0_rk - exp( -self%gamma * par / self%rmax )
         production = self%rmax * func * nut_lim_tot
      !   _SET_DIAGNOSTIC_(self%id_rate, production )
         !do i = 1,NUM_NUTRIENT
         !   exudation(i) =  ( uptake(i) - minval( uptake ) ) * chem_stoichiometry(i) ! Add DOX as a dependency
         !end do

         ! All losses: respiration + sinking
         respiration = self%resp                     ! C loss for DIN-uptake
         sinking     = self%s0                       ! TODO - move to FABM's set_settling
         mort        = self%qmort * phytoplankton_C  ! density-dependent mortality (e.g., virus)
         
         _SET_DIAGNOSTIC_(self%id_nut, mort )

         ! combine to get the temporal derivative for phytoplankton C
         rhs_phy = (production  - sinking - respiration - mort) * phytoplankton_C

         ! Nutrient dynamics
         ! TODO generalize, detect partitioning based on chem2nut (tame_types)
         ! NO3+NH4 partitioning in DIN uptake 
         !   here equal above critical threshold, then down-regulated
         part = 1.0_rk
         do i = 1,2 ! ToDO replace by index related to chem2nut(NUM_CHEM) = (/    1,    1,    2/)
           part(i) = 1.0_rk - exp(-dix_chemical(i)/nut_minval(chem2nut(i)))
         end do
         part_safe = 1.0_rk - exp(-part/0.1_rk)
         part(1:2) = part(1:2)/(sum(part(1:2))+small)

         avg_rhs_phyC = rhs_phy
        !  set quota either as flexible or constant (Redfield)
         if (self%FlexStoich) then
         ! Flexible regulation of non-Redfield stoichiometry (C:N:P)
         !   based on MAECS output
         ! call response functions for intracellular N:C and P:C quotas
            do i = 1, NUM_NUTRIENT           
               j  = nut2othernut(i) ! index of complementary, co-limiting nutrient
               ie = nut2elem(i)
               quota(ie)  = calc_quota(nutrient(i),nutrient(j), par, temp, i,j)  ! calc quota from two nutrient conc. and env. conditions
            end do
           ! _SET_DIAGNOSTIC_(self%id_nut2, nutrient_lim(1))

            if (doy0 .ge. 0 .AND. doy0 .lt. 367) then 
             !  avg_rhs_phyC = self%num_dum*old_rhs_phyC+(1.0_rk-self%num_dum)*rhs_phy
             !  avg_rhs_phyC = 0.5_rk*old_rhs_phyC + 0.5_rk*rhs_phy
               do i = 1,NUM_ELEM
                  if (ElementList(i:i) .NE. 'C') then
                     _GET_(self%id_Phy_X_old(i), phy_X_old(i))
                     phy_X = phytoplankton_C*quota(i)
                     ! compare with critical phy_X change per time-step and correct  
                     ! TODO: separate Q and B (new memory for phy_C)
                     dtime = (doy-doy0+0.001_rk*days_per_sec)
                     dphyXdt  = (phy_X - phy_X_old(i))/dtime
                     ! freeze at too abrupt phy_X changes
                     ! smoothed truncation of temporal change in phyX 
                     rdphy = dphyXdt/dphyXdt_crit(i)
                     quota_2 = quota(i)
                     if (abs(rdphy) .gt. 0.0001_rk ) then
                        ! reset calculated quota to avoid too abrupt shifts
                        phy_X_change(i) = dphyXdt_crit(i) * ( -1._rk + 2._rk/(1._rk + exp(-2*rdphy)) )
                        quota(i) = (phy_X_old(i) + phy_X_change(i)*dtime) / (phytoplankton_C+small)
                        !write (*,'(A10,4F9.2) ') 't/phy_X=',doy0,dtime*1E3,phy_X , phy_X_old(i)
                     !  if (abs(rdphy) .gt. 0.1_rk .AND. i==-3) write (*,'(F7.3,A2,I2,A8,5F9.3,A5,2F9.3,A3,2F9.4) ') doy,ElementList(i:i),i,'r/N/dPdt',rdphy,nutrient(elem2nut(i)),dphyXdt,phy_X_change(i),phytoplankton_C,'rhs:',rhs_phy,old_rhs_phyC,' Q=',quota_2*1E3,quota(i)*1E3
                        !quota(i) = (phy_X_old(i) + dphy * func) / (phytoplankton_C+small)
                     else
                        phy_X_change(i) = dphyXdt
                     endif
                     !if (i == 3) _SET_DIAGNOSTIC_(self%id_nut, (quota(i)-quota_2)*1E3 )

                  else
                     quota(i) = fixed_stoichiometry(i)
                  endif
               end do
            else ! initial period with memory TODO: refine for better mass balance
               phy_X_change = 0.0_rk
            endif
         else
            do i = 1,NUM_ELEM !
               quota(i) = fixed_stoichiometry(i)
               phy_X_change(i) = 0.0_rk
            end do
         endif

         ! get change in nutrients for calculating feed-back: nutrient demand by quota changes
         ! change in chemical: uptake related to (1) new production and (2) quota changes
         do i = 1,NUM_CHEM
            j = chem2elem(i) ! index of element for each chemical - TODO generalize for molecules of >1 resolved element
            ni = chem2nut(i)
            ! ---------- Nutrient sink due to uptake/release by phyto ----------
            dphyt = phy_X_change(j)- avg_rhs_phyC*quota(j) !  rhs_phy*quota(j)
            chem_change = -part(i)*(production * quota(j)* phytoplankton_C + dphyt) 
         ! TODO: check for too strong sraw-down to negative values
            _ADD_SOURCE_(self%id_var(i), chem_change * days_per_sec) 
         end do
      !   if (nutrient(2)<nut_minval(2)*0.2_rk) write (*,'(A4,4F9.4,A1,2F9.4) ') 'PO4 ',nutrient(2),nut_lim_tot,chem_change,phy_X_change(j),'-',rhs_phy*quota(j),(phy_X_change(j)-rhs_phy*quota(j))

         ! here the final change in phy C 
          _ADD_SOURCE_(self%id_phytoplankton_C, rhs_phy  * days_per_sec)

         ! Exudation to DOM (proportional to C-respiration)         

         loss = respiration * phytoplankton_C
         do i = 1,NUM_ELEM
            if (ElementList(i:i) .NE. 'C') then
               _ADD_SOURCE_(self%id_var(dom_index(i)), loss*quota(i) * days_per_sec) ! Nutrients sink TODO check dom_index
               _SET_DIAGNOSTIC_(self%id_Q(i), quota(i))
               _SET_DIAGNOSTIC_(self%id_dQ_dt(i), phy_X_change(i) ) !*secs_per_day
               _SET_DIAGNOSTIC_(self%id_phy_elem(i), quota(i) * phytoplankton_C )
               !write (*,'(2A4,F9.4) ') 'phy',trim(ElementList(i:i)), quota(i) * phytoplankton_C
            else
               _SET_DIAGNOSTIC_(self%id_phy_elem(i), doy )
            endif
         end do
         _SET_DIAGNOSTIC_(self%id_rhs_phyC, rhs_phy)
         _SET_DIAGNOSTIC_(self%id_nut2, phy_X_change(3)/50)

         ! sinking and mortality to POM
         loss = (sinking+mort) * phytoplankton_C

         do i = 1,NUM_ELEM  ! C, N, P (Si, Fe)
            _ADD_SOURCE_(self%id_var(det_index(i)), loss*quota(i) * days_per_sec) !
             !write (*,'(A3,I3,F6.3) ') ElementList(i:i),i,loss*quota(i)
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
end module tame_phyto
