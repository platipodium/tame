! SPDX-License-Identifier: CC0-1.0
! SPDX-FileCopyrightText: 2024-2025 Helmholtz-Zentrum hereon GmbH
! SPDX-FileContributor Kai Wirtz <kai.wirtz@hereon.de>
! SPDX-FileContributor Thomas Imbert <thomas.imbert@hereon.de>
! TAME phytoplankton module with flexible stoichiometry
!   and partly based on MAECS (Kai Wirtz et al)

#include "fabm_driver.h"
module tame_phyto

use fabm_types
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
      type (type_dependency_id)          :: id_par, id_temp,id_nut_change(NUM_CHEM),id_Q_old(NUM_ELEM)  ! PAR light
      type (type_diagnostic_variable_id) :: id_nut,id_nut2, id_din,id_rate,id_Q(NUM_ELEM),id_dQ_dt(NUM_ELEM),id_phy_elem(NUM_ELEM)
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
      logical :: FlexStoich

   contains
      procedure :: initialize
      procedure :: do
   end type

   integer :: det_index(NUM_ELEM), dom_index(NUM_ELEM) ! TODO: global setting!

contains
   subroutine initialize(self, configunit)
      class (type_tame_phyto), intent(inout), target :: self
      integer,                        intent(in)            :: configunit
      integer :: i ! Indice dummy
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
      call self%get_parameter(self%FlexStoich,  'FlexStoich',  '',    'Is FlexStoich?',        default=.true.)
      
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
      do i = 1,NUM_CHEM !
         call self%register_state_dependency(self%id_var(i), chemicals(i),'mmol m-3',chemicals(i))
         if (self%FlexStoich) call self%register_dependency(self%id_nut_change(i),'RHS_'//chemicals(i), 'mmol m-3 d-1','dummy')
      end do
      !  retrieve OM variables for each element
      do i = 1,NUM_ELEM ! e.g., C, N, P (Si, Fe)
         det_index(i) = NUM_CHEM+2*i-1  ! TODO: move to tame_types
         dom_index(i) = NUM_CHEM+2*i
         elem = ElementList(i:i)
         call self%register_state_dependency(self%id_var(det_index(i)), 'det_' // elem,'mmol-' // elem // ' m-3','Detritus ' // trim(ElementName(i)))
         call self%register_state_dependency(self%id_var(dom_index(i)), 'dom_' // elem,'mmol-' // elem // ' m-3','Dissolved Organic ' // trim(ElementName(i)))
  ! print *,det_index(i), ElementList(i:i),dom_index(i)

         if (elem .NE. 'C') then  ! here only non-carbon elements as Q_C=1 and phytoplankton biomass assumed to be in carbon units  
            call self%register_diagnostic_variable(self%id_Q(i), 'Q_' // elem,'mol-' // elem // ' mol-C-1', elem // ':C-quota')
            call self%register_diagnostic_variable(self%id_dQ_dt(i), 'dQ_dt_' // elem,'mol-' // elem // ' mol-C-1 d-1', 'change in ' // elem // ':C-quota')
            call self%register_diagnostic_variable(self%id_phy_elem(i), 'phy_' // elem,'mol-' // elem // ' m-3', 'phytoplankton ' // elem)
!            call self%register_dependency(self%id_Q_old(i), 'old_Q_' // elem,'mol-' // elem // ' mol-C-1', elem // ':C-quota')
!         if (elem == 'P') call self%register_dependency(self%id_Q_old(3), int_change_in_phosphorus)

            call self%add_to_aggregate_variable(type_universal_standard_variable(name='total_' // trim(ElementName(i)), units='mmol-' // elem // ' m-3', &
              aggregate_variable=.true., conserved=.true.), self%id_phy_elem(i))
!            call self%add_to_aggregate_variable(type_universal_standard_variable(name='total_' // trim(ElementName(i)), units='mmol-' // elem // ' m-3', &
!              aggregate_variable=.true., conserved=.true.), self%id_phytoplankton_C, scale_factor = fixed_stoichiometry(i))
         endif
      end do
      call self%register_diagnostic_variable(self%id_nut, 'nut1','mmol-? m-3', 'nutrient related')
      call self%register_diagnostic_variable(self%id_nut2,'nut2','mmol-? m-3', 'nutrient related')
      call self%register_diagnostic_variable(self%id_rate,'rate','mmol-? m-3', 'dummy')
      call self%register_diagnostic_variable(self%id_din,'din', 'mmol-N m-3','dummy')

   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_tame_phyto), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_
      real(rk)            :: phytoplankton_C, par, temp, din, nut, func,rhs_phy
      real(rk)            :: production, respiration, sinking, new, loss,sgn, nut_lim_tot
      real(rk)            :: chem_change, delta_q,ANutC(NUM_NUTRIENT,NUM_NUTRIENT),iANutC(NUM_NUTRIENT,NUM_NUTRIENT)
      real(rk)            :: nutrient_lim(NUM_NUTRIENT), exudation(NUM_NUTRIENT)
      real(rk)            :: nutrient(NUM_NUTRIENT),rhs_nut(NUM_NUTRIENT),nut_change(NUM_NUTRIENT)
      real(rk)            :: quota(NUM_ELEM), quota_old(NUM_ELEM), quota_change(NUM_ELEM)
      real(rk)            :: dNut_dt0(NUM_NUTRIENT), quota2, dNut, dQ_dNut(NUM_ELEM,NUM_NUTRIENT)
      real(rk)            :: dix_chemical(NUM_CHEM),total_elem(NUM_ELEM), part(NUM_CHEM),rhs_chem(NUM_CHEM)
      logical             :: ncrit(NUM_CHEM), IsPhosporus
      integer             :: i,j,ie ! Indices

      ! Enter spatial loops (if any)
      !print *,stoichiometry
      _LOOP_BEGIN_
         ! Retrieve current (local) state variable values.
         _GET_(self%id_phytoplankton_C, phytoplankton_C)     ! phytoplankton carbon

         do i = 1,NUM_CHEM ! e.g., CO2, NO3, NH4 (PO4)
            _GET_(self%id_var(i), dix_chemical(i))  ! Dissolved Inorganic Nutrient DIX in mmol-X/m**3
         end do
!          _GET_(self%id_Q_old(3), total_elem(3)) 
!         write (*,'(F8.5) ') total_elem(3)

         ! Retrieve current environmental conditions.
         _GET_(self%id_par,par)          ! local photosynthetically active radiation
         _GET_(self%id_par,temp)         ! water temperature

         ! get change in nutrients for calculating feed-back: nutrient demand by quota changes
         if (self%FlexStoich) then
            do i = 1,NUM_CHEM !
               _GET_(self%id_nut_change(i),rhs_chem(i))
            end do         
!            do i = 1,NUM_ELEM 
!               if (ElementList(i:i) .NE. 'C') _GET_(self%id_Q_old(i),quota_old(i))
!            end do
            rhs_nut  = 0.0_rk
         endif
         ! tODO replace by CHEM -> NUT matrix operation,
         ! TODO replace by TransIndex_DOMDIX, TransIndex2_DOMDIX, which should be set globally
         nutrient = 0.0_rk
         do i = 1,NUM_CHEM ! e.g., CO2, NO3, NH4 (PO4)
            nutrient(chem2nut(i)) = nutrient(chem2nut(i))+ dix_chemical(i)
            if (self%FlexStoich) rhs_nut(chem2nut(i))  = rhs_nut(chem2nut(i)) + rhs_chem(i)
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
!         production = self%rmax * light_absorb(self%rmax, self%gamma, par) !* minval( nutrient_lim )
      !   _SET_DIAGNOSTIC_(self%id_rate, production )
         !_SET_DIAGNOSTIC_(self%id_rate,  new)
         !do i = 1,NUM_NUTRIENT
         !   exudation(i) =  ( uptake(i) - minval( uptake ) ) * chem_stoichiometry(i) ! Add DOX as a dependency
         !end do

         ! Losses
         respiration = self%resp  ! C loss for DIN-uptake
         sinking = self%s0

         ! Nutrient dynamics
         ! TODO generalize, detect partitioning based on chem2nut (tame_types)
         ! NO3+NH4 partitioning in DIN uptake 
         !   here equal above critical threshold, then down-regulated
         part = 1.0_rk
         do i = 1,2 ! ToDO replace by index related to chem2nut(NUM_CHEM) = (/    1,    1,    2/)
           part(i) = 1.0_rk - exp(-dix_chemical(i)/nut_minval(chem2nut(i)))
         end do
         part(1:2) = part(1:2)/(sum(part(1:2))+small)

        !  set quota either as flexible or constant (Redfield)
         if (self%FlexStoich) then
         ! Flexible regulation of non-Redfield stoichiometry (C:N:P)
         !   based on MAECS output
         ! call response functions for intracellular N:C and P:C quotas
            dQ_dNut = 0.0_rk
            do i = 1, NUM_NUTRIENT           
               j  = nut2othernut(i) ! index of complementary, co-limiting nutrient
               ie = nut2elem(i)
               quota(ie)      = calc_quota(nutrient(i),nutrient(j), par, temp, i,j)  ! calc quota from two nutrient conc. and env. conditions
               ! nutrient change = remineralisation (from tame/bgc) - uptake = new production times stoichiometry 
               !   --- uptake by quota change times stock to be evaluated below
               dNut_dt0(i)    = rhs_nut(i) - production * phytoplankton_C * quota(ie)
               ! calculate derivatives of quota on all nutrients, here diagonal entries
               dNut           = sign(nut_minval(i)*0.001_rk,dNut_dt0(i)) ! small change in first nutrient
               quota2         = calc_quota(nutrient(i) + dNut,nutrient(j), par, temp, i,j)  ! calc quota with varied nutrient conc.
               dQ_dNut(ie,i) = (quota2 - quota(ie))/dNut ! derivative = simple numerical difference
            end do
            ! calculate derivatives of quota on all nutrients, non-diagonal entries !TODO generalize to > 2 nutrients
            do i = 1, NUM_NUTRIENT           
               j  = nut2othernut(i) ! index of 2nd, complementary, co-limiting nutrient
               ie = nut2elem(i)
               dNut           = sign(nut_minval(j)*0.001_rk,dNut_dt0(j)) ! small change in 2nd nutrient
               quota2         = calc_quota(nutrient(i),nutrient(j)+ dNut, par, temp, i,j)  ! calc quota with varied nutrient conc.
               dQ_dNut(ie,j) = (quota2 - quota(ie))/dNut ! again, derivative = simple numerical difference
            end do
            
            ! set feed-back in nutrient changes : nutrient demand by nutrient-induced quota changes
            ! solve inverse linear matrix problem A*X = X0  ->  X= A-1 * dNi_dt0 with X = dNi_dt  
            !    with inverse matrix A-1 = 1/(a11*a22-a12*a21)*(a22 -a12 / -a21 a11)
            do i = 1, NUM_NUTRIENT
               ie = nut2elem(i)
               do j = 1, NUM_NUTRIENT
                  sgn = merge(1._rk,0._rk,i==j)
                  ANutC(i,j) = sgn + dQ_dNut(ie,j)*phytoplankton_C
               end do
            end do
            if (NUM_NUTRIENT .eq. 2) then
               iANutC = invmatrix_2(ANutC)
            else
              write(*,*) 'NUM_NUTRIENT differs from 2. Please implement another inverse matrix function'
              stop
               ! https://fortranwiki.org/fortran/show/Matrix+inversion
            end if
            nut_change = matmul(iANutC,dNut_dt0)
!               new =  0.0_rk*dQ_dNut(ie,j)* dNut_dt0(j) * phytoplankton_C ! TODO recheck!
!                 _SET_DIAGNOSTIC_(self%id_nut, dQ_dNut(ie,j) )
            _SET_DIAGNOSTIC_(self%id_nut2, nut_change(2))
!            _SET_DIAGNOSTIC_(self%id_nut, dQ_dNut(3,2) )
!            _SET_DIAGNOSTIC_(self%id_nut, dQ_dNut(3,2) )
            !iANutC(2,2)
            ! re-estimated temporal change in quota based on partial derivaties dQ_dNut and changes in nutrients
            quota_change = matmul(dQ_dNut,nut_change)
               
               ! TODO: complete with ALL derivatives (complicated :-(
               !        check for accuracy first
   ! get change in nutrients for calculating feed-back: nutrient demand by quota changes
         else
            do i = 1,NUM_ELEM !
               quota(i) = fixed_stoichiometry(i)
               quota_change(i) = 0.0_rk
            end do
         endif
         ! change in chemical: uptake related to (1) new production and (2) quota changes
         do i = 1,NUM_CHEM
            j = chem2elem(i) ! index of element for each chemical - TODO generalize for molecules of >1 resolved element
            ! ---------- Nutrient sink due to uptake/release by phyto ----------
            chem_change = -part(i)*( production * quota(j) + quota_change(j)) * phytoplankton_C

            ! if sum is negative: sink of DIX
          !!  if (chem_change .lt. 0.0_rk) then
               _ADD_SOURCE_(self%id_var(i), chem_change * days_per_sec) 
          !!  else    ! if total sum is positive: exudation to DOX
          !!     _ADD_SOURCE_(self%id_var(dom_index(chem2elem(i))), chem_change * days_per_sec)
          !!  end if
         end do
            _SET_DIAGNOSTIC_(self%id_nut, rhs_chem(3)+chem_change)

         ! temporal derivative for phytoplankton C
         rhs_phy = (production  - sinking - respiration) * (phytoplankton_C + self%p0)
         _ADD_SOURCE_(self%id_phytoplankton_C, rhs_phy  * days_per_sec)

         ! Exudation to DOM (proportional to C-respiration)
         ! TODO: check how 2nd _ADD_SOURCE_ works -> join?
         loss = respiration * phytoplankton_C
         do i = 1,NUM_ELEM
            if (ElementList(i:i) .NE. 'C') then
               _ADD_SOURCE_(self%id_var(dom_index(i)), loss*quota(i) * days_per_sec) ! Nutrients sink TODO check dom_index
               _SET_DIAGNOSTIC_(self%id_Q(i), quota(i) )
               _SET_DIAGNOSTIC_(self%id_dQ_dt(i), quota_change(i) )
               _SET_DIAGNOSTIC_(self%id_phy_elem(i), quota(i) * phytoplankton_C )
            endif
         end do
!         _SET_DIAGNOSTIC_(self%id_nut, din )

         ! sinking to POM
         loss = sinking * phytoplankton_C
         do i = 1,NUM_ELEM  ! C, N, P (Si, Fe)
            _ADD_SOURCE_(self%id_var(det_index(i)), loss*quota(i) * days_per_sec) !
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
