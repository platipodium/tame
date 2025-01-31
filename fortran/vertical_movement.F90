
!> @brief handles vertical movement for depth-varying movement rates
!> @details phyto sinking rate depends on the nutritional state, so for each node:
!! \n \f$ phy\%relQ \f$ obtained by calling calc_internal_states(self,phy,det,dom,zoo)
!! \n then \f$ phyQstat=phy\%relQ\%N * phy\%relQ\%P \f$
!! \n finally, vs_phy = tame_functions::sinking(self\%vS_phy, phyQstat, vs_phy)
subroutine tame_get_vertical_movement(self,_ARGUMENTS_GET_VERTICAL_MOVEMENT_)

use tame_functions
use tame_types

implicit none
!
! !INPUT PARAMETERS:
 class (type_hzg_tame),intent(in) :: self
_DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_
 !   REALTYPE, intent(in)              ::vstokes
type (type_tame_om):: det, dom 
logical  :: IsCritical = .false. ! phyC and phyN below reasonable range ?

!
! !LOCAL VARIABLES:
REALTYPE    :: phyQstat, ef, vs_phy,vs_det, phyEner, minPigm,min_Cmass, minc, zmax, par, vs_zoo
!REALTYPE    :: aggf, agge=16.d0
REALTYPE, parameter :: secs_pr_day = 86400.d0
!EOP
!-----------------------------------------------------------------------
!BOC

!define _REPLNAN_(X) X !changes back to original code
#define _REPLNAN_(X) nan_num(X)

_FABM_LOOP_BEGIN_

   ! Retrieve phtoplankton state

   !Retrieve the 'phyQstat' directly as a diagnostic variable: does not work yet.
   !fabm_get_bulk_diagnostic_data(self%id_phyqstat,phyQstatD) !where, phyQstat=relQ%N*relQ%P
   !_GET_(self%id_phyqstat,phyQstatD)

   !Calculate manually
   _GET_(self%id_phyC, phy%C)  ! Phytplankton Carbon in mmol-C/m**3
   _GET_(self%id_phyN, phy%N)  ! Phytplankton Nitrogen in mmol-N/m**3
   _GET_(self%id_detC, det%C)  ! Detritus Nitrogen in mmol-N/m**3
!   _GET_(self%id_detN, det%N)  ! Detritus Nitrogen in mmol-N/m**3
   _GET_(self%id_domC, dom%C)  ! DONitrogen in mmol-N/m**3
    if (self%PhosphorusOn) then
!      _GET_(self%id_detP, det%P)  ! Detritus Phosphorus in mmol-P/m**3
      _GET_(self%id_phyP, phy%P)  ! Phytplankton Phosphorus in mmol-P/m**3
    endif
 !
 ! ****** SINKING AS A FUNCTION OF INTERNAL STATES *************
 !
  !!  call min_mass(self,phy, min_Cmass, IsCritical, method=2)
  !!  call calc_internal_states(self,phy,det,dom,zoo)
    call calc_internal_states(self,det,dom)
    !write (*,'(A,2(F10.3))') 'phy%relQ%N, phy%relQ%P=', phy%relQ%N, phy%relQ%P

    phyQstat = 1.0_rk
! infected cells sink faster to save population
    if (self%vir_loss .gt. self%small_finite .or. self%VirusOn ) then
      _GET_(self%id_vir, phy%vir)  ! Virus C density in cells in -
      phyQstat = phyQstat * 1.0_rk/(1.0_rk + exp(-self%vir_infect*(1.0_rk-phy%vir/(phy%C+self%small_finite))))
! threshold virus with multi-stage replication
    endif
! non-linear response of sinking speed to nutrient limitation
    if(self%sink_phys .gt. 0) then ! weak non-linearity
      phyQstat = phyQstat * phy%relQ%N**2/(HALFQ**2+phy%relQ%N**2)
      if (self%PhosphorusOn) then
        phyQstat =phyQstat * phy%relQ%P**2/(HALFQ**2+phy%relQ%P**2)
      end if
     ! Calculate sinking speed
      vs_phy = -self%vS_phy * exp( -self%sink_phys * phyQstat)
    else   ! strong non-linearity
      phyQstat = phyQstat * 1.0_rk/(1.0_rk + exp(self%sink_phys*(phy%relQ%N-HALFQ)))
      if (self%PhosphorusOn) then
        phyQstat =phyQstat * 1.0_rk/(1.0_rk + exp(self%sink_phys*(phy%relQ%P-HALFQ)))
      end if
      vs_phy = -self%vS_phy * (1.0_rk-phyQstat)
   endif
   ! nutrient limitation ; TODO check product rule and add other elements such as Si
   ! phyQstat = phy%relQ%N * phy%relQ%P

! ascending of top-conditioned cells
   if(self%genMeth .gt. 0) then
     vs_phy = vs_phy + self%vS_phy * exp(-3.0d0+self%genMeth*0.2d0)
   elseif(self%genMeth .lt. 0) then
     vs_phy = vs_phy - self%vS_phy * exp(-3.0d0-self%genMeth*0.2d0)
   endif
   !if (self%RateDiagOn) then
   ! write (*,'(A)',advance='no') '' ! Silly Fix to 'NETCDF: Numeric conversion not representable' problem ??
   !  _SET_DIAGNOSTIC_(self%id_vsinkr, _REPLNAN_(vs_phy)) !average Relative Sinking Velocity
   !end if
   !CONSTANT SINKING
   !vs_phy = self%vS_phy

   _SET_DIAGNOSTIC_(self%id_vsinkr, _REPLNAN_(vs_phy)) !average Relative_Sinking_Rate_

  _GET_HORIZONTAL_(self%id_zmax, zmax)  ! max depth
! increase vDet in shallow water:co-agulation with lithogenic particles
!   vs_det = vs_det * (1.0_rk+ 4*exp(-(zmax/20.0_rk)**2))
   vs_phy = vs_phy / secs_pr_day
   !write (*,'(A,2(F10.3))') 'phyQstat, vs_phy=', phyQstat, vs_phy
!   vs_det = -self%vS_det*aggf/secs_pr_day
   vs_det = -self%vS_det / secs_pr_day

   vs_det = vs_det * (0.01*det%C)**0.38_rk

! slowing down of vertical velocities at high and very low concentration to smooth numerical problems in shallow, pesitional boxes
   ef     = 20_rk/(1+zmax)
   vs_det = vs_det * 1.0_rk/(1.0_rk+((0.003*det%C +0.02*det%N + 100*self%small_finite/(det%C+self%small_finite))*ef )**4 )
! additional slowdown in very shallow waters
!   if (vs_det*secs_pr_day .gt. 1.0_rk .and. zmax .lt. 8.0_rk) then
!     ef = exp(2*(zmax-5.0_rk))
!     vs_det = vs_det * (1.0_rk+ef)/(3.0_rk+ef)
!   endif
   vs_phy = vs_phy * 1.0_rk/(1.0_rk+((0.002*phy%C + 100*self%small_finite/(phy%C+self%small_finite))*ef )**4 )

  !set the rates
   _SET_VERTICAL_MOVEMENT_(self%id_detC,vs_det)
   _SET_VERTICAL_MOVEMENT_(self%id_detN,vs_det)
   _SET_VERTICAL_MOVEMENT_(self%id_phyN,vs_phy)
   _SET_VERTICAL_MOVEMENT_(self%id_phyC,vs_phy)
!   if (ef .lt. 1.0_rk) then
!    _GET_(self%id_par, par)  ! light photosynthetically active radiatio
!    vs_zoo = 22.0_rk * (1.0_rk-2.0_rk/(1.0_rk + exp(1.0d0-par)))
!    _SET_VERTICAL_MOVEMENT_(self%id_zooC,vs_zoo / secs_pr_day)
!   else
!    _SET_VERTICAL_MOVEMENT_(self%id_zooC,0.0_rk)
!   endif
   if (self%PhosphorusOn) then
      _SET_VERTICAL_MOVEMENT_(self%id_phyP,vs_phy)
      _SET_VERTICAL_MOVEMENT_(self%id_detP,vs_det)
   end if
   if (self%SiliconOn) then
      _SET_VERTICAL_MOVEMENT_(self%id_phyS,vs_phy)
      _SET_VERTICAL_MOVEMENT_(self%id_detS,vs_det)
   end if
   if (self%PhotoacclimOn) then
      _SET_VERTICAL_MOVEMENT_(self%id_chl, vs_phy)
      _SET_VERTICAL_MOVEMENT_(self%id_Rub, vs_phy)
   end if
   if (self%VirusOn) then
      _SET_VERTICAL_MOVEMENT_(self%id_vir, vs_phy)
   end if

_FABM_LOOP_END_

end subroutine tame_get_vertical_movement

subroutine tame_do_surface(self,_ARGUMENTS_DO_SURFACE_)
   use tame_functions

   class (type_hzg_tame), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_

   real(rk) :: tot_vi_C,tot_vi_N,tot_vi_P,tot_vi_S, O2flux,O2airbl,oxy,tot_vi_GPPR,tot_vi_Denitr

!define _REPLNAN_(X) X !changes back to original code
#define _REPLNAN_(X) nan_num(X)

   _HORIZONTAL_LOOP_BEGIN_

      if (self%Budget2DDiagOn) then
      _GET_HORIZONTAL_(self%id_totN_vertint,tot_vi_N)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_totN_vertint_diag,_REPLNAN_(tot_vi_N))
      _GET_HORIZONTAL_(self%id_totC_vertint,tot_vi_C)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_totC_vertint_diag,_REPLNAN_(tot_vi_C))
      if (self%PhosphorusOn) then
         _GET_HORIZONTAL_(self%id_totP_vertint,tot_vi_P)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_totP_vertint_diag,_REPLNAN_(tot_vi_P))
      end if
      if (self%SiliconOn) then
         _GET_HORIZONTAL_(self%id_totS_vertint,tot_vi_S)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_totS_vertint_diag,_REPLNAN_(tot_vi_S))
      end if
      end if

! --- wet and dry deposition of NO3
      _SET_SURFACE_EXCHANGE_(self%id_nutN, self%N_depo UNIT)
! --- atmospheric deposition of PO4
      if (self%PhosphorusOn) then
         _SET_SURFACE_EXCHANGE_(self%id_nutP, self%P_depo UNIT)
      end if

! --- oxygen flux between sea water and air -----
      if (self%BioOxyOn) then
! O2 flux across the boundary layer
! O2airbl is the saturation concentration of O2
! airsea_ex is the average diffusivity coefficient (m2/sec) divided by the thickness of the boundary layer.
! for O2 in mmol m-3, the rate of exchange in mmol m-2 s-1).
! Positive values imply a flux into the water, negative: out of the water.
         O2airbl = self%O2_sat
!        _GET_HORIZONTAL_(self%id_O2airbl, O2airbl)! boundary layer dissolved oxygen in mmolO2/m**3
        _GET_(self%id_oxy, oxy)   ! sea water dissolved oxygen in mmolO2/m**3

        O2flux  = self%ex_airsea * (O2airbl - oxy)!
        _SET_SURFACE_EXCHANGE_(self%id_oxy, O2flux )
        if (self%BGC2DDiagOn) then
          _SET_HORIZONTAL_DIAGNOSTIC_(self%id_O2flux_diag, _REPLNAN_(O2flux)) ! converts mmol/m2.s to mmol/m2.d
        end if
      endif
   _HORIZONTAL_LOOP_END_
end subroutine tame_do_surface