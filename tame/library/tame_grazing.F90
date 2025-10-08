! SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
! SPDX-FileContributor: Kai Wirtz <kai.wirtz@hereon.de>
! SPDX-License-Identifier: Apache-2.0

#include "fabm_driver.h"
module tame_grazing

! !USES:
   use fabm_types
   use tame_types
!   use tame_functions
   private
   public   grazing, grazing_excretion
 contains
 !> @brief  calculates grazing rate
!> @details
!> rate= @f$ I_{max} * F^2/(K^2+F^2) @f$
!> (Holling-III response function)

real(rk) function grazing(preyC,HalfSatC)
  implicit none
  real(rk), intent(in)       :: HalfSatC,PreyC

  grazing   = PreyC**2 /(HalfSatC**2+PreyC**2)          ! [d^{-1}]
  !      ivlev     = 1.0_rk - exp(-prey / saturation)
  ! TODO: implement further functional responses, also adaptive down-regulation
  end function grazing

!---------------------------------------------------------
!> @brief  C-N-P excretion rates of grazers to variable prey stoichiometry
!> @details
!> assumes a constant C:N:P unit feeding on variable C:N:P food

subroutine grazing_excretion(prey_Q,ingestion,resp_hetero,excretion)
  implicit none
  !type (type_tame_switch), intent(in)  :: mswitch
  real(rk), intent(in)            :: prey_Q(NUM_ELEM)
  real(rk), intent(in)            :: ingestion,resp_hetero
  real(rk), intent(out)           :: excretion(NUM_ELEM)
  real(rk) :: Excess_C_upt(NUM_ELEM)
  integer  :: i

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

! -------------- loss by floppy feeding and egestion to detritus
!----------------------------------------------------------------
! -------------- homeostasis for zooplankton
!     default rates, given that N:C and P:C of prey match 'const_NC_zoo' and 'const_PC_zoo'

! -------------- basal respiration  ; TODO: add activity respiration
! -------------- assuming homeostasis for zooplankton relaxation towards Redfield ratio
!                 nitrogen excretion (urea or ammonia)
!      if (mswitch%isTotIng) then  ! compensate with unused prey components
!         if (lossZDet%P .lt. 0.0d0 ) then ! compensate by respiration & exudation
!  if (mswitch%) then ! Si-release
! neglecting silicification in zooplankton such as in choanoflagellates, radiolarians
!    lossZNut%Si   = zoo%yield * zoo%feeding * Q_prey%Si !(mmolSi m^{-3} d^{-1})
! -------------- loss by floppy feeding to detritus
!  end if ! mswitch%isSi

  end subroutine grazing_excretion
end module tame_grazing
