! SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
! SPDX-License-Identifier: Apache-2.0

!> @file tame_stoich_functions.f90
!> @author kai wirtz
!---------------------------------------------------------
! !module: tame_functions
!> @brief stoichiometry related functions (derived from MAECS output) called by tame
#include "fabm_driver.h"

module tame_phy_stoich
use fabm_types
use tame_types
use tame_stoich_functions
 implicit none
! !uses:
 public flexstoich 
 contains
 subroutine flexstoich(phy_C,phy_X_old,nutrient,par, temp, dtime,quota,phy_X_change)
 real(rk), intent(in)   :: phy_X_old(NUM_ELEM),nutrient(NUM_NUTRIENT)
 real(rk), intent(in)   :: phy_C,par, temp, dtime
 real(rk), intent(out)  :: quota(NUM_ELEM),phy_X_change(NUM_ELEM)
 real(rk) :: dphyX_dt,phy_X,rdphy
 integer  :: i,j,ie

   quota(1) = 1._rk
   do i = 1, NUM_NUTRIENT           
      j  = nut2othernut(i) ! index of complementary, co-limiting nutrient
      ie = nut2elem(i)
      quota(ie) = calc_quota(nutrient(i), nutrient(j), par, temp, i, j)
      !if (i==2) write(*,'(3I3,A,3F8.4,A,3F8.4)') i,j,ie,'Q=',quota(ie),nutrient(i), nutrient(j),'P/T/PQo=', par, temp,phy_X_old(2)
   end do
   !write(*,'(A,1F8.4)') 'Q=',quota(3)

   if (abs(dtime) .lt. 1._rk) then 
      do i = 2, NUM_ELEM
         phy_X   = phy_C*quota(i)
         dphyX_dt = (phy_X - phy_X_old(i))/dtime
         rdphy   = dphyX_dt/dphyXdt_crit(i)
         ! checks for too abrupt changes, apply smoothing 
         if (abs(rdphy) .gt. 0.0001_rk) then
            phy_X_change(i) = dphyXdt_crit(i) * (-1._rk + 2._rk/(1._rk + exp(-2*rdphy)))
            quota(i) = (phy_X_old(i) + phy_X_change(i)*dtime) / (phy_C + small)
         else
            phy_X_change(i) = dphyX_dt
         endif
         !if (i==3) write(*,'(A,4F8.4)') 'Q=',quota(3),phy_X_old(2),phy_X_change(2),dphyX_dt

      end do
   else
      phy_X_change = 0._rk
   endif

end subroutine flexstoich

end module
