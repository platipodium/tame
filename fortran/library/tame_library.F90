! SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
!
! SPDX-License-Identifier: Apache-2.0

!> @file tame_library.F90
!> @brief tame_library module
! SPDX-FileContributor Carsten Lemmen <carsten.lemmen@hereon.de
! SPDX-FileContributor Ovidio

#include "fabm_driver.h"

!> @brief  Library functions types used in fabm_tame are defined here
!> @details
!! Parser could do it just like in the initialize subroutine.
! by adding:
! \\describepar{model_par, symbol, : some description}
! lines before each type
module tame_library

use fabm_types
use tame_types
use fabm_expressions

public

real(rk), PARAMETER :: Pi = 3.1415927_rk
real(rk), PARAMETER :: d_per_s = 1.0_rk/86400.0_rk
real(rk), PARAMETER :: small_const = 1.0E-9_rk

! colimitation gfunction
! Written for 2 inputs, if you apply it successively, it could
! deal with n inputs
! @todo there should be a function gfunc_infinite that calls this one.
! It is a surrogate function developed by Kai for the queuing function
!@param x  a limiting quantity (or aggregate of some limiting quantity)
!@param y  another limiting quantity (or aggregate of some limiting quantity)
!@param n  shape parameter
real(rk) function gfunc(x,y,n)

    implicit none
    real(rk), intent(in) :: x,y,n
    real(rk) :: r

    r = (x+small_const)/(y+small_const)
    r = r*(1.0_rk-r**n)/(1.0_rk-r**(n+1.0_rk))
    gfunc = y*r*(1.0_rk+y*x/n+log(4.0_rk**(-1.0_rk/n) + 1.0_rk/(2.0_rk*n)))
end function gfunc

end module
