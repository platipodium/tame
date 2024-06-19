!> @file tame_library.F90
!> @brief tame_library module
! SPDX-FileContributor Carsten Lemmen <carsten.lemmen@hereon.de

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

real(rk), PARAMETER ::  Pi = 3.1415927_rk


end module
