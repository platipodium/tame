#include "fabm_driver.h"
# define NUM_ELEM 3
# define NUM_CHEM 3
!----------------------------------------
!	tame/bgc
!
!> @brief ** BioGeoChemical Equations ** 
!> This is the main BGC routine where right-hand-sides are calculated
!! for organic material (Det, DOM) with arbitrary element units such carbon, nitrogen, & phosphorus
!!  and arbitrary dissolved chemical species such as NO3, NH4, or CO2 as state variables.
!----------------------------------------
! !INTERFACE:
module tame_bgc

use fabm_types
!use tame_types
!use tame_functions 

implicit none

 private
 ! !PUBLIC_DERIVED_TYPES:
 type,extends(type_base_model),public :: type_tame_bgc
  type (type_dependency_id) :: id_temp

contains
	procedure :: initialize
	procedure :: do
! procedure :: do_bottom
! procedure :: get_sinking_rate
end type
!----------------------------------------
contains
!----------------------------------------!BOP
!
! !IROUTINE: Initialise the tame/bgc model
! Reading namelist and registration of variables with FABM
!
subroutine initialize(self,configunit)
 class (type_tame_bgc), intent(inout), target :: self
 integer,		intent(in)		:: configunit

 call self%register_dependency(self%id_temp, standard_variables%temperature)

end subroutine initialize
!----------------------------------------
!

! !IROUTINE: Right hand sides of tame/bgc model
!
! !INTERFACE:
 subroutine do(self,_ARGUMENTS_DO_)
! !LOCAL VARIABLES:
 class (type_tame_bgc),intent(in) :: self
	real(rk) :: par, temp, ddix, dummy(10)

_LOOP_BEGIN_

!---------- first get ambient conditions ----------
!_GET_(self%id_temp, env%temp)  ! water temperature 
_GET_(self%id_temp, ddix)  ! water temperaturedummy

_LOOP_END_

end subroutine do
