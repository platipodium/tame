#include "fabm_driver.h"

module tame_nutrient
   use fabm_types

   implicit none

   private

   type, extends(type_base_model),public :: type_tame_nutrient
      ! Variable identifiers
      type (type_state_variable_id) :: id_n
   contains
      procedure :: initialize
   end type type_tame_nutrient

contains

   subroutine initialize(self, configunit)
      class (type_tame_nutrient), intent(inout), target :: self
      integer,                        intent(in)            :: configunit

      ! Register state variables
      call self%register_state_variable(self%id_n, 'c', 'mmol m-3', 'concentration', 1.0_rk, minimum=0.0_rk, no_river_dilution=.true.)

      ! Register contribution of state to global aggregate variables.
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_n)
   end subroutine initialize

end module tame_nutrient
