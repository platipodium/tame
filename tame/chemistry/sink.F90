#include "fabm_driver.h"

module tame_chemistry_sink
   use fabm_types

   implicit none
   private

   type, extends(type_base_model),public :: type_tame_chemistry_sink

      type (type_state_variable_id)     :: id_sink

      contains
      procedure :: initialize
      procedure :: do
   end type type_tame_chemistry_sink

contains

   subroutine initialize(self, configunit)
      class (type_tame_chemistry_sink), intent(inout), target :: self
      integer,                        intent(in)              :: configunit

      call self%register_state_variable(self%id_sink, 'sink','mmol m-3',  'sink concentration', 0.0_rk, &
            minimum=0.0_rk)
      
   end subroutine initialize

subroutine do(self,_ARGUMENTS_DO_)
   class(type_tame_chemistry_sink), INTENT(IN) :: self
   _DECLARE_ARGUMENTS_DO_

   real(rk) :: sink
   _LOOP_BEGIN_

      ! this is not really needed
      _GET_(self%id_sink, sink)
      _ADD_SOURCE_(self%id_sink,0.0_rk)

   _LOOP_END_

   end subroutine do

end module tame_chemistry_sink
