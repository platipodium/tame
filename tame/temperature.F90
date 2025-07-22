#include "fabm_driver.h"

!module examples_npzd_phy ! TAME temperature module
module tame_temperature

use fabm_types
use fabm_expressions

! for demonstration only

   implicit none
   private
   !type, extends(type_base_model), public :: type_examples_npzd_phy
   type, extends(type_base_model), public :: type_tame_temperature

      type (type_dependency_id)                  :: id_temp
      type (type_dependency_id)                  :: id_temp_tempmean
      type (type_global_dependency_id)           :: id_doy
      type (type_diagnostic_variable_id)         :: id_temp_tempmean_diag
      type (type_diagnostic_variable_id)         :: id_temp_change, id_doy_diag

   contains
      procedure :: initialize
      procedure :: do
   end type


contains
    subroutine initialize(self, configunit)
      class (type_tame_temperature), intent(inout), target :: self
      integer,                    intent(in)            :: configunit

      call self%register_dependency(self%id_temp, standard_variables%temperature)
      !@todo check for parametrizing the period and resolution
      call self%register_dependency(self%id_temp_tempmean, temporal_mean(self%id_temp, period=1._rk*3600._rk, resolution=3600._rk))
      call self%register_dependency(self%id_doy, standard_variables%number_of_days_since_start_of_the_year)

      call self%register_diagnostic_variable(self%id_temp_tempmean_diag, 'temp_mean',  'degree_C', '1-hour running mean of temperature')
      call self%register_diagnostic_variable(self%id_temp_change, 'temp_change',  'degree_C h-1', 'hourly change of temperature')
      call self%register_diagnostic_variable(self%id_doy_diag, 'doy',  '', 'Day of the year')
   end subroutine initialize

 subroutine do(self,_ARGUMENTS_DO_)
      class (type_tame_temperature), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: temp_mean, temp, doy

      _GET_GLOBAL_(self%id_doy, doy)

      _LOOP_BEGIN_
         _GET_(self%id_temp_tempmean, temp_mean)
         _GET_(self%id_temp, temp)
         _SET_DIAGNOSTIC_(self%id_temp_tempmean_diag ,temp_mean)
         _SET_DIAGNOSTIC_(self%id_temp_tempmean_diag ,temp_mean)
         _SET_DIAGNOSTIC_(self%id_temp_change ,2*(temp - temp_mean))
         _SET_DIAGNOSTIC_(self%id_doy_diag ,doy)
         ! 2 may be 1 + period/resolution
      _LOOP_END_
   end subroutine do

end module tame_temperature
