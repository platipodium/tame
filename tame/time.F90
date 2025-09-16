! SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
! SPDX-License-Identifier: CC0-1.0
! SPDX-FileContributor Carsten Lemmen <carsten.lemmen@hereon.de>
!
! The `tame_time` model is a demonstration model to showcase the capability of
! FABM models using diagnostics and self-dependencies to keep track of
!   1. the time step since the last invocation of FABM.
!   2. the past value of a 3D variable and its change to the current time step.
!
! Attention: the time is diluted, I have not figured out how to make this a global diagnostic

#include "fabm_driver.h"

module tame_time

   use fabm_types
   use fabm_expressions

   implicit none
   private

   type, extends(type_base_model), public :: type_tame_time

      ! A time state variable
      type(type_state_variable_id)        :: id_time
      type(type_dependency_id)            :: id_previous_time
      type (type_diagnostic_variable_id)  :: id_current_time

      ! A concentration state variable that changes over time
      type(type_state_variable_id)        :: id_conc
      type(type_dependency_id)            :: id_previous
      type (type_diagnostic_variable_id)  :: id_current

   contains
      procedure :: initialize
      procedure :: do
   end type

contains

subroutine initialize(self,configunit)

   class (type_tame_time), intent(inout), target :: self
   integer,            intent(in)            :: configunit

   ! These should *not* be 3D variables. They are just scalars.
   call self%register_state_variable(self%id_time, 'time', 's', 'Time', 0.0_rk, minimum=0.0_rk)
   call self%register_diagnostic_variable(self%id_current_time,'current_time','s','Current time')
   call self%register_dependency(self%id_previous_time, 'previous_time', 's','Previous time')

   call self%register_state_variable(self%id_conc, 'concentration', 'mM', 'Concentration', 0.0_rk, minimum=0.0_rk)
   call self%register_diagnostic_variable(self%id_current,'current','mM','Current concentration')
   call self%register_dependency(self%id_previous, 'previous', 'mM','Previous concentration')

   end subroutine initialize

subroutine do(self,_ARGUMENTS_DO_)
   class(type_tame_time), INTENT(IN) :: self
   _DECLARE_ARGUMENTS_DO_

   real(rk) :: conc, previous
   real(rk) :: time, previous_time

   _LOOP_BEGIN_

      _GET_(self%id_conc,conc)
      _GET_(self%id_previous, previous)
      _GET_(self%id_time,time)
      _GET_(self%id_previous_time, previous_time)

      if (previous < 0) previous = conc
      if (previous_time < 0) previous_time = time


      _ADD_SOURCE_(self%id_conc,0.1*cos(0.1_rk*time))
      _ADD_SOURCE_(self%id_time,1.0_rk)
      _SET_DIAGNOSTIC_(self%id_current,conc)
      _SET_DIAGNOSTIC_(self%id_current_time,time)

      write(*,'(A,I7.0,X,A,I3.0,X,A,F5.2,X,A,F6.2)') &
        't=', int(time), 'dt=', int(time - previous_time), 'C=', conc,  'dC= ', conc - previous
   _LOOP_END_

   end subroutine do
end module tame_time
