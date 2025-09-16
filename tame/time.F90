! SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
! SPDX-License-Identifier: CC0-1.0
! SPDX-FileContributor Carsten Lemmen <carsten.lemmen@hereon.de>
!
! The `tame_time` model is a demonstration model to showcase the capability of
! FABM models using diagnostics and self-dependencies to keep track of
!   1. the time step since the last invocation of FABM.
!   2. the past value of a 3D variable and its change to the current time step.
!

#include "fabm_driver.h"

module tame_time

   use fabm_types
   use fabm_expressions

   implicit none
   private

   integer, parameter :: NUM_C = 2
   type, extends(type_base_model), public :: type_tame_time
      ! A time state variable
      type(type_state_variable_id)        :: id_time
      type(type_dependency_id)            :: id_previous_time
      type (type_diagnostic_variable_id)  :: id_current_time

      ! A concentration state variable that changes over time
      type(type_state_variable_id)        :: id_conc(NUM_C)
      type(type_dependency_id)            :: id_previous(NUM_C)
      type (type_diagnostic_variable_id)  :: id_current(NUM_C)

   contains
      procedure :: initialize
      procedure :: do
   end type

contains

subroutine initialize(self,configunit)

   class (type_tame_time), intent(inout), target :: self
   integer,            intent(in)            :: configunit
   integer  :: i ! concentration index

   ! These should *not* be 3D variables. They are just scalars.
   call self%register_state_variable(self%id_time, 'time', 's', 'Time', 0.0_rk, minimum=0.0_rk)
   call self%register_diagnostic_variable(self%id_current_time,'current_time','s','Current time')
   call self%register_dependency(self%id_previous_time, 'previous_time', 's','Previous time')
   do i = 1,NUM_C       
      call self%register_state_variable(self%id_conc(i), 'concentration_' // char(48+i), 'mM', 'Concentration ', 0.0_rk, minimum=0.0_rk)
      call self%register_diagnostic_variable(self%id_current(i),'current_' // char(48+i),'mM','Current concentration ' // char(48+i))
      call self%register_dependency(self%id_previous(i), 'previous_' // char(48+i), 'mM','Previous concentration ' // char(48+i))
   end do
   end subroutine initialize

subroutine do(self,_ARGUMENTS_DO_)
   class(type_tame_time), INTENT(IN) :: self
   _DECLARE_ARGUMENTS_DO_

   real(rk) :: conc(NUM_C), previous(NUM_C)
   real(rk) :: time, previous_time
   integer  :: i ! concentration index

   _LOOP_BEGIN_

      do i = 1,NUM_C     
         _GET_(self%id_conc(i),conc(i))
         _GET_(self%id_previous(i), previous(i))
         if (previous(i) < 0) previous(i) = conc(i)
      end do

      _GET_(self%id_time,time)
      _GET_(self%id_previous_time, previous_time)

      if (previous_time < 0) previous_time = time

      do i = 1,NUM_C     
         _ADD_SOURCE_(self%id_conc(i),0.1*cos(0.1_rk*i*time))
         _SET_DIAGNOSTIC_(self%id_current(i),conc(i))
      end do
      _SET_DIAGNOSTIC_(self%id_current_time,time)
      
      _ADD_SOURCE_(self%id_time,1.0_rk)

      write(*,'(A,I7.0,X,A,I3.0,X,A,2F5.2,X,A,2F6.2)') &
        't=', int(time), 'dt=', int(time - previous_time), 'C=', conc,  'dC= ', conc - previous

   end subroutine do
end module tame_time
