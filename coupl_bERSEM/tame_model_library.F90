! SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
!
! SPDX-License-Identifier: Apache-2.0

module tame_model_library

   use fabm_types, only: type_base_model_factory,type_base_model
   use tame_npzdep
      ! Add new tame modules here

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory), save, target, public :: tame_model_factory

contains

   subroutine create(self,name,model)

      class (type_factory), intent(in) :: self
      character(*),         intent(in) :: name
      class (type_base_model), pointer :: model

      select case (name)
         case ('npzdep');   allocate(type_tame_npzdep::model)
         ! Add new tame models here
      end select

   end subroutine

end module
