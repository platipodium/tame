module tame_model_library

   use fabm_types, only: type_base_model_factory, type_base_model

   implicit none

   private

   type, extends(type_base_model_factory) :: type_factory
   contains
      procedure :: create
   end type

   type (type_factory), save, target, public :: tame_model_factory

contains

   subroutine create(self, name, model)

!      use tame_light

      class (type_factory), intent(in) :: self
      character(*),         intent(in) :: name
      class (type_base_model), pointer :: model

      select case (name)
;         case ('light');  allocate(type_tame_light::model)
      end select

   end subroutine

end module
