module tame_model_library

   use fabm_types, only: type_base_model_factory,type_base_model

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory), save, target, public :: tame_model_factory

contains

   subroutine create(self,name,model)

   ! Add new tame modules here
   use tame_bgc
   use chemistry
   use tame_phytoplankton
   use tame_detritus
   use tame_nutrient
   use tame_phyto

   !use tame_npzdep

      class (type_factory), intent(in) :: self
      character(*),         intent(in) :: name
      class (type_base_model), pointer :: model

      select case (name)
 !      case ('npzdep');   allocate(type_tame_npzdep::model)
      case ('bgc');   allocate(type_tame_bgc::model)
      case ('chlorophyll');   allocate(type_chlorophyll::model)
      case ('phytoplankton');   allocate(type_tame_phytoplankton::model)
      case ('detritus');   allocate(type_tame_detritus::model)
      case ('nutrient');   allocate(type_tame_nutrient::model)
         case ('phyto');   allocate(type_tame_phyto::model)
         ! Add new tame models here
      end select

   end subroutine

end module
