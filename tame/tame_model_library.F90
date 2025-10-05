! SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
!
! SPDX-License-Identifier: Apache-2.0

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
   use tame_life
   use chemistry
   !use tame_chemistry_phytoplankton
   use tame_chemistry_detritus
   use tame_chemistry_nutrient
   use tame_chemistry_sink
   !use tame_phyto
   ! use tame_zooplankton
   !use tame_temperature
   use tame_time

   !use tame_npzdep

      class (type_factory), intent(in) :: self
      character(*),         intent(in) :: name
      class (type_base_model), pointer :: model

      select case (name)
 !      case ('npzdep');   allocate(type_tame_npzdep::model)
      case ('bgc');   allocate(type_tame_bgc::model)
      case ('chlorophyll');   allocate(type_chlorophyll::model)
      !case ('chemistry_phytoplankton');   allocate(type_tame_chemistry_phytoplankton::model)
      case ('chemistry_detritus');   allocate(type_tame_chemistry_detritus::model)
      case ('chemistry_nutrient');   allocate(type_tame_chemistry_nutrient::model)
      case ('chemistry_sink');   allocate(type_tame_chemistry_sink::model)
         !case ('chlorophyll');   allocate(type_chlorophyll::model)
!         case ('phytoplankton');   allocate(type_tame_phytoplankton::model)
!         case ('phyto');   allocate(type_tame_phyto::model)
         case ('life');   allocate(type_tame_life::model)
         case ('time');   allocate(type_tame_time::model)
!         case ('zooplankton');   allocate(type_tame_zooplankton::model)
!         case ('temperature');   allocate(type_tame_temperature::model)

         ! Add new tame models here
      end select

   end subroutine

end module
