! SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
! SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de>
! SPDX-License-Identifier: Apache-2.0
#include "fabm_driver.h"

module tame_chemistry_nutrient
   use fabm_types
   use chemistry_types

   implicit none

   private

   type, extends(type_base_model),public :: type_tame_chemistry_nutrient
      type (type_state_variable_id), allocatable     :: id_nutrient(:)
      type (type_group), allocatable :: state_variables(:)

   contains
      procedure :: initialize
   end type type_tame_chemistry_nutrient

contains

   subroutine initialize(self, configunit)
      class (type_tame_chemistry_nutrient), intent(inout), target :: self
      integer,                        intent(in)            :: configunit

      integer :: i, n

      call global_molecule_table%register('Nitrate','NO3')
      call global_molecule_table%register('Nitrite','NO3')
      call global_molecule_table%register('Phosphate','PO4')

      n = 2
      allocate(self%state_variables(n))
      call self%state_variables(1)%create('DIN',(/'Nitrite','Nitrate'/))
      call self%state_variables(2)%create('DIP',(/'Phosphate'/))

      allocate(self%id_nutrient(n))


      do i=1,n
        call self%register_state_variable(self%id_nutrient(i), trim(self%state_variables(i)%name), &
           'mmol m-3',  'concentration', 4.5_rk, minimum=0.0_rk)
     enddo

     do i=1,n
        ! @todo
        !if self%state_variables(i)%contains('N') call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_d)
        !if self%state_variables(i)%contains('P') call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_d)
      enddo


   end subroutine initialize

end module tame_chemistry_nutrient
