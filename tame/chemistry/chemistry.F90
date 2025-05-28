! SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
! SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de>
! SPDX-FileContributor: Kai Wirtz <kai.wirt@hereon.de>
! SPDX-License-Identifier: Apache-2.0
!
#include "fabm_driver.h"

module chemistry

   use fabm_types
   use chemistry_types, only : global_molecule_table
   implicit none
   private

   type, extends(type_base_model), public :: type_chlorophyll

      type(type_state_variable_id) :: id_chlorophyll

   contains
      procedure :: initialize
   end type

contains

subroutine initialize(self,configunit)

   class (type_chlorophyll), intent(inout), target :: self
   integer,            intent(in)            :: configunit
   
   call global_molecule_table%register('Methane','CH4')
   call global_molecule_table%register('Chlorophyll-a','C55H72MgN4O5')
   
   !C55H72MgN4O5
   !Molar mass	893.509 g·mol−1

   ! Register state variables
   call self%register_state_variable(self%id_chlorophyll,'chl','', '')


   end subroutine initialize

end module chemistry