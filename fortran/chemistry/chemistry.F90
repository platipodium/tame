! SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
! SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de>
! SPDX-FileContributor: Kai Wirtz <kai.wirt@hereon.de>
! SPDX-License-Identifier: Apache-2.0
!
#include "fabm_driver.h"

module chemistry

   use fabm_types
   use chemistry_types, only : global_molecule_table, type_group
   use iso_fortran_env, only : stdout => output_unit, stderr => error_unit
   implicit none
   private

   type, extends(type_base_model), public :: type_chlorophyll

      type(type_state_variable_id) :: id_chlorophyll
      type(type_group) :: din, nutrient

   contains
      procedure :: initialize
   end type

contains

subroutine initialize(self,configunit)

   class (type_chlorophyll), intent(inout), target :: self
   integer,            intent(in)            :: configunit

   character(len=3) :: din_molecules(3) = (/'NO3','NH4','NO2'/)

   call global_molecule_table%register('Methane','CH4')
   call global_molecule_table%register('Chlorophyll-a','C55H72MgN4O5')
   call global_molecule_table%register('Nitrate','NO3')
   call global_molecule_table%register('Nitrite','NO2')
   call global_molecule_table%register('Ammonium','NH4')
   call global_molecule_table%register('Phosphate','PO4')
   
   call self%din%create('DIN',din_molecules)
   call self%nutrient%create('Nutrient',(/'NO3','NH4','NO2','PO4'/))

   !C55H72MgN4O5
   !Molar mass	893.509 g·mol−1

   ! Register state variables
   call self%register_state_variable(self%id_chlorophyll,'chl','', '')

   !write(stderr, *) global_molecule_table%get('Chlorophyll-a')

end subroutine initialize

end module chemistry