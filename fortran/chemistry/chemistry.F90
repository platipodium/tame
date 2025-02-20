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

      type(type_state_variable_id), allocatable :: id_var(:)
      type(type_group), allocatable :: state_variables(:)

   contains
      procedure :: initialize
      procedure :: do
   end type

contains

subroutine initialize(self,configunit)

   class (type_chlorophyll), intent(inout), target :: self
   integer,            intent(in)            :: configunit

   integer :: i, n

   call global_molecule_table%register('Methane','CH4')
   call global_molecule_table%register('Chlorophyll-a','C55H72MgN4O5')
   call global_molecule_table%register('Nitrate','NO3')
   call global_molecule_table%register('Nitrite','NO2')
   call global_molecule_table%register('Ammonium','NH4')
   call global_molecule_table%register('Phosphate','PO4')

   ! Define here the (number of) groups that you would like to
   ! integrate as interior state variables.  Groups are
   ! single molecules or molecule
   ! groups consisting of molecules registered in the global_molecule_table
   allocate(self%state_variables(2))
   call self%state_variables(1)%create('DIN',(/'NO3','NH4','NO2'/))
   call self%state_variables(2)%create('PO4',(/'PO4'/))

   ! Define here the list of elements that should be considered as state
   ! variables
   !allocate(elements(3))
   !elements = (/'C','N','P'/)

   n = ubound(self%state_variables,1)

   allocate(self%id_var(n))
   do i = 1,n
      call self%register_state_variable(self%id_var(i), &
        self%state_variables(i)%name,'dummy unit','dummy long name')
   enddo

end subroutine initialize

subroutine do(self, _ARGUMENTS_DO_)
   _DECLARE_ARGUMENTS_DO_
   class (type_chlorophyll), intent(in) :: self

   integer       :: i, n
   real(kind=rk), allocatable :: rhs(:)
   real(kind=rk), allocatable :: state_variables(:)

   n = ubound(self%id_var,1)
   if (n<1) return

   allocate(rhs(n))
   allocate(state_variables(n))

   _LOOP_BEGIN_

      ! Get state variables into local variable
      do i = 1, n
         _GET_(self%id_var(i), state_variables(i))
      enddo

      !do i = 1, num_elements ! e.g., N  ( C, Si, Fe, P)
      !_GET_(self%id_var(det_index(i)), det_element(i))  ! Detritus Organics in mmol-C/m**3
      !_GET_(self%id_var(dom_index(i)), dom_element(i))  ! Dissolved Organics in mmol-C/m**3
      !end do

   ! here, nutrients are only remineralised (e.g., uptake in tame_phy)
   !do i = 1,num_chemicals
   !   rhs(dix_index(i)) = remin_chemical(i) !+ nut_prod(i)
   !end do

   ! Communicate right hand sides to FABM
   do i = 1, n
      _ADD_SOURCE_(self%id_var(i), rhs(i))
   end do

   deallocate(rhs)
   _LOOP_END_
end subroutine do

end module chemistry
