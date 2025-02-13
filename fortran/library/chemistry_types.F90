! SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
! SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de>
! SPDX-FileContributor: Kai Wirtz <kai.wirt@hereon.de>
! SPDX-License-Identifier: Apache-2.0
!
module chemistry_types

  implicit none

  type type_molecule
    type(type_element), allocatable :: elements(:)
    integer, allocatable :: stoichiometry(:)
    real(kind=4) :: weight
    contains
    procedure :: calculate_weight
  end type type_molecule

  type type_element
    character(len=20) :: name
    character(len=2)  :: symbol
    real(kind=4) :: weight, electronegativity
    integer(kind=2) :: number
    contains
    procedure :: register
  end type type_element

  type type_element_table
    class(type_element), allocatable :: elements(:)
    contains
    procedure :: load
  end type type_element_table

  contains

  function calculate_weight(self) result(w)

    class(type_molecule), intent(inout) :: self
    real(kind=4) :: w
    integer :: i

  w = 0.0

  do i=lbound(self%elements,1), ubound(self%elements,1)
    w = w + self%stoichiometry(i) * self%elements(i)%weight
  enddo

end function calculate_weight

subroutine register(self, number, name, symbol, weight, electronegativity)
    class(type_element), intent(inout) :: self
    character(len=*), intent(in) :: name, symbol
    real(kind=4), intent(in) :: weight, electronegativity
    integer(kind=2), intent(in) :: number

    self%number = number
    self%name = name
    self%symbol = symbol
    self%weight = weight
    self%electronegativity = electronegativity
end subroutine register

subroutine load(self)
  use yaml_types, only : type_node
  use yaml_settings, only : type_settings
  use yaml, only : parse
  use iso_fortran_env, only : output_unit

  class(type_element_table), intent(inout) :: self

  ! Later we want to read from a yaml
  ! class (type_node), pointer :: yaml_node
  ! class (type_settings), pointer :: yaml_settings, child

  !call yaml_settings%load('elements.yaml', unit=100)
  ! child => yaml_settings%get_child('elements')
   !do while (associated(child))
   !enddo
    !call child%dump(unit=output_unit)
   !call yaml_settings%finalize()
   !deallocate(yaml_settings)

   ! Do this statically here for now
   allocate(self%elements(2))
   !self%elements(1)%register('Carbon','C',12.011,2.55)
   !self%elements(2)%register(7,'Nitrogen','N',14.007,3.04)

end subroutine load

end module chemistry_types
