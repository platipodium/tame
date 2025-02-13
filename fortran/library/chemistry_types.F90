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
    character(len=20) :: name
    contains
    procedure :: calculate_weight
    procedure :: create
    procedure :: decompose
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
    class(type_element), pointer :: elements(:)
    contains
    procedure :: load
    procedure :: clear
    procedure :: get

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
    integer(kind=4), intent(in) :: number

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
   call self%elements(1)%register(1,'Hydrogen','H',1.008,2.20)
   call self%elements(2)%register(6,'Carbon','C',12.011,2.55)
   call self%elements(3)%register(7,'Nitrogen','N',14.007,3.04)
   call self%elements(4)%register(8,'Oxygen','O',15.999,3.44)
   call self%elements(5)%register(14,'Silicon','Si',28.085,1.90)
   call self%elements(6)%register(15,'Phosphorus','P',30.974,2.19)

end subroutine load

subroutine clear(self)
  class(type_element_table), intent(inout) :: self
  if (associated(self%elements)) deallocate(self%elements)
end subroutine clear


function get(self, symbol) result(elem)
  class(type_element_table), intent(inout) :: self
  character(len=2), intent(in) :: symbol
  class(type_element), pointer :: elem
  integer :: i

  do i=1, ubound(self%elements,1)
    if (len_trim(self%elements(i)%symbol) /= len_trim(symbol)) cycle
    elem => self%elements(i)
  enddo
end function get

subroutine create(self, name, composition)

  class(type_molecule), intent(inout) :: self
  character(len=*) :: name, composition

  self%name = name
  call self%decompose(composition)

end subroutine create

subroutine decompose(self, composition)
  class(type_molecule), intent(inout) :: self
  character(len=*) :: composition
  integer :: i,j,k,n,length

  character(len=26), parameter :: majuscules = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  character(len=26), parameter :: minuscules = 'abcdefghijklmnopqrstuvwxyz'
  character(len=10), parameter :: numbers = '0123456789'
  character(len=20) :: number

  ! Count the number of uppercase letters to determine number of elements
  length = len_trim(composition)
  n = 0
  do i = 1, length
    if (index(majuscules,composition(i:i)) < 1) cycle
    n = n + 1
  enddo

  if (.not.allocated(self%elements)) allocate(self%elements(n))
  if (.not.allocated(self%stoichiometry)) allocate(self%stoichiometry(n))

  n = 0
  do i = 1, length
    if (index(majuscules,composition(i:i)) < 1) cycle 
    n = n + 1
    ! If there is a minuscule letter following a capital letter, we have a two-letter symbol
    if (i < length .and. index(minuscules,composition(i+1:i+1)) > 0) then 
      !self%elements(n) = composition(i:i+1)
    else 
      !self%elements(n) = composition(i:i)
    endif

    ! Check whether there are any numbers before the next majuscule
    j = index(majuscules,composition(i:length))
    k = index(numbers,composition(i:length))
    if (k > 0) then 
      number = composition(i+k:length)
      if (j > 0 .and. k < j - 1) number = composition(i+k:i+j-1)
      read(number, '(I3)') self%stoichiometry(n)
    else 
      self%stoichiometry(n) = 1
    endif
  end do

end subroutine decompose

end module chemistry_types
