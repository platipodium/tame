! SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
! SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de>
! SPDX-License-Identifier: Apache-2.0
!
module chemistry_types

  use iso_fortran_env, only : stdout => output_unit, stderr => error_unit
  implicit none
  private

  type type_element
    character(len=20) :: name
    character(len=2)  :: symbol
    real(kind=4)      :: weight
    real(kind=4)      :: electronegativity ! What scale is this on?
    integer(kind=2)   :: number
  contains
    procedure :: element_equals
    generic :: operator(==) => element_equals
    procedure :: element_assign
    generic :: assignment(=) => element_assign
    procedure :: element_write_formatted
    generic :: write(formatted) => element_write_formatted
    procedure :: element_write_unformatted
    generic :: write(unformatted) => element_write_unformatted

  end type type_element

  type type_elementPtr
    type(type_element), pointer :: element
  end type type_elementPtr

  ! Type molecule represents a chemical molecule.  A molecule has a name, a
  ! list of constituent chemical elements and the stoichiometry for these
  ! elements.  To prevent storing multiple copies of elements, the elements
  ! are represented by an array of pointers to chemical elements residing in
  ! the global element table.
  type type_molecule
    type(type_elementPtr), allocatable :: elementPtr(:)
    integer, allocatable :: stoichiometry(:)
    real(kind=4) :: weight
    character(len=20) :: name
    !character(len=20) :: formula
    contains
    procedure :: calculate_weight
    !procedure :: calculate_formula
    procedure :: decompose
    procedure :: assign => molecule_assign
    generic :: assignment(=) => assign
    procedure :: clear => molecule_clear
    procedure :: molecule_write_unformatted
    generic :: write(unformatted) => molecule_write_unformatted
    procedure :: molecule_write_formatted
    generic :: write(formatted) => molecule_write_formatted
 end type type_molecule

  ! The type element table contains an pointer to an array of elements.
  ! We choose a pointer instead of an allocatable array to make this array
  ! a target for pointing to it outside the element table object.
  type type_element_table
    type(type_element), pointer :: elements(:)
    contains
    procedure :: load
    procedure :: clear
    procedure :: get
    procedure :: register => table_register_element
    procedure :: contains => table_contains_symbol
    procedure :: size => element_table_size
    procedure :: dump => element_table_dump

  end type type_element_table

  ! The type molecule table serves as a storage for molecules,
  ! making sure that molecules only live once in memory
  type type_molecule_table
    class(type_molecule), allocatable :: molecules(:)
    contains
    procedure :: register => table_register_molecule
    procedure :: contains => molecule_table_contains_name
    procedure :: size => molecule_table_size
    procedure :: dump => molecule_table_dump

  end type type_molecule_table

  ! The global element table is an instance of the element table
  ! holding the storage for all elements
  type(type_element_table), save, public, target :: global_element_table
  type(type_molecule_table), save, public, target :: global_molecule_table

contains

  logical function element_equals(a, b) result(is_equal)
    class(type_element), intent(in) :: a, b
    is_equal = (trim(a%name) == trim(b%name)) .and. &
             (trim(a%symbol) == trim(b%symbol)) .and. &
             (a%number == b%number)
  end function element_equals

  subroutine element_assign(target, source)
    class(type_element), intent(out) :: target
    class(type_element), intent(in)  :: source

    target%name = source%name
    target%symbol = source%symbol
    target%weight = source%weight
    target%electronegativity = source%electronegativity
    target%number = source%number
  end subroutine element_assign

  subroutine molecule_assign(target, source)
    class(type_molecule), intent(out) :: target
    class(type_molecule), intent(in)  :: source

    integer :: i, size

    target%name =  source%name
    size = ubound(source%stoichiometry,1)
    allocate(target%stoichiometry(size))
    allocate(target%elementPtr(size))
    target%stoichiometry =source%stoichiometry

    do i=1,size
      allocate(target%elementPtr(i)%element)
      ! target%elementPtr(j)%element%name = table%molecules(i)%elementPtr(j)%element%name
      ! target%elementPtr(j)%element%symbol = table%molecules(i)%elementPtr(j)%element%symbol
      ! target%elementPtr(j)%element%electronegativity = table%molecules(i)%elementPtr(j)%element%electronegativity
      ! target%elementPtr(j)%element%weight = table%molecules(i)%elementPtr(j)%element%weight
      target%elementPtr(i)%element = source%elementPtr(i)%element
    enddo
  end subroutine molecule_assign

  subroutine molecule_clear(molecule)
    class(type_molecule), intent(inout) :: molecule

    integer :: i, size

    size = ubound(molecule%stoichiometry,1)
    do i=1,size
      deallocate(molecule%elementPtr(i)%element)
    enddo
    deallocate(molecule%elementPtr)
    deallocate(molecule%stoichiometry)
  end subroutine molecule_clear

  subroutine element_write_formatted(element, unit, iotype, vlist, iostat, iomsg)
    class(type_element), intent(in) :: element
    integer, intent(in) :: unit
    character(len=*), intent(in) :: iotype
    integer, intent(in) :: vlist(:)
    integer, intent(out) :: iostat
    character(len=*), intent(inout) :: iomsg

    iostat = 0

    write(unit, '(A,1x,A,1x,I3,1x,F8.4,1x,F6.2)') trim(element%name), trim(element%symbol), &
    element%number, element%weight, element%electronegativity
end subroutine element_write_formatted

subroutine element_write_unformatted(element, unit, iostat, iomsg)
  class(type_element), intent(in) :: element
  integer, intent(in) :: unit
  integer, intent(out) :: iostat
  character(len=*), intent(inout) :: iomsg

  iostat = 0
  write(unit, *) trim(element%name), trim(element%symbol), element%number, element%weight, element%electronegativity
end subroutine element_write_unformatted


  function element_table_size(table) result(size)
    class(type_element_table), intent(in) :: table

    integer :: size

    size = 0
    if (associated(table%elements)) size = ubound(table%elements,1)
  end function element_table_size

function table_contains_number(table, number) result(contains)
  class(type_element_table), intent(in) :: table
  integer(kind=4), intent(in) :: number

  logical :: contains
  integer :: i

  contains = .false.
  do i=1,table%size()
    if (table%elements(i)%number /= number) cycle
    contains = .true.
  enddo
end function table_contains_number

function table_contains_symbol(table, symbol) result(contains)
    class(type_element_table), intent(in) :: table
    character(len=*), intent(in) :: symbol

    logical :: contains
    integer :: i

    contains = .false.
    do i=1,table%size()
      if (trim(table%elements(i)%symbol) /= symbol) cycle
      contains = .true.
    enddo
  end function table_contains_symbol

! This subroutine registers an element within the global element table
! only if an entry with the same symbol does not exist. It also creates
! the instance of this element.
subroutine table_register_element(table, number, name, symbol, weight, electronegativity)

  class(type_element_table), intent(inout) :: table
  character(len=*), intent(in) :: name, symbol
  real(kind=4), intent(in) :: weight, electronegativity
  integer(kind=4), intent(in) :: number

  type(type_element_table) :: temporary
  integer :: i, n

  if (table%contains(symbol)) return

  n = table%size() + 1

  allocate(temporary%elements(n))
  do i=1,n-1
    temporary%elements(i) = table%elements(i)
  enddo

  temporary%elements(n)%number =  number
  temporary%elements(n)%name = name
  temporary%elements(n)%symbol =  symbol
  temporary%elements(n)%weight =  weight
  temporary%elements(n)%electronegativity = electronegativity

  if (associated(table%elements)) deallocate(table%elements)
  allocate(table%elements(n))
  do i=1,n
    table%elements(i) =  temporary%elements(i)
  enddo
  deallocate(temporary%elements)

end subroutine table_register_element

subroutine element_table_dump(table)

  class(type_element_table), intent(in) :: table
  integer :: i, size

  size = table%size()
  if (size < 1) return

  write(stdout,*) '  .. global element table contains ',size,' entries'

  do i = 1, size
    write(stdout, *) '  ..',table%elements(i)
  enddo

end subroutine element_table_dump


function molecule_table_size(table) result(size)
  class(type_molecule_table), intent(in) :: table

  integer :: size

  size = 0
  if (allocated(table%molecules)) size = ubound(table%molecules,1)
end function molecule_table_size

function molecule_table_contains_name(table, name) result(contains)
  class(type_molecule_table), intent(in) :: table
  character(len=20), intent(in) :: name

  logical :: contains
  integer :: i

  contains = .false.
  do i=1,table%size()
    if (trim(table%molecules(i)%name) /= trim(name)) cycle
    contains = .true.
  enddo
end function molecule_table_contains_name


subroutine molecule_table_dump(table)

  class(type_molecule_table), intent(in) :: table
  integer :: i, size

  size = table%size()
  if (size < 1) return

  write(stdout,*) '  .. molecule element table contains ',size,' entries'

  do i = 1, size
    write(stdout, *) '  ..',table%molecules(i)
  enddo

end subroutine molecule_table_dump


! This subroutine registers a molecule within the global molecule table
! only if an entry with the same symbol does not exist. It also creates
! the instance of this molecule.
subroutine table_register_molecule(table, name, composition)

  class(type_molecule_table), intent(inout) :: table
  character(len=*), intent(in) :: name, composition

  type(type_molecule_table) :: temporary
  integer :: i, j, n, size

  if (global_element_table%size() < 1) then
    call global_element_table%load()
    call global_element_table%dump()
  endif

  write(stdout,*) '  .. searching for name "'//trim(name)//'" in global molecule table'
  if (table%contains(name)) return

  n = table%size() + 1
  write(stdout,*) '  .. item "'//trim(name)//'" not found in table of size ', n-1

  allocate(temporary%molecules(n))
  do i=1,n-1
    !write(stdout,*) '  .. saving "'//trim(table%molecules(i)%name)//'" to temporary table'
    temporary%molecules(i) = table%molecules(i)
  enddo
  temporary%molecules(n)%name = name
  call temporary%molecules(n)%decompose(composition)

  if (allocated(table%molecules)) then
    do i=1,n-1
      call table%molecules(i)%clear()
    enddo
    deallocate(table%molecules)
  endif

  allocate(table%molecules(n))
  do i=1,n
    !write(stdout,*) '  .. restoring molecule "'//trim(temporary%molecules(i)%name)//'" from temporary table'
    table%molecules(i) = temporary%molecules(i)
  enddo

  do i=1,n-1
    call temporary%molecules(i)%clear()
  enddo
  deallocate(temporary%molecules)

  call table%dump()

end subroutine table_register_molecule

subroutine load(table)
  !use yaml_types, only : type_node
  !use yaml_settings, only : type_settings
  !use yaml, only : parse
  use iso_fortran_env, only : output_unit

  class(type_element_table), intent(inout) :: table

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
   call table%register(1,'Hydrogen','H',1.008,2.20)
   call table%register(6,'Carbon','C',12.011,2.55)
   call table%register(7,'Nitrogen','N',14.007,3.04)
   call table%register(8,'Oxygen','O',15.999,3.44)
   call table%register(14,'Silicon','Si',28.085,1.90)
   call table%register(15,'Phosphorus','P',30.974,2.19)
   call table%register(12,'Magnesium','Mg',24.305,1.31)

end subroutine load

  function calculate_weight(molecule) result(w)

    class(type_molecule), intent(inout) :: molecule
    real(kind=4) :: w
    integer :: i

    w = 0.0

    do i=lbound(molecule%elementPtr,1), ubound(molecule%elementPtr,1)
      w = w + molecule%stoichiometry(i) * molecule%elementPtr(i)%element%weight
    enddo

end function calculate_weight

subroutine register(element, number, name, symbol, weight, electronegativity)
    class(type_element), intent(inout) :: element
    character(len=*), intent(in) :: name, symbol
    real(kind=4), intent(in) :: weight, electronegativity
    integer(kind=4), intent(in) :: number

    element%number = number
    element%name = name
    element%symbol = symbol
    element%weight = weight
    element%electronegativity = electronegativity
end subroutine register

subroutine clear(self)
  class(type_element_table), intent(inout) :: self
  if (associated(self%elements)) deallocate(self%elements)
end subroutine clear


function get(table, symbol) result(elem)
  class(type_element_table), intent(in) :: table
  character(len=2), intent(in) :: symbol
  type(type_element), pointer :: elem
  integer :: i

  nullify(elem)
  do i=1, ubound(table%elements,1)
    if (trim(table%elements(i)%symbol) /= trim(symbol)) cycle
    elem => table%elements(i)
  enddo
end function get

subroutine molecule_write_formatted(molecule, unit, iotype, vlist, iostat, iomsg)
  class(type_molecule), intent(in) :: molecule
  integer, intent(in) :: unit
  character(len=*), intent(in) :: iotype
  integer, intent(in) :: vlist(:)
  integer, intent(out) :: iostat
  character(len=*), intent(inout) :: iomsg

  integer :: i, size
  character(len=40) :: format

  iostat = 0
  size = ubound(molecule%stoichiometry,1)
  write(format, '(A,I2,A)') '(A,',size,'(X,A,I2))'
  write(unit, format) trim(molecule%name), (trim(molecule%elementPtr(i)%element%symbol), molecule%stoichiometry(i), i= 1,size)
end subroutine molecule_write_formatted

subroutine molecule_write_unformatted(molecule, unit, iostat, iomsg)
  class(type_molecule), intent(in) :: molecule
  integer, intent(in) :: unit
  integer, intent(out) :: iostat
  character(len=*), intent(inout) :: iomsg
  integer :: i, size

  iostat = 0
  size = ubound(molecule%stoichiometry,1)
  write(unit, *) molecule%name, (molecule%elementPtr(i)%element%symbol, molecule%stoichiometry(i), i= 1,size)
end subroutine molecule_write_unformatted

subroutine decompose(molecule, composition)
  class(type_molecule), intent(inout) :: molecule
  character(len=*) :: composition
  integer :: i,j,k,n,length

  character(len=26), parameter :: majuscules = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  character(len=26), parameter :: minuscules = 'abcdefghijklmnopqrstuvwxyz'
  character(len=10), parameter :: numbers = '0123456789'
  character(len=2) :: symbol
  type(type_element), pointer :: elementPtr

  ! Count the number of uppercase letters to determine number of elements
  length = len_trim(composition)
  n = 0
  do i = 1, length
    if (index(majuscules,composition(i:i)) < 1) cycle
    n = n + 1
  enddo

  if (.not.allocated(molecule%elementPtr)) allocate(molecule%elementPtr(n))
  if (.not.allocated(molecule%stoichiometry)) allocate(molecule%stoichiometry(n))

  write(stdout,'(A,I2,A)') '   .. molecule '//trim(molecule%name)//' contains ',n,' elements'

  n = 0
  do i = 1, length
    ! Cycle unless we detect a new element by its capital letter
    if (index(majuscules,composition(i:i)) < 1) cycle

    n = n + 1
    ! If there is a minuscule letter following a capital letter, we have a two-letter symbol
    if (i < length .and. index(minuscules,composition(i+1:i+1)) > 0) then
      symbol = composition(i:i+1)
    else
      symbol = composition(i:i)
    endif

    if (.not.global_element_table%contains(symbol)) then
      write(stderr, *) 'Global element table does not contain element with symbol '//trim(symbol)
      stop
    endif

    elementPtr => global_element_table%get(symbol)
    allocate(molecule%elementPtr(n)%element)
    molecule%elementPtr(n)%element = elementPtr
    !write(stderr, *) '  .. assigned pointer to  "'//trim(molecule%elementPtr(n)%element%name)//'" in global element table '
    ! Check whether there are any numbers up next
    j = i + 1

    do while (index(numbers,composition(j:j)) > 0)
      j = j + 1
    enddo

    if (i+1 > j - 1) then
      molecule%stoichiometry(n) = 1
    else
      read(composition(i+1:j-1), '(I3)') molecule%stoichiometry(n)
    endif
    !write(stderr, *) '  .. assigned stoichiometry '//trim(symbol), molecule%stoichiometry(n)

  end do
  write(stdout,*) '  .. decomposed ',molecule

end subroutine decompose

end module chemistry_types
