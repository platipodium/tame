! SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
!
! SPDX-License-Identifier: Apache-2.0

module ecosystem_model
   use fabm_pelagic
   use organism_module
   implicit none

   type, extends(type_pelagic) :: type_ecosystem
      real :: phyto_biomass   ! Phytoplankton biomass
      real, allocatable :: phyto_stoichiometry(:)  ! Phytoplankton stoichiometry

      real :: zoo_biomass     ! Zooplankton biomass
      real, allocatable :: zoo_stoichiometry(:)    ! Zooplankton stoichiometry

      real, allocatable :: nutrients(:)            ! Nutrients

   contains
      procedure :: initialize => initialize_ecosystem
      procedure :: update => update_ecosystem
   end type type_ecosystem

contains

   subroutine initialize_ecosystem(self, model, config_file)
      class(type_ecosystem), intent(inout) :: self
      type(type_model), intent(in) :: model
      character(len=*), intent(in) :: config_file

      ! Number of stoichiometric components
      integer, parameter :: num_components = 10

      ! Allocate arrays
      allocate (self%phyto_stoichiometry(num_components))
      allocate (self%zoo_stoichiometry(num_components))
      allocate (self%nutrients(num_components))

      ! Initialize state variables
      self%phyto_biomass = 1.0
      self%phyto_stoichiometry = [0.5, 0.1, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

      self%zoo_biomass = 0.5
      self%zoo_stoichiometry = [0.25, 0.05, 0.005, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

      self%nutrients = [1.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

      call self%initialize_pelagic(model, config_file)
   end subroutine initialize_ecosystem

   subroutine update_ecosystem(self, model, dt)
      class(type_ecosystem), intent(inout) :: self
      type(type_model), intent(in) :: model
      real, intent(in) :: dt

      ! Update processes
      call phytoplankton_growth(self, dt)
      call zooplankton_grazing(self, dt)
      call nutrient_uptake(self, dt)

      call self%update_pelagic(model, dt)
   end subroutine update_ecosystem

   subroutine phytoplankton_growth(self, dt)
      class(type_ecosystem), intent(inout) :: self
      real, intent(in) :: dt

      ! Example: Simple linear growth
      real :: growth_rate
      growth_rate = 0.1
      self%phyto_biomass = self%phyto_biomass + growth_rate*self%phyto_biomass*dt

      ! Update stoichiometry
      self%phyto_stoichiometry = self%phyto_stoichiometry + 0.01*growth_rate*self%phyto_biomass*dt
   end subroutine phytoplankton_growth

   subroutine zooplankton_grazing(self, dt)
      class(type_ecosystem), intent(inout) :: self
      real, intent(in) :: dt

      ! Example: Simple grazing on phytoplankton
      real :: grazing_rate
      grazing_rate = 0.05
      self%zoo_biomass = self%zoo_biomass + grazing_rate*self%phyto_biomass*dt
      self%phyto_biomass = self%phyto_biomass - grazing_rate*self%phyto_biomass*dt

      ! Update stoichiometry
      self%zoo_stoichiometry = self%zoo_stoichiometry + 0.01*grazing_rate*self%phyto_stoichiometry*dt
   end subroutine zooplankton_grazing

   subroutine nutrient_uptake(self, dt)
      class(type_ecosystem), intent(inout) :: self
      real, intent(in) :: dt

      ! Example: Nutrient uptake by phytoplankton
      real :: uptake_rate
      integer :: i
      uptake_rate = 0.05

      do i = 1, size(self%nutrients)
         self%nutrients(i) = self%nutrients(i) - uptake_rate*self%phyto_biomass*dt
         self%phyto_stoichiometry(i) = self%phyto_stoichiometry(i) + uptake_rate*self%phyto_biomass*dt
      end do
   end subroutine nutrient_uptake

end module ecosystem_model
