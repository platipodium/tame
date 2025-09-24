! SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
! SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de>
! SPDX-FileContributor: Kai Wirtz <kai.wirtz@hereon.de>
! SPDX-License-Identifier: Apache-2.0
#include "fabm_driver.h"

module tame_chemistry_detritus
   use fabm_types
   use chemistry_types
   use tame_functions

   implicit none

   private

   type, extends(type_base_model),public :: type_tame_chemistry_detritus
      ! Variable identifiers
      type (type_state_variable_id), allocatable     :: id_detritus(:)
      type (type_state_variable_id), allocatable     :: id_target(:)

      type (type_group), allocatable :: state_variables(:)
      type (type_dependency_id) :: id_temp, id_par

      ! Model parameters
      real(rk) :: remineral,hydrolysis,alloc_N,Nqual,CNref,DenitKNO3,denit,T_ref,rq10,dil, tlim
 
      real(rk) :: rdn
   contains
      procedure :: initialize
      procedure :: do
   end type type_tame_chemistry_detritus

contains

   subroutine initialize(self, configunit)
      class (type_tame_chemistry_detritus), intent(inout), target :: self
      integer,                        intent(in)            :: configunit

      real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
      real(rk)            :: w_d, kc
      integer             :: i,n

      call global_molecule_table%register('det_fast','C55H62N4')
      call global_molecule_table%register('det_slow','C55H12PO4')
      call global_molecule_table%register('nitrate','NO3')
      call global_molecule_table%register('nitrite','NO2')
      call global_molecule_table%register('ammonia','NH4')
      call global_molecule_table%register('phosphate','PO4')
      
      ! chemicals(NUM_CHEM)  = (/'NO3','NH4','PO4'/)

      ! Define here the (number of) groups that you would like to
      ! integrate as interior state variables.  Groups are
      ! single molecules or molecule
      ! groups consisting of molecules registered in the global_molecule_table
      n = 4
      allocate(self%state_variables(n))
      call self%state_variables(1)%create('fast',(/'det_fast'/))
      call self%state_variables(2)%create('slow',(/'det_slow'/))
      call self%state_variables(3)%create('DIN',(/'nitrate','nitrite','ammonia'/))
      call self%state_variables(4)%create('PO4',(/'phosphate'/))
      

      allocate(self%id_detritus(n))
      allocate(self%id_target(n))

      ! Store parameter values in our own derived type
      ! NB: all rates must be provided in values per day and are converted here to values per second.
      call self%get_parameter(w_d,      'w_d', 'm d-1',     'vertical velocity (<0 for sinking)', default=-5.0_rk, scale_factor=d_per_s)
      call self%get_parameter(kc,       'kc',  'm2 mmol-1', 'specific light extinction',          default=0.03_rk)
      call self%get_parameter(self%rdn, 'rdn', 'd-1',       'remineralization rate',              default=0.003_rk, scale_factor=d_per_s)
      call self%get_parameter(self%remineral, 'remineral','d-1','DOM remineralisation rate', default=0.1_rk )
      call self%get_parameter(self%hydrolysis, 'hydrolysis','d-1','detritus hydrolysis rate', default=0.05_rk )
      call self%get_parameter(self%alloc_N, 'alloc_N','-','nh4 - NO3 product ratio remineralisation', default=0.5_rk )
      call self%get_parameter(self%Nqual, 'Nqual','-','OM fraction w quality prop to N:Cratio ', default=1.0_rk )
      call self%get_parameter(self%CNref, 'CNref','Redfield','POM quality relative to carbon : nitrogen ratio (mol C/mol N)', default=6.625_rk )
      call self%get_parameter(self%DenitKNO3, 'DenitKNO3','mmol N/m3','half-saturation NO3 denitrification', default=1.0_rk )
      call self%get_parameter(self%denit, 'denit','1/d','pelagic denitrification rate', default=0.01_rk )
      call self%get_parameter(self%T_ref, 'T_ref','Kelvin','reference temperature', default=293.0_rk )
      call self%get_parameter(self%rq10, 'rq10','-','temperature dependence Q10', default=0.175_rk )
      !call self%get_parameter(self%tlim, 'tlim','0: none, 1: flagellate-style, 2: cyanobacteria-style','temperature limitation of growth', default=0 )

      do i=1,n
         call self%register_state_variable(self%id_detritus(i), trim(self%state_variables(i)%name),'mmol m-3',  'concentration', 4.5_rk, &
            minimum=0.0_rk, vertical_movement=w_d, specific_light_extinction=kc)
      enddo

      ! Register contribution of state to global aggregate variables.
      do i=1,n
        ! @todo
        !if self%state_variables(i)%contains('N') call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_d)
      enddo

      call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
      call self%register_dependency(self%id_temp, standard_variables%temperature)

      do i=1,n
         call self%register_state_dependency(self%id_target(i), trim(self%state_variables(i)%name)//'_target','mmol m-3', 'sink for remineralized matter')
      enddo
   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_tame_chemistry_detritus), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: d
      real(rk), allocatable :: concentration(:)
      integer :: i,n

      n = ubound(self%state_variables,1)
      allocate(concentration(n))

      _LOOP_BEGIN_

         do i=1,n
            _GET_(self%id_detritus(i), concentration(i)) ! detritus

            _ADD_SOURCE_(self%id_detritus(i), -self%rdn*concentration(i))
            _ADD_SOURCE_(self%id_target(i), self%rdn*concentration(i))
         enddo

      ! Leave spatial loops (if any)
      _LOOP_END_

   end subroutine do

end module tame_chemistry_detritus
