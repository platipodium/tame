! SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
!
! SPDX-License-Identifier: Apache-2.0

!> @file tame_types.F90
!> @brief tame_types module
#include "fabm_driver.h"

!> @brief  Data types used in fabm_hzg_tame are defined here
module tame_types
use fabm_types
use fabm_expressions
 implicit none

public type_tame_sensitivities, type_tame_om, type_tame_chemical
public type_tame_env, type_tame_elem
public secs_per_day, days_per_sec, small
public NUM_ELEM, NUM_CHEM, NUM_NUTRIENT, chemicals, ElementList, ElementName, fixed_stoichiometry, zoo_stoichiometry,dphyXdt_crit
public chem2elem,chem2nut,nutrient_name,nut2elem,elem2nut,nut_minval,nut2othernut,det_index, dom_index,num_chem_of_nut, share_nut_chemindex,tame_index_set

real(rk),parameter :: small = 1.E-6_rk
integer, parameter :: NUM_ELEM = 3
integer, parameter :: NUM_CHEM = 3
integer, parameter :: NUM_NUTRIENT = 2
integer, parameter :: NUM_PHY = 2
integer, parameter :: NUM_PREDATOR = 2
character(len=3)   :: chemicals(NUM_CHEM) = (/'NO3','NH4','PO4'/)
integer, parameter :: chem2elem(NUM_CHEM) = (/    2,    2,    3/) ! index of element for each chemical TODO: generalize for molecules of >1 resolved element
integer, parameter :: chem2nut(NUM_CHEM)     = (/ 1,    1,    2/)
character(len=3)   :: nutrient_name(NUM_NUTRIENT) = (/'DIN','PO4'/)
integer, parameter :: nut2elem(NUM_NUTRIENT)    = (/ 2,  3/)
integer, parameter :: elem2nut(NUM_ELEM)        = (/-1,  1,  2/)
integer, parameter :: nut2othernut(NUM_NUTRIENT)= (/ 2,  1/) ! complementary nutrient N -> P, P -> N ! TODO resolve more than two nutrients
real(rk), parameter:: nut_minval(NUM_NUTRIENT)=  (/0.6, 0.02/)
character(len=3) ::  ElementList= 'CNP'!SF'
character(len=10) ::  ElementName(NUM_ELEM)= (/'carbon    ','nitrogen  ','phosphorus'/)
real(rk), parameter :: fixed_stoichiometry(NUM_ELEM) = (/ 1._rk, 1._rk/16_rk, 1._rk/106_rk /)! Redfield ratio C-based
real(rk), parameter :: zoo_stoichiometry(NUM_ELEM) = (/ 1._rk, 1._rk/16_rk, 1._rk/106_rk /)! Redfield ratio C-based
real(rk), parameter :: dphyC_dt = 66.67_rk ! max Carbon change per day * secs_per_day/200.0_rk
real(rk) :: dphyXdt_crit(NUM_ELEM) = dphyC_dt*fixed_stoichiometry

! converts biological unit d-1 into physical FABM/driver unit s-1 for RHS
real(rk),parameter :: secs_per_day = 86400.0_rk
real(rk),parameter :: days_per_sec = 1.0_rk/secs_per_day

! C-based stoichiometry of all chemicals, so NO3-N to C, NH4-N to C, PO4-P to C
real(rk)            :: chem_stoichiometry(NUM_CHEM)=(/0.0625_rk, 0.0625_rk, 0.0094_rk/) ! Redfield TODO

integer :: det_index(NUM_ELEM), dom_index(NUM_ELEM), dix_index(NUM_CHEM)
integer :: num_chem_of_nut(NUM_NUTRIENT), share_nut_chemindex(NUM_NUTRIENT,NUM_CHEM)

! standard fabm model types
type,extends(type_base_model),public :: type_tame_base_model

!type (type_global_dependency_id)     :: id_doy
!type (type_horizontal_dependency_id) :: id_lat, id_lon
!type (type_dependency_id)            :: id_totC, id_totN, id_totP
!type (type_horizontal_dependency_id) :: id_totC_vertint, id_totN_vertint, id_totP_vertint
!type (type_horizontal_diagnostic_variable_id)  :: id_totC_vertint_diag,id_totN_vertint_diag,id_totP_vertint_diag
!type (type_horizontal_dependency_id) :: id_zmax, id_o2flux, id_oduflux, id_tke_bot
!logical  ::  detritus_no_river_dilution, plankton_no_river_dilution, nutrient_no_river_dilution
end type type_tame_base_model

!
!!-------------------------------------------------------------------
type type_tame_env
 real(rk) :: temp,par,doy !RNit, nh3, oxy, odu, ,CO2,attpar vphys_dep,GPPR_dep,GPPR_vertint,Denitr_vertint
end type

type type_tame_sensitivities
   real(rk) :: f_T       ! temperature dependency of metabolic rates
!   type (type_tame_om) :: upt_pot ! potential uptake rates
				   ! depending on ambient concentration incl. light limitation [dimensionless]
end type type_tame_sensitivities

type type_tame_chemical_index
   integer  :: NO3,NH4,PO4,DIN !,CO2,O2,SiO2,FeS,DIP,DISi,DIC
end type
! new meta structure for pointing/looping over chemicals (DIX)
type type_tame_chemical
   real(rk),pointer :: NO3,NH4,PO4,DIN !,CO2,O2,SiO2,FeS,DIP,DISi,DIC
   real(rk) :: chemical(10)
   type(type_tame_chemical_index) :: index
!   integer  :: index(10)
end type

type type_tame_elem_index
   integer  :: C,N,P !,Si,Fe
end type

type type_tame_elem
    real(rk), pointer :: C,N,P !,Si,Fe
    type(type_tame_elem_index) :: index
end type

type,extends(type_tame_elem) :: type_tame_om
   logical :: IsParticulate
!   real(rk) :: elem(10)
end type

contains

subroutine tame_index_set()
implicit none
integer :: i, j

! ==== partitioning of chemical-nutrient uptake by phy  ======
num_chem_of_nut = 0 
do i = 1, NUM_CHEM
   j = chem2nut(i)
   num_chem_of_nut(j) = num_chem_of_nut(j) + 1 
   share_nut_chemindex(j,num_chem_of_nut(j)) = i 
end do

end subroutine tame_index_set

end module tame_types