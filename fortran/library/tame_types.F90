!> @file tame_types.F90
!> @brief tame_types module
#include "fabm_driver.h"
! converts biological unit d-1 into physical FABM/driver unit s-1 for RHS
#define UNIT *1.1574074074E-5_rk  

!> @brief  Data types used in fabm_hzg_tame are defined here
module tame_types
use fabm_types
use fabm_expressions

public type_tame_sensitivities, type_tame_om, type_tame_elem, small

real(rk),parameter :: small = 1.E-4_rk
integer, parameter :: NUM_ELEM = 3
integer, parameter :: NUM_CHEM = 3
integer, parameter :: NUM_NUTRIENT = 2
character(len=3) :: chemicals(NUM_CHEM) = (/'NO3','NH4','PO4'/)
character(len=3) :: uptake_chemicals(NUM_NUTRIENT) = (/'DIN','PO4'/)
character(len=3) ::  ElementList= 'CNP'!SF'
character(len=10) ::  ElementName(NUM_ELEM)= (/'Carbon    ','Nitrogen  ','Phosphorus'/)
real(rk), parameter :: stoichiometry(NUM_ELEM) = (/ 1._rk, 1._rk/16_rk, 1._rk/106_rk /)! Redfield ratio C-based
real(rk),parameter :: secs_per_day = 86400.0_rk
real(rk),parameter :: days_per_sec = 1.0_rk/secs_per_day

! standard fabm model types
type,extends(type_base_model),public :: type_tame_base_model

! standard fabm model types
!type (type_global_dependency_id)     :: id_doy
!type (type_horizontal_dependency_id) :: id_lat, id_lon
!type (type_dependency_id)            :: id_totC, id_totN, id_totP
!type (type_horizontal_dependency_id) :: id_totC_vertint, id_totN_vertint, id_totP_vertint
!type (type_horizontal_diagnostic_variable_id)  :: id_totC_vertint_diag,id_totN_vertint_diag,id_totP_vertint_diag
!type (type_horizontal_dependency_id) :: id_zmax, id_o2flux, id_oduflux, id_tke_bot
!real(rk) :: remineral,hydrolysis,alloc_N,Nqual,CNref,DenitKno3,denit,T_ref,rq10
!logical  ::  detritus_no_river_dilution, plankton_no_river_dilution, nutrient_no_river_dilution
!integer  ::  tlim
end type type_tame_base_model

!
!!-------------------------------------------------------------------
type type_tame_env
 real(rk) :: temp,par,doy !RNit, nh3, oxy, odu, ,CO2,attpar vphys_dep,GPPR_dep,GPPR_vertint,GPPR_vertint_diag,Denitr_dep,Denitr_vertint,Denitr_vertint_diag,zmax,O2flux_diag
end type

type type_tame_sensitivities
   real(rk) :: f_T       ! temperature dependency of metabolic rates
!   type (type_tame_om) :: upt_pot ! potential uptake rates
				   ! depending on ambient concentration incl. light limitation [dimensionless]
end type type_tame_sensitivities

type type_tame_chemical_index
   integer  :: NO3,NH4,PO4,CO2,O2,SiO2,FeS,DIN,DIP,DISi,DIC
end type
! new meta structure for pointing/looping over chemicals (DIX)
type type_tame_chemical
   real(rk),pointer :: NO3,NH4,PO4,CO2,O2,SiO2,FeS,DIN,DIP,DISi,DIC
   real(rk) :: chemical(10)
   type(type_tame_chemical_index) :: index
!   integer  :: index(10)
end type

type type_tame_elem_index
   integer  :: C,N,P,Si,Fe
end type

type type_tame_elem
    real(rk), pointer :: C,N,P,Si,Fe
!    real(rk) :: element(5)
    type(type_tame_elem_index) :: index
end type

type,extends(type_tame_elem) :: type_tame_om
   logical :: IsParticulate 
!   real(rk) :: elem(10)
end type

type,extends(type_tame_om) :: type_tame_life
      type (type_tame_om)  :: Q              ! quotaformer QN,QP,QPN
      real(rk)   :: odu                   ! degradation rate
end type type_tame_life

end module
