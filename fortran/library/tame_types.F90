!> @file tame_types.F90
!> @brief tame_types module
#include "fabm_driver.h"
!> @brief  Data types used in fabm_hzg_tame are defined here
module tame_types
use fabm_types
use fabm_expressions
public type_tame_sensitivities, type_tame_om, type_tame_elem

!  meta structure for pointing/looping over elements

type,extends(type_state_variable_id) :: type_tame_elem_id
   real(rk) :: C,N,P,Si,Fe
end type

! standard fabm model types
type,extends(type_base_model),public :: type_tame_base_model

! standard fabm model types
!type (type_state_variable_id)        :: id_nutN,id_nutP,id_nutS,id_detC,id_detN,id_detP,id_detS,id_domC,id_domN,id_domP,id_RNit,id_nh3,id_oxy,id_odu
!type (type_global_dependency_id)     :: id_doy
!type (type_horizontal_dependency_id) :: id_lat, id_lon
!!real(kind=rk), allocatable, target ::  id_dix(:)
!!type (type_tame_elem_index)  :: Index_Det,Index_DOM

!type (type_dependency_id)            :: id_totC, id_totN, id_totP
!type (type_horizontal_dependency_id) :: id_totC_vertint, id_totN_vertint, id_totP_vertint
!type (type_horizontal_diagnostic_variable_id)  :: id_totC_vertint_diag,id_totN_vertint_diag,id_totP_vertint_diag
!type (type_horizontal_dependency_id) :: id_zmax, id_o2flux, id_oduflux, id_tke_bot
!real(rk) :: remineral,hydrolysis,alloc_N,Nqual,CNref,DenitKno3,denit,T_ref,rq10
!real(rk) ::  nutN_initial, nutP_initial, nutS_initial, detC_initial, detN_initial, detP_initial, detS_initial, domC_initial, domN_initial, domP_initial, RNit_initial, nh3_initial, oxy_initial, odu_initial
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

! new meta structure for pointing/looping over chemicals (DIX)
type type_tame_chemical
   real(rk),pointer :: no3,nh4,po4,co2,o2,sio4,FeS,din,dip,dis,dic
   real(rk) :: chemical(10)
end type

type type_tame_elem_index
   integer  :: C,N,P,Si,Fe
end type

type type_tame_elem
   real(rk),pointer :: C,N,P,Si,Fe
   real(rk) :: element(5)
   type (type_tame_elem_index)  ::  index  
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
