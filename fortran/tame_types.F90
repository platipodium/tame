!> @file tame_types.F90
!> @brief tame_types module
#include "fabm_driver.h"
!> @brief  Data types used in fabm_hzg_tame are defined here
module tame_types
use fabm_types
use fabm_expressions
public type_tame_sensitivities, type_tame_om

! standard fabm model types
type,extends(type_base_model),public :: type_tame_base_model
type (type_tame_nutindex) :: nutind

! standard fabm model types
type (type_state_variable_id)        :: id_nutN,id_nutP,id_nutS,id_detC,id_detN,id_detP,id_detS,id_domC,id_domN,id_domP,id_RNit,id_nh3,id_oxy,id_odu
type (type_dependency_id)            :: id_temp, id_sal,id_par,id_CO2,id_attpar
type (type_global_dependency_id)     :: id_doy
type (type_horizontal_dependency_id) :: id_lat, id_lon

type (type_dependency_id)            :: id_totC, id_totN, id_totP
type (type_horizontal_dependency_id) :: id_totC_vertint, id_totN_vertint, id_totP_vertint
type (type_horizontal_diagnostic_variable_id)  :: id_totC_vertint_diag,id_totN_vertint_diag,id_totP_vertint_diag
type (type_horizontal_dependency_id) :: id_zmax, id_o2flux, id_oduflux, id_tke_bot
type (type_horizontal_diagnostic_variable_id)            :: id_O2flux_diag
type (type_diagnostic_variable_id)   :: id_pPads, id_datt, id_vphys, id_GPPR, id_Denitr, id_dPAR, id_chl2C, id_Theta, id_fracR, id_fracT, id_fracNU, id_DNP, id_QNP, id_QN, id_QP, id_QSi, id_aVN, id_aVP, id_aVSi, id_faN, id_faP, id_faSi, id_rQN, id_rQP, id_rQSi, id_tmp, id_fac4, id_fac5, id_fac3, id_fac1, id_fac2, id_phyUR, id_phyRER, id_phyELR, id_phyALR, id_phyVLR, id_phyGLR, id_vsinkr, id_zoomort, id_qualPOM, id_qualDOM, id_no3, id_bgatt
real(rk) ::  nutN_initial, nutP_initial, nutS_initial, detC_initial, detN_initial, detP_initial, detS_initial, domC_initial, domN_initial, domP_initial, RNit_initial, nh3_initial, oxy_initial, odu_initial
real(rk) ::  a_water, a_minfr, a_spm, a_phyc, a_doc, a_fz, a_chl, rel_co2, frac_PAR, small, maxVal, dil, ex_airsea, O2_sat, N_depo, P_depo
real(rk) ::  rq10 
logical  ::  PhosphorusOn, SiliconOn, GrazingOn, BioOxyOn, DebugDiagOn, Budget0DDiagOn, Budget2DDiagOn, BGC0DDiagOn, BGC2DDiagOn, ChemostatOn, SwitchOn, NResOn
logical  ::  detritus_no_river_dilution, plankton_no_river_dilution, nutrient_no_river_dilution
integer  ::  genMeth
end type type_tame_base_model

type type_tame_nutindex
   integer :: iN, iP, iSi
   integer :: nutnum, nhi
   integer :: nfV, nSRN
end type

!
type type_tame_env
 real(rk) :: RNit, nh3, oxy, odu, temp,par,CO2,attpar,doy !vphys_dep,GPPR_dep,GPPR_vertint,GPPR_vertint_diag,Denitr_dep,Denitr_vertint,Denitr_vertint_diag,zmax,O2flux_diag
end type
type type_tame_rhs
 real(rk) :: nutC,nutN,nutP,nutS,detC,detN,detP,detS,domC,domN,domP,RNit
end type
! standard fabm model types
!#SP#

!!-------------------------------------------------------------------
type type_tame_switch
   logical :: isP, isSi
end type

type type_tame_elem
   real(rk) :: C,N,P,Si,Fe
end type

type,extends(type_tame_elem) type_tame_om
   logical :: IsParticulate 
   real(rk) :: elem(10)
end type

type,extends(type_tame_om) :: type_tame_life
      type (type_tame_om)  :: Q              ! quotaformer QN,QP,QPN
      real(rk)   :: odu                   ! degradation rate
end type type_tame_life

type type_tame_sensitivities
   real(rk) :: f_T       ! temperature dependency of metabolic rates
!   type (type_tame_om) :: upt_pot ! potential uptake rates
				   ! depending on ambient concentration incl. light limitation [dimensionless]
end type type_tame_sensitivities

! new meta structure for pointing/looping over elements

end module
