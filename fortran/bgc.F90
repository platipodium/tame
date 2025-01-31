!#include "fabm_driver.h"
!> @brief This is the main BGC routine where right-hand-sides are calculated
!> ** BioGeoChemical Equations **
!! for organic material (Det, DOM) with arbitrary element units such carbon, nitrogen, & phosphorus
!!  and arbitrary dissolved chemical species such as NO3, NH4, or CO2 as state variables.

subroutine bgc(self,_ARGUMENTS_DO_)

use fabm_types
use tame_types
use tame_functions 

! !INPUT PARAMETERS:
 class (type_hzg_tame),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
type (type_tame_rhs)    :: rhsv
type (type_tame_om)     :: dom, det, nut
type (type_tame_env)    :: env
type (type_tame_switch) :: mswitch
type (type_tame_sensitivities) :: sens
!type (stoich_pointer), dimension(5)::elem ! struct-pointer addressing elements wthin loops

! --- LOCAL MODEL VARIABLES:
integer  :: i,  ihour, iloop
real(rk) :: remineralisation , hydrolysis       ! Temp dependent remineralisation and hydrolysis rates
! --- QUOTA and FRACTIONS

! --- AGGREGATION
real(rk) :: aggreg_rate ! particle aggregation 
logical  :: out = .true.
!   if(36000.eq.secondsofday .and. mod(julianday,1).eq.0 .and. outn) out=.true.
real(rk) :: no3
real(rk),parameter :: relaxO2 = 0.04_rk
#define _DEBUG_ 0
#define UNIT *1.1574074074E-5_rk ! 1/86400
 _LOOP_BEGIN_

! First retrieve current (local) state  variable values
!#S_GET
!---------- GET for each state variable ----------
do i = 1,num_stoich ! e.g., CO2, NO3, NH4 (PO4)
!  if (_AVAILABLE_(self%id_dix(i))) then
   _GET_(self%id_dix(i), dix%stoich(i))  ! Dissolved Inorganic Nutrient DIX in mmol-X/m**3
!  end if
end do
do i = 1,num_elements ! e.g., N  ( C, Si, Fe, P)
  _GET_(self%id_det(i), det%elem(i))  ! Detritus Organics in mmol-C/m**3
  _GET_(self%id_dom(i), dom%elem(i))  ! Dissolved Organics in mmol-C/m**3
end do
!do i = 1,num_phyclass ! e.g., DIA, FLA, CYA  (or 1,5,30,100 µm)
!  _GET_(self%id_phy(i), phy%class(i))  ! Detritus Organics in mmol-C/m**3
!end do
!---------- GET ambient physico-stoich conditions ----------
_GET_(self%id_temp, env%temp)  ! water temperature
! _GET_(self%id_par, env%par)    ! light photosynthetically active radiation

!_SET_DIAGNOSTIC_(self%id_PAR_diag,env%par)         !average Photosynthetically_Active_Radiation_

call calc_sensitivities(self,sens,env)
! call photosynthesis(self,sens,phy,nut,uptake,exud,acclim)

!___________________________________________________________________
!  ---  POM&DOM quality, relative to Refield ?
!  Nqual = 1 full N:C dependency   0: only fresh material
qualDet   = (1.0d0-self%Nqual) + self%Nqual * det%N /(det%C + self%small) * self%CNref
qualDOM   = (1.0d0-self%Nqual) + self%Nqual * dom%N /(dom%C + self%small) * self%CNref
! distribute  preferential degradation rate to elements (N:fast; C:quality dep; others: intermediate)
if (associated(self%Index_DetN))  qualDetv(self%Index_DetN) = 1.0d0
if (associated(self%Index_DON))   qualDOMv(self%Index_DON)  = 1.0d0
if (associated(self%Index_DetC))  qualDetv(self%Index_DetC) = qualDet
if (associated(self%Index_DOC))   qualDOMv(self%Index_DOC)  = qualDOM
if (associated(self%Index_Det_No_NorC)) qualDetv(self%Index_Det_No_NorC) = (1.0d0 + qualDet)/2
if (associated(self%Index_DOX_No_NorC)) qualDOMv(self%Index_DOX_No_NorC) = (1.0d0 + qualDOM)/2
!_____________________________________________________________________________
!     denitrification 
! pelagic N-loss by denitrification, emulating benthic pool and suboxic micro-environments
! check for the carbon index 
! calculate substrate (all OM)
!if (associated(self%dom%C)) sum_OM = sum_OM + dom%C* qualDOM ! 

if (associated(self%det%C)) then ! POC Glud LO 2015 (suboxic spots in particles) 
  ! calculate oxidant (NO3)
  if (associated(self%dix%no3)) then ! TODO: add nitrite NO2
    nitrate = dix%no3
  elseif (associated(self%dix%N)) then
    nitrate = dix%N    ! lower denitrication self%denit
  else
    nitrate = 0.1_rk ! TODO: replace by SMALL
  endif
  denitrate = self%denit * sens%f_T * det%C * qualDet * (1.0d0-exp(-nitrate/self%DenitKNO3)) 
elseif
  denitrate = 0.0_rk
endif

!  ---  hydrolysis & remineralisation rate (temp dependent)
hydrolysis  = self%hydrolysis * sens%f_T 
remineral   = self%remineral  * sens%f_T 

!________________________________________________________________________________
!
!  --- DETRITUS C
! _GET_(self%id_, phy%mort, zoo%mort)
det_prod(i)  =   !(floppZ%C + zoo_mort) * zoo%C + aggreg_rate * phy%C  
dom_prod(i)  =   ! exud%C * phy%C 

do i = 1,num_elements ! e.g., N  ( C, Si, Fe, P)
  hydrol = hydrolysis * qualDetv(i) * det%elem(i)
  remin  = remineral  * qualDOMv(i) * dom%elem(i)
  rhsv%det%elem(i)   = det_prod(i) - hydrol
  rhsv%dom%elem(i)   = dom_prod(i) + hydrol - remin
  
  ! transfer matrix of remineralised DOX to DIX
  j = self%TransIndex_DOMDIX(i)
  if (j .gt. 0) then
     remin_stoich(j) = remin
  elseif (j .lt. 0) then ! potential partitioning between NO3 and NH4
     remin_stoich(-j) = remin * self%alloc_N
     remin_stoich(self%TransIndex_DON) = remin * (1.0_rk - self%alloc_N)
  endif
end do
! add denitrification of POC(!) Glud et al LO 2015 (suboxic spots in particles)
if(associated(self%Index_DetC)) rhsv%det%elem(self%Index_DetC) = rhsv%det%elem(self%Index_DetC) + denitrate 

!Index_DetN Index_DON Index_DetC Index_DOC Index_Det_No_NorC Index_DOX_No_NorC
nut_prod(i)  =   !  -uptake%N * phy%C + lossZ%N * zoo%C 

do i = 1,num_stoich 
  rhsv%nut%stoich(i) = remin_stoich(i) + nut_prod(i) 
end do
if(associated(self%Index_NO3)) rhsv%nut%stoich(self%Index_NO3) = rhsv%nut%stoich(self%Index_NO3) - 0.8_rk*denitrate 

!  chemostat mode 
if (self%dil .gt. 0.0d0) then
  do i = 1,num_elements ! e.g., N  ( C, Si, Fe, P)
    rhsv%det%elem(i) = rhsv%det%elem(i) - self%dil * det%elem(i)
    rhsv%dom%elem(i) = rhsv%dom%elem(i) - self%dil * dom%elem(i)
  end do
  do i = 1,num_stoich 
    rhsv%nut%stoich(i) = rhsv%nut%stoich(i) - self%dil * (self%dix%stoich0(i) - self%dix%stoich(i))
  end do
endif

! tell FABM about right hand sides ....
do i = 1,num_elements 
  _SET_ODE_(self%id_det(i), rhsv%det%elem(i) UNIT)
  _SET_ODE_(self%id_dom(i), rhsv%dom%elem(i) UNIT)
end do
do i = 1,num_stoich 
  _SET_ODE_(self%id_dix(i), rhsv%nut%stoich(i) UNIT)
end do

!_SET_DIAGNOSTIC_(self%id_vphys, exp(-self%sink_phys*phy%relQ%N * phy%relQ%P))       !average
! experimental formulation for emulating P-adsorption at particles in the water column and at the bottom interface
!_GET_HORIZONTAL_(self%id_o2flux, flO2)!_GET_HORIZONTAL_(self%id_oduflux, flODU)!_GET_HORIZONTAL_(self%id_zmax, zmax)  ! max depth
!aPO4 = (flODU-flO2)/(zmax+self%small)
!_SET_DIAGNOSTIC_(self%id_vphys, aPO4)       !average Temporary_diagnostic_
!________________________________________________________________________________

!define _REPLNAN_(X) X !changes back to original code
#define _REPLNAN_(X) nan_num(X)

if (self%BGC0DDiagOn) then
  _SET_DIAGNOSTIC_(self%id_qualDet, _REPLNAN_(qualDet))      !average Quality_of_POM_
  _SET_DIAGNOSTIC_(self%id_qualDOM, _REPLNAN_(qualDOM))      !average Quality_of_DOM_
end if

  _LOOP_END_

end subroutine tame_do

