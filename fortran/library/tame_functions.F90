! SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
!
! SPDX-License-Identifier: Apache-2.0

!> @file tame_functions.F90
!> @author Kai Wirtz
#include "fabm_driver.h"
!---------------------------------------------------------
! !MODULE: tame_functions --- more to come
!> @brief  functions called by tame
   module tame_functions

! !USES:
   use fabm_types
   use tame_types
   public   queuefunc1, queuederiv , smooth_small, sinking, min_mass,  &
            calc_sensitivities, calc_internal_states, nan_num,  &
            set_pointer, set_chem_pointer
   private  upper

 contains
 !------------------------------------------------------
!> @brief calculate the internal states
!> @details
!> @todo the theta-related calculations should obviously be related to \latexonly Eq. \ref{eq:ftheta} \endlatexonly but I get lost. See the Q's therein
subroutine calc_internal_states(tame,det,dom) !,phy,zoo

implicit none
class (type_tame_base_model),intent(in)     :: tame
!type (type_tame_phy), intent(inout) :: phy
type (type_tame_elem), intent(inout) :: det
type (type_tame_elem), intent(inout) :: dom
!type (type_tame_zoo), intent(inout) :: zoo

!> @fn tame_functions::calc_internal_states()
!> 1. Calculate elemental absolute and relative quotas (Q and relQ):
!>   - phy\%Q\%X = phy\%X / phy\%C where x=N,P,Si
!>   - phy\%relQ\%X= (phy\%Q\%X - tame\%qN_phy_0) / tame\%iK_QN where x=N,P,Si
!!  phy%Q%N    = phy%reg%N / phy%reg%C
! added for mixing effects in estuaries kw Jul, 15 2013
!write (*,'(A,2(E19.5))') 'QN-1:',phy%Q%N,tame%small_finite* tame%QN_phy_max
!! phy%Q%N  = smooth_small(phy%Q%N, tame%QN_phy_0 + tame%small_finite * tame%QN_phy_max)
! fraction of free (biochemically available) intracellular nitrogen
!! phy%relQ%N  = (phy%Q%N - tame%QN_phy_0) * tame%iK_QN
!! phy%relQ%N  = smooth_small(phy%relQ%N,tame%small)
end subroutine

! set indices of elements vectors
function upper(text) result(up_text)
implicit none
character(len=*), intent(in) :: text
character(len=len(text)) :: up_text
integer :: i, ichar_val

up_text = text
do i = 1, len(text)
   ichar_val = ichar(text(i:i))

   ! Convert lowercase letters (a-z) to uppercase (A-Z)
   if (ichar_val >= ichar('a') .and. ichar_val <= ichar('z')) then
      up_text(i:i) = char(ichar_val - 32)
   end if
end do
end function upper

! set indices of elements vectors
! set pointer to indexed elements vector
subroutine set_pointer(object, vector, name, value)
    type(type_tame_elem), intent(inout) :: object
    real(rk), target, intent(in)        :: vector(:)
    character(*), intent(in) :: name
    integer, intent(in) :: value

    ! Add bounds checking
    if (value < 1 .or. value > 5) then
        print *, "Error: value out of bounds"
        return
    end if
    select case (upper(trim(adjustl(name))))  ! Make case matching more robust
    case ('C')
        object%C => vector(value)
        object%index%C = value
    case ('N')
        object%N => vector(value)
        object%index%N = value
    case ('P')
        object%P => vector(value)
        object%index%P = value
   !  case ('S') !Si
   !       object%Si => vector(value)
   !       object%index%Si = value
   !  case ('F') !Fe
   !       object%Fe => vector(value)
   !       object%index%Fe = value
    ! ... add other elements ...
    case default
        print *, "Error: unknown element name"
    end select
end subroutine set_pointer

subroutine set_chem_pointer(object, vector, name, value)
    type(type_tame_chemical), intent(inout) :: object
    real(rk), target, intent(in)          :: vector(:)
    character(len=3), intent(in)              :: name
    integer, intent(in)                   :: value
    !  bounds checking
    if (value < 1 .or. value > 10) then
        print *, "Error: value out of bounds"
        return
    end if
    select case (upper(trim(adjustl(name))))  ! Make case matching more robust
    ! NO3,NH4,PO4,CO2,O2,SiO2,FeS,DIN,DIP,DISi,DIC
    case ('NO3')
        object%NO3 => vector(value)
        object%index%NO3 = value
    case ('NH4')
        object%NH4 => vector(value)
        object%index%NH4 = value
    case ('PO4')
        object%PO4 => vector(value)
        object%index%PO4 = value
   !  case ('CO2')
   !       object%CO2 => vector(value)
   !       object%index%CO2 = value
   !  case ('FeS')
   !       object%FeS => vector(value)
   !       object%index%FeS = value
   !  case ('SiO') !
   !       object%SiO2 => vector(value)
   !       object%index%SiO2 = value
     case ('DIN') !
          object%DIN => vector(value)
          object%index%DIN = value
   !  case ('DIP') !
   !       object%DIP => vector(value)
   !       object%index%DIP = value
   !  case ('DIC') !
   !       object%DIC => vector(value)
   !       object%index%DIC = value
   !  case ('DIS') !
   !       object%DISi => vector(value)
   !       object%index%DISi = value
    ! ... add other elements ...
    case default
        print *, "Error: unknown chemical name ",name
    end select
end subroutine set_chem_pointer

!------------------------------------------------------
!> @brief calculate sensitivities
!> @details Details:
!> - sens\%f\_T \latexonly see eq. \ref{eq:arrhenius} \endlatexonly
subroutine calc_sensitivities(sens,env,q10,tref)

implicit none
!class (type_tame_base_model), intent(in) :: tame
real(rk), intent(in)          :: q10, tref
type (type_tame_sensitivities), intent(out) :: sens
type (type_tame_env),intent(in) :: env
real(rk) :: par, T_Kelv

!> @fn tame_functions::calc_sensitivities()
!> 1. calculate (sens\%) f\_T, P\_max\_T, a\_light, upt\_pot\%C
! par          = tame%frac_PAR * env%par ! use active  fraction frac_PAR of available radiation
T_Kelv       = env%Temp + 273.d0 ! temperature in Kelvin
! ----------------------------------------------------------------------------
! +++ determine rates for photoautotophic growth ++++++++++++++++++++++++++++++++++++++++++++++++++
! --- temperature dependence of metabolic rates (with T_ref given in units [Kelvin]) --------------
!standard Q10 RULE
sens%f_T     = q10**((T_Kelv-tref)/10.0)

end subroutine

!-----------------------------------------------------------------------
!> @brief the queuing function
!> @details
!> provides both the queuing function and it's derivative
!> with the parameter n->inf :liebig and n~1:product
!> \latexonly see: Section \ref{sec:colim} \endlatexonly \n
!> @todo: add equations
real function  queuefunc0(syn,x)
   implicit none
   real(rk), intent(in)          :: x, syn
   real(rk)                      :: px, dn
! synchrony: inf :Blackman/linear, 2:Ivlev  1:Michaelis-Menten/Holling-II

   if(abs(1.0_rk-x) .lt. 1E-2) then
      queuefunc0 = syn/(syn+1.0_rk)
   else
      px    = x**(syn+1.0_rk)
      dn    = 1.0_rk / (1.0_rk-px)
      queuefunc0 =  (x-px) * dn
   endif
end function queuefunc0

!-----------------------------------------------
!> @brief numerical approximation of the queue function
!> @details
!> n->inf :liebig  n~1:product\n
!> Here is an example of adding a snip of code:
!> @snippet tame_functions.F90 queuefunc1_snippet
subroutine queuefunc1(n,x,qfunc,dq_dx)
   implicit none
   real(rk), intent(in)       :: x, n
   real(rk), intent(out)      :: qfunc, dq_dx
   real(rk)                   :: nn, hh,x0,en
   ! [queuefunc1_snippet]
   nn = n+1.
   hh = 1./nn
   x0 = (log(exp(nn)-1))*hh
   en = exp(-nn*(x/(1+hh*x) - x0))
   qfunc =  1. - hh*log(1.+ en)
   dq_dx = 1. / ( (1. + hh * x)**2 * (1.+1./en))
   ! [queuefunc1_snippet]

end subroutine queuefunc1

!-------------------------------------------------------------
!> @brief derivative of the queue function ??
!> @details
!> @return queuederiv
!> @todo: add equations
real(rk) function queuederiv(n,x)
!
   implicit none
   real(rk), intent(in)          :: x, n
   real(rk)                      :: nqueue_aa
   real(rk)                      :: nqueue_bb = 1.3863d0
!-----------------------------------------------------------------------
   nqueue_aa = 6.9315d-1*(1.d0+1.0d0/(n))
   queuederiv = exp(-nqueue_aa*(x)**nqueue_bb)
end function queuederiv

!------------------------------------------------------
!> @brief minimum boundary for numerical stability in 3D simulations
   !---------------------------------------------------------
!> @brief  continuous smoothing function by kw Apr 2012
!> @details
!! smoothly converges x to **eps/2** for x<eps
!! \f[ x=eps+(x-eps)*e^{x/eps}/(1+e^{x/eps}) \f]
pure real(rk) function smooth_small(x, eps)
   implicit none
   real(rk), intent(in)          :: x, eps
   real(rk)                      :: arg, larger, larger2
!   integer                       :: na = 2
   integer                       :: nb
!--------------------------------------------------------------
   nb      = 8
   if (x .lt. nb*eps) then
     arg     = x/(eps+1E-7)
     larger  = exp(2*arg)
     larger2 = exp(arg/2)
     smooth_small  = ((nb+larger2)*eps + x*larger)/(nb+larger)
   else
    smooth_small  = x
   endif
end function smooth_small

end module tame_functions
!------------------------------------------------------
