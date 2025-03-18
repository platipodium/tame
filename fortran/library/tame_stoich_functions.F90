!> @file tame_stoich_functions.f90
!> @author kai wirtz
!---------------------------------------------------------
! !module: tame_functions 
!> @brief stoichiometry related functions (derived from MAECS output) called by tame 
module tame_stoich_functions
  implicit none

! !uses:
   public   quota_response,quota_params,queuefunc
 contains  
 !
 ! function for quota dependence on ambient nutrient
real(8) function quota_response(param, nut)
  real(8), intent(in) :: param(4), nut
  real(8) :: arg, syn, q0  
  ! minimal quota  - clipping at zero
  q0  = max(0.0, param(1))
  ! synchrony: inf:blackman/linear, 2:ivlev 1:michaelis-menten/holling-ii
  syn = min(param(4), 3.0) ! TODO: remove clipping at 3 ?
  !arg = nut / max(nut * 1.0e-2, param(3))
  quota_response = q0 + (param(2) - q0) * queuefunc(syn, nut / param(3))
end function quota_response

real function  queuefunc(syn,x)
   implicit none
   real(8), intent(in)          :: x, syn
   real(8)                      :: px, dn
! synchrony: inf :Blackman/linear, 2:Ivlev  1:Michaelis-Menten/Holling-II

   if(abs(1.0-x) .lt. 1E-2) then
      queuefunc = syn/(syn+1.0)
   else
      px    = x**(syn+1.0)
      dn    = 1.0 / (1.0-px)
      queuefunc =  (x-px) * dn
   endif
end function queuefunc

! function for array of trade-offs in the primer parameters of the Q_X(DIX) relation
subroutine quota_params(dix, par, temp, i, q_param)
implicit none
!  integer, intent(in) :: pi    ! index for response parameter 
integer, intent(in) :: i     ! nutrient index 1:N 2:P
real(8), intent(in) :: dix, par, temp ! ambient 2nd nutrient, light and temperature
!real(8), intent(in) :: trade_param(4) ! trade-off parameter trade_param_f=[20 20 0.04 0.2]
real(8), intent(out) :: q_param(4)
real(8) :: param(2) ! input response parameter
real(8) :: env_exp, xtemp, xpar, xpart, kx, nut, q0, expar
real(8) :: exparv(4)= (/ 0.5, 0.20, -0.25, 1.0 /)
real(8) :: trade_param(4)  = (/ 20., 20., 0.04, 0.2 /)
real(8) :: ekx, xpf, epf, sgn, soff, pf
integer :: j,pi
real(8) :: pp(2, 4, 2) = RESHAPE([ &
    ! pp{1,1:4}
    0.0578, 4.8138, &
    0.0556, 1.7993, &
    6.2234, 0.9400, &
    0.0481, 20.5375, &
    ! pp{2,1:4}
    0.7885, 1.7478, &
    5.4712, 1.0775, &
    0.0413, 0.0, &
    3.0538, 0.0 ], [2, 4, 2], ORDER=[3, 2, 1])
  
!do i = 1,2 ! ! nutrient index 1:N 2:P
  do pi = 1,4 ! ! 
    param=pp(i, pi, 1:2)

    ! general dependencies on PAR and Temp
    ! exponent of scaled PAR and Temp dependencies; functional sensitivity
    env_exp = 0.5 * (1.0 + merge(1.0,0.0, pi==1 .or. pi==3) * merge(1.0,0.0, i==1)+merge(1.0,0.0, pi<=3)*merge(1.0,0.0, i==2))
    ! scaled Temp dependency
    xtemp = (temp/trade_param(2))**(1.0-merge(0.5,0.0,pi>=3)) 
    ! exponent of scaled PAR  dependency
    if (i == 2) then ! DIP -> Q_P
        expar = exparv(pi)
    else             ! DIN -> Q_N
        expar = 0.5
    end if
    ! scaled PAR  dependency
    xpar = (1.0/(1.0+par/trade_param(1)))**expar
    ! combined scaled PAR & Temp dependency
    xpart = (xpar*xtemp)**env_exp
    
    ! exponent of scaled Temp dependency AND 2nd Nutrient factor
    ekx = 0.5*(1.0+merge(1.0,0.0,i==1 .and. pi==2))
    ! scaling constant for ambient nutrient concentration
    kx = trade_param(2+i) * (1.+2.*ekx*(xtemp**ekx)) 
    ! scaled (and log) complementary nutrient (N, P, ..)
    nut = log(dix/kx)
    
    if (i == 1) then  ! DIN -> Q_N
        if (pi > 1) then
            xpf = nut/(1. + 0.5 * ekx * nut)
            epf = exp(-0.5*xpf**2)
        end if
        
        select case(pi)
        case(1)
            xpf = (0.5 - 1.0/(1.0+exp(-1.0*nut)))/(1.0+nut*nut+0.5*nut)
            q_param(pi) = param(1)*(1.0+param(2)*xpart*xpf)
        case(2)
            soff = sqrt(xpar) + sqrt(xtemp)
            pf = soff + param(2)*xpart*epf
            q_param(pi) = param(1)*pf
        case default
            sgn = 1.0-merge(2.0,0.0,pi==4)
            q_param(pi) = param(1)*(1.0+sgn*param(2)*(sqrt(xpart)*(0.5*(1.0-epf) - 1.0)))
        end select   
    else ! DIP -> Q_P
        xpf = nut/(1.0+nut)
        select case(pi)
        case(1)
            q_param(pi) = param(1) * exp(-param(2)*xpart)
        
        case(2)
            pf = 1.0 + param(2)*xpart
            q_param(pi) = param(1)*pf*xpf      
        case(3)
            q_param(pi) = param(1)*(xpart*xpf)
        
        case(4)
            q_param(pi) = param(1)*(xpart*xpf)
            if (q_param(pi) > 2.0) then
                q_param(pi) = 2.0
            end if
        end select
    end if
  end do
end subroutine quota_params
end module 