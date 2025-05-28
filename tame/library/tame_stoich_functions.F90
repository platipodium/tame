! SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
!
! SPDX-License-Identifier: Apache-2.0

!> @file tame_stoich_functions.f90
!> @author kai wirtz
!---------------------------------------------------------
! !module: tame_functions
!> @brief stoichiometry related functions (derived from MAECS output) called by tame
module tame_stoich_functions
   use tame_types

   implicit none
! #ifndef fabm
!#define rk real(8)
!#else
!use fabm, only :: rk
!#endif
! !uses:
   public quota_response, calc_quota, quota_params, queuefunc, invmatrix_2
contains
   !
   ! linear quota equation
   ! function for quota dependence on ambient nutrient
   real(8) function quota_response(param, nut, IsP)
      real(8), intent(in) :: param(4), nut
      logical, intent(in) :: IsP
      real(8) :: arg, syn, q0
      ! minimal quota  - clipping at zero
      q0 = max(0.0, param(1))
      ! synchrony: inf:blackman/linear, 2:ivlev 1:michaelis-menten/holling-ii
      syn = min(param(4), 3.0) ! TODO: remove clipping at 3 ?
      !arg = nut / max(nut * 1.0e-2, param(3))
      quota_response = q0 + (param(2) - q0)*queuefunc(syn, nut/param(3))
      if (IsP) quota_response = quota_response*1.E-3  ! back-transformation of P-quota
      !if (IsP)  write (*,'(5F12.7) ') nut,q0,param(2),syn,nut / param(3)

   end function quota_response

! function for quota dependence on ambient nutrient
   ! TODO: check for accuracy
   ! TODO: complete with ALL derivatives (complicated :-(
   real(8) function quota_nut_deriv(param, nut, IsP)
      real(8), intent(in) :: param(4), nut
      logical, intent(in) :: IsP
      real(8) :: arg, syn, q0
      ! minimal quota  - clipping at zero
      q0 = max(0.0, param(1))
      ! synchrony: inf:blackman/linear, 2:ivlev 1:michaelis-menten/holling-ii
      syn = min(param(4), 3.0) ! TODO: remove clipping at 3 ?
      !arg = nut / max(nut * 1.0e-2, param(3))
      quota_nut_deriv = (param(2) - q0)*qfunc_deriv(syn, nut/param(3))/param(3)
      if (IsP) quota_nut_deriv = quota_nut_deriv*1.E-3  ! back-transformation of P-quota
   end function quota_nut_deriv

!-----------------------------------------------------------------------
!> @brief the queue function
!> @details
!> provides both the queuing function and it's derivative
!> with the parameter n->inf :liebig and n~1:product
   real(8) function queuefunc(syn, x)
      implicit none
      real(8), intent(in) :: x, syn
      real(8) :: px, dn
! synchrony: inf :Blackman/linear, 2:Ivlev  1:Michaelis-Menten/Holling-II
      if (abs(1.0 - x) .LT. 1.E-4) then
         queuefunc = syn/(syn + 1.0)
      else
         px = x**(syn + 1.0)
         dn = 1.0/(1.0 - px)
         queuefunc = (x - px)*dn
      end if
   end function queuefunc

   real(8) function qfunc_deriv(syn, x)
      implicit none
      real(8), intent(in) :: x, syn
      real(8) :: px, dn
      if (abs(1.0 - x) .LT. 1.E-5) then
         qfunc_deriv = 0.5*syn/(syn + 1.0)
      else
         px = x**(syn + 1.0)
         dn = 1.0/(1.0 - px)
         qfunc_deriv = (1.0 - (syn + 1.0)*x**syn + syn*px)*dn*dn
      end if
   end function qfunc_deriv

! ----------------------------
! retrieve parameters of linear quota equation
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
      real(8) :: exparv(4) = (/0.5, 0.20, -0.25, 1.0/)
      real(8) :: trade_param(4) = (/20., 20., 0.04, 0.2/)
      real(8) :: ekx, xpf, epf, sgn, soff, pf
      integer :: j, pi
      real(8) :: pp(2, 4, 2) = reshape([ &
                                       ! pp{1,1:4}
                                       0.0578, 4.8138, &
                                       0.0556, 1.7993, &
                                       6.2234, 0.9400, &
                                       0.0481, 20.5375, &
                                       ! pp{2,1:4}
                                       0.7969, 1.7666, &
                                       5.4002, 1.0735, &
                                       0.0414, 0.0, &
                                       3.0491, 0.0], [2, 4, 2], ORDER=[3, 2, 1])

!do i = 1,2 ! ! nutrient index 1:N 2:P
      do pi = 1, 4 ! !
         param = pp(i, pi, 1:2)

         ! general dependencies on PAR and Temp
         ! exponent of scaled PAR and Temp dependencies; functional sensitivity
         env_exp = 0.5 * (1.0 + merge(1.0,0.0, pi==1 .or. pi==3) * merge(1.0,0.0, i==1)+merge(1.0,0.0, pi<=3)*merge(1.0,0.0, i==2))
         ! scaled Temp dependency
         xtemp = (temp/trade_param(2))**(1.0 - merge(0.5, 0.0, pi >= 3))
         ! exponent of scaled PAR  dependency
         if (i == 2) then ! DIP -> Q_P
            expar = exparv(pi)
         else             ! DIN -> Q_N
            expar = 0.5
         end if
         ! scaled PAR  dependency
         xpar = (1.0/(1.0 + par/trade_param(1)))**expar
         ! combined scaled PAR & Temp dependency
         xpart = (xpar*xtemp)**env_exp

         ! exponent of scaled Temp dependency AND 2nd Nutrient factor
         ekx = 0.5*(1.0 + merge(1.0, 0.0, i == 1 .AND. pi == 2))
         ! scaling constant for ambient nutrient concentration
         kx = trade_param(2 + i)*(1.+2.*ekx*(xtemp**ekx))
         ! scaled (and log) complementary nutrient (N, P, ..)
         nut = log(dix/kx)

         if (i == 1) then  ! DIN -> Q_N
            if (pi > 1) then
               xpf = nut/(1.+0.5*ekx*nut)
               epf = exp(-0.5*xpf**2)
            end if

            select case (pi)
            case (1)
               xpf = (0.5 - 1.0/(1.0 + exp(-1.0*nut)))/(1.0 + nut*nut + 0.5*nut)
               q_param(pi) = param(1)*(1.0 + param(2)*xpart*xpf)
            case (2)
               soff = sqrt(xpar) + sqrt(xtemp)
               pf = soff + param(2)*xpart*epf
               q_param(pi) = param(1)*pf
            case default
               sgn = 1.0 - merge(2.0, 0.0, pi == 4)
               q_param(pi) = param(1)*(1.0 + sgn*param(2)*(sqrt(xpart)*(0.5*(1.0 - epf) - 1.0)))
            end select
         else ! DIP -> Q_P
            xpf = nut/(1.0 + nut)
            select case (pi)
            case (1)
               q_param(pi) = max(0.3, param(1)*exp(-param(2)*xpart))
            case (2)
               pf = 1.0 + param(2)*xpart
               q_param(pi) = max(4.0, param(1)*pf*xpf)
            case (3)
               q_param(pi) = param(1)*(xpart*xpf)

            case (4)
               q_param(pi) = param(1)*(xpart*xpf)
               if (q_param(pi) > 2.0) then
                  q_param(pi) = 2.0
               end if
            end select
         end if
      end do
   end subroutine quota_params

   real(8) function calc_quota(nut_i, nut_j, par, temp, i, j)
      real(8), intent(in) :: nut_i, nut_j, par, temp
      integer, intent(in) :: i, j     ! nutrient index 1:N 2:P
      logical :: IsPhosporus
      real(8) :: q_param(4), arg, syn, q0
      ! clip for too low values for response function
      ! TODO : same for PAR :: unrealistic values at night
      call quota_params(max(nut_j, nut_minval(j)), par, temp, i, q_param)  ! retrieve parameters of linear quota equation

      IsPhosporus = (nutrient_name(i) == 'PO4')
      calc_quota = quota_response(q_param, max(nut_i, nut_minval(i)), IsPhosporus) ! linear quota equation
      !  if (quota(ie) .lt. fixed_stoichiometry(ie)/8)  then
      !    write (*,'(I3,9F10.4) ') ie,nutrient(i),nut,quota(ie),q_param(4),q_param(1)*100,q_param(2),nutrient(i)/q_param(3),1E-3*(q_param(2))*queuefunc(q_param(4), nut / q_param(3))
      !stop
      !  end if
      !   if (i==2)  write (*,'(6F12.8) ') max(nut_i,nut_minval(i)),calc_quota,q_param

   end function calc_quota

! calculates the inverse of a 2×2 matrix.
   function invmatrix_2(A) result(iA)
      real(8), dimension(2, 2), intent(in) :: A  !! input matrix
      real(8), dimension(2, 2) :: iA   !! inverse matrix
      real(8) :: detinv

      ! inverse determinant of the matrix
      detinv = 1./(A(1, 1)*A(2, 2) - A(1, 2)*A(2, 1))

      ! Calculate the inverse of the matrix
      iA(1, 1) = +detinv*A(2, 2)
      iA(2, 1) = -detinv*A(2, 1)
      iA(1, 2) = -detinv*A(1, 2)
      iA(2, 2) = +detinv*A(1, 1)
   end function

end module
