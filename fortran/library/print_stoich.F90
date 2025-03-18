PROGRAM print_stoich
  use tame_stoich_functions
  ! Ausgabe des Arrays zur Überprüfung
    INTEGER :: i, in, j, k
    real(8) :: din(4)  = (/ 2., 8., 14., 32./)
    real(8) :: dip(4)  = (/ 0.05, 0.2, 0.4, 1./)
    real(8) :: par(4)  = (/ 10., 20., 40., 120./)
    real(8) :: q_param(4),q(2)
    real(8) :: temp=24

    ! loop over env conditions
    i=2 ! N, P
    DO in = 1, 4
     DO j = 1, 4
      DO k = 1, 4
       call quota_params(din(j), par(k), temp, i, q_param)
       !write (*,'(4F7.4) ') q_param
       q(i)=quota_response(q_param,dip(in))
      write (*,'(A, F3.0,A,F4.2,A,F4.0,A, F6.3)') 'din=',din(in),' dip=',dip(j),' par=',par(k),'-> QP=',q(i)
      END DO
     END DO
    END DO  
END PROGRAM print_stoich
