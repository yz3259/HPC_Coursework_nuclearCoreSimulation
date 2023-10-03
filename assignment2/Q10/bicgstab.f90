subroutine bicgstab(M1,M2,S,F1,F2, FU1,FU2, sigma, maxit_inner, tau, u1,u2,k,P1,P2,Y1,Y2,R1,R2,Z1,Z2)
  !====================================
  ! Bicgstab algorithm
  ! It overwrites the right hand side FU1,FU2
  ! with the initial residual
  !====================================
  use header
  implicit none

  integer, intent(in) :: maxit_inner
  type(Matrix), intent(inout) :: M1, M2, S, F1, F2
  type(Vector), intent(inout) :: u1, u2, FU1, FU2
  real(kind=rk), intent(in) :: tau, sigma

  ! Internal variables
  integer,intent(out):: k
  integer :: n, ibeg, iend
  type(Vector),intent(inout) :: P1, P2, Y1, Y2, R1, R2, Z1, Z2
  real(kind=rk) :: alpha, omega, beta, rnrm, rnrm0
  ! External function
  real(kind=rk) :: vec_dot

  n = u1%n
  ibeg = u1%ibeg
  iend = u1%iend
  
! P1%n = n; P1%ibeg = ibeg; P1%iend = iend
!  P2%qn = n; P2%ibeg = ibeg; P2%iend = iend
!  Y1%n = n; Y1%ibeg = ibeg; Y1%iend = iend
!  Y2%n = n; Y2%ibeg = ibeg; Y2%iend = iend
!  R1%n = n; R1%ibeg = ibeg; R1%iend = iend
!  R2%n = n; R2%ibeg = ibeg; R2%iend = iend
!  Z1%n = n; Z1%ibeg = ibeg; Z1%iend = iend
!  Z2%n = n; Z2%ibeg = ibeg; Z2%iend = iend
!  allocate(P1%xx(n), P2%xx(n), R1%xx(n), R2%xx(n), Y1%xx(n), Y2%xx(n), Z1%xx(n), Z2%xx(n))

  ! Store initial residual r0 in FU1,FU2
  ! since we need it in all steps of BICGSTAB
  ! u1,u2 contain the initial guess
  ! subtract (A-sigma*F)*U = [M1*u1 - sigma*F1*u1 - sigma*F2*u2]
  !                          [M2*u2 - S*u1                     ]
  ! from the right hand side
  ! We have a tuned mat_mult (similar to dgemv) for that
  call mat_mult(-1.0_rk, M1, u1, 1.0_rk, FU1)
  call mat_mult(-1.0_rk, M2, u2, 1.0_rk, FU2)
  call mat_mult(1.0_rk,   S, u1, 1.0_rk, FU2)
  call mat_mult(sigma,   F1, u1, 1.0_rk, FU1)
  call mat_mult(sigma,   F2, u2, 1.0_rk, FU1)

  ! Copy initial residual to working vectors R1, R2
  call dcopy(iend-ibeg+1, FU1%xx(ibeg), 1, R1%xx(ibeg), 1)
  call dcopy(iend-ibeg+1, FU2%xx(ibeg), 1, R2%xx(ibeg), 1)    
  ! Copy residual to P1,P2
  call dcopy(iend-ibeg+1, R1%xx(ibeg), 1, P1%xx(ibeg), 1)
  call dcopy(n, R2%xx(ibeg), 1, P2%xx(ibeg), 1)
  
  ! Backup the norm of the initial residual
  rnrm0 = sqrt(vec_dot(R1,R2, R1,R2))
  rnrm = rnrm0

  ! Initialise extra vectors we use in mat_mult
  Y1%xx = 0.0_rk
  Y2%xx = 0.0_rk
  Z1%xx = 0.0_rk
  Z2%xx = 0.0_rk  
  
  do k=0,maxit_inner
     ! check convergence
     rnrm = sqrt(vec_dot(R1,R2, R1,R2))
     if (rnrm<tau*rnrm0) exit

     ! Y = A*P -- this must be backed up, since we need it later
     ! A*P = [M1  0]*[P1] - sigma*[F1 F2]*[P1]
     !       [-S M2] [P2]         [ 0  0] [P2]
     call mat_mult(1.0_rk, M1, P1, 0.0_rk, Y1)
     call mat_mult(1.0_rk, M2, P2, 0.0_rk, Y2)
     call mat_mult(-1.0_rk, S, P1, 1.0_rk, Y2)
     call mat_mult(-sigma, F1, P1, 1.0_rk, Y1)
     call mat_mult(-sigma, F2, P2, 1.0_rk, Y1)

     ! backup r0'*r since we need it in beta later
     beta = vec_dot(FU1,FU2, R1,R2)
     ! compute alpha -- the current update step
     alpha = beta/vec_dot(FU1,FU2, Y1,Y2)
     ! Correct solution and residual
     call vec_axpy(alpha, P1,P2, u1,u2)
     call vec_axpy(-alpha, Y1,Y2, R1,R2)
     ! check convergence
     rnrm = sqrt(vec_dot(R1,R2, R1,R2))
     if (rnrm<tau*rnrm0) exit

     ! The second update step
     ! Z = A*R, stored in Z1, Z2
     call mat_mult(1.0_rk, M1, R1, 0.0_rk, Z1)
     call mat_mult(1.0_rk, M2, R2, 0.0_rk, Z2)
     call mat_mult(-1.0_rk, S, R1, 1.0_rk, Z2)
     call mat_mult(-sigma, F1, R1, 1.0_rk, Z1)
     call mat_mult(-sigma, F2, R2, 1.0_rk, Z1)
     
     omega = vec_dot(R1,R2, Z1,Z2)/vec_dot(Z1,Z2, Z1,Z2)
     ! Correct solution and residual
     call vec_axpy(omega, R1,R2, u1,u2)
     call vec_axpy(-omega, Z1,Z2, R1,R2)

     ! Correct the direction P
     beta = (vec_dot(FU1,FU2, R1,R2)/beta)*(alpha/omega)
     call dscal(iend-ibeg+1, beta, P1%xx(ibeg), 1)
     call dscal(iend-ibeg+1, beta, P2%xx(ibeg), 1)
     call vec_axpy(1.0_rk, R1,R2, P1,P2)
     call vec_axpy(-beta*omega, Y1,Y2, P1,P2)
  end do

  write(*, '(4X,A,I5,A,es10.3)'), 'bicgstab conducted ', k, ' iterations up to residual = ', rnrm/rnrm0
  k = k
 ! deallocate(P1%xx, P2%xx, R1%xx, R2%xx, Y1%xx, Y2%xx, Z1%xx, Z2%xx)
end subroutine bicgstab
