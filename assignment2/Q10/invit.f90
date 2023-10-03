subroutine invit(M1, M2, S, F1, F2, maxit, eps, sigma, maxit_inner, tau, U1, U2, rho)
  !========================
  ! Inverse Power iteration
  !========================
  use header
  implicit none
  include "mpif.h"

  integer, intent(inout) :: maxit, maxit_inner
  type(Matrix), intent(inout) :: M1, M2, S, F1, F2
  type(Vector), intent(inout) :: U1, U2
  type(Vector) :: P1, P2, Y1, Y2, R1, R2, Z1, Z2
  real(kind=rk), intent(inout) :: eps, tau, sigma
  real(kind=rk), intent(inout) :: rho
  real(kind=rk)::Vec_Dot  
  real(kind=rk)::residual,nrm2
  real(kind=rk)::t_start,t_finish, bictime,totalt,bicite 
  integer::its!,myid
  integer::k
  type(Vector)::FU1,FU2,AU1,AU2,R11,R21
  ! Those definitions are for bicgstab()

  !===========================
  ! Implement Algorithm 1 here
  !===========================
  ! initial value of U1 and U2 are taken as input
  
  allocate(FU1%xx(U1%n),FU2%xx(U2%n),R21%xx(U2%n))
  allocate(AU1%xx(U1%n),AU2%xx(U2%n),R11%xx(U1%n))
  allocate(P1%xx(U1%n), P2%xx(U1%n), R1%xx(U1%n), R2%xx(U1%n),&
                  Y1%xx(U1%n), Y2%xx(U1%n), Z1%xx(U1%n), Z2%xx(U1%n))

  FU1%n    = U1%n; FU1%ibeg    = U1%ibeg; FU1%iend    = U1%iend
  FU2%n    = U2%n; FU2%ibeg    = U2%ibeg; FU2%iend    = U2%iend

  AU1%n    = U1%n; AU1%ibeg    = U1%ibeg; AU1%iend    = U1%iend
  AU2%n    = U2%n; AU2%ibeg    = U2%ibeg; AU2%iend    = U2%iend

  R11%n    = U1%n; R11%ibeg    = U1%ibeg; R11%iend    = U1%iend
  R21%n    = U1%n; R21%ibeg    = U1%ibeg; R21%iend    = U1%iend 

         P1%n =U1%n; P1%ibeg = U1%ibeg; P1%iend =U1%iend
         P2%n =U1%n; P2%ibeg = U1%ibeg; P2%iend = U1%iend
         Y1%n =U1%n; Y1%ibeg = U1%ibeg; Y1%iend = U1%iend
         Y2%n =U1%n; Y2%ibeg = U1%ibeg; Y2%iend = U1%iend
         R1%n =U1%n; R1%ibeg = U1%ibeg; R1%iend = U1%iend
         R2%n =U1%n; R2%ibeg = U1%ibeg; R2%iend = U1%iend
         Z1%n =U1%n; Z1%ibeg = U1%ibeg; Z1%iend = U1%iend
         Z2%n =U1%n; Z2%ibeg = U1%ibeg; Z2%iend = U1%iend
        


  totalt = 0.0_8 ! This is the initialize of total bic iteration time
  bicite = 0.0_8!initialize the counter for bic loop
  do its = 0,maxit
!     print*,'The processor number:',myid ! these two lines are for testing of progress
!     print*,'its', its
     ! start the inexact inverse calculation on each processor
     
        !=====================================================
        ! A*U = [ M1 0 ]*[U1]     F*U = [F1 F2]*[U1] 
        !       [ -S M2] [U2]           [0  0 ] [U2]
        !      Calculate products with A and F:
        !=====================================================

        !  A*U1 = M1*U1 
        ! Mat_Mult(alpha,A,u,beta,b): b = alpha*A*u+beta*b
        AU1%xx(U1%ibeg:U1%iend) = 0.0_8!initialize AU1
        call Mat_Mult(1.0_8,M1,U1,0.0_8,AU1)

        ! A*U2 = M2*U2 - S*U1
        AU2%xx(U1%ibeg:U1%iend) = 0.0_8 ! initialize AU2
        call Mat_Mult(1.0_8,S, U1,0.0_8,AU2)
        call Mat_Mult(1.0_8,M2,U2,-1.0_8,AU2)

        ! F*U1 = F1*U1+F2*U2; F*U2 = 0
        FU1%xx(U1%ibeg:U1%iend) = 0.0_8
        call Mat_Mult(1.0_8,F1,U1,0.0_8,FU1)
        call Mat_Mult(1.0_8,F2,U2,1.0_8,FU1)
        !initialise FU, FU2 is always 0
        FU2%xx(U1%ibeg:U1%iend) = 0.0_8                   
        ! Rayleigh quitient     
        rho = Vec_Dot(U1,U2,AU1,AU2)/&
              Vec_Dot(U1,U2,FU1,FU2)
        
        ! calculate the residual and residual 2norm
        ! The Vec_AXPY() only gives result on one processor
        ! [R1] = [AU1] - [FU1]*rho 
        ! [R2]   [AU2]   [FU2]
        R11%xx = 0.0_8;R21%xx = 0.0_8
        call Vec_AXPY(-rho,FU1,FU2,R11,R21)
        call Vec_AXPY(1.0_8,AU1,AU2,R11,R21)
        ! calculate the norm 
        residual =  sqrt(Vec_Dot(R11,R21,R11,R21))                          

        if (residual < eps) exit
        !print*,'residual',residual ! This part is for testing
        ! impleneting bicgstab()to solve equations 
         !        U1%xx= U1%xx/(rho-sigma)
         call dscal(U1%iend-U1%ibeg+1,1/(rho-sigma),U1%xx(U1%ibeg),1)
         !        U2%xx= U2%xx/(rho-sigma)
         call dscal(U2%iend-U2%ibeg+1,1/(rho-sigma),U2%xx(U2%ibeg),1)

        
         t_start = MPI_Wtime()         
         call bicgstab (M1,M2,S,F1,F2,FU1,FU2,sigma,&
                       maxit_inner,tau,U1,U2,k,P1,P2,Y1,Y2,R1,R2,Z1,Z2)
         t_finish = MPI_Wtime() 
         ! This records the time used for k bic iterations
         bictime = t_finish-t_start        
         totalt = totalt+bictime        
        ! This counter is counting how many single 
        ! loops of bicgstab was run. Updated by every invit iteration
         bicite = bicite +k
       
        ! normalise Unew = Unew/U_2norm
        nrm2 = sqrt(vec_dot(U1,U2,U1,U2))
        call dscal(U1%iend-U1%ibeg+1,1/nrm2,U1%xx(U1%ibeg),1)
        call dscal(U1%iend-U1%ibeg+1,1/nrm2,U2%xx(U2%ibeg),1)
!        print*,'the 2norm of U',sqrt(Vec_Dot(U1,U2,U1,U2))  
  end do
 rho = rho
 print*, 'the norm of resicual R:',residual,'/'
 print*, 'The total number of bic iterations:',bicite,'/'
 print*, 'The total time spent on bic iterations:',totalt,'/'
 print*, 'The time for one Bic iteration:', totalt/(bicite*1.0_8)
 deallocate(FU1%xx,FU2%xx,R21%xx)
 deallocate(AU1%xx,AU2%xx,R11%xx)
  deallocate(P1%xx, P2%xx, R1%xx, R2%xx, Y1%xx, Y2%xx, Z1%xx, Z2%xx)  
end subroutine invit

