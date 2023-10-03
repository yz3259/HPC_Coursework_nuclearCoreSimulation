program neutron
  !===============================
  ! Sequential program for solving
  ! 2D neutron diffision equations
  !===============================
  use header
  implicit none
  include "mpif.h"
  
  integer       :: n, m, maxit, maxit_inner
  integer       :: ierr
  integer       :: myid,nprocs,nrows,ibeg,iend
  type(Matrix)  :: M1, M2, S, F1, F2
  type(Vector)  :: U1, U2 ! solution
  real(kind=rk) :: eps, tau, rho, sigma, t_start, t_finish
  
  
  !=====================================================
  !Beginning of program - Initialisation of MPI context
  !=====================================================

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

  ! Read in parameters from the input.dat on processor 0
  ! Print out the parameters on processor 0
  if (myid == 0) then
      print*,'The input parameters are:'

      open(unit=2,file="input.dat")
      read(2,*) m; read(2,*) eps; read(2,*) maxit
      read(2,*) sigma; read(2,*) tau; read(2,*) maxit_inner

      write(*,'(A,I5)'),'Value of m =',m 
      write(*,'(A,es10.3)'),'Value of eps =',eps 
      write(*,'(A,I5)'),'Value of maxit  =',maxit 
      write(*,'(A,es10.3)'),'Value of sigma =',sigma 
      write(*,'(A,I5)'),'Value of maxit_inner =',maxit_inner 
      write(*,'(A,es10.3)'),'Value of tau =',tau 
      close(2)
  end if

  !=======================================================
  !  Broadcast the 6 paramenters to other processes
  !=======================================================
  call MPI_Bcast(m,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)  
  call MPI_Bcast(eps,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)  
  call MPI_Bcast(maxit,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)  
  call MPI_Bcast(sigma,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)  
  call MPI_Bcast(maxit_inner,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)  
  call MPI_Bcast(tau,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)  

   ! the following was for testing whether each processor
   ! got the same correct parameters
!      print*,'The processor number is:',myid
!      write(*,'(A,I5)'),'Value of m =',m 
!      write(*,'(A,es10.3)'),'Value of eps =',eps 
!      write(*,'(A,I5)'),'Value of maxit  =',maxit 
!      write(*,'(A,es10.3)'),'Value of sigma =',sigma 
!      write(*,'(A,I5)'),'Value of maxit_inner =',maxit_inner 
!      write(*,'(A,es10.3)'),'Value of tau =',tau 
     
  !=======================================================
  !  Calculate the start and end indices of the rows 
  !  to be held locally
  !=======================================================
  n     = m*m
  nrows = n/nprocs
  ibeg  = myid*nrows +1  
  iend = (myid+1)*nrows
  
  ! Assemble matrices
  call create_matrices(m, M1, M2, S, F1, F2, ibeg,iend)

  ! Allocate vectors
  U1%n = nrows ;U2%n = U1%n
  allocate(U1%xx(U1%n), U2%xx(U2%n))

  U2%ibeg = ibeg; U2%iend = iend

  U1%ibeg = ibeg; U1%iend = iend 
  ! Initial guess
  U1%xx = sqrt(1.0_rk/(2*U1%n))
  U2%xx = U1%xx
  
  ! Solve the eigenvalue problem
  t_start = MPI_Wtime()
  call invit(M1, M2, S, F1, F2, maxit, eps, sigma, maxit_inner, tau, U1, U2, rho)
  t_finish = MPI_Wtime()  

  ! Save solutions to a file on processor 0

     print *, 'cpu time elapsed = ', t_finish - t_start
     print *, 'multiplication constant lambda_0 = ', 1.0_rk/rho
  ! we need to collect all the truncated U1 and U2 back to processor 0 before saving them

     call MPI_Gather(U1%xx(U1%ibeg),iend-ibeg+1,MPI_DOUBLE_PRECISION,&
        U1%xx,iend-ibeg+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     call MPI_Gather(U2%xx(U2%ibeg),iend-ibeg+1,MPI_DOUBLE_PRECISION,&
        U2%xx,iend-ibeg+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  if(myid ==0)then  
     call save_fields(U1,U2,'solution.dat')
  end if
!  call save_fields(U1,U2,'solution.dat')
  deallocate(M1%ii, M2%ii, S%ii, F1%ii, F2%ii)
  deallocate(M1%aa, M2%aa, S%aa, F1%aa, F2%aa)
  deallocate(M1%jj, M2%jj, S%jj, F1%jj, F2%jj)
  deallocate(U1%xx, U2%xx)

  call MPI_Finalize(ierr)

end program neutron

