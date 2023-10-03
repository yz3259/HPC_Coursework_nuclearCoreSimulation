!-----------------------------------------------------------------------
! This is for testing the dotproduct and vecaxpy
! -----------------------------------------------------------------------

program testq1

  use header
  implicit none
  include "mpif.h"

  type(Vector) :: u1,u2, v1,v2,x1,x2
 
  real (kind=8) :: alpha,d,d_real
  integer       :: n=8
  integer       :: myid, nprocs, nblock, ibeg, iend
  integer       :: ierr

  real(kind=8)  ::  Vec_Dot


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Beginning of program - Initialisation of MPI context
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calculate the start and end indices of the rows to be held locally
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nblock = n/nprocs
  ibeg   = myid*nblock +1
  if (myid<nprocs -1)then
     iend = (myid+1)* nblock
  else
     iend = n
     nblock = n - myid*nblock
  end if 

 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Allocate memory for u1, u2, v1,v2 and set dimensions
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  allocate(u1%xx(n),u2%xx(n),v1%xx(n),v2%xx(n),x2%xx(n),x1%xx(n))

  u1%n    = n
  u1%ibeg = ibeg
  u1%iend = iend

  u2%n    = n
  u2%ibeg = ibeg
  u2%iend = iend

  v1%n    = n
  v1%ibeg = ibeg
  v1%iend = iend

  v2%n    = n
  v2%ibeg = ibeg
  v2%iend = iend
 
  x1%n    = n
  x1%ibeg = ibeg
  x1%iend = iend

  x2%n    = n
  x2%ibeg = ibeg
  x2%iend = iend

  call random_number(u1%xx)
  call random_number(u2%xx)
  call random_number(v1%xx)
  call random_number(v2%xx)

!Testing the result from dot_product(u1,v1)+dot_product(u2,v2) and the Vec_Dot
  x1%xx = 0.0_8
  x2%xx = 1.0_8
  alpha = 1.0_8
  d = Vec_Dot(u1,u2,v1,v2)
  
!  if (myid == 0) then
   
    
    d_real=dot_product(u1%xx,v1%xx)+dot_product(u2%xx,v2%xx)
    print*, "err of dot_product:",d_real-d, 'This is done by processor:',myid
    
 ! end if

    call Vec_AXPY(alpha,u1,u1,x1,x2)
    call MPI_Allgather(x1%xx(x1%ibeg),iend-ibeg+1,MPI_DOUBLE_PRECISION,&
        x1%xx,iend-ibeg+1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)    

    call MPI_Allgather(x2%xx(x2%ibeg),iend-ibeg+1,MPI_DOUBLE_PRECISION,&
        x2%xx,iend-ibeg+1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)    


!  if (myid ==0) then 
    print*,'processor number',myid, "err of Vec_Axpy 1",(x1%xx-u1%xx)
    print*,'processor number',myid,"err of Vec_Axpy2",x2%xx-u1%xx
 
 ! end if

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Deallocate memory
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  deallocate(u1%xx,u2%xx,v1%xx,v2%xx,x1%xx,x2%xx)

call MPI_Finalize(ierr)

end program testq1




