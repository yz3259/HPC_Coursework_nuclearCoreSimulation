function Vec_Dot(u1,u2,v1,v2) result (d)
!=====================================================
!dot production of vectors of blocks
!    d = dot_product(u1(ibeg:iend),v1(ibeg:iend))&
!     +dot_product(u2(ibeg:iend),v2(ibeg:iend))
!=====================================================
  use header

  implicit none

  include "mpif.h"

  type(Vector),intent(inout)::u1,u2,v1,v2
  real(kind=8)::d

  real(kind=8)::local_sum
  integer:: ierr,i

! do the local dot_product (u1(ibeg:iend),v1(ibeg:iend))
! and sums the second part form the vectors u2 and v2     
! I can possibly use ddot(N,X,INCX,Y,INCY)
     local_sum = 0.0_8
  do i = u1%ibeg,u1%iend

     local_sum = local_sum + u1%xx(i)*v1%xx(i) + u2%xx(i)*v2%xx(i)

  end do  

!collect the local sum from all processors

  call MPI_Allreduce(local_sum,d,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

end function Vec_Dot
