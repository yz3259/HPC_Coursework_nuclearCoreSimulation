!------------------------------------------------------------------------
!
!  Subroutine to multiply a matrix and a vector in compressed row storage 
!  format, i.e to calculate b = alpha*A*u + beta*b (parallel version).
!
!------------------------------------------------------------------------

subroutine Mat_Mult(alpha,A,u,beta,b)

  use header
  implicit none
  include "mpif.h" !this is added for calling MPI subroutines

  real(kind=rk), intent(in) :: alpha ! Multiplier of A*u
  type(Matrix), intent(inout)  :: A      ! Matrix to multiply with
  type(Vector), intent(inout)  :: u      ! Input vector u
  real(kind=rk), intent(in) :: beta  ! Multiplier of b
  type(Vector), intent(inout) :: b   ! Output vector b = alpha*A*u + beta*b

  integer :: i,i_j,j
  
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !    Gather entire vector u on each processor
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  call sparsegather(u,A%BW)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calculate each component of b by taking the scalar product of the
!     i-th row of A and the vector u.
!
!       Note that in the compressed row storage format the nonzero 
!       entries of row i are stored in 
!
!         A%aa(A%ii(i)), A%aa(A%ii(i)+1), ..., A%aa(A%ii(i+1)-1)
!
!       the according (global) column numbers are stored in
!
!         A%jj(A%ii(i)), A%jj(A%ii(i)+1), ..., A%jj(A%ii(i+1)-1)   
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  !  b%xx = beta * b%xx
  call dscal(b%iend-b%ibeg+1,beta,b%xx(b%ibeg),1)
  
  do i =b%ibeg,b%iend
     do i_j = A%ii(i), A%ii(i+1)-1
        j = A%jj(i_j)
        b%xx(i) = b%xx(i) + alpha * A%aa(i_j) * u%xx(j)
     end do
  end do

end subroutine Mat_Mult
