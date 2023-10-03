subroutine Vec_AXPY(alpha, x1,x2, y1,y2)
  !=====================================================
  ! AXPY operation [y1] = a*[x1] + [y1] of block vectors
  !                [y2]     [x2]   [y2]
  ! With only local part for each processor
  !=====================================================
  use header
  implicit none
  
  real(kind=8), intent(in) :: alpha
  type(Vector), intent(inout) :: x1,x2
  type(Vector), intent(inout) :: y1,y2
    
!  print*,y1%ibeg
!  print*,y1%iend
  call daxpy(x1%iend-x1%ibeg+1, alpha, x1%xx(x1%ibeg), 1, y1%xx(y1%ibeg), 1)
  call daxpy(x1%iend-x1%ibeg+1, alpha, x2%xx(x2%ibeg), 1, y2%xx(y2%ibeg), 1)
end subroutine Vec_AXPY
