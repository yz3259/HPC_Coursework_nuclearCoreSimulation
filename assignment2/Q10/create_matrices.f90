subroutine create_matrices(m, M1, M2, S, F1, F2,ibeg,iend)
  !====================================================
  ! A subroutine for assembling discretization matrices
  ! in the compressed row storage format
  !====================================================
  use header
  implicit none

  integer, intent(inout) :: m
  type(Matrix), intent(inout) :: M1, M2, S, F1, F2
  
  integer :: i,j, irow, inzd, inzm, n
  real(kind=rk) :: h, xl, xr, yt, yb
  integer, intent(inout) :: ibeg,iend
  
  ! External coefficient functions
  real(kind=rk) :: K1, K2, Sigma_a1, Sigma_a2, Sigma_s, Sigma_f1, Sigma_f2

  h = 1.0_rk/m
  n = m*m

  M1%n = n; M1%ibeg = ibeg; M1%iend = iend; M1%BW = m
  M2%n = n; M2%ibeg = ibeg; M2%iend = iend; M2%BW = m
  S%n  = n; S%ibeg  = ibeg; S%iend  = iend;! S%BW  = m
  F1%n = n; F1%ibeg = ibeg; F1%iend = iend;! F1%BW = m
  F2%n = n; F2%ibeg = ibeg; F2%iend = iend;! F2%BW = m
  S%BW = 0; F1%BW=0; F2%BW = 0; ! This is one optimisation. Due to the shape of
! those diagonal matrices S,F1,F2.
  allocate(M1%ii(n+1), M2%ii(n+1), S%ii(n+1), F1%ii(n+1), F2%ii(n+1))
  allocate(M1%aa(5*(iend-ibeg+1)), M2%aa(5*(iend-ibeg+1)),&
          S%aa(iend-ibeg+1), F1%aa(iend-ibeg+1), F2%aa(iend-ibeg+1))
  allocate(M1%jj(5*(iend-ibeg+1)), M2%jj(5*(iend-ibeg+1)),&
          S%jj(iend-ibeg+1), F1%jj(iend-ibeg+1), F2%jj(iend-ibeg+1))

  inzm = 1 ! index of current nonzero element in M1 and M2
  inzd = 1 ! index of current nonzero element in S, F1, F2
  do irow=ibeg,iend     ! this part is changed for the purpose of parallel computation
     ! Calculate cartesian index splitting
     j = (irow-1)/m + 1
     i = irow - (j-1)*m
     ! left, right, top, bottom midpoints
     xl = h*(i-0.5_rk)
     xr = xl + h
     yb = h*(j-0.5_rk)
     yt = yb + h

     ! Init current row position
     M1%ii(irow) = inzm
     M2%ii(irow) = inzm
     S%ii(irow) = inzd
     F1%ii(irow) = inzd
     F2%ii(irow) = inzd

     ! Diagonal elements of all matrices
     S%aa(inzd) = (Sigma_s(xl,yb)+Sigma_s(xr,yb)+Sigma_s(xl,yt)+Sigma_s(xr,yt))*0.25_rk
     S%jj(inzd) = irow
     F1%aa(inzd) = (Sigma_f1(xl,yb)+Sigma_f1(xr,yb)+Sigma_f1(xl,yt)+Sigma_f1(xr,yt))*0.25_rk
     F1%jj(inzd) = irow
     F2%aa(inzd) = (Sigma_f2(xl,yb)+Sigma_f2(xr,yb)+Sigma_f2(xl,yt)+Sigma_f2(xr,yt))*0.25_rk
     F2%jj(inzd) = irow
     
     M1%aa(inzm) = (K1(xl,yb)+K1(xr,yb)+K1(xl,yt)+K1(xr,yt))/(h**2) &
                 + (Sigma_a1(xl,yb)+Sigma_a1(xr,yb)+Sigma_a1(xl,yt)+Sigma_a1(xr,yt))*0.25_rk &
                 + S%aa(inzd)
     M1%jj(inzm) = irow
     M2%aa(inzm) = (K2(xl,yb)+K2(xr,yb)+K2(xl,yt)+K2(xr,yt))/(h**2) &
                 + (Sigma_a2(xl,yb)+Sigma_a2(xr,yb)+Sigma_a2(xl,yt)+Sigma_a2(xr,yt))*0.25_rk
     M2%jj(inzm) = irow

     ! Done with diagonal elements, shift the counter
     inzd = inzd + 1
     inzm = inzm + 1

     ! Off-diagonal elements in M1 and M2
     ! left
     if (i>1) then
        M1%aa(inzm) = -(K1(xl,yb)+K1(xl,yt))*0.5_rk/(h**2)
        M2%aa(inzm) = -(K2(xl,yb)+K2(xl,yt))*0.5_rk/(h**2)
        ! Neumann boundary conditions at x==1
        if (i==m) then
           M1%aa(inzm) = M1%aa(inzm) - (K1(xr,yb)+K1(xr,yt))*0.5_rk/(h**2)
           M2%aa(inzm) = M2%aa(inzm) - (K2(xr,yb)+K2(xr,yt))*0.5_rk/(h**2)
        end if
        M1%jj(inzm) = i-1 + (j-1)*m ! column index: "left" is (i-1,j)
        M2%jj(inzm) = i-1 + (j-1)*m ! column index
        ! shift the nnz counter
        inzm = inzm + 1
     end if
     
     ! right
     if (i<m) then
        M1%aa(inzm) = -(K1(xr,yb)+K1(xr,yt))*0.5_rk/(h**2)
        M2%aa(inzm) = -(K2(xr,yb)+K2(xr,yt))*0.5_rk/(h**2)
        M1%jj(inzm) = i+1 + (j-1)*m ! column index: "right" is (i+1,j)
        M2%jj(inzm) = i+1 + (j-1)*m ! column index
        ! shift the nnz counter
        inzm = inzm + 1
     end if
     
     ! bottom
     if (j>1) then
        M1%aa(inzm) = -(K1(xl,yb)+K1(xr,yb))*0.5_rk/(h**2)
        M2%aa(inzm) = -(K2(xl,yb)+K2(xr,yb))*0.5_rk/(h**2)
        ! Neumann boundary conditions at y==1
        if (j==m) then
           M1%aa(inzm) = M1%aa(inzm) - (K1(xl,yt)+K1(xr,yt))*0.5_rk/(h**2)
           M2%aa(inzm) = M2%aa(inzm) - (K2(xl,yt)+K2(xr,yt))*0.5_rk/(h**2)
        end if
        M1%jj(inzm) = i + (j-2)*m ! column index: "bottom" is (i,j-1)
        M2%jj(inzm) = i + (j-2)*m ! column index
        ! shift the nnz counter
        inzm = inzm + 1
     end if
     
     ! top
     if (j<m) then
        M1%aa(inzm) = -(K1(xl,yt)+K1(xr,yt))*0.5_rk/(h**2)
        M2%aa(inzm) = -(K2(xl,yt)+K2(xr,yt))*0.5_rk/(h**2)
        M1%jj(inzm) = i + j*m ! column index: "top" is (i,j+1)
        M2%jj(inzm) = i + j*m ! column index
        ! shift the nnz counter
        inzm = inzm + 1
     end if
  end do

  ! Finalise row positions
  M1%ii(irow) = inzm
  M2%ii(irow) = inzm
  S%ii(irow) = inzd
  F1%ii(irow) = inzd
  F2%ii(irow) = inzd

  ! Record the exact number of nonzeros
  M1%nnz = inzm-1
  M2%nnz = inzm-1
  S%nnz = inzd-1
  F1%nnz = inzd-1
  F2%nnz = inzd-1

  
end subroutine create_matrices
