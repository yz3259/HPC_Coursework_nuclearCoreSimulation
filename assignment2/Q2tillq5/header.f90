!----------------------------------------------------------------- 
!
!  Module with general definitions and data structures
!
!----------------------------------------------------------------- 

module header

  ! Size of floating point number
  integer, parameter :: rk=8
  ! Numerical value of pi
  real(kind=rk) :: PI = 3.14159265358979323846264338327950288_rk

  ! Declare BLAS routines as external
  external daxpy, dcopy, dscal
  real(kind=rk), external :: dnrm2, ddot

  ! Data type for matrix in compressed row storage (CRS)
  type Matrix
    real (kind=rk), allocatable, dimension(:) :: aa
    integer, allocatable, dimension(:) :: ii
    integer, allocatable, dimension(:) :: jj
    integer :: n
    integer :: nnz
    integer :: ibeg 
    integer :: iend
    integer :: BW
  end type Matrix

  ! Data type for vector
  type Vector
    real (kind=rk), allocatable, dimension(:) :: xx
    integer :: n
    integer :: ibeg 
    integer :: iend
  end type Vector

end module header
