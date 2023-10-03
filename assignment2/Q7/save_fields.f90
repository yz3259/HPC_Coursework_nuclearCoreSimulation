subroutine save_fields(U1,U2,filename)
  ! ===========================================================
  ! Save solutions to disk.
  ! Writes a text file which contains the values of the vectors
  ! u_1 and u_2 at the nodal points.
  ! ===========================================================

  use header
  implicit none
  ! Velocity field $u$
  type(Vector), intent(in) :: U1, U2
  ! Name of file to write to
  character(len=*), intent(in) :: filename
  ! Loop variable
  integer :: i
  integer :: file_id

  file_id = 99 ! avoid conflicts if file of parameters is not closed

  ! Write first component
  open(unit=file_id,file=trim(filename))
  do i=1, U1%n
     write(file_id,'(E20.8e3," ")',advance='no') U1%xx(i)
  end do
  write(file_id,'("")') ! delimiter
  ! Write second component
  do i=1, U2%n
     write(file_id,'(E20.8e3," ")',advance='no') U2%xx(i)
  end do
  write(file_id,'("")')
  close(file_id)

end subroutine save_fields
