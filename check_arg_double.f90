double precision function check_arg_double(arg,iarg,narg,length,ierr)
  !! Defining parameters and local variables
  implicit none
  integer, intent(in):: narg, length
  integer, intent(inout):: iarg
  character(len=length), intent(in):: arg
  character(len=length) :: argtmp
  integer :: ierr, itmp

  ierr = 0
  itmp = scan(arg(3:),"-.1234567890") + 2
  if (itmp > 2) then
     print *, arg
     read(arg(itmp:),*) check_arg_double
  else if (itmp == 2 .and. iarg /= narg) then
     call get_command_argument(iarg + 1,argtmp)
     print *, argtmp
     read(argtmp,*) check_arg_double
     iarg = iarg + 1
  else
     ierr = 1
  end if
end function check_arg_double
