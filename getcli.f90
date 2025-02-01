subroutine getcli(n, ntsteps, idebug, tau, x0, p0, alpha, mode, output)
  !-----------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------------
  ! Description: 
  !------------------------------------------------------------------------------
  ! (What is the purpose of this subroutine? State the theory and the
  ! algorithm used in a few sentences. Refer to literature or
  ! documentation for detailed derivations, proofs, or pseudocode.
  !------------------------------------------------------------------------------
  !
  implicit none
  !------------------------------------------------------------------------------
  ! Modules and Global Variables
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  ! External Functions
  !------------------------------------------------------------------------------
  double precision, external :: check_arg_double
  integer, external :: check_arg_int
  character(len=40), external :: check_arg_char
  !------------------------------------------------------------------------------
  ! Input Parameters
  !------------------------------------------------------------------------------
  
  !------------------------------------------------------------------------------
  ! Output Parameters
  !------------------------------------------------------------------------------
  integer, intent(out):: n, ntsteps, idebug
  double precision, intent(out):: tau, x0, p0, alpha
  character(len=20), intent(out):: mode
  character(len=80), intent(out):: output
  !------------------------------------------------------------------------------
  ! Input/Output Parameters
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !  Local Variables
  !------------------------------------------------------------------------------
  integer :: num_args, ix, itmp, ierr
  character(len=40) :: arg
  !------------------------------------------------------------------------------
  !  Local Constants 
  !------------------------------------------------------------------------------
  integer, parameter :: length=40
  ! Allocatables
  num_args = command_argument_count()
  if (num_args == 0) then
     print*, "HELPPPPP function!"
     stop
  end if

  ! Set default values
  n = 100
  ntsteps = 10
  idebug = 0
  tau = 0.1d0
  x0 = 0.5d0
  p0 = 0.5d0
  alpha = 10d0
  mode = "exact" ! if mode = "help" ...
  output = ""

  ! Read command line arguments
  do ix = 1, num_args
     call get_command_argument(ix,arg)
     if (arg(1:7) == "--alpha" .or. arg(1:2) == "-a") then
        alpha = check_arg_double(arg,ix,num_args,length,ierr)
     else if (arg(1:12) == "--gridpoints" .or. arg(1:2) == "-n") then
        n = check_arg_int(arg,ix,num_args,length,ierr)
     else if (arg(1:9) == "--ntsteps") then
        ntsteps = check_arg_int(arg,ix,num_args,length,ierr)
     else if (arg(1:7) == "--debug" .or. arg(1:2) == "-d") then
        idebug = 1
        idebug = check_arg_int(arg,ix,num_args,length,ierr)
     else if (arg(1:5) == "--tau" .or. arg(1:10) == "--timestep") then
        tau = check_arg_double(arg,ix,num_args,length,ierr)
     else if (arg(1:17) == "--initialposition" .or. arg(1:2) == "-x") then
        x0 = check_arg_double(arg,ix,num_args,length,ierr)
     else if (arg(1:17) == "--initialmomentum" .or. arg(1:2) == "-p") then
        p0 = check_arg_double(arg,ix,num_args,length,ierr)
     else if (arg(1:12) == "--mode" .or. arg(1:2) == "-m") then
        mode = check_arg_char(arg,ix,num_args,length,ierr)
     else if (arg(1:8) == "--output" .or. arg(1:2) == "-o") then
        output = check_arg_char(arg,ix,num_args,length,ierr)
     else
        print*, "HELPPPPP!"
        stop
     end if
  end do

  
end subroutine getcli
