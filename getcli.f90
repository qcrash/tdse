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
  integer :: num_args, ix, itmp
  character(len=40) :: arg
  !------------------------------------------------------------------------------
  !  Local Constants 
  !------------------------------------------------------------------------------
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
        itmp = scan(arg(3:),"-.1234567890") + 2
        if (itmp > 0) then
           read(arg(itmp:),*) alpha
        end if
     else if (arg(1:12) == "--gridpoints" .or. arg(1:2) == "-n") then
        itmp = scan(arg(3:),"-.1234567890") + 2
        if (itmp > 0) then
           read(arg(itmp:),*) n
        end if
     else if (arg(1:9) == "--ntsteps") then
        itmp = scan(arg(3:),"-.1234567890") + 2
        if (itmp > 0) then
           read(arg(itmp:),*) ntsteps
        end if
     else if (arg(1:7) == "--debug" .or. arg(1:2) == "-d") then
        itmp = scan(arg(3:),"-.1234567890") + 2
        if (itmp > 0) then
           read(arg(itmp:),*) idebug
        end if
     else if (arg(1:5) == "--tau" .or. arg(1:10) == "--timestep") then
        itmp = scan(arg(3:),"-.1234567890") + 2
        if (itmp > 0) then
           read(arg(itmp:),*) tau
        end if
     else if (arg(1:17) == "--initialposition" .or. arg(1:2) == "-x") then
        itmp = scan(arg(3:),"-.1234567890") + 2
        if (itmp > 0) then
           read(arg(itmp:),*) x0
        end if
     else if (arg(1:17) == "--initialmomentum" .or. arg(1:2) == "-p") then
        itmp = scan(arg(3:),"-.1234567890") + 2
        if (itmp > 0) then
           read(arg(itmp:),*) p0
        end if
     else if (arg(1:12) == "--mode" .or. arg(1:2) == "-m") then
        itmp = scan(arg(3:),'"') + 2
        if (itmp > 0) then
           read(arg(itmp:),*) mode
        end if
     else if (arg(1:8) == "--output" .or. arg(1:2) == "-o") then
        itmp = scan(arg(3:),'"') + 2
        if (itmp > 0) then
           read(arg(itmp:),*) output
        end if
     else
        print*, "HELPPPPP!"
        stop
     end if
  end do

  
end subroutine getcli
