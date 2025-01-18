subroutine getcli(n, ntsteps, idebug, tau, x0, p0, alpha, mode)
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

  ! Set default values
  n = 100
  ntsteps = 10
  idebug = 0
  tau = 0.1d0
  x0 = 0.5d0
  p0 = 0.5d0
  alpha = 10d0
  mode = "exact" ! if mode = "help" ...

  ! Read command line arguments
  do ix = 1, num_args
     call get_command_argument(ix,arg)
     select case (arg)
     case ("--alpha", "-a")
        itmp = scan(arg(3:),"-.1234567890")
        if (itmp > 0) then
           read(arg(itmp:),*) alpha
        end if
     case ("--help", "-h")
        print*, "HELPPPPP!"
     end select
  end do

  
end subroutine getcli
