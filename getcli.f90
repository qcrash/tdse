subroutine getcli(n, ntsteps, idebug, tau, x0, p0, t0, alpha, mode, output, repr)
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
  double precision, intent(out):: tau, x0, p0, t0, alpha
  character(len=20), intent(out):: mode, repr
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
  arg = ""
  
  ! Set default values
  n = 100
  ntsteps = 10
  idebug = 0
  tau = 0.1d0
  x0 = 0.5d0
  p0 = 0.5d0
  t0 = 0d0
  alpha = 10d0
  mode = "exact" ! if mode = "help" ...
  output = ""
  repr = "psi"

  ! Read command line arguments
  ix = 0
  do while (ix .le. num_args)
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
     else if (arg(1:13) == "--initialtime" .or. arg(1:2) == "-t") then
        t0 = check_arg_double(arg,ix,num_args,length,ierr)
     else if (arg(1:12) == "--mode" .or. arg(1:2) == "-m") then
        mode = check_arg_char(arg,ix,num_args,length,ierr)
     else if (arg(1:8) == "--output" .or. arg(1:2) == "-o") then
        output = check_arg_char(arg,ix,num_args,length,ierr)
     else if (arg(1:16) == "--representation" .or. arg(1:2) == "-r") then
        repr = check_arg_char(arg,ix,num_args,length,ierr)
     else
        print '(A)', "Usage: tdse [OPTION] [VALUE]" // new_line('A') // &
             & "" // new_line('A') // &
             & "Solve time-dependent Schroedinger equation for a " // &
             & "particle in a 1D box" // new_line('A') // &
             & "All quantities are in atomic (Hartree) units" // &
             & new_line('A') // &
             & "" // new_line('A') // &
             & "Options:" // new_line('A') // &
             & "" // new_line('A') // &
             & "" // new_line('A') // &             
             & "--alpha, -a                     " // &
             & "set the reciprocal wavepacket width" // new_line('A') // &
             & "                                default: 10d0" // &
             & new_line('A') // &
             & "" // new_line('A') // &             
             & "--debug, -d                     " // &
             & "set debug level" // new_line('A') // &
             & "                                default: 0" // &
             & new_line('A') // &
             & "" // new_line('A') // &             
             & "--gridpoints, -n                " // &
             & "set the number of points to use in discretizing" // &
             & new_line('A') // &
             & "                                real space (must be an" // &
             & " integer)"// new_line('A') // &
             & "                                default: 100" // &
             & new_line('A') // &
             & "" // new_line('A') // &             
             & "--ntsteps                       " // &
             & "set the number of time steps (must be an integer)" // &
             & new_line('A') // &
             & "                                default: 10" // &
             & "" // new_line('A') // &                          
             & new_line('A') // "--timestep, --tau               " // &
             & "set the size of the time step when discretizing time" // &
             & new_line('A') // &
             & "                                default: 0.1d0" // &
             & "" // new_line('A') // &             
             & new_line('A') // &
             & "--initialposition, -x           " // &
             & "set the wavepacket initial position" // new_line('A') // &
             & "                                must be within" // &
             & "[-1d0,+1d0]" // &
             & new_line('A') // &
             & "                                default: 0.5d0" // &
             & "" // new_line('A') // &             
             & new_line('A') // &
             & "--initialmomentum, -p           " // &
             & "set the wavepacket initial momentum" // new_line('A') // &
             & "                                default: 0.5d0" // &
             & "" // new_line('A') // &             
             & new_line('A') // &
             & "--initialtime, -t               " // &
             & "set the start time" // new_line('A') // &
             & "                                default: 0d0" // &
             & "" // new_line('A') // &             
             & new_line('A') // &
             & "--mode, -m                      " // &
             & "set the propagation scheme used" // new_line('A') // &
             & "                                options include: " // &
             & "'central' 'ab' 'fft' 'green' 'trap' " // &
             & new_line('A') // &
             & "                                default: exact" // &
             & "" // new_line('A') // &             
             & new_line('A') // &
             & "--output, -o                    " // &
             & "set the name of the output file" // new_line('A') // &
             & "                                default: standard output" // &
             & "" // new_line('A') // &             
             & new_line('A') // &
             & "--representation, -r            " // &
             & "choose type of output data " // &
             & new_line('A') // &
             & "                                options include: " // &
             & "'psi' 'sb' 'wigner' 'observables'" // &
             & "                                    default: psi"    
        stop
     end if
     ix = ix + 1
  end do

  
end subroutine getcli
