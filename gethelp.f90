  subroutine gethelp
  !-----------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------------
  ! Description: 
  !------------------------------------------------------------------------------
    ! Print help text to stdout
    !
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
  !------------------------------------------------------------------------------
  ! Input/Output Parameters
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !  Local Variables
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  !  Local Constants 
  !------------------------------------------------------------------------------

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

end subroutine gethelp
