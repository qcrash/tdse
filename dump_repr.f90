subroutine dump_repr(n,h,psi,t,tau,chi,output,repr)
!------------------------------------------------------------------------------
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
  character(len=20), intent(in):: output, repr
  integer, intent(in):: n
  double precision, intent(in):: h, tau, t
  double complex, intent(in):: psi(n,2), chi(n,2)
!------------------------------------------------------------------------------
! Output Parameters
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Input/Output Parameters
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!  Local Variables
!------------------------------------------------------------------------------
  integer :: iunit
  double precision, external :: xval, variance, pval, pvar
!------------------------------------------------------------------------------
!  Local Constants 
!------------------------------------------------------------------------------
  if (output /= "") then
     iunit = 69
     open(iunit, file = output, position = "append")
  else
     iunit = 6
  end if

  write(iunit,*) "t = ", t
  select case (repr)
  case ("psisq")
     call dump_psi(n,h,iunit,psi)
  case ("sb")
     call segal_bargmann(n,h,iunit,psi,tau,chi)
  case ("wigner")
     call wigner_distr(n,h,iunit,psi,tau,chi)
  case ("observables")
     write(iunit,*) 'Wavepacket position =',xval(n,h,psi),xval(n,h,psi(1,2)) ! intial position expectation value
     write(iunit,*) 'Wavepacket position uncertainty =',variance(n,h,psi)&
          &,variance(n,h,psi(1,2)) ! initial position uncertainty
!!$     write(iunit,*) 'Wavepacket momentum =',pval(n,h,1,psi) ! initial momentum expectation value
     write(iunit,*) 'Wavepacket momentum =',pval(n,h,0,psi),pval(n,h,0,psi(1,2)) ! initial momentum expectation value
     write(iunit,*) 'Wavepacket momentum uncertainty =',pvar(n,h,psi),pvar(n&
          &,h,psi(1,2)) ! initial momentum uncertainty
     write(iunit,*) 'Wavepacket kinetic energy =',(pvar(n,h,psi)**2 + pval(n,h,0&
          &,psi)**2)/2d0,(pvar(n,h,psi(1,2))**2 + pval(n,h,0,psi(1,2))**2)/2d0
     write(iunit,*) 'Uncertainty product =',variance(n,h,psi)*pvar(n,h,psi)&
          &,variance(n,h,psi(1,2))*pvar(n,h,psi(1,2)) ! must be greater
     ! than 1/2
     write(iunit,*)
  end select
  
  if (iunit /= 6) close(iunit)
end subroutine dump_repr
