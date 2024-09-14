subroutine propagate_central(n, h, psi, psi_old, tau, chi)
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------------
  ! Description: e^(-iHt)=1-iHt propagate the wavefunction at time t on psi
  ! on t + delta_t (same as propagate subroutine but using central difference)
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
  integer, intent(in):: n
  double precision, intent(in):: h, tau
  double complex, intent(in):: psi(n), psi_old(n)
  !------------------------------------------------------------------------------
  ! Output Parameters
  !------------------------------------------------------------------------------
  double complex, intent(out):: chi(n)
  !------------------------------------------------------------------------------
  ! Input/Output Parameters
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !  Local Variables
  !------------------------------------------------------------------------------
  integer :: igrid
  !------------------------------------------------------------------------------
  !  Local Constants 
  !------------------------------------------------------------------------------
  

  ! second derivative loop
  do igrid = 2, n-1
     chi(igrid) = -psi(igrid-1) + 2d0*psi(igrid) - psi(igrid+1)
  end do

  ! first and last points special case
  chi(1) = 2d0*psi(1) - psi(2)
  chi(n) = -psi(n-1) + 2d0*psi(n)

  ! combine psi: psi(t+tau) = (2h/i)*Hpsi(t) + psi(t-tau)
  chi = psi_old - dcmplx(0d0,tau/(h*h))*chi

end subroutine propagate_central
