subroutine propagate_ab(n, h, psi, psi_old, tau, chi)
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------------
  ! Description: e^(-iHt)=1-iHt propagate the wavefunction at time t on psi
  ! on t + delta_t (same as propagate subroutine but using two-step
  ! Adams-Bashforth)
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
  double complex :: tmp1, tmp2, tmp3
  !------------------------------------------------------------------------------
  !  Local Constants 
  !------------------------------------------------------------------------------
  
  ! second derivative loop
  do igrid = 2, n-1
     tmp1 = 3d0*psi(igrid) - psi_old(igrid)
     tmp2 = 3d0*psi(igrid-1) - psi_old(igrid-1)
     tmp3 = 3d0*psi(igrid+1) - psi_old(igrid+1)
     chi(igrid) = -tmp2 + 2d0*tmp1 - tmp3
  end do

  ! first and last points special case
  chi(1) = 6d0*psi(1) - 3d0*psi(2) - 2d0*psi_old(1) + psi_old(2)
  chi(n) = -3d0*psi(n-1) + 6d0*psi(n) + psi_old(n-1) - 2d0*psi_old(n)

  ! combine psi: psi(t+tau) = (2h/i)*Hpsi(t) + psi(t-tau)
  chi = psi + dcmplx(0d0,-tau*0.5d0/(h*h))*chi
end subroutine propagate_ab
