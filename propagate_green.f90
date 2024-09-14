subroutine propagate_green(n, h, psi, tau, chi)
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------------
  ! Description: e^(-iHt)=1-iHt propagate the wavefunction at time t on psi
  ! on t + delta_t (using the Green's function)
  ! psi(t + tau) = psi(t) - i*tau/2 H[psi(t) + psi(t + tau)] OR
  ! [1 + i*tau/2 H]psi(t + tau) = [1 - i*tau/2 H]psi(t)
  !------------------------------------------------------------------------------
  !
  implicit none
  !------------------------------------------------------------------------------
  ! Input Parameters
  !------------------------------------------------------------------------------
  integer, intent(in):: n
  double precision, intent(in):: h, tau
  double complex, intent(in):: psi(n)
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
  integer :: i, j
!!$  double complex, allocatable ::
!!$  integer, allocatable :: 
  double complex :: tmp, tmp2, xtemp, xptemp, prefactor
  double complex, external :: scalar
  !------------------------------------------------------------------------------
  !  Local Constants 
  !------------------------------------------------------------------------------
  double precision, parameter :: pi = 4d0*atan(1d0)
  
!!$  allocate()
  
  do i=1,n
     xtemp = dble(i)*h - 1d0
     tmp2 = dcmplx(0d0,0d0)
     do j=1,n
        xptemp = dble(j)*h - 1d0
        tmp = exp(-(xtemp-xptemp)**2/dcmplx(0d0,2d0*tau))
        tmp2 = tmp2 + tmp*psi(j)
     end do
     chi(i) = tmp2 
  end do
  prefactor = 1/sqrt(dcmplx(0d0,2*pi*tau))
  chi = prefactor*h*chi
  
!!$  deallocate()

end subroutine propagate_green
