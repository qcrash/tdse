subroutine propagate_convert(n, h, psi, tau, chi)
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------------
  ! Description: 
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
  integer :: i, j, istatus
  double precision :: tmp
  double precision, parameter :: pi = 4d0*atan(1d0)
  double complex, allocatable :: tmppsi(:)
  !------------------------------------------------------------------------------
  !  Local Constants 
  !------------------------------------------------------------------------------

  allocate(tmppsi(n))

  ! Applying eigenvector matrix of momentum operator (transform to
  ! momentum basis)
  chi = 0d0
  do i = 1,n
     do j = 1,i-1
        tmp = sin(dble(i*j)*pi/dble(n+1))
        chi(i) = chi(i) + tmp*psi(j)
        chi(j) = chi(j) + tmp*psi(i)
     end do
     chi(i) = chi(i) + sin(dble(i*i)*pi/dble(n+1))*psi(i)
  end do
  
  do i = 1,n
     tmp = cos(pi*dble(i)/dble(n+1))/h ! momentum eigenvalue
     tmppsi(i) = chi(i)*exp(cmplx(0d0,-tmp*tmp*tau*0.5d0)) ! propagated in
     ! momentum basis
  end do

  ! Applying eigenvector matrix of momentum operator (transform to
  ! position basis)
  chi = 0d0
  do i = 1,n
     do j = 1,i-1
        tmp = sin(dble(i*j)*pi/dble(n+1))
        chi(i) = chi(i) + tmp*tmppsi(j)
        chi(j) = chi(j) + tmp*tmppsi(i)
     end do
     chi(i) = chi(i) + sin(dble(i*i)*pi/dble(n+1))*tmppsi(i)
  end do
  
  tmp = 2d0/dble(n+1) ! multiply by prefactor of eigenvector matrix
  ! element twice
  chi = chi*tmp
  
  deallocate(tmppsi)
end subroutine propagate_convert
