subroutine propagate_trap(n, h, psi, tau, chi, Hpsi)
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------------
  ! Description: e^(-iHt)=1-iHt propagate the wavefunction at time t on psi
  ! on t + delta_t (using the trapezoidal rule)
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
  double complex, intent(in):: psi(n), Hpsi(n)
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
  integer :: igrid, ierr, i, j
  double complex, allocatable :: psitmp(:), sub(:), main(:), super(:), &
       & toeplitz(:,:), rhs(:)
  integer, allocatable :: ipiv(:)
  double complex :: tmp
  double complex, external :: scalar
  !------------------------------------------------------------------------------
  !  Local Constants 
  !------------------------------------------------------------------------------
  allocate(psitmp(n), sub(n), main(n), super(n), toeplitz(7,n), ipiv(n),&
       & rhs(n))
  
  ! compute [1 - i*tau/2 H]psi(t) and store on psitmp
  psitmp = psi + dcmplx(0d0,-0.5d0*tau)*Hpsi

  ! storing diagonals of [1 - tau/2i H] = [1 + i*tau/2 H]
  main = dcmplx(1d0, 0.125d0*tau/(h*h))
  sub = dcmplx(0d0, -0.0625d0*tau/(h*h))
  super = sub
  sub(1) = dcmplx(0d0, 0d0)
  sub(2) = dcmplx(0d0, 0d0)
  super(n-1) = dcmplx(0d0, 0d0)
  super(n) = dcmplx(0d0, 0d0)
  main(1) = 0.5d0*main(1) ! boundary conditions
  main(n) = 0.5d0*main(n)

  ! define toeplitz [1 - tau/2i H] = [1 + i*tau/2 H]
  toeplitz = dcmplx(0d0,0d0)
  do j = 1,n
     toeplitz(3,j) = super(j) ! superdiagonal
     toeplitz(5,j) = main(j) ! main diagonal
     toeplitz(7,j) = sub(j) ! subdiagonal
  end do
!  rhs = psitmp
  call zgbsv(n, 2, 2, 1, toeplitz, 7, ipiv, psitmp, n, ierr)
  print *, 'Info ZGBSV =', ierr
  chi = psitmp

  
!  print *, 'Residual norm =', dble(sqrt(scalar(n,h,chi,chi)))

  
  ! tridiagonal matrix algorithm
!!$  do igrid = 2, n
!!$     tmp = sub(igrid)/main(igrid-1)
!!$     main(igrid) = main(igrid) - tmp*super(igrid-1)
!!$     psitmp(igrid) = psitmp(igrid) - tmp*psitmp(igrid-1)
!!$  end do
!!$  chi(n) = psitmp(n)/main(n)
!!$  do igrid = n-1, 1, -1
!!$     chi(igrid) = (psitmp(igrid) - super(igrid)*chi(igrid+1))/main(igrid)
!!$  end do
  
  deallocate(psitmp, sub, main, super, toeplitz, ipiv, rhs)

end subroutine propagate_trap
