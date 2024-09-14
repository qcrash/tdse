subroutine wigner_distr(n, h, iunit, psi, tau, chi)
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------------
  ! Description: 
  !------------------------------------------------------------------------------
  ! Calculates the Wigner distribution using a triple nested loop: one
  ! over x, k, and xi where...
  ! W(x,k_n) = (2/pi)*sum[psi*(x+xi_m)psi(x-xi_m)exp(2pi*xi_m*k_n*i)]
  ! from 1 to min(x,n-x)
  implicit none
  !------------------------------------------------------------------------------
  ! Input Parameters
  !------------------------------------------------------------------------------
  integer, intent(in):: n, iunit
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
  integer :: i, j, k, istatus, kmax
  double precision :: tmp
  double precision, parameter :: pi = 4d0*atan(1d0)
  double precision, allocatable :: w(:)
  character :: fmt*40, fmt0*40
  double precision :: xi, eta
  !------------------------------------------------------------------------------
  !  Local Constants 
  !------------------------------------------------------------------------------

  allocate(w(n))
  
  write(fmt,*) "(f16.10,",n,"(2x, f16.10))"
  write(fmt0,*) "(i10,",n,"(2x, f16.10))"
  
  write(70,fmt0) n+1, (h*dble(j) - 1d0, j =1, n)
  do i = 1,n ! x loop
     do j = 1,n ! k loop
        tmp = 0d0
        eta = -1d0 + dble(2*j)/dble(n+1)

        kmax = min(i, n+1-i) - 1

        do k = 1,kmax
           xi = dble(2*k)/dble(n+1)
           tmp = tmp + real(conjg(psi(i+k))*psi(i-k)*exp(dcmplx(0d0,2d0*pi&
                &*xi*eta)))
        end do
        
        tmp = 2d0*tmp
        tmp = tmp + real(psi(i))**2 + aimag(psi(i))**2
        
           
        tmp = tmp*h*2d0/pi
        w(j) = tmp
     end do
     write(70,fmt) h*dble(i) - 1d0,(w(j), j = 1,n)
  end do
  
  deallocate(w)
end subroutine wigner_distr
