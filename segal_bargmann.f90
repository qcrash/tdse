subroutine segal_bargmann(n, h, iunit, psi, tau, chi)
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------------
  ! Description: 
  !------------------------------------------------------------------------------
  ! Calculates the Segal Bargmann trasnform using a triple nested loop: one
  ! over x, k, and xi where...
  ! psi_bar(x,z) = e^(0.25*z^2)*sum_m[(psi(xi_m)exp(-0.5*(z - xi_m)^2)]
  ! from 1 to n
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
  integer :: i, j, k, istatus
  double precision, parameter :: pi = 4d0*atan(1d0), sqpinv = 1d0&
       &/sqrt(pi), sqrt2 = sqrt(2d0)
  double complex, allocatable :: psi_bar(:,:)
  character :: fmt*60, fmt0*40, fmt1*40
  double precision :: xi, eta, upsilon
  double complex :: z, tmp
  !------------------------------------------------------------------------------
  !  Local Constants 
  !------------------------------------------------------------------------------

  allocate(psi_bar(n,n))
  
!!$  write(fmt,*) "(f16.10,",n,"(2x, f16.10, 1x, f16.10))"
  write(fmt1,*) "(f16.10,",n,"(2x, f16.10))"
  write(fmt0,*) "(i10,",n,"(2x, f16.10))"
  
  write(70,fmt0) n+1, (h*dble(j) - 1d0, j =1, n)
  do i = 1,n ! x loop
     upsilon = -1d0 + dble(2*i)/dble(n+1) ! real part of z
     do j = 1,n ! k loop
        tmp = dcmplx(0d0,0d0)
        eta = -1d0 + dble(2*j)/dble(n+1) ! imaginary part of z

        z = sqrt2*dcmplx(upsilon, eta)

        do k = 1, n
           xi = dble(2*k)/dble(n+1)
           tmp = tmp + psi(k)*exp(-0.5d0*(z - dcmplx(xi,0d0))**2)
        end do
        
        tmp = exp(0.5d0*(0.5d0*z**2-abs(z)**2))*tmp
        
           
        tmp = tmp*h*sqpinv
        psi_bar(j,i) = tmp
     end do
!!$     write(70,fmt) h*dble(i) - 1d0,(real(psi_bar(j,i)),
     !!aimag(psi_bar(j,i)), j = 1,n)
     write(70, fmt1) h*dble(i) - 1d0, (abs(psi_bar(j,i))**2, j = 1,n)
  end do
  
  deallocate(psi_bar)
end subroutine segal_bargmann
