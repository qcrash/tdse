program tdse
  implicit none
  !! Declaring variables
  integer, parameter :: idebug = 1
  integer :: info, i, j, n, ntsteps, itsteps ! LAPACK status, loop vars, vector space dim
  double precision :: h, psinorm, tau, t ! grid spacing, time step, time param
  double precision :: hinv, imz0sq, xtemp, exptemp ! temp vars
  double precision, parameter :: x0 = -0.8d0, p0 = 0.3d0, alpha = 1000d0, &
       & pi = 4d0*atan(1d0)
  double complex :: temp
  double complex, allocatable :: psi(:), psi0(:), chi(:), &
       & psi_old(:) ! psi_old is wavefunction in real orthonormal basis
  double precision, allocatable :: ham(:,:), tkin(:,:), &
       & vpot(:), work(:)
  double precision, external :: dznrm2
  double complex, external :: scalar
  double precision, external :: xval, variance, pval, pvar
  
  !! User-defined parameters
  print *, 'Number of grid points ='
  read (*,*) n
  ! Include user prompts for x0, p0, alpha
  h = 2d0/dble(n-1)

  !! Initializing time variables
  print *, 'Number of time steps ='
  read (*,*) ntsteps
  print *, 'Size of time step ='
  read (*,*) tau
  t = 0d0

  !! Allocating higher order tensors
  allocate (psi(n),psi0(n),psi_old(n),ham(2,n),tkin(2,n),vpot(n), &
       & work(3*n),chi(n))
  
  !! Discretizing initial wavepacket
  imz0sq = 0.25d0*(p0/alpha)**2
  psi0(1) = 0d0
  psi0(n) = 0d0 ! boundary conditions
  psi0_loop: do i = 2, n-1
     xtemp = dble(i-1)*h - 1d0 - x0
     exptemp = exp(-alpha*(xtemp*xtemp - imz0sq))
     psi0(i) = dcmplx(exptemp*cos(p0*xtemp), exptemp*sin(p0*xtemp))
     ! psi0(i) = dcmplx(1d0,0d0) ! debug wavepacket
  end do psi0_loop
  
  !! Normalizing initial wavepacket
  psinorm = dble(sqrt(scalar(n,h,psi0,psi0)))
  call zdscal(n,1d0/psinorm,psi0,1)
  if (idebug > 0) then
     print *, 'Initial wavepacket norm =',psinorm
     psinorm = dble(sqrt(scalar(n,h,psi0,psi0)))
     print *, 'Normalized wavepacket norm =',psinorm
  end if

  !! Testing initial position and variance
  print *, 'Initial wavepacket position =',xval(n,h,psi0) ! intial position expectation value
  print *, 'Initial wavepacket position uncertainty =',variance(n,h,psi0) ! initial position uncertainty
  print *, 'Initial wavepacket momentum =',pval(n,h,psi0) ! initial momentum expectation value
  print *, 'Initial wavepacket momentum uncertainty =',pvar(n,h,psi0) ! initial momentum uncertainty
  print *, 'Initial uncertainty product =',variance(n,h,psi0)*pvar(n,h,psi0) ! must be greater than 1/2
     
  !! Defining Hamiltonian operator
  vpot = 0d0 ! local potential at grid points
  hinv = 1d0/(h*h)
  t_loop_diagonal: do i = 1,n
     tkin(1,i) = hinv ! diagonal
     tkin(2,i) = -0.5d0*hinv ! offdiagonal
     ! tkin(2,i) = 0d0 ! debug
  end do t_loop_diagonal
  h_loop: do i = 1,n
     ham(1,i) = tkin(1,i)+vpot(i)
     ham(2,i) = tkin(2,i)
  end do h_loop

  !! Graphing eigenstates
  ! print *, "Static Hamiltonian eigenvalues :", omega
  ! do i=1,n
  !    write (88,*) dble(i-1)*h - 1d0, statevec(i,n)
  ! end do

  psi = psi0 ! save initial wavefunction
  ! propagate backward by 1 time step for psi(t-tau) or psi_old
  call propagate(n, h, psi, -tau, psi_old)

  !! Propagating wavefunction by one time step  
  do itsteps = 1, ntsteps
     !     call propagate(n, h, psi, tau, chi)
     call propagate_ab(n, h, psi, psi_old, tau, chi)
     ! save old wavefunction
     psi_old = psi
     psi = chi ! use result  as new input in next iteration

     
!     do i = 1, n
!       print *, 'Probability amplitude =', abs(psi(i))**2
!       psi(i) = dcmplx(cos(omega(i)*tau),-sin(omega(i)*tau))*psi(i) ! propagation
!     end do

     !! Testing initial position and variance
     print *, 'Absolute time =', t + tau*itsteps ! current time
     
     psinorm = dble(sqrt(scalar(n,h,psi,psi)))
     print *, 'norm before normalized: ', psinorm
     call zdscal(n,1d0/psinorm,psi,1)
     
     
     print *, 'Wavepacket position =',xval(n,h,psi) ! intial position expectation value
     print *, 'Wavepacket position uncertainty =',variance(n,h,psi) ! initial position uncertainty
     print *, 'Wavepacket momentum =',pval(n,h,psi) ! initial momentum expectation value
     print *, 'Wavepacket momentum uncertainty =',pvar(n,h,psi) ! initial momentum uncertainty
     print *, 'Uncertainty product =',variance(n,h,psi)*pvar(n,h,psi) ! must be greater than 1/2
     print *

     
     
  end do
  
  !! Deallocation of higher order tensors
  deallocate (psi,psi0,psi_old,ham,tkin,vpot,work,chi)
end program tdse

double complex function scalar(n,h,psi1,psi2)
  ! Calculate scalar product <psi1|psi2> on the grid
  ! using trapezoidal rule and the boundary conditions
  ! <psi1|psi2> = int(psi1*(x)psi2(x)dx)
  
  !! Defining parameters and local variables
  implicit none
  integer, intent(in):: n
  double precision, intent(in) :: h
  double complex, intent(in) :: psi1(n), psi2(n)
  integer :: igrid

  !! Calculating scalar product
  scalar = 0.5d0*(conjg(psi1(1))*psi2(1) &
       & + conjg(psi1(n))*psi2(n)) ! First and last point
  do igrid = 2, n-1
     scalar = scalar + conjg(psi1(igrid))*psi2(igrid)
  end do 
  scalar = scalar*h
end function scalar

double precision function xval(n,h,psi)
  implicit none
  integer, intent(in):: n
  double precision, intent(in) :: h
  double complex, intent(in) :: psi(n)
  integer :: igrid

  xval = 0.5d0*(-conjg(psi(1))*psi(1) &
       & + conjg(psi(n))*psi(n)) ! First and last point
  do igrid = 2, n-1
     xval = xval + conjg(psi(igrid))*psi(igrid)*(dble(igrid-1)*h - 1d0)
  end do 
  xval = xval*h
end function xval

double precision function variance(n,h,psi)
  implicit none
  integer, intent(in):: n
  double precision, intent(in) :: h
  double complex, intent(in) :: psi(n)
  double precision :: tmp
  double precision, external :: xval
  integer :: igrid
  
  ! var^2 = expt[(x-expt(x))^2]
  tmp = xval(n,h,psi)
  variance = conjg(psi(1))*psi(1)*(-0.5d0 - tmp)**2 &
       & + conjg(psi(n))*psi(n)*(0.5d0 - tmp)**2 ! First and last point
  do igrid = 2, n-1
     variance = variance + conjg(psi(igrid))*psi(igrid)*(dble(igrid-1)*h - 1d0 - tmp)**2
  end do
  variance = sqrt(variance*h)
end function variance

double precision function pval(n,h,psi)
  implicit none
  integer, intent(in):: n
  double precision, intent(in) :: h
  double complex, intent(in) :: psi(n)
  integer :: igrid

  pval = 0.5d0*(aimag(conjg(psi(1))*psi(2))/(2d0*h) &
       & - aimag(conjg(psi(n))*psi(n-1))/(2d0*h)) ! First and last point
  do igrid = 2, n-1
     pval = pval + aimag(conjg(psi(igrid))*(psi(igrid+1)-psi(igrid-1))/(2d0*h))
  end do 
  pval = pval*h
end function pval

double precision function pvar(n,h,psi)
  implicit none
  integer, intent(in):: n
  double precision, intent(in) :: h
  double complex, intent(in) :: psi(n)
  double precision :: tmp
  double precision, external :: pval
  integer :: igrid
  
  ! var^2 = expt[(p-expt(p))^2] = (abs(p*psi - expt(p)*psi))^2
  tmp = pval(n,h,psi)
  pvar = 0.5d0*((abs((psi(2))/(2d0*h) - tmp*psi(1)))**2 &
       & + (abs((-psi(n-1))/(2d0*h) - tmp*psi(n)))**2) ! First and last point
  do igrid = 2, n-1
     pvar = pvar + (abs(((psi(igrid+1)-psi(igrid-1)))/(2d0*h) - tmp*psi(igrid)))**2
  end do
  pvar = sqrt(pvar*h)
end function pvar


subroutine propagate(n, h, psi, tau, chi)
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------------
  ! Description: e^(-iHt)=1-iHt propagate the wavefunction at time t on psi
  ! on t + delta_t
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

  ! combine psi: (1-iHt)*psi = psi - iHpsit
  chi = psi - dcmplx(0d0,tau/(2d0*h*h))*chi

end subroutine



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

end subroutine
