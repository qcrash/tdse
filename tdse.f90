program tdse
  implicit none  
  !! Declaring variables
  integer, parameter :: idebug = 1
  integer :: info, i, j, n, ntsteps, itsteps ! LAPACK status, loop vars, vector space dim
  double precision :: h, psinorm, tau, t ! grid spacing, time step, time param
  double precision :: hinv, imz0sq, xtemp, exptemp, x0, p0, alpha ! temp vars
  double precision, parameter :: pi = 4d0*atan(1d0)
  double complex :: temp
  double complex, allocatable :: psi(:,:), psi0(:,:), chi(:,:), &
       & psi_old(:,:), work(:) ! psi_old is wavefunction in real orthonormal basis
  double precision, allocatable :: ham(:,:), tkin(:,:), &
       & vpot(:)
  double precision, external :: dznrm2
  double complex, external :: scalar
  double precision, external :: xval, variance, pval, pvar
  character*15 :: filename
  
  !! User-defined parameters
  print *, 'Number of grid points ='
  read (*,*) n ! Read number of grid points BETWEEN -1 and 1 (not including -1 and 1)
  ! Include user prompts for x0, p0, alpha in the future
  h = 2d0/dble(n+1)

  !! Initializing time variables
  print *, 'Number of time steps ='
  read (*,*) ntsteps
  print *, 'Size of time step ='
  read (*,*) tau
  t = 0d0

  !! Initialize wavepacket parameters
  print *, 'Initial wavepacket position ='
  read (*,*) x0
  print *, 'Initial wavepacket momentum ='
  read (*,*) p0
  print *, 'Initial wavepacket width ='
  read (*,*) alpha
    
  !! Allocating higher order tensors
  allocate (psi(n,2),psi0(n,2),psi_old(n,2),ham(2,n),tkin(2,n),vpot(n), &
       & work(n),chi(n,2))
  ! Remove ham, tkin, and vpot in the future
  ! Rename psi_old to hpsi in the future
  
  !! Discretizing initial wavepacket
  ! Comment in the explicit formula for psi0 in the future
  imz0sq = 0.25d0*(p0/alpha)**2
  print *, 'imz0sq =',imz0sq
  psi0_loop: do i = 1, n
     xtemp = dble(i)*h - 1d0 - x0
     exptemp = exp(-alpha*(xtemp*xtemp - imz0sq))
     print *, 'exptemp =', exptemp
     psi0(i,1) = dcmplx(exptemp*cos(p0*xtemp), exptemp*sin(p0*xtemp))
     psi0(i,2) = dcmplx(exptemp*cos(p0*xtemp), exptemp*sin(-p0*xtemp))
     ! psi0(i,1) = dcmplx(1d0,0d0) ! debug wavepacket
  end do psi0_loop
  
  !! Normalizing and saving initial wavepacket with forward propagation
  psinorm = dble(sqrt(scalar(n,h,psi0,psi0)))
  call zdscal(n,1d0/psinorm,psi0,1)
  if (idebug > 0) then
     print *, 'Initial forward wavepacket norm =',psinorm
     psinorm = dble(sqrt(scalar(n,h,psi0,psi0)))
     print *, 'Normalized forward wavepacket norm =',psinorm ! should be 1
  end if
  
  !! Normalizing and saving initial wavepacket with backward propagation
  psinorm = dble(sqrt(scalar(n,h,psi0(1,2),psi0(1,2))))
  call zdscal(n,1d0/psinorm,psi0(1,2),1)
  if (idebug > 0) then
     print *, 'Initial backward wavepacket norm =',psinorm
     psinorm = dble(sqrt(scalar(n,h,psi0(1,2),psi0(1,2))))
     print *, 'Normalized backward wavepacket norm =',psinorm ! should be 1
  end if
  psi = psi0 ! save initial wavefunction
  
  !! Testing initial position and variance
  print *, 'Initial wavepacket position =',xval(n,h,psi0),xval(n,h&
       &,psi0(1,2)) ! intial position expectation value
  print *, 'Initial wavepacket position uncertainty =',variance(n,h&
       &,psi0),variance(n,h,psi0(1,2)) ! initial position uncertainty
!!$  print *, 'Initial wavepacket momentum =',pval(n,h,0,psi0) ! initial momentum expectation value
  print *, 'Initial wavepacket momentum =',pval(n,h,0,psi0),pval(n,h,0,psi0(1,2)) ! initial momentum expectation value
  print *, 'Initial wavepacket momentum uncertainty =',pvar(n,h,psi0)&
       &,pvar(n,h,psi0(1,2)) ! initial momentum uncertainty
  print *, 'Initial wavepacket kinetic energy =',(pvar(n,h,psi0)**2&
       & + pval(n,h,0,psi0)**2)/2d0,(pvar(n,h,psi0(1,2))**2 + pval(n,h,0,psi0(1,2))**2)/2d0
  print *, 'Initial uncertainty product =',variance(n,h,psi0)*pvar(n,h&
       &,psi0),variance(n,h,psi0(1,2))*pvar(n,h,psi0(1,2)) ! must be greater than 1/2
  print *

  write (filename, "(a4,i0)") "psi_", 1
  open(69, file = filename)
  call dump_psi(n, h, 69, psi0)
  close(69)

  write(filename, "(a7,i0)") "wigner_", 1
  open(70, file = filename)
  call wigner_distr(n,h,70,psi,tau,chi)
  close(70)
  
!!$  !! Defining Hamiltonian operator
!!$  vpot = 0d0 ! local potential at grid points
!!$  hinv = 1d0/(h*h)
!!$  t_loop_diagonal: do i = 1,n
!!$     tkin(1,i) = hinv ! diagonal
!!$     tkin(2,i) = -0.5d0*hinv ! offdiagonal
!!$     ! tkin(2,i) = 0d0 ! debug
!!$  end do t_loop_diagonal
!!$  h_loop: do i = 1,n
!!$     ham(1,i) = tkin(1,i)+vpot(i)
!!$     ham(2,i) = tkin(2,i)
!!$  end do h_loop

!!$  !! Graphing eigenstates
!!$  print *, "Static Hamiltonian eigenvalues :", omega
!!$  do i=1,n
!!$     write (88,*) dble(i-1)*h - 1d0, statevec(i,n)
!!$  end do

  !! Propagating wavefunction by one time step
!!$  call propagate(n, h, psi, -tau, psi_old)
!!$  propagate backward by 1 time step for psi(t-tau) aka psi_old
  do itsteps = 1, ntsteps
!!$     call propagate(n, h, psi, tau, chi)
!!$     call propagate_ab(n, h, psi, psi_old, tau, chi)     

!!$     psinorm = dble(sqrt(scalar(n,h,psi(1,2),psi(1,2))))
!!$     print *, 'Norm before backward propagation =', psinorm
!!$     call ham_psi(n, h, psi, psi_old)
!!$     call ham_psi(n, h, psi(1,2), psi_old(1,2))
!!$     call propagate_trap(n, h, psi, tau, chi, psi_old)
!!$     call propagate_trap(n, h, psi(1,2), -tau, chi(1,2), psi_old(1,2))
     call propagate_convert(n, h, psi0, tau*itsteps, chi)
     call propagate_convert(n, h, psi0(1,2), tau*itsteps, chi(1,2))
!!$     do i = 1, n/2
!!$        work(i) = temp
!!$        work(i) = work(n-i+1) 
!!$        work(n-i+1) = temp
!!$     end do
!!$     chi = 0.5d0*(conjg(work)+chi)
!!$     psi_old = psi ! save old wavefunction
     
!!$     psinorm = dble(sqrt(scalar(n,h,chi,chi))) ! renormalization
     psi(:,1) = chi(:,1)!/psinorm ! use result as new input in next iteration
!!$     psinorm = dble(sqrt(scalar(n,h,chi(1,2),chi(1,2)))) ! renormalization
     psi(:,2) = chi(:,2)!/psinorm ! use result as new input in next iteration
     
!!$     do i = 1, n
!!$        print *, 'Probability amplitude =', abs(psi(i))**2
!!$        psi(i) = dcmplx(cos(omega(i)*tau),-sin(omega(i)*tau))*psi(i) ! propagation
!!$     end do
     ! Someone remind me the purpose of these lines bc I forgot ~Toby

     !! Testing position and variance after propagation
     print *, 'Absolute time =', t + tau*itsteps ! current time
     print *, 'Current time step =', itsteps
     psinorm = dble(sqrt(scalar(n,h,psi,psi))) ! norm
     print *, 'Norm =', psinorm
!!$     psinorm = dble(sqrt(scalar(n,h,conjg(psi(:,2)),psi)))     
!!$     print *, '<psi*(-t)|psi(t)> after propagation =', psinorm ! Commented out bc we forced renormalization
     print *, 'Wavepacket position =',xval(n,h,psi),xval(n,h,psi(1,2)) ! intial position expectation value
     print *, 'Wavepacket position uncertainty =',variance(n,h,psi)&
          &,variance(n,h,psi(1,2)) ! initial position uncertainty
!!$     print *, 'Wavepacket momentum =',pval(n,h,1,psi) ! initial momentum expectation value
     print *, 'Wavepacket momentum =',pval(n,h,0,psi),pval(n,h,0,psi(1,2)) ! initial momentum expectation value
     print *, 'Wavepacket momentum uncertainty =',pvar(n,h,psi),pvar(n&
          &,h,psi(1,2)) ! initial momentum uncertainty
     print *, 'Wavepacket kinetic energy =',(pvar(n,h,psi)**2 + pval(n,h,0&
          &,psi)**2)/2d0,(pvar(n,h,psi(1,2))**2 + pval(n,h,0,psi(1,2))**2)/2d0
     print *, 'Uncertainty product =',variance(n,h,psi)*pvar(n,h,psi)&
          &,variance(n,h,psi(1,2))*pvar(n,h,psi(1,2)) ! must be greater than 1/2
     print *
     
     write (filename, "(a4,i0)") "psi_", itsteps+1
     open(69, file = filename)
     call dump_psi(n, h, 69, psi)
     close(69)

     write(filename, "(a7,i0)") "wigner_", itsteps+1
     open(70, file = filename)
     call wigner_distr(n,h,70,psi,tau,chi)
     close(70)
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
  ! scalar = 0.5d0*(conjg(psi1(1))*psi2(1) &
  !     & + conjg(psi1(n))*psi2(n)) ! First and last point
  scalar = 0d0
  do igrid = 1, n
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

!!$  xval = 0.5d0*(-conjg(psi(1))*psi(1) &
!!$       & + conjg(psi(n))*psi(n)) ! First and last point
  xval = 0d0
  do igrid = 1, n
     xval = xval + conjg(psi(igrid))*psi(igrid)*(dble(igrid)*h - 1d0)
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
!!$  variance = conjg(psi(1))*psi(1)*(-0.5d0 - tmp)**2 &
!!$       & + conjg(psi(n))*psi(n)*(0.5d0 - tmp)**2 ! First and last point
  variance = 0d0
  do igrid = 1, n
     variance = variance + conjg(psi(igrid))*psi(igrid)*(dble(igrid)*h - 1d0 - tmp)**2
  end do
  variance = sqrt(variance*h)
end function variance

double precision function pval(n,h,idebug,psi)
  implicit none
  integer, intent(in):: n, idebug
  double precision, intent(in) :: h
  double complex, intent(in) :: psi(n)
  double complex :: tmp
  integer :: igrid

  pval = aimag(conjg(psi(1))*psi(2))/(2d0*h) &
       & - aimag(conjg(psi(n))*psi(n-1))/(2d0*h) ! First and last point
  do igrid = 2, n-1
     pval = pval + aimag(conjg(psi(igrid))*(psi(igrid+1)-psi(igrid-1))/(2d0*h))
  end do 
  pval = pval*h
  
  if (idebug == 1) then
     tmp = conjg(psi(1))*psi(2)/(2d0*h) &
          & - conjg(psi(n))*psi(n-1)/(2d0*h) ! First and last point
     do igrid = 2, n-1
        tmp = tmp + conjg(psi(igrid))*(psi(igrid+1)-psi(igrid-1))/(2d0*h)
     end do
     print *, 'Brute force momentum =',tmp*h
  end if
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
  tmp = pval(n,h,0,psi)
  pvar = (abs((psi(2))/(2d0*h) - tmp*psi(1)))**2 &
       & + (abs((-psi(n-1))/(2d0*h) - tmp*psi(n)))**2 ! First and last point
  do igrid = 2, n-1
     pvar = pvar + (abs(((psi(igrid+1)-psi(igrid-1)))/(2d0*h) - tmp*psi(igrid)))**2
  end do
  pvar = sqrt(pvar*h)
end function pvar


subroutine propagate(n, h, psi, tau, chi)
  !-----------------------------------------------------------------------------
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
end subroutine propagate_ab

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

subroutine ham_psi(n, h, psi, Hpsi)
!------------------------------------------------------------------------------
! Description:
!------------------------------------------------------------------------------
! Quickly compute H applied to psi given n and h.
!------------------------------------------------------------------------------
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
  integer, intent(in) :: n
  double precision, intent(in) :: h
  double complex, intent(in) :: psi(n)
!------------------------------------------------------------------------------
! Output Parameters
!------------------------------------------------------------------------------
  double complex, intent(out) :: Hpsi(n)
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
  do igrid = 3, n-2
     Hpsi(igrid) = -(psi(igrid-2) - 2d0*psi(igrid) + psi(igrid+2))
     ! -d2/dx^2 psi => Hpsi
  end do

  ! first and last points special case
  Hpsi(1) = -(-psi(1) + psi(3))
  Hpsi(2) = -(-2d0*psi(2) + psi(4))
  Hpsi(n-1) = -(psi(n-3) - 2d0*psi(n-1))
  Hpsi(n) = -(psi(n-2) - psi(n))

  ! multiply -d2/dx^2 by 1/(2h^2) to get Hpsi
  Hpsi = dcmplx(0.125d0/(h*h),0d0)*Hpsi ! 1/2 from p^2/2 and 1/4 from
  ! central difference formula
end subroutine ham_psi


subroutine dump_psi(n, h, i_unit, psi)
!------------------------------------------------------------------------------
! Description: 
!------------------------------------------------------------------------------
! dump psi on unit i unit.
!------------------------------------------------------------------------------
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
  integer, intent(in) :: n, i_unit
  double complex, intent(in) :: psi(n)
  double precision, intent(in) :: h
!------------------------------------------------------------------------------
! Output Parameters
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Input/Output Parameters
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!  Local Variables
!------------------------------------------------------------------------------
  integer :: igrid
  double precision :: tmp1, tmp2
!------------------------------------------------------------------------------
!  Local Constants 
!------------------------------------------------------------------------------
  ! second derivative loop
  do igrid = 1, n
     tmp1 = real(psi(igrid))
     tmp2 = aimag(psi(igrid))
     if (abs(tmp1) < 1d-10) tmp1 = 0d0
     if (abs(tmp2) < 1d-10) tmp2 = 0d0
     write(i_unit,"(f16.10,2(2x, g24.17))") (h*dble(igrid) - 1d0),&
          & tmp1, tmp2
  end do
end subroutine dump_psi

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

!!$subroutine propagate_fft(n, h, psi, tau, chi)
!!$  !------------------------------------------------------------------------------
!!$  !------------------------------------------------------------------------------
!!$  !
!!$  !------------------------------------------------------------------------------
!!$  ! Description: 
!!$  !------------------------------------------------------------------------------
!!$  !
!!$  use mkl_dfti
!!$  implicit none
!!$  !------------------------------------------------------------------------------
!!$  ! Input Parameters
!!$  !------------------------------------------------------------------------------
!!$  integer, intent(in):: n
!!$  double precision, intent(in):: h, tau
!!$  double complex, intent(in):: psi(n)
!!$  !------------------------------------------------------------------------------
!!$  ! Output Parameters
!!$  !------------------------------------------------------------------------------
!!$  double complex, intent(out):: chi(n)
!!$  !------------------------------------------------------------------------------
!!$  ! Input/Output Parameters
!!$  !------------------------------------------------------------------------------
!!$
!!$  !------------------------------------------------------------------------------
!!$  !  Local Variables
!!$  !------------------------------------------------------------------------------
!!$  integer :: i, j, istatus
!!$  type(dfti_descriptor), pointer :: My_Desc1_Handle, My_Desc2_Handle
!!$  double precision :: tmp, tmp_expt
!!$  double precision, parameter :: pi = 4d0*atan(1d0)
!!$  double complex, allocatable :: tmpfft(:)
!!$  !------------------------------------------------------------------------------
!!$  !  Local Constants 
!!$  !------------------------------------------------------------------------------
!!$
!!$  allocate(tmpfft(n))
!!$  
!!$  istatus = DftiCreateDescriptor(My_Desc1_Handle, DFTI_DOUBLE, DFTI_COMPLEX, 1, n)
!!$  istatus = DftiSetValue(My_Desc1_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
!!$  istatus = DftiCommitDescriptor(My_Desc1_Handle)
!!$  istatus = DftiComputeForward(My_Desc1_Handle, psi, tmpfft) ! forward FFT
!!$
!!$  tmp_expt = 0d0
!!$!  
!!$  do i = 1,n
!!$     tmp = sin(2d0*pi*dble(i-1)/dble(n))/h ! momentum eigenvalue
!!$     tmp_expt = tmp_expt + abs(tmpfft(i))**2*tmp ! momentum expectation
!!$     ! value in energy basis
!!$     tmpfft(i) = tmpfft(i)/dble(n)*exp(cmplx(0d0,-tmp*tmp*tau*0.5d0)) ! propagated in energy basis
!!$  end do
!!$
!!$  
!!$  
!!$  istatus = DftiComputeBackward(My_Desc1_Handle, tmpfft, chi) ! backward FFT
!!$  istatus = DftiFreeDescriptor(My_Desc1_Handle)
!!$! result is given by {X_out(1),X_out(2),...,X_out(32)}
!!$  print *, 'momentum expectation value in energy basis =', h*tmp_expt/dble(n)
!!$  deallocate(tmpfft)
!!$end subroutine propagate_fft

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
  integer :: i, j, k, istatus
  double complex :: tmp
  double precision, parameter :: pi = 4d0*atan(1d0)
  double complex, allocatable :: w(:)
  character :: fmt*40
  double precision :: xi, eta
  !------------------------------------------------------------------------------
  !  Local Constants 
  !------------------------------------------------------------------------------

  allocate(w(n))

  write(fmt,*) "(f16.10,",2*n,"(2x, g24.17))"
  print *,fmt
  do i = 1,n ! x loop
     do j = 1,n ! k loop
        tmp = dcmplx(0d0,0d0)
        eta = -1d0 + dble(2*j)/dble(n+1)

        do k = i + 1, n - i !xi loop
           xi = -1d0 + dble(2*k)/dble(n+1)
           tmp = tmp + conjg(psi(i+k))*psi(i-k)*exp(dcmplx(0d0,2d0*pi&
                &*xi*eta))
        end do

        do k = n - i, i + 1 !xi loop
           xi = -1d0 + dble(2*k)/dble(n+1)
           tmp = tmp + conjg(psi(i+k))*psi(i-k)*exp(dcmplx(0d0,2d0*pi&
                &*xi*eta))
        end do
           
        tmp = tmp*dble(h*2d0/pi)
        w(j) = tmp
     end do
     write(70,fmt) h*dble(i) - 1d0&
          &,(real(w(j)), aimag(w(j)), j = 1,n)
  end do
  
  deallocate(w)
end subroutine wigner_distr
