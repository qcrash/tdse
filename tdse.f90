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

  write(filename, "(a3,i0)") "sb_", 1
  open(70, file = filename)
  call segal_bargmann(n,h,70,psi,tau,chi)
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

     write(filename, "(a3,i0)") "sb_", itsteps+1
     open(70, file = filename)
     call segal_bargmann(n,h,70,psi,tau,chi)
     close(70)
  end do
  
  !! Deallocation of higher order tensors
  deallocate (psi,psi0,psi_old,ham,tkin,vpot,work,chi)
end program tdse
