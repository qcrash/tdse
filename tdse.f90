program tdse
  implicit none  
  !! Declaring variables
  integer :: info, i, j, n, ntsteps, itsteps, idebug, iunit ! LAPACK status, loop vars, vector space dim
  double precision :: h, psinorm, tau, t0 ! grid spacing, time step, time param
  double precision :: hinv, imz0sq, xtemp, exptemp, x0, p0, alpha ! temp vars
  double precision, parameter :: pi = 4d0*atan(1d0)
  double complex :: temp
  double complex, allocatable :: psi(:,:), psi0(:,:), chi(:,:), &
       & psi_old(:,:), work(:) ! psi_old is wavefunction in real orthonormal basis
  double precision, external :: dznrm2
  double complex, external :: scalar
  double precision, external :: xval, variance, pval, pvar
  character(len=15) :: filename
  character(len=20) :: mode, repr
  character(len=80) :: output

  call getcli(n, ntsteps, idebug, tau, x0, p0, t0, alpha, mode, output, repr)
  
  h = 2d0/dble(n+1)
    
  !! Allocating higher order tensors
  allocate (psi(n,2),psi0(n,2),psi_old(n,2),work(n),chi(n,2))
  
  !! Discretizing initial wavepacket
  ! Comment in the explicit formula for psi0 in the future
  imz0sq = 0.25d0*(p0/alpha)**2
  psi0_loop: do i = 1, n
     xtemp = dble(i)*h - 1d0 - x0
     exptemp = exp(-alpha*(xtemp*xtemp - imz0sq))
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
  
  call dump_repr(n,h,psi0,t0,tau,chi,output,repr)

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
     if (idebug > 0) then
        print *, 'Absolute time =', t0 + tau*itsteps ! current time
        print *, 'Current time step =', itsteps
        psinorm = dble(sqrt(scalar(n,h,psi,psi))) ! norm
        print *, 'Norm =', psinorm
     end if
!!$     psinorm = dble(sqrt(scalar(n,h,conjg(psi(:,2)),psi)))     
!!$     print *, '<psi*(-t)|psi(t)> after propagation =', psinorm ! Commented out bc we forced renormalization
          
     call dump_repr(n,h,psi,t0 + tau*itsteps,tau,chi,output,repr)
  end do
  
  !! Deallocation of higher order tensors
  deallocate (psi,psi0,psi_old,work,chi)
end program tdse
