program tdse
  implicit none
  !! Declaring variables
  integer, parameter :: idebug = 1
  integer :: info, i, j, n ! LAPACK status, loop vars, vector space dim
  double precision :: h, psinorm, tau, t ! grid spacing, time step, time param
  double precision :: hinv, imz0sq, xtemp, exptemp ! temp vars
  double precision, parameter :: x0 = -0.8d0, p0 = 0.1d0, alpha = 1000d0, &
       & pi = 4d0*atan(1d0)
  double complex :: temp
  double complex, allocatable :: psi(:), psi0(:), transmat(:,:), &
       & psitilda(:) ! psitilda is wavefunction in real orthonormal basis
  double precision, allocatable :: ham(:,:), tkin(:,:), &
       & vpot(:), statevec(:,:), work(:), omega(:)
  double precision, external :: dznrm2
  double complex, external :: scalar
  double precision, external :: xval, variance, pval
  
  !! User-defined parameters
  print *, 'Number of grid points ='
  read (*,*) n
  ! Include user prompts for x0, p0, alpha
  h = 2d0/dble(n-1)

  !! Initializing time variables
  t = 0d0
  tau = 1d0

  !! Allocating higher order tensors
  allocate (psi(n),psi0(n),ham(2,n),tkin(2,n),vpot(n), &
       & work(3*n),omega(n),statevec(n,n),transmat(n,n),psitilda(n))
  
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
  
  !! Transforming to real orthonormal basis
  temp = dcmplx(sqrt(h),0d0)
  do i = 1,n
     psitilda(i) = psi0(i)*temp ! rescale to make basis orthonormal
  end do
  psitilda(1) = psitilda(1)/dcmplx(sqrt(2d0),0d0)
  psitilda(n) = psitilda(n)/dcmplx(sqrt(2d0),0d0)
  if (idebug > 0) then
     print *, 'Norm in real orthonormal basis =', dznrm2(n,psitilda,1)
  end if
  
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
  call dscal(2*n,1d0/h,ham,1)
  ham(1,1) = 2d0*ham(1,1)
  ham(2,1) = sqrt(2d0)*ham(2,1)
  ham(1,n) = 2d0*ham(1,n)
  ham(2,n) = sqrt(2d0)*ham(2,n)
  
  !! Diagonalizing Hamiltonian operator
  omega = 0d0
  statevec = 0d0
  call dsbev('v','l',n,1,ham,2,omega,statevec,n,work,info)
  if (info /=  0) then
     print *, "Info = ", info
     stop 'error in dsbev'
  end if
  do i = 1, n
     do j = 1, n
        transmat(i,j) = dcmplx(statevec(i,j),0d0) ! make statevec dcmplx
     end do
  end do

  !! Graphing eigenstates
  ! print *, "Static Hamiltonian eigenvalues :", omega
  ! do i=1,n
  !    write (88,*) dble(i-1)*h - 1d0, statevec(i,n)
  ! end do

  !! Transforming wavefunction to energy eigenbasis
  do i = 1, n
     psi(i) = dcmplx(0d0,0d0)
     do j = 1, n
        psi(i) = psi(i)+transmat(j,i)*psitilda(j) ! transformation
     end do
  end do
  ! call zgemv('t',n,n,dcmplx(1d0,0d0),transmat,n,psi0,1,dcmplx(0d0,0d0),psi,1)
  if (idebug > 0) then
     print *, 'Norm in energy basis =', dznrm2(n,psi,1)
  end if
     
  !! Propagating momentum wavefunction by one time step
  do i = 1, n
     psi(i) = exp(-dcmplx(0d0,1d0)*omega(i)*tau)*psi(i) ! propagation
  end do
  if (idebug > 0) then
     print *, 'Norm after time step in energy basis =', dznrm2(n,psi,1)
  end if
  
  !! Transforming wavefunction back to real orthonormal basis
  do i = 1, n
     psitilda(i) = dcmplx(0d0,0d0)
     do j = 1, n
        psitilda(i) = psitilda(i)+transmat(j,i)*psi(j)
     end do
  end do
  ! call zgemv('n',n,n,dcmplx(1d0,0d0),transmat,n,psi,1,dcmplx(0d0,0d0),psi0,1)
  if (idebug > 0) then
     print *, 'Norm after time step in real orthonormal basis =', dznrm2(n,psitilda,1)
  end if
  
  !! Deallocation of higher order tensors
  deallocate (psi,psi0,ham,tkin,vpot,work,omega,statevec,transmat,psitilda)
end program tdse

double complex function scalar(n,h,psi1,psi2)
  ! Calculate scalar product <psi1|psi2> on the grid
  ! using trapezoidal rule and the boundary conditions
  ! <psi1|psi2> = int(psi1*(x)psi2(x)dx)
  
  !! Defining parameters and local variables
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
  integer, intent(in):: n
  double precision, intent(in) :: h
  double complex, intent(in) :: psi(n)
  double precision :: tmp
  double precision, external :: xval
  
  ! var^2 = exp[(x-exp(x))^2]
  tmp = xval(n,h,psi)
  variance = conjg(psi(1))*psi(1)*(-0.5d0 - tmp)**2 &
       & + conjg(psi(n))*psi(n)*(0.5d0 - tmp)**2 ! First and last point
  do igrid = 2, n-1
     variance = variance + conjg(psi(igrid))*psi(igrid)*(dble(igrid-1)*h - 1d0 - tmp)**2
  end do
  variance = sqrt(variance*h)
end function variance

double precision function pval(n,h,psi)
  integer, intent(in):: n
  double precision, intent(in) :: h
  double complex, intent(in) :: psi(n)
  integer :: igrid

  pval = 0.5d0*(aimag(conjg(psi(1))*psi(2))/(2d0*h) &
       & - aimag(conjg(psi(n))*psi(n-1))/(2d0*h)) ! First and last point
  do igrid = 2, n-1
     pval = pval + aimag(conjg(psi(igrid))*(psi(igrid+1)-psi(igrid-1)))/(2d0*h)
  end do 
  pval = pval*h
end function pval
