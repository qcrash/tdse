program tdse

  ! Local variables
  integer :: n ! vector space dimension
  integer :: i, j ! loop variables
  integer :: info ! LAPACK return status
  double precision :: tau, t ! time step, time variable
  ! Local floating point intermediates
  double precision :: h, hinv, sigmainv, psinorm, imz0sq, x, xtemp, exptemp
  ! Local parameters
  double precision, parameter :: x0 = -0.8d0, p0 = 0.1d0, alpha = 1000d0, sigma = 0.1d0, pi = 4d0*atan(1d0)
  ! External functions
  double precision, external :: dznrm2
  ! Allocatable arrays
  double complex, allocatable :: psi(:), psi0(:) 
  double precision, allocatable :: work(:), omega(:)
  double precision, allocatable :: ham(:,:), tkin(:,:), vpot(:,:), u(:,:)

  ! Assign values to scalars
  read (*,*) n ! prompts user-defined input for number of grid points
  h = 2d0/dble(n-1)
  !!!!! Before, h = 2d0/dble(n+1), but shouldn't h = 2d0/dble(n-1) since there are n-1 number of trapezoids? I changed h here, but you guys can revert it if h was correct before. ~Toby 10/25/22

  ! Allocation of higher order tensors
  allocate (psi(n),psi0(n),ham(n,n),tkin(n,n),vpot(n,n),u(n,n),work(3*n),omega(n))
  
  ! Discretizing the initial wavepacket
  ! sigmainv = 1d0 / (sigma*sigma)
  !!!!! sigma and sigmainv are no longer in use, so I commented it out; feel free to delete the line entirely. ~Toby 10/25/22
  imz0sq = 0.25d0*(p0/alpha)**2
  psi0(1) = 0d0
  psi0(n) = 0d0 ! setting the first and last grid point of the initial wavepacket to be zero in order to satisfy the boundary conditions
  psi0_loop: do i = 2, n-1 ! skipped i=1 and i=n since they're zero due to the boundary condition
     !!!!! psi0(i) = cmplx(exp(-(-1d0 + dble(i)*h - x0)**2 * signmainv),0d0) I kept our old wavepacket if you guys still wanted to reference them ~Toby 10/25/22
     x = dble(i-1)*h - 1d0
     xtemp = x+x0
     exptemp = exp(-alpha*(xtemp*xtemp - imz0sq))
     psi0(i) = cmplx(exptemp*cos(p0*xtemp), exptemp*sin(p0*xtemp))
     ! psi0(i) = cmplx(1d0,0d0) ! debug wavepacket
  end do psi0_loop
  
  ! Normalizing the initial wavepacket
  psinorm = sqrt(h)*dznrm2(n,psi0,1)
  print *, 'Actual unscaled norm =',psinorm
  call zdscal(n,1d0/psinorm,psi0,1)
  psinorm = sqrt(h)*dznrm2(n,psi0,1)
  print *, 'Actual rescaled norm =',psinorm
  print *, '2-h =',2d0-h !!!!! Can someone remind me why we care about 2-h again? T-T ~Toby 10/25/22
  
  ! Loop to set vpot to 0
  v_loop: do i = 1,n
     do j = 1,n
        vpot(j,i) = 0d0 ! Read the elements of this matrix with j first since in this loop j is defined, then go to outer loop
     end do
  end do v_loop

  ! Loop to define tkin matrix
  hinv = 1d0 / (h*h)
  tkin = 0d0
  t_loop_diagonal: do i = 1,n-1
     tkin(i,i+1) = -0.5d0 * hinv
     tkin(i,i) = 1d0 * hinv
     tkin(i+1,i) = -0.5d0 * hinv
  end do t_loop_diagonal
  tkin(n,n) = 1d0 / h**2

  ! Loop to define ham matrix
  h_loop: do i = 1,n
     do j = 1,n
        ham(j,i) = tkin(j,i) + vpot(j,i)
     end do
     ! This is a matrix addition of the tkin and vpot matrices
  end do h_loop

  call ssyev('v','u',n,ham,n,omega,work,3*n,info)
  if (info /=  0) stop 'error in ssyev'
  
  ! Deallocation of higher order tensors
  deallocate (psi,psi0,ham,tkin,vpot,u,work,omega)
end program tdse

double complex function scalar(n,psi1,psi2)
  ! Calculate scalar product <psi1|psi2> on the grid
  ! using trapezoidal rule and
  ! zero boundary conditions

  
  ! Input variables
  integer, intent(in):: n
  double complex, intent(in) :: psi1(n), psi2(n)

  ! Local variables
  
  integer :: igrid


  
  
end function scalar
