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
  integer, parameter :: idebug = 0
  ! External functions
  double precision, external :: dznrm2
  double complex, external :: scalar
  ! Allocatable arrays
  double complex, allocatable :: psi(:), psi0(:) 
  double precision, allocatable :: work(:), omega(:)
  double precision, allocatable :: ham(:,:), tkin(:,:), vpot(:), u(:,:)

  ! Assign values to scalars
  read (*,*) n ! prompts user-defined input for number of grid points
  h = 2d0/dble(n-1)
  !!!!! Before, h = 2d0/dble(n+1), but shouldn't h = 2d0/dble(n-1) since there are n-1 number of trapezoids? I changed h here, but you guys can revert it if h was correct before. ~Toby 10/25/22

  ! Allocation of higher order tensors
  allocate (psi(n),psi0(n),ham(n,2),tkin(n,2),vpot(n),u(n,n),work(3*n),omega(n))
  
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
  psinorm=dble(sqrt(scalar(n,h,psi0,psi0)))
  call zdscal(n,1d0/psinorm,psi0,1)
  if (idebug > 0) then
     psinorm = dble(sqrt(scalar(n,h,psi0,psi0)))
     print *, 'Actual rescaled norm =',psinorm
     print *, '2-h =',2d0-h !!!!! Can someone remind me why we care about 2-h again? T-T ~Toby 10/25/22
  end if

  ! vpot is vector containing values of local potential
  ! at grid points 
  ! Loop to set vpot to 0
  vpot = 0d0 ! array initialization

  ! Loop to define tkin matrix
  hinv = 1d0 / (h*h)
  t_loop_diagonal: do i = 1,n-1
     tkin(i,2) = -0.5d0*hinv
     tkin(i,1) = hinv
  end do t_loop_diagonal
  tkin(n,1) = 1d0 / h**2
  tkin(n,2) = 0d0
  
  ! Loop to define ham matrix
  h_loop: do i = 1,n
        ham(i,1) = tkin(i,1)+vpot(i)
        ham(i,2) = tkin(i,2)
     end do
  end do h_loop

  call ssyev('v','u',n,ham,n,omega,work,3*n,info)
  if (info /=  0) stop 'error in ssyev'
  
  ! Deallocation of higher order tensors
  deallocate (psi,psi0,ham,tkin,vpot,u,work,omega)
end program tdse

double complex function scalar(n,h,psi1,psi2)
  ! Calculate scalar product <psi1|psi2> on the grid
  ! using trapezoidal rule and
  ! the boundary conditions
  
  ! Input variables
  integer, intent(in):: n
  double precision, intent(in) :: h
  double complex, intent(in) :: psi1(n), psi2(n)

  ! Local variables
  
  integer :: igrid

  scalar = (0d0,0d0)
  do igrid = 2, n-1
     scalar = scalar + conjg(psi1(igrid))*psi2(igrid)
  end do
  ! First and last point
  scalar = scalar+0.5d0*(conjg(psi1(1))*psi2(1)&
       & +conjg(psi1(n))*psi2(n)) 
  scalar = scalar*h
  
end function scalar
