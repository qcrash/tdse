program tdse

  ! Local variables
  integer :: n ! vector space dimension
  integer :: i, j ! loop variables
  integer :: info ! LAPACK return status
  double precision :: tau, t ! time step, time variables
  ! Local floating point intermediates
  double precision :: h, hinv, sigmainv, psinorm
  ! Local parameters
  double precision, parameter :: x0 = -0.8d0, sigma = 0.1d0
  ! External functions
  double precision, external :: dznrm2
  ! Allocatable arrays
  double complex, allocatable :: psi(:), psi0(:) 
  double precision, allocatable :: work(:), omega(:)
  double precision, allocatable :: ham(:,:), tkin(:,:), vpot(:,:), u(:,:)

  !Variable initialization
  !Assign values to scalars
  read (*,*) n
  h = 2d0 / dble (n+1)

  !Allocation of variables
  allocate (psi(n),psi0(n),ham(n,n),tkin(n,n),vpot(n,n),u(n,n),work(3*n),omega(n))
  
  !Initialization
  sigmainv = 1d0 / (sigma*sigma)
  psi0_loop: do i = 1, n
     !     psi0(i) = cmplx(exp(-(-1d0 + dble(i)*h - x0)**2 *
     !     signmainv),0d0)
     psi0(i) = cmplx(1d0,0d0)  ! debug
  end do psi0_loop
  ! Normalize
  psinorm = sqrt(h)*dznrm2(n,psi0,1)
  print *, 'Psinorm^2 = ', psinorm*psinorm
  print *, '2-h     = ', 2d0-h
  call zdscal(n,cmplx(1d0/psinorm,0d0),psi0,1)
  
  
  !Loop to set vpot to 0
  v_loop: do i = 1,n
     do j = 1,n
        vpot(j,i) = 0d0 !Read the elements of this matrix with j first since in this loop j is defined, then go to outer loop
     end do
  end do v_loop

  !Loop to define tkin matrix
  hinv = 1d0 / (h*h)
  tkin = 0d0
  t_loop_diagonal: do i = 1,n-1
     tkin(i,i+1) = -0.5d0 * hinv
     tkin(i,i) = 1d0 * hinv
     tkin(i+1,i) = -0.5d0 * hinv
  end do t_loop_diagonal
  tkin(n,n) = 1d0 / h**2

  !Loop to define ham matrix
  h_loop: do i = 1,n
     ham(j,i) = tkin(j,i) + vpot(j,i)
     !This is a matrix addition of the tkin and vpot matrices
  end do h_loop

  call ssyev('v','u',n,ham,n,omega,work,3*n,info)
  if (info /=  0) stop 'error in ssyev'
  
  !Deallocation of variables
  deallocate (psi,psi0,ham,tkin,vpot,u,work,omega)
end program tdse
  
