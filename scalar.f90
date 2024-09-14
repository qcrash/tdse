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
