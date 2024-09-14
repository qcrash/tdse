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
