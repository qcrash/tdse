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
