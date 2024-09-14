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
