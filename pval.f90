double precision function pval(n,h,idebug,psi)
  implicit none
  integer, intent(in):: n, idebug
  double precision, intent(in) :: h
  double complex, intent(in) :: psi(n)
  double complex :: tmp
  integer :: igrid

  pval = aimag(conjg(psi(1))*psi(2))/(2d0*h) &
       & - aimag(conjg(psi(n))*psi(n-1))/(2d0*h) ! First and last point
  do igrid = 2, n-1
     pval = pval + aimag(conjg(psi(igrid))*(psi(igrid+1)-psi(igrid-1))/(2d0*h))
  end do 
  pval = pval*h
  
  if (idebug == 1) then
     tmp = conjg(psi(1))*psi(2)/(2d0*h) &
          & - conjg(psi(n))*psi(n-1)/(2d0*h) ! First and last point
     do igrid = 2, n-1
        tmp = tmp + conjg(psi(igrid))*(psi(igrid+1)-psi(igrid-1))/(2d0*h)
     end do
     print *, 'Brute force momentum =',tmp*h
  end if
end function pval
