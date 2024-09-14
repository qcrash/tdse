subroutine propagate_fft(n, h, psi, tau, chi)
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------------
  ! Description: 
  !------------------------------------------------------------------------------
  !
  use mkl_dfti
  implicit none
  !------------------------------------------------------------------------------
  ! Input Parameters
  !------------------------------------------------------------------------------
  integer, intent(in):: n
  double precision, intent(in):: h, tau
  double complex, intent(in):: psi(n)
  !------------------------------------------------------------------------------
  ! Output Parameters
  !------------------------------------------------------------------------------
  double complex, intent(out):: chi(n)
  !------------------------------------------------------------------------------
  ! Input/Output Parameters
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !  Local Variables
  !------------------------------------------------------------------------------
  integer :: i, j, istatus
  type(dfti_descriptor), pointer :: My_Desc1_Handle, My_Desc2_Handle
  double precision :: tmp, tmp_expt
  double precision, parameter :: pi = 4d0*atan(1d0)
  double complex, allocatable :: tmpfft(:)
  !------------------------------------------------------------------------------
  !  Local Constants 
  !------------------------------------------------------------------------------

  allocate(tmpfft(n))
  
  istatus = DftiCreateDescriptor(My_Desc1_Handle, DFTI_DOUBLE, DFTI_COMPLEX, 1, n)
  istatus = DftiSetValue(My_Desc1_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
  istatus = DftiCommitDescriptor(My_Desc1_Handle)
  istatus = DftiComputeForward(My_Desc1_Handle, psi, tmpfft) ! forward FFT

  tmp_expt = 0d0
!  
  do i = 1,n
     tmp = sin(2d0*pi*dble(i-1)/dble(n))/h ! momentum eigenvalue
     tmp_expt = tmp_expt + abs(tmpfft(i))**2*tmp ! momentum expectation
     ! value in energy basis
     tmpfft(i) = tmpfft(i)/dble(n)*exp(cmplx(0d0,-tmp*tmp*tau*0.5d0)) ! propagated in energy basis
  end do

  
  
  istatus = DftiComputeBackward(My_Desc1_Handle, tmpfft, chi) ! backward FFT
  istatus = DftiFreeDescriptor(My_Desc1_Handle)
! result is given by {X_out(1),X_out(2),...,X_out(32)}
  print *, 'momentum expectation value in energy basis =', h*tmp_expt/dble(n)
  deallocate(tmpfft)
end subroutine propagate_fft
