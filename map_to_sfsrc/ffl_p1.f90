!
! This routine performs a P1 transform, given all (x, y, z), and returns
! (h .ge.0, all k, all l). The real and imaginary parts are adjacent in h.
!
subroutine ffl_p1f(rho)
  use mapdefs
  type (map) rho

  integer k,j, n1, n2, n3
  real temp1(0:rho%mxyz(1)+1, 0:rho%mxyz(2)-1)
  real temp3(0:rho%mxyz(1)+1, 0:rho%mxyz(3)-1)		! work arrays.
  real, pointer	:: final(:,:,:)

  print *,"Now in subroutine ffl_p1f (ffl_p1.f90)"
!
! Save some typing
!
  final => rho%sfs
  n1 = rho%mxyz(1)
  n2 = rho%mxyz(2)
  n3 = rho%mxyz(3)
!
! 2-D transform each section, and reconstruct all n3.
!
  do k = 0,n3-1
!
! Copy input section to temporary array
!
    temp1(0:n1-1, : ) = rho%rho(:, :, k)
!
!!!call fftreal_recip(test, testout, n1*n2*n3, n1,n2,n3)
!
    call ffl_p1p(temp1, n1, n2)
    final(:, :, k) = temp1
  enddo

  print *,"Done 2D transform on each section"

!
! Now the n3 transform 
!
  do j = 0,n2-1
    temp3 = final(:, j, :)
    call ffl_gcf(temp3, n1+2, n3)
!
! take the complex conjugate, as all the processing is done exp(-),
! whereas + is needed.
!
    temp3( 1::2, : ) = - temp3(1::2, :)
    final(:,j,:) = temp3
  enddo

  print *,"Done n3 transform"
  print *,"Returning from subroutine ffl_p1f (ffl_p1.f90)"

  return
end
