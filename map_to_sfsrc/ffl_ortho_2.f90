! RKB 13-JUN-1996 11:29:50
! This routines performs a forward transform on orthorhombic cells
! with a z-2fold.
! So far, works for
! P2_1 2_1 2 (spacegroup 18)	Thk-z = (-)^(h+k)Th-kz
! P222		16		Thk-z = Th-kz
! Given all x, half y, half z, returns h,k,l all .ge. 0.
! The real and imaginary parts are on ajacent l-planes.
! some care is needed on the special planes and lines, which contain
! duplicates within the supplied volume.

SUBROUTINE ffl_ortho_2f(rho)
  USE mapdefs
  TYPE (map) rho

!  map	:: in(0:n1-1, 0:n2/2, 0:n3/2)
!  sfs	:: final(0:n1/2, 0:n2/2, 0:n3+1)

INTEGER k,j, n1, n2, n3
REAL temp1(0:rho%mxyz(1)/2, 0:rho%mxyz(2)+1)
REAL temp3(0:rho%mxyz(1)/2, 0:rho%mxyz(3)+1)		! work arrays.
REAL, POINTER	:: final(:,:,:)

final => rho%sfs	! save retyping.
!!WRITE(6,*) ' ffl_g18 -- Bound of final',lbound(final),ubound(final)

n1 = rho%mxyz(1)
n2 = rho%mxyz(2)
n3 = rho%mxyz(3)

! 2-D transform each section using symmetry, and reconstruct all n3
! I222 is similar, but only to z/4; the sign change is similar, but
! reflected in z/2 instead of z.
DO k = 0,n3/2

! Copy input section to temporary array, with even/odd ordering.
   DO j = 0,n2/2
     temp1(0:n1/2-1, 2*j)   = rho%rho(0:n1-1:2, j, k)
     temp1(0:n1/2-1, 2*j+1) = rho%rho(1:n1-1:2, j, k)
   END DO

!! and apply symmetry on first and last cols (note that nf is even)
!!mapout(nf2-1:(nf2+2)/2:-1, (/ 0, ns /))   = mapout(1:(nf2-1)/2, (/ 0, ns /))
!!mapout(nf2-1:(nf2+1)/2:-1, (/ 1, ns+1 /)) = mapout(0:(nf2-2)/2, (/ 1, ns+1 /))
! symmetery is checked in ffl_p2p

   CALL ffl_p2p(temp1, n1, n2)
   final(:,:,k) = temp1(:,0:n2/2)
   IF( MOD(k, n3/2) .eq. 0 ) CYCLE
   final(:, 1:n2/2-1, n3-k) = temp1(:, n2-1:n2/2+1:-1)
   final(:, (/ 0, n2/2 /), n3-k) = temp1(:, (/ 0, n2/2 /))
! For P222 and I222 omit the sign change
! For P21212 change signs for -z if n1+n2 odd
   IF ( rho%lspgrp .eq. 18 ) THEN
     final(0:n1/2:2, 1:n2/2:2, n3-k) = - final(0:n1/2:2, 1:n2/2:2, n3-k) 
     final(1:n1/2:2, 0:n2/2:2, n3-k) = - final(1:n1/2:2, 0:n2/2:2, n3-k) 
   END IF
END DO

!! Now the n3 transform 
! I222 uses centering here. P21212 and P222 are the same
DO j = 0,n2/2
  temp3(:,0:n3-1) = final(:,j,0:n3-1)
  CALL ffl_g1f(temp3, n1/2+1, n3)
  final(:,j,:) = temp3
END DO

RETURN
END

!!! .... and the inverse transform
!!! one could improve -- don't use transpose, and final IN only
!!subroutine ffl_g18b(final, n1, n2, n3, out)
!!  integer, intent(in)	:: n1, n2, n3
!!  real, intent(inout)	:: final(0:n1/2,0:n2/2,0:n3+1)
!!  real, intent(out)	:: out(0:n2/2,0:n1+1,0:n3/2)
!!
!!integer k,j
!!real temp3(0:n1/2,0:n3+1), in2(0:n2-1,0:n1/2)		! work arrays.
!!
!!do j = 0,n2/2
!!  temp3 = final(:,j,:)
!!  call ffl_g1r(temp3, n1/2+1, n3)
!!  final(:,j,0:n3-1) = temp3(:,0:n3-1)
!!end do
!!do k = 0,n3/2
!!   in2(0:n2/2, 0:n1/2) = transpose(final(:,:,k))
!!   in2(n2-1:n2/2:-1, 0:n1/2) = transpose(final( :, 1:, modulo(n3-k,n3)))
!!   in2(n2-1:n2/2:-2, 0:n1/2:2) = - in2(n2-1:n2/2:-2, 0:n1/2:2) 
!!   in2(n2-2:n2/2:-2, 1:n1/2:2) = - in2(n2-2:n2/2:-2, 1:n1/2:2) 
!!   call ffl_p2p(in2, n2, n1, out(0,0,k))
!!end do
!!return
!!end
