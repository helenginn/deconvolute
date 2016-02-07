! This subroutine performs the Fourier transforms of a 2-D planes with
! a perpendicular 2-fold through the origin.
! The input map is real, and dimension Nfast x Nslow/2 + 1 x Ng
! 2-fold symmetry means there is redundancy on the first and last Nslow cols -
!   the second halves are filles in from the first.
! First FT real -> cconjg on Fast ( Nfast/2+1 x Nslow/2+1) x Ng
! then FT cconjg -> real on slow  ( Nfast/2+1 x Nslow) x Ng
! However, this leave the storage out of normal order (because of the
! different directions), as well as silly dimensioning problems.
! So the input array is ( Nfast/2+1 x Nslow+2) x Ng, and with
! the data as x00 x20 ...., x10, x30..., etc.

subroutine ffl_p2p(mapout, nf, ns)

! we must use explicitly-dimensioned arrays, as assumed size
! arrays cannot have specific elements sent as actual arguments
! to dummies (in realft, etc) which are themselves arrays.

integer, intent(in)	:: nf, ns
real, intent(inout)	:: mapout(0:nf/2,0:ns+1)

integer d(5)	! Ten Eyck dimensioning array.
integer nf2, ns2

if ( modulo(nf,2) .ne. 0 .or. modulo(ns,2) .ne. 0 ) &
write(6,*) ' Dimensions must be even - actually ', nf, ns

nf2 = nf/2
ns2 = ns/2

!!! check the symmetry of first and last cols.
!!if( any(mapout(1:nf2-1, 0) .ne. mapout(nf2-1:1:-1, 0)))	&
!!  Write(6,*) ' Symmetry error on col 0 - even'
!!if( any(mapout(0:nf2-1, 1) .ne. mapout(nf2-1:0:-1, 1)))	&
!!  Write(6,*) ' Symmetry error on col 0 - odd'
!!if( any(mapout(1:nf2-1, ns) .ne. mapout(nf2-1:1:-1, ns)))	&
!!  Write(6,*) ' Symmetry error on last col - even'
!!if( any(mapout(0:nf2-1, ns+1) .ne. mapout(nf2-1:0:-1, ns+1)))	&
!!  Write(6,*) ' Symmetry error on last col - odd'

! Perform first FT
d(1) = (nf2+1) * (ns+2)
d(2) = 1	
d(3) = nf+2
d(4) = 1
d(5) = 1
call realft(mapout, mapout(0,1), nf2, d)

!write(6,'(/(3g14.6))') mapout
! Put the Ns/2 real parts in the 0 Imag parts
mapout(0:nf2,1) = mapout(0:nf2,ns)
!write(6,'(/(3g14.6))') mapout

! Perform second FT
d(2) = nf + 2
d(3) = (nf + 2) * (ns2 + 1)
d(4) = nf2 + 1
call hermft(mapout, mapout(0,1), ns2, d)
!write(6,'(/(3g14.6))') mapout

! and that's it.
return
end

subroutine ffl_g1f(map,nh,nk)
! forward real -> cc group trans on second index of 2-d array 
! with positive in the exponent.
! nh,nk are the real space dimensions	
! the input map should be nh x 0:nk+1.
integer, intent(in)	:: nh, nk
real, intent(inout)	:: map(nh, 0:nk+1)

integer index(5)

if(mod(nk,2) .ne. 0) &
          stop 'FFL _g1f -- transform dimension must be even'

index(1) = nh*(nk+2)
index(2) = 2*nh
index(3) = index(1)
index(4) = nh
index(5) = 1
call realft(map(1,0), map(1,1), nk/2, index)
! reverse sign on Imaginaries
map(:,1::2) = -map(:,1::2)
return
end

SUBROUTINE FFL_g1r(map,nh,nk)
! reverse trans with negative in the exponent.
! nh,nk are the real space dimensions	
integer, intent(in)	:: nh, nk
real, intent(inout)	:: map(nh, 0:nk+1)

integer index(5)

if(mod(nk,2) .ne. 0) &
          stop 'FFL_g1r -- transform dimension must be even'

!store real(f(nk/2)) in im(f(0))
map(:,1) = map(:,nk)

index(1) = nh*(nk+2)
index(2) = 2*nh
index(3) = index(1)
index(4) = nh
index(5) = 1

call hermft(map(1,0), map(1,1), nk/2, index)
return
end
