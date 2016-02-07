! This subroutine performs the Fourier transforms of a 2-D planes with
! a perpendicular 2-fold through the origin.
! The input map is real, and dimension Nfast x Nslow/2 + 1 x Ng
! 2-fold symmetry means there is redundancy on the first and last Nslow cols -
!   the second halves are filles in from the first.
! First FT real -> cconjg on Fast ( Nfast/2+1 x Nslow/2+1) x Ng
! then FT cconjg -> real on slow  ( Nfast/2+1 x Nslow) x Ng
! However, this leave the storage out of normal order (because of the
! different directions), as well as silly dimensioning problems.
! So the output array is ( Nfast/2+1 x Nslow+2) x Ng, and initially
! the data is copied as x00 x20 ...., x10, x30..., etc.

subroutine ffl_p2p(mapin, nf, ns, mapout)

! we must use explicitly-dimensioned arrays, as assumed size
! arrays cannot have specific elements sent as actual arguments
! to dummies (in realft, etc) which are themselves arrays.

integer, intent(in)	:: nf, ns
real, intent(in)	:: mapin(0:nf-1,0:ns/2)
real, intent(out)	:: mapout(0:nf/2,0:ns+1)

real temp(0:nf/2,0:ns+1)

integer d(5)	! Ten Eyck dimensioning array.
integer j, nf2, ns2

print *,'Made it to ffl_p2p'

if ( modulo(nf,2) .ne. 0 .or. modulo(ns,2) .ne. 0 ) &
write(6,*) ' Dimensions must be even - actually ', nf, ns
!write(6,'(/(4g14.6))') mapin

nf2 = nf/2
ns2 = ns/2

print *,'Position 1'


! Copy input array to output.
do j = 0,ns2
  mapout(0:nf2-1, 2*j)   = mapin(0:nf-1:2, j)
  mapout(0:nf2-1, 2*j+1) = mapin(1:nf-1:2, j)
end do

print *,'Position 2'

! and apply symmetry on first and last cols (note that nf is even)
mapout(nf2-1:(nf2+2)/2:-1, (/ 0, ns /))   = mapout(1:(nf2-1)/2, (/ 0, ns /))
mapout(nf2-1:(nf2+1)/2:-1, (/ 1, ns+1 /)) = mapout(0:(nf2-2)/2, (/ 1, ns+1 /))

print *,'Beyond problem area'

!write(6,'(/(3g14.6))') mapout
! check the symmetry of first and last cols.
if( any(mapout(1:nf2-1, 0) .ne. mapout(nf2-1:1:-1, 0)))	&
  Write(6,*) ' Symmetry error on col 0 - even'
if( any(mapout(0:nf2-1, 1) .ne. mapout(nf2-1:0:-1, 1)))	&
  Write(6,*) ' Symmetry error on col 0 - odd'
if( any(mapout(1:nf2-1, ns) .ne. mapout(nf2-1:1:-1, ns)))	&
  Write(6,*) ' Symmetry error on last col - even'
if( any(mapout(0:nf2-1, ns+1) .ne. mapout(nf2-1:0:-1, ns+1)))	&
  Write(6,*) ' Symmetry error on last col - odd'

temp = mapout

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
! dxml test

return
end

subroutine ffl_g1f(map,nh,nk)
! forward real -> cc group trans on second index of 2-d array 
! with negative in the exponent.
! nh,nk are the real space dimensions	
! the input map should be nh x 0:nk+1.
integer, intent(in)	:: nh, nk
real, intent(inout)	:: map(nh, 0:nk+1)
integer status, sfft_grp
if(mod(nk,2) .ne. 0) &
          stop 'FFL _g1f -- transform dimension must be even'

status = sfft_grp('r','c','f',map,map,nk,nh,nh,1,1)
if( status .ne. 0 ) write(6,*) 'ffl_g1f status ', status

return
end

SUBROUTINE FFL_g1r(map,nh,nk)
! reverse trans with positive in the exponent.
! nh,nk are the real space dimensions	
integer, intent(in)	:: nh, nk
real, intent(inout)	:: map(nh, 0:nk+1)

integer status, sfft_grp

if(mod(nk,2) .ne. 0) &
          stop 'FFL_g1r -- transform dimension must be even'

status = sfft_grp('c','r','b',map,map,nk,nh,nh,1,1)
if( status .ne. 0 ) write(6,*) 'ffl_g1r status ', status
map = map*nk	! the dxml backward transforms divides by n.
return
end
