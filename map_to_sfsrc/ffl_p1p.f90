! This routine performs an in-place real -> cconjg transform
! on a 2D array.
! This array must be nf+2 x ns, for a nf x ns real input,
! going to  nf/2+1 x ns complex output.
! This transform has a negative in the exponent.
subroutine ffl_p1p(map, nf, ns)

integer, intent(in)	:: nf, ns
real, intent(inout)	:: map(0:nf+1, 0:ns-1)

integer d(5)	! Ten Eyck dimensioning array.
integer nf2

if(modulo(nf,2) .ne. 0) &
write(6,*) ' FFl_p1p: first dimension must be even - actually ', nf, ns

nf2 = nf/2

d(1) = ( nf + 2 ) * ns
d(2) = 2
d(3) = nf + 2
d(4) = 2
d(5) = 2
call realft(map, map(1,0), nf/2, d )
d(2) = nf + 2
d(3) = (nf + 2) * ns
d(4) = nf + 2
call cmplft(map, map(1,0), ns, d )

return
end


subroutine ffl_gcf(map,nh,nk)
! complex -> complex group trans on second index of 2-d array 
! with negative in the exponent.
! nh,nk are the real space dimensions	
! The complex numbers are conventionally adjacent elements in
! the first dimension.
integer, intent(in)	:: nh, nk
real, intent(inout)	:: map(nh, nk)

integer index(5)

index(1) = nh*nk
index(2) = nh
index(3) = index(1)
index(4) = nh
index(5) = 2
call cmplft(map(1,1), map(2,1), nk, index)
return
end
