	subroutine fftreal_recip(map,trans,ntot,nh,nk,nl)
C forward trans with negative in the exponent.
C nh,nk,nl are the real space dimensions	
	implicit none
	integer ntot,nh,nk,nl
	real map(ntot),trans((nh+2)*nk*nl)
	integer nhp2,npr,index(5)
c the complex form is 1+nh/2,nk,nl

	if(mod(nh,2) .ne. 0) 
     1          stop 'FFTreal_recip -- first dimension must be even'
	nhp2 = nh+2
	npr = nhp2*nk*nl
c shift each row to make room and replicate.
	call shift_up(map,trans,nh,nk*nl)

	index(5)=2
	index(1)=nhp2*nk*nl
	index(2)=index(5)
	index(3)=nhp2
	index(4)=index(5)
	call realft(trans(1),trans(2),nh/2,index)
c	write(2,213)' After real ft',(trans(i),i=1,npr)
	index(2)=index(3)
	index(3)=index(3)*nk
	index(4)=index(2)
	call cmplft(trans(1),trans(2),nk,index)
c	write(2,213)' After complex y ft',(trans(i),i=1,npr)
	index(2)=index(3)
	index(3)=index(3)*nl
	index(4)=index(2)
	call cmplft(trans(1),trans(2),nl,index)
c	write(2,213)' After complex z ft',(trans(i),i=1,npr)
	return
	end
	SUBROUTINE FFTrecip_real(trans,map,ntot,nh,nk,nl)
C reverse trans with positive in the exponent.
C nh,nk,nl are the real space dimensions	
	implicit none
	integer ntot,nh,nk,nl
	real map(ntot),trans((nh+2)*nk*nl)
	integer nhp2,npr,index(5),i
c the complex form is 1+nh/2,nk,nl

	if(mod(nh,2) .ne. 0) 
     1          stop 'FFTrecip_real -- first dimension must be even'

	nhp2 = nh+2
	npr = nhp2*nk*nl
c reverse sign on imaginaries
	index(5)=2
	index(1)=nhp2*nk*nl
	do i=2,index(1),2
		trans(i)=-trans(i)
	end do
c	write(2,213)' After im sign rev.',(trans(i),i=1,npr)
c transform in z and y
	index(2)=nhp2*nk
	index(3)=index(2)*nl
	index(4)=index(2)
	call cmplft(trans(1),trans(2),nl,index)
c	write(2,213)' After z trans',(trans(i),i=1,npr)
	index(2)=nhp2
	index(3)=index(2)*nk
	index(4)=index(2)
	call cmplft(trans(1),trans(2),nk,index)
c	write(2,213)' After y trans.',(trans(i),i=1,npr)
c store real(f(nh/2)) in im(f(0))
	call t_to_t(trans,nh+2,nk*nl)
	index(2)=index(5)
	index(3)=nhp2
	index(4)=index(5)
	call hermft(trans(1),trans(2),nh/2,index)
c	write(2,213)' After herm trans.',(trans(i),i=1,npr)
c compress rows
	call shift_down(map,trans,nh,nk*nl)
c	write(2,212)' After compression.',(trans(i),i=1,npr)
	return
	end
