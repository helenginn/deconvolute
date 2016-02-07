	subroutine div_oddh_i(a,hklist,nhk,nl)

        integer nhk,nl,ihk,l
	integer hklist(2,nhk)
	complex a(nhk,nl)

	do ihk=1,nhk
		if(mod(hklist(1,ihk),2).eq.1) then
			do l=1,nl
			  a(ihk,l) = cmplx(aimag(a(ihk,l)),-real(a(ihk,l)))
			end do
		end if
	end do
	return
	end

	subroutine mult_oddh_i(a,hklist,nhk,nl)

        integer nhk,nl,ihk,l
	complex a(nhk,nl)
	integer hklist(2,nhk)

	do ihk=1,nhk
		if(mod(hklist(1,ihk),2).eq.1) then
			do l=1,nl
			  a(ihk,l) = cmplx(-aimag(a(ihk,l)),real(a(ihk,l)))
			end do
		end if
	end do
	return
	end

	subroutine shift_up(a,c,nh,nr)
	real a(1),c(1)
	integer h,nh,nr,l

	do  l = nr ,1,-1
		do h = nh,1,-1
			c(h+(nh+2)*(l-1)) = a(h+nh*(l-1))
		end do
	end do
	return
	end

	subroutine shift_down(a,c,nh,nr)
	real a(1),c(1)
	integer h,nh,nr,l

	do  l = 1 , nr
		do h = 1 , nh
			a(h+nh*(l-1)) = c(h+(nh+2)*(l-1))
		end do
	end do
	return
	end
	subroutine t_to_t(a,nhp2,nr)
c store re(n/2) in im(1)
        integer nhp2,nr,l
	real a(nhp2,nr)

	do l = 1,nr
		a(2,l) = a(nhp2-1,l)
	end do
	return
	end

	subroutine expand19(a,nhd,nh,nk,nlp1)
c the plane l=nl+1 is stored on the plane l=1 as
c   x= n/2-1... y=0..n/2-1   contains x=n/2...n-1 y=0..n/2-1
c        ""       n/2..n-1              0...n/2-1     " "
c appropriate bits are copied
c then x-screw on l=0 and y-screw on l=nl+1
        integer nhd,nh,nk,nlp1,nho2,k1,k2
	real a(nhd,nk,nlp1)
	integer h1,h2,h1a,h2a

	nho2 = nh/2
	do k1 = 1,nk/2
		k2 = k1+nk/2
		do h1 = 1,nho2
			h2 = h1+nho2
			h1a = mod(nh+1-h1,nh)+1
			h2a = mod(nh+1-h2,nh)+1

			a(h1,k1,nlp1) = a(h2,k2,1)
			a(h1a,k2,nlp1) = a(h2,k2,1)	! y screw
			a(h2,k1,nlp1) = a(h2,k1,1)
			a(h2a,k2,nlp1) = a(h2,k1,1)
		end do
	end do
c x screw on l=0
	do k1 = 1,nk
		k2 = mod(nk+nk/2+1-k1,nk)+1
		do h1 = 1,nho2
			h2 = h1+nho2
			a(h2,k2,1) = a(h1,k1,1)
		end do
	end do
	return
	end


	subroutine contract19(a,nhd,nh,nk,nlp1)
c the plane l=nl+1 is stored on the plane l=1 as
c   x= n/2-1... y=0..n/2-1   contains x=n/2...n-1 y=0..n/2-1
c        ""       n/2..n-1              0...n/2-1     " "
        integer nhd,nh,nk,nlp1,nho2,k1,k2
	real a(nhd,nk,nlp1)
	integer h1,h2

	nho2 = nh/2
	do k1 = 1,nk/2
		k2 = k1+nk/2
		do h1 = 1,nho2
			h2 = h1+nho2
			a(h2,k2,1) = a(h1,k1,nlp1)
			a(h2,k1,1) = a(h2,k1,nlp1)
		end do
	end do
	return
	end
