! Oh, dear. I don't really want to use the ccp4 parser, but...
	subroutine read_options(title, smin, smax, fullcheck, pretend_sg)

	character*80 title
	real Smin,Smax	! these are in A^-2
	logical fullcheck
	integer pretend_sg

      integer maxtok
* Max number of tokens for 'parser'
      parameter (maxtok = 40)
      
      integer dummy, ntok
      integer ibeg(maxtok), iend(maxtok), ityp(maxtok), idec(maxtok)
      real fvalue(maxtok)
      character str*500
      character*4 key, cvalue(maxtok)
	integer itok
	real reSmin,reSmax

      ntok = maxtok
	title = ' '
	smin = 1.0
	smax = 1.0
	fullcheck = .false.
	pretend_sg = 0

* Parse user input
	do while(.true.)
	read (*,'(a)',end=500) str

	call parser (key, str, ibeg, iend, ityp, fvalue, cvalue,
     +	 	idec, ntok, dummy, .true.)


	select case (key)

	case('TITL')
		title = str(ibeg(2):iend(ntok))
	case('RESO')
		ITOK = 2
		CALL RDRESO(ITOK,ITYP,FVALUE,NTOK,ResMin,ResMax,Smin,Smax)
!	print *,ResMin,ResMax,Smin,Smax
	case ('CHEC')
		fullcheck = .true.

	case ('END', 'GO')
		exit

	case ('PRET')
		pretend_sg = fvalue(2)

	case default
		write(6,*) ' Keyword '//key//' unrecognised and ignored.'
!!	        call ccperr (1, 'Bad keyword: '//key)
	end select

	end do

!*     LABOUT keyword
!      else if (key .eq. 'LABO') then
!         k = 0
!         if (ntok.gt.maxlab) call ccperr (1, 'Too many labels')
!         do 20 i = 2, ntok
!            k = k+1
!            lsprgo(k) = str(ibeg(i):iend(i))
! 20      continue
!         nlprgo = (ntok-1)
!         write (*,*) 'Number of columns to be read in:', nlprgo
!         labout = .true.
!         goto 10


 500  continue

	return

	end ! subroutine read_options
