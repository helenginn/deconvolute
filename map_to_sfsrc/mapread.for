C----------------------------------------------------------------------
C
C Read a map either in CCP4 format, or (Mapexchange) formatted
C
C  Input files:
C     MAPIN (stream 1)
C               input map, formatted (ascii) or binary CCP4 format. The program
C               reads the first line of the map to determine automatically
C               whether the file is formatted or not
C     SYMOP (stream 4)
C               symmetry library, required if the output file is CCP4 format
C
C Argument is a map structures, defined in mapdefs
	subroutine mapread(rho)	! and other args
C     ==================

	use mapdefs
	type (map) rho
	logical lform
	integer nsec,nw1,nw2,nu1,nu2,nv1,nv2
	real scale,plus
	real, allocatable	:: temp(:,:)

	integer allostat

C Open input file & find out if it is formatted or not
C LFORM is true if input file in formatted
	call opnfil(lform)
	print *,'Result from opnfil',lform
C Read header information
	call rdhead

C Check mode
	if (rho%lmode .ne. 2) then
	   write(6,*) '****** Only mode 2 maps supported ******'
	   call ccperr(1,' stop mode wrong')
	endif

C Check array size
!	nsize = (nu2-nu1+1)*(nv2-nv1+1)
	allocate(rho%rho(rho%limits(1,1):rho%limits(1,2),
     +	 	rho%limits(2,1):rho%limits(2,2),
     +	 	rho%limits(3,1):rho%limits(3,2)),
     +	 	stat=allostat)
	if ( allostat .ne. 0 ) then
		write(6,*) ' status returned from allocate ',allostat
!		call lib$stop(allostat)
	end if
	allocate(temp(rho%extlim(1,1):rho%extlim(1,2),
     +	 	rho%extlim(2,1):rho%extlim(2,2)))

c	write(6,*) ' Shape of input rho ',shape(rho%rho)
c	write(6,*)   rho%title, 
c	1	rho%limits, 
c	2	rho%extlim, rho%iuvw,rho%mxyz,nsec,
c	3	nw1,nw2,nu1,nu2,nv1,nv2, 
c	4	rho%lspgrp,rho%lmode, 
c	5	rho%cell,rho%rhmin,rho%rhmax,rho%rhmean,rho%rhrms

C Copy sections
	call readsections
	return

	contains

	subroutine rdhead
C     ==================================================================

C Read header information
C
	if (lform) then
C  Formatted, read keyworded header
	   call rdfhdr(1)
	else
C  Read Map format (CCP4) header
	   call mrdhdr(1,'MAPIN',rho%title,nsec,rho%iuvw,rho%mxyz,
     +	 	nw1,nu1,nu2,nv1,nv2,rho%cell,
     +	 	rho%lspgrp,rho%lmode,rho%rhmin,rho%rhmax,rho%rhmean,rho%rhrms)
C  Calculate scale & offset to put densities in range -998 to +9998 
C	 (for safety)
	   scale = (+9998. + 998.)/(rho%rhmax-rho%rhmin)
	   plus  = -998. - scale*rho%rhmin
	nw2 = nw1 + nsec - 1
C (sigh) this is a bit pedestrian...	
	rho%extlim = reshape( (/nu1,nv1,nw1,nu2,nv2,nw2 /), (/3,2/) )

	rho%Limits(rho%iuvw(:),:) = rho%extlim
	endif

	return
	end subroutine rdhead

	subroutine readsections

C read all sections into rho
	integer k,nu,nv,ier

	nu = nu2-nu1+1
	nv = nv2-nv1+1

	write(6,*) 'Bounds of rho ', lbound(rho%rho),ubound(rho%rho)
	do k = nw1,nw1+nsec-1

	   if (lform) then
		call rdfsec(1,ier,temp,nu,nv)
		if (ier .ne. k) then
		   write(6,*) '****** Wrong section read :',ier
		   call ccperr(1,' stop section error')
		endif
	   else
		call mgulp(1,temp,ier)
	   endif

! copy current section to location in main array
! this copies from fast med slow to x y z, which we do in a mundane way
! Notice that yzx and zyx ording give particularly bad indexing.
	select case(rho%iuvw(1))
	case (1)			! x fast
		select case(rho%iuvw(2))
		case (2)		! y medium
			rho%rho(:,:,k) = temp
		case (3)		! z medium
			rho%rho(:,k,:) = temp
		end select
	case (2)			! y fast
		select case(rho%iuvw(2))
		case (1)		! x medium
			rho%rho(:,:,k) = transpose(temp)
		case (3)		! z medium
			rho%rho(k,:,:) = temp
		end select
	case (3)			! z fast
		select case(rho%iuvw(2))
		case (1)		! x medium
			rho%rho(:,k,:) = transpose(temp)
		case (2)		! y medium
			rho%rho(k,:,:) = transpose(temp)
		end select
	end select

!	select case(rho%iuvw(1)*3+rho%iuvw(2) )
!	case(5)
!			rho%rho(:,:,k) = temp	! x fast, y medium
!	case(6)
!			rho%rho(:,k,:) = temp	! x fast, z medium
!	case(7)
!			rho%rho(:,:,k) = transpose(temp)  ! y fast, x medium
!	case(9)
!			rho%rho(k,:,:) = temp	! y fast, z medium
!	case(10)
!			rho%rho(:,k,:) = transpose(temp)  ! z fast, x medium
!	case(11)
!			rho%rho(k,:,:) = transpose(temp)  ! z fast, y medium
!	case default
!		print *,' Mess on axis ordering'
!	end select
	end do

	return
	end subroutine readsections


	subroutine rdfsec(lun,isec,rho,nu,nv)
C     ================================================

C Read formatted section from unit lun

	integer lun,isec,nu,nv
	real rho(nu,nv)

	integer j,line(nu)

	read (lun,'(/7X,I8/)') isec

	do j=1,nv
	   read(lun,'(20I4)') line
		rho(:,j) = (float(line(:))-plus)/scale
	end do
	return
	end subroutine rdfsec

	subroutine rdfhdr(lun)
C     ===============================================================

C  Formatted, read keyworded header, fixed format
	integer lun
	integer templim(2,3)
	integer kdummy,kfail
	character*1 axes(3),xyz*3
	data xyz/'XYZ'/
	KDUMMY = 80
	KFAIL = 0
	call ccpdpn(lun,'MAPIN','READONLY','F',KDUMMY,KFAIL)

	read(lun,6000) 
 6000 format('MAPEXCHANGE HEADER')
	write(6,6100) 
 6100 format(' MAPEXCHANGE HEADER')
	read(lun,6010) rho%title
 6010 format(/A)
	write(6,6110) rho%title
 6110 format(' TITLE'/1X,A)
	read(lun,6020) axes
 6020 format('AXIS    ',3(7X,A1))
	write(6,6120) axes
 6120 format(' AXIS    ',3(7X,A1))
	read(lun,6030) rho%mxyz
 6030 format('GRID    ',3I8)
	write(6,6130) rho%mxyz
 6130 format(' GRID    ',3I8)
	read(lun,6040) templim
 6040 format('XYZLIM  ',6I8)
	write(6,6140) templim
	rho%limits = transpose(templim)
 6140 format(' XYZLIM  ',6I8)
	read(lun,6050) rho%lspgrp
 6050 format('SPACEGROUP',6X,I8)
	write(6,6150) rho%lspgrp
 6150 format(' SPACEGROUP',6X,I8)
	read(lun,6060) rho%lmode
 6060 format('MODE    ',I8)
	write(6,6160) rho%lmode
 6160 format(' MODE    ',I8)
	read(lun,6070) rho%cell
 6070 format('CELL    ',6F10.3)
	write(6,6170) rho%cell
 6170 format(' CELL    ',6F10.3)
	read(lun,6080) rho%rhmin,rho%rhmax,rho%rhmean,rho%rhrms
 6080 format('RHOLIM  ',4G16.6)
	write(6,6180) rho%rhmin,rho%rhmax,rho%rhmean,rho%rhrms
 6180 format(' RHOLIM  ',4G16.6)
	read(lun,6090) scale,plus
 6090 format('PACK    ',2G16.6)
	write(6,6190) scale,plus
 6190 format(' PACK    ',2G16.6)
	read(lun,6095) 
 6095 format('END HEADER')
	write(6,6195) 
 6195 format(' END HEADER')

	rho%iuvw = index(xyz,axes)	! array assignment
	rho%extlim = rho%Limits(rho%iuvw,:)
      nu1 = rho%limits(rho%iuvw(1),1)
      nu2 = rho%limits(rho%iuvw(1),2)
      nv1 = rho%limits(rho%iuvw(2),1)
      nv2 = rho%limits(rho%iuvw(2),2)
      nw1 = rho%limits(rho%iuvw(3),1)
      nw2 = rho%limits(rho%iuvw(3),2)
	nsec= rho%extlim(3,2) - rho%extlim(3,1) + 1
	return
	end subroutine rdfhdr
	end 


	subroutine opnfil(lform)
C     ========================

C Find out if input file MAPIN is formatted or not
C Leave file closed

C LFORM set true if input file in formatted
	logical lform
	integer index
	character*80 line,filename
	integer kdummy,kfail

	lform = .false.
	KDUMMY = 80
	KFAIL = 0

c	call ccpdpn(1,'MAPIN','READONLY','F',KDUMMY,KFAIL)

	call ugtenv('MAPIN',filename)
	open(1,file=filename,status='OLD')

c	inquire(1,opened=lform,name=line)
c	if(lform)then
c	  print *,'Unit 1 is connected to file ',line
c	else
c	  print *,'Unit 1 is not connected'
c	endif

	lform=.false.


C Read 1st line as is formatted, should contain string MAPEXCHANGE if so
	read (1,6001,err=100) line
 6001 format(A)

	lform = index(line,'MAPEXCHANGE') .gt. 0

C Close file
 100  close (unit=1)
	return 
	end

