program map_to_sf
! this is intended to perform space-group-specific transforms of maps
! to structure factors.
! Note:- The sign convention for Fourier transforms is
!    map  -> structure factors 
! uses  exp +i  
! Ten Eyck routines are all exp -i
! DXML use exp -i for "forward" transforms, exp +i for 'inverse',
!   and also scale the inverse by 1/n.
!
use space_group
use mapdefs
type (map), target	::	rho, rhotemp
type (sg_element), pointer	:: group

interface
  subroutine ffl_ortho_2f(rho)
  use mapdefs
  type (map) rho
  end subroutine ffl_ortho_2f
end interface
interface
  subroutine mapread(rho)
  use mapdefs
  type (map) rho
  end subroutine mapread
end interface
interface
  subroutine write_mtz(rho, smin, smax)
  use mapdefs
  type (map) rho
  real smin, smax
  end subroutine write_mtz
end interface

character*80 title
real smin, smax

integer allocstat, pretend_sg
logical fullcheck

!call lib$init_timer
!call printstats('Starting')
! I suppose we'd better make this look like a ccp4 program.
call ccpfyp
CALL CCPRCS(6,'map_to_mtz','31-MAY-1996 12:26:11')

! and I suppose we ought to read 'input cards' here.
! resolution limits, title. (Column labels?...)

call read_options(title, smin, smax, fullcheck, pretend_sg)

call mapread(rho)
!call printstats('After reading map')
rho%true_spgrp = rho%lspgrp

!write(6,*) ' Shape of input rho and total size',shape(rho%rho),size(rho%rho),&
!lbound(rho%rho),ubound(rho%rho)
!write(6,*)   rho%title, &
!rho%limits, &
!rho%extlim, &
!rho%iuvw,rho%mxyz, &
!rho%lspgrp,rho%lmode, &
!rho%cell,rho%rhmin,rho%rhmax,rho%rhmean,rho%rhrms

!write(6,'((1x8f10.2))') rho%rho(:,rho%limits(2,1),rho%limits(2,2))

call cell_metric(rho%cell, rho%basis, rho%metric, rho%recip_metric)

call space_group_setup(rho)

write(6,*) 'Rho limits ',rho%limits
write(6,*) 'Rho au_box ',rho%au_box
write(6,*) 'Rho sf_box ',rho%sf_box

write(6,*) 'Space group generators'

call print_sg(rho%generators)

write(6,*) ' Generating space-group operators '

call gen_to_group(rho%generators, group)

call group_to_grid(rho%mxyz, group, rho%grid_ops)
write (6,*) ' Grid operators created' 

if ( fullcheck) call print_grid_ops(rho%grid_ops)

call check_symmetry(rho, 3)
!call printstats('After face symmetry check')
if ( fullcheck) call check_symmetry(rho, 4)

!if ( fullcheck) call printstats('After symmetry check')

if ( pretend_sg .ne. 0 ) then	! the true spacegroup is unaltered.
	call change_sg(rho, rhotemp, pretend_sg)
	if ( fullcheck) CALL print_grid_ops(rhotemp%grid_ops)
	if ( .not. associated( rho%rho, rhotemp%rho ) ) deallocate(rho%rho)
	rho = rhotemp		! but we should deallocate everything!
!	call printstats('After extraction')
	if ( fullcheck) call check_symmetry(rho, 4)
!	if ( fullcheck) call printstats('After symmetry check')
end if

if( any( rho%Limits .ne. rho%au_box ) ) then

  write(6,'(a40,6i5)') ' Grid limits are non-canonical: ', rho%limits

  rhotemp = rho		! copy map to temp (including pointers)
  nullify ( rho%rho )
  rho%Limits = rho%au_box 
  write(6,'(a40,6i5)') ' Extracting to new grid: ', rho%limits

  allocate(rho%rho(rho%limits(1,1):rho%limits(1,2),	&
	rho%limits(2,1):rho%limits(2,2),	&
	rho%limits(3,1):rho%limits(3,2)),	&
	stat=allocstat)
	if ( allocstat .ne. 0 ) then
		write(6,*) ' status returned from allocate ',allocstat
!		call lib$stop(allocstat)
	end if

  call extract( rhotemp, rho )
!  call printstats('After extracting map')

  deallocate ( rhotemp%rho )
!  call printstats('After deallocate old map')
end if

! allocate space for the structure factors

allocate(rho%sfs(rho%sf_box(1,1):rho%sf_box(1,2),	&
	rho%sf_box(2,1):rho%sf_box(2,2),	&
	rho%sf_box(3,1):rho%sf_box(3,2)),	&
	stat=allocstat)

if ( allocstat .ne. 0 ) then
	write(6,*) ' status returned from allocate ',allocstat
!	call lib$stop(allocstat)
end if

!call printstats('After alloc for Fs')

select case (rho%lspgrp)

case (16, 18 ) 	! P222, P2_12_12, I222	!   23 soon

  if(fullcheck)call check_symmetry(rho,3)
!  call printstats('After face symmetry check')
! this is difficult to do efficiently, except in a s-g-specific way.

! perform the transform
  call ffl_ortho_2f(rho)

case (1)
print *, 'Performing P1 transform.'
print *, rho%mxyz(1)
print *, rho%mxyz(2)
print *, rho%mxyz(3)

  call ffl_p1f(rho)
   
case default
  stop ' This spacegroup not yet implemented '
end select

! write out the results.
!call printstats('After transform')

call write_mtz(rho, smin, smax)

!call printstats('Final stats')
stop

end program map_to_sf

!subroutine printstats(string)
!
!character*(*) string
!
!integer status, frecnt, cpu, faults, pagfilcnt
!integer*8 quadtime
!character*11 elapsed
!include '($jpidef)'
!integer lib$getjpi, lib$stat_timer, lib$sys_asctim
!external lib$getjpi, lib$stat_timer, lib$sys_asctim
!
!status = lib$getjpi(jpi$_freptecnt,,, frecnt)
!	if ( .not. status ) call lib$signal(%val(status))
!status = lib$getjpi(jpi$_pagfilcnt,,, pagfilcnt)
!	if ( .not. status ) call lib$signal(%val(status))
!status = lib$stat_timer(1,quadtime)		! in system time format.
!	if ( .not. status ) call lib$signal(%val(status))
!status = lib$sys_asctim(,elapsed,quadtime,1)
!	if ( .not. status ) call lib$signal(%val(status))
!status = lib$stat_timer(2,cpu)
!	if ( .not. status ) call lib$signal(%val(status))
!status = lib$stat_timer(5,faults)
!	if ( .not. status ) call lib$signal(%val(status))
!
!write(6,'(1x,a/1x,a,a11,f10.2,i10/1x,a,2i10)') string,      &
!'Cumul. elapsed time, cpu and faults ',elapsed,0.01*cpu,faults, &
!'Remaining vm pages and pagefilequota',frecnt,pagfilcnt
!
!!write(6,'(1x,2a,2i10)') string,      &
!!', remaining vm pages and pagefilequota',frecnt,pagfilcnt
!!call lib$show_timer()
!end subroutine printstats
