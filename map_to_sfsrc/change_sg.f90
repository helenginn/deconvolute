! This creates a new map structure with a different space-group
! to the original. The generators, canonical boxes, etc., of the
! new map have to be constructed, and then a new map extracted.
! It might be the case that the original limits satisfy the
! new box limits, in which case a pointer assignment rather
! than an 'extract' is all that is necessary.
SUBROUTINE change_sg(rhoold, rhonew, new_sg)

USE space_group
USE mapdefs

TYPE (map)	::	rhoold, rhonew
INTEGER		::	new_sg

TYPE (sg_element), POINTER	::	group
integer allocstat


rhonew		= rhoold	! copy everything
rhonew%lspgrp	= new_sg	! and make the changes ...
NULLIFY(rhonew%generators, rhonew%grid_ops)

WRITE(6,*) ' Setting up new spacegroup',new_sg
CALL space_group_setup(rhonew)

WRITE(6,'(1x,a12,6i5)') 'Rho limits ',rhonew%limits
WRITE(6,'(1x,a12,6i5)') 'Rho au_box ',rhonew%au_box
WRITE(6,'(1x,a12,6i5)') 'Rho sf_box ',rhonew%sf_box
     
!! WRITE(6,*) 'Space group generators'

!! CALL print_sg(rhonew%generators)

WRITE(6,*) ' Generating space-group operators '

CALL gen_to_group(rhonew%generators, group)
    
CALL group_to_grid(rhonew%mxyz, group, rhonew%grid_ops)
WRITE (6,*) ' Grid operators created'
    
IF ( .NOT. ( rhonew%grid_ops .subgroupof. rhoold%grid_ops ) )		&
	WRITE(6,'(a,i4,a,i4)') ' Warning -- spacegroup ',rhonew%lspgrp,	&
	' is not a subgroup of ',rhoold%lspgrp

rhonew%limits = rhonew%au_box	! set the canonical box

IF ( ALL( rhonew%limits .eq. rhoold%limits ) ) THEN
	rhonew%rho => rhoold%rho
	WRITE(6,*) ' Original box ok, pointer copied.'
ELSE
	WRITE(6,'(a40,6i5)') ' Extracting to new grid: ', rhonew%limits

	ALLOCATE(rhonew%rho(rhonew%limits(1,1):rhonew%limits(1,2),	&
		rhonew%limits(2,1):rhonew%limits(2,2),	&
		rhonew%limits(3,1):rhonew%limits(3,2)),	&
		STAT=allocstat)
	IF ( allocstat .ne. 0 ) THEN
		WRITE(6,*) ' status returned from allocate ',allocstat
!		CALL lib$stop(allocstat)
	END IF
	CALL extract( rhoold, rhonew )
END IF

RETURN
END SUBROUTINE change_sg
