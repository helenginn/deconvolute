! This subroutine sets various space-group specific stuff:-
! 	canonical bounding box
!	canonical structure factor box
!	space group generators
! Which means lots of case statements.

SUBROUTINE space_group_setup(rho)
USE space_group
USE mapdefs
TYPE (map)			::	rho

TYPE (sg_operator), target	::	&
  xscr =  sg_operator(reshape((/1,0,0,0,-1,0,0,0,-1/),(/3,3/)), (/12,12,0/)), &
  ztfold  = sg_operator(reshape((/-1,0,0,0,-1,0,0,0,1/),(/3,3/)), (/0,0,0/)), &
  xtfold  = sg_operator(reshape((/1,0,0,0,-1,0,0,0,-1/),(/3,3/)), (/0,0,0/)), &
  Icent = sg_operator(unit_rot , (/12,12,12/)), &
  permxyz = sg_operator(reshape((/0,1,0,0,0,1,1,0,0/),(/3,3/)), (/0,0,0/))

TYPE (sg_element), pointer	::	last
TYPE (sg_element), target 	::	first

NULLIFY(first%next)
last => first

SELECT CASE ( rho%lspgrp )	! Space group operators

CASE (1)		! P1
  CALL add_op(last, unit_op )

CASE ( 18 )		! P2_1 2_1 2

  CALL add_op(last, xscr)
  CALL add_op(last, ztfold)

CASE ( 16, 23 )		! P2 2 2, I222

  CALL add_op(last, xtfold)
  CALL add_op(last, ztfold)
  if ( rho%lspgrp .eq. 23 )  CALL add_op(last, Icent)

CASE ( 197 )		! I23

  CALL add_op(last, xtfold)
  CALL add_op(last, ztfold)
  CALL add_op(last, Icent)
  CALL add_op(last, permxyz)

CASE DEFAULT
  STOP 'This space group not yet implemented.'

END SELECT

rho%generators => first%next

SELECT CASE ( rho%lspgrp )	! Bounding boxes

CASE (1)

  rho%au_box(:,1) = (/0,0,0/)
  rho%au_box(:,2) = (/ rho%mxyz(1)-1, rho%mxyz(2)-1, rho%mxyz(3)-1 /)

  rho%sf_box(:,1) =  (/0,0,0/)
  rho%sf_box(:,2) =  (/ rho%mxyz(1)+1, rho%mxyz(2)-1, rho%mxyz(3)-1 /)

CASE ( 16, 18 )		! Most orthorhombic .... (?)

  rho%au_box(:,1) = (/0,0,0/)
  rho%au_box(:,2) = (/ rho%mxyz(1)-1, rho%mxyz(2)/2, rho%mxyz(3)/2 /)

  rho%sf_box(:,1) =  (/0,0,0/)
  rho%sf_box(:,2) =  (/ rho%mxyz(1)/2, rho%mxyz(2)/2, rho%mxyz(3)+1 /)

CASE ( 23, 197 )	! conventionlly, I23 stored as I222

  rho%au_box(:,1) = (/0,0,0/)
  rho%au_box(:,2) = (/ rho%mxyz(1)-1, rho%mxyz(2)/2, rho%mxyz(3)/4 /)

  rho%sf_box(:,1) =  (/0,0,0/)
  rho%sf_box(:,2) =  (/ rho%mxyz(1)/2, rho%mxyz(2)/2, rho%mxyz(3)/2+1 /)

CASE DEFAULT
  STOP 'This space group not yet implemented.'

END SELECT

RETURN

CONTAINS

SUBROUTINE add_op( current, oper )
TYPE (sg_operator) 		::	oper
TYPE (sg_element), POINTER	::	current

  ALLOCATE(current%next)
  current => current%next
  current%op = oper		! operator is actually copied
  NULLIFY(current%next)

RETURN
END SUBROUTINE add_op

END SUBROUTINE space_group_setup
