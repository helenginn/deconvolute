! This routine checks the symmetry of the map, by applying 
! each grid operator in turn.
! The check levels available are:-
! 0	no check
! 1	corners
! 2	edges
! 3	faces
! 4	entire box
! Obviously each level includes all the lower levels.
! Note that the new scheme for boxes (3,2) is used.
! The performance should be checked; particularly for level = 4,
! it should be compared with a fixed stride of 1 in the loops.
SUBROUTINE check_symmetry(rho, level)
USE space_group
USE mapdefs
TYPE (map), INTENT(IN)	:: rho
INTEGER, INTENT(IN)	:: level

TYPE (grid_element), POINTER	:: curgrid
INTEGER n1,n2,n3, n1step, n2step, n3step
LOGICAL n2ends, n3ends
INTEGER, PARAMETER	:: maxprint = 100
INTEGER npoints, ntest, nfail, nprint
INTEGER related(3), stride(3)
REAL diff_limit, max_fail, difference
INTEGER maxbas(3), maxrel(3)

diff_limit = rho%rhrms * 1e-4

npoints = 0
ntest = 0
nfail = 0
max_fail = 0.0
nprint = maxprint

! Cycle over space_group operators.
curgrid => rho%grid_ops
DO
  IF ( level .GE. 2 ) THEN
    n3step = 1
  ELSE
    n3step = rho%limits(3,2) - rho%limits(3,1)
  END IF
  DO n3 = rho%limits(3,1), rho%limits(3,2), n3step
    n3ends = n3 .EQ. rho%limits(3,1) .OR. n3 .EQ. rho%limits(3,2)
    IF ( level .GE. 3 .OR. ( level .EQ. 2 .AND. n3ends ) ) THEN
      n2step = 1
    ELSE
      n2step = rho%limits(2,2) - rho%limits(2,1)
    END IF
    DO n2 = rho%limits(2,1), rho%limits(2,2), n2step
      n2ends = n2 .EQ. rho%limits(2,1) .OR. n2 .EQ. rho%limits(2,2)
      IF ( level .GE. 4 .OR. &
	 ( level .GE. 3 .AND.( n2ends .OR. n3ends ) ) .OR.&
	 ( level .GE. 2 .AND. n2ends .AND. n3ends ) ) THEN
	n1step = 1
      ELSE
	n1step = rho%limits(1,2) - rho%limits(1,1)
      END IF
      stride = MATMUL( curgrid%op%rot, (/ n1step,0,0 /) )
      related = curgrid%op * (/ rho%limits(1,1)-n1step, n2, n3 /)
      DO n1 = rho%limits(1,1), rho%limits(1,2), n1step
	npoints = npoints + 1
	related = MODULO( related + stride, rho%mxyz )
!!! calculation check (not for production)
!!IF ( any ( related .ne. curgrid%op * (/ n1,n2,n3 /) ) ) THEN
!!print *, ' Check symm operator error at',n1,n2,n3
!!print *, ' Operator', curgrid%op%rot, curgrid%op%trans, curgrid%op%grid
!!print *, ' Related point',related
!!print *, ' n1 stride', stride
!!print *, ' Directly calculated point ',curgrid%op * (/ n1,n2,n3 /)
!!call print_grid_ops(rho%grid_ops)
!!stop
!!END IF
	IF ( ANY( related .LT. rho%limits(:,1) ) .OR.	&
	     ANY( related .GT. rho%limits(:,2) ) ) CYCLE
	ntest = ntest + 1
	difference = ABS( rho%rho(n1, n2, n3 ) -		&
	      rho%rho( related(1), related(2), related(3) ) )
	IF (difference .GT. diff_limit ) nfail = nfail + 1
	IF (difference .GT. 100.0*diff_limit .AND. nprint .GE. 0) THEN
	  IF ( nprint .EQ. maxprint ) PRINT '(a,i5,a)',' First',	&
		maxprint,' differences more than rhorms/100:-'
	  PRINT '(1x,g14.6,3i5,a4,g14.6,3i5)', rho%rho(n1,n2,n3),	&
	    n1,n2,n3, ' <=>', &
	    rho%rho( related(1), related(2), related(3) ), related
	  nprint = nprint - 1
	END IF
	IF ( difference .GT. max_fail ) THEN
	  max_fail = difference
	  maxbas = (/ n1, n2, n3 /)
	  maxrel = related
	END IF
      END DO
    END DO
  END DO
  IF ( .NOT. ASSOCIATED( curgrid%next ) ) EXIT
  curgrid => curgrid%next
END DO

PRINT '(/a/a40,f10.6/3(a40,i10/),a40,2g14.6/a40,3i5,a4,3i5)', &
' Symmetry test statistics....',	&
' Checking to accuracy of rhorms times:', diff_limit/rho%rhrms, &
' Number of map points:', npoints,	&
' Number of tests:', ntest,		&
' Number of failures:', nfail,		&
' Maximum failure and rhorms:', max_fail, rho%rhrms,	&
' Point of max failure:', maxbas, ' <=>', maxrel

RETURN
END
