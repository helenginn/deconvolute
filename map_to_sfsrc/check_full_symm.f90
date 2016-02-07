! This routines checks any symmetry of the map that is defined by the
! grid operators of rho, by applying each symmetry operator
! to each point in the box.
subroutine check_full_symm(rho)
use space_group
use mapdefs
type (map), intent(in)	:: rho

type (grid_element), pointer	:: curgrid
type (grid_element), target	:: startgrid
integer n1,n2,n3
integer npoints, ntest, nfail, nprint
integer base(3), related(3)
real diff_limit, max_fail, difference
integer maxbas(3), maxrel(3)

startgrid%next => rho%grid_ops
diff_limit = rho%rhrms * 1e-4

npoints = 0
ntest = 0
nfail = 0
max_fail = 0.0
nprint = 100


!!    curgrid => startgrid
!!    do
!!      if ( .not. associated( curgrid%next ) ) exit
!!      curgrid => curgrid%next
!!do n3 = rho%limits(1,3), rho%limits(2,3)
!!do n2 = rho%limits(1,2), rho%limits(2,2)
!!do n1 = rho%limits(1,1), rho%limits(2,1)
!!  npoints = npoints + 1
!!
!!  base = (/ n1, n2, n3 /)
!!      related = curgrid%op * base
!!      if ( .not. (related .within. rho%limits) ) cycle
!!      ntest = ntest + 1
!!
!!	difference = abs( rho%rho(n1,n2,n3) -	&
!!		rho%rho( related(1), related(2), related(3) ) )
!!	if( difference .gt. diff_limit ) nfail = nfail + 1
!!	if( difference .gt. 100.0*diff_limit .and. nprint .gt. 0 ) then
!!	if ( nprint .eq. 100 ) write(6,*) &
!!		' First 100 differences more than rhorms/100:-'
!!	write(6,'(1x,g14.6,3i5,a4,g14.6,3i5)') rho%rho(n1,n2,n3), base, ' <=>',&
!!	rho%rho( related(1), related(2), related(3) ), related
!!	nprint = nprint - 1
!!	end if
!!	if ( difference .gt. max_fail ) then
!!		max_fail = max(max_fail, difference)
!!		maxbas = base
!!		maxrel = related
!!	end if
!!end do
!!end do
!!end do
!!    end do
!!
do n3 = rho%limits(3,1), rho%limits(3,2)
do n2 = rho%limits(2,1), rho%limits(2,2)
do n1 = rho%limits(1,1), rho%limits(1,2)
  npoints = npoints + 1

  base = (/ n1, n2, n3 /)
    curgrid => startgrid
    do
      if ( .not. associated( curgrid%next ) ) exit
      curgrid => curgrid%next
      related = curgrid%op * base
      if ( .not. (related .within. rho%limits) ) cycle
      ntest = ntest + 1

	difference = abs( rho%rho(n1,n2,n3) -	&
		rho%rho( related(1), related(2), related(3) ) )
	if( difference .gt. diff_limit ) nfail = nfail + 1
	if( difference .gt. 100.0*diff_limit .and. nprint .gt. 0 ) then
	if ( nprint .eq. 100 ) write(6,*) &
		' First 100 differences more than rhorms/100:-'
	write(6,'(1x,g14.6,3i5,a4,g14.6,3i5)') rho%rho(n1,n2,n3), base, ' <=>',&
	rho%rho( related(1), related(2), related(3) ), related
	nprint = nprint - 1
	end if
	if ( difference .gt. max_fail ) then
		max_fail = max(max_fail, difference)
		maxbas = base
		maxrel = related
	end if

    end do

end do
end do
end do

write(6,'(/a/a40,f10.6/3(a40,i10/),a40,2g14.6/a40,3i5,a4,3i5)') &
' Symmetry test statistics....', &
' Checking to accuracy of rhorms times:', diff_limit/rho%rhrms,	&
' Number of map points:', npoints,	&
' Number of tests:', ntest,		&
' Number of failures:', nfail,		&
' Maximum failure and rhorms:', max_fail, rho%rhrms,	&
' Point of max failure:', maxbas, ' <=>',maxrel

return
end
