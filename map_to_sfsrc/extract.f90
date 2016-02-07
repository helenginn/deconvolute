! Extracts rho2 from rho1, with matching cells but different bounding boxes.
subroutine extract(rho1, rho2)
use space_group
use mapdefs
type (map), intent(in)	:: rho1
type (map), intent(out)	:: rho2

! cycle over all points of output map.
! transform output point to input points, using grid operators of input map,
! until a point within the bounding box of the input map is found.
! retain operator for next go.

type (grid_element), pointer	:: curgrid, prevgrid
type (grid_element), target	:: startgrid
integer n1,n2,n3
integer npoints, nsearch, ntest, nignore
integer source(3), stride(3)

startgrid%next => rho1%grid_ops
curgrid => rho1%grid_ops

npoints = 0
nsearch = 0
ntest = 0
nignore = 0

do n3 = rho2%limits(3,1), rho2%limits(3,2)
do n2 = rho2%limits(2,1), rho2%limits(2,2)
source = curgrid%op * (/ rho2%limits(1,1)-1, n2, n3 /)
stride = MATMUL( curgrid%op%rot , (/ 1, 0, 0 /) )
do n1 = rho2%limits(1,1), rho2%limits(1,2)
  npoints = npoints + 1

!!  source = curgrid%op * (/ n1, n2, n3 /)
  source = modulo( source + stride, rho1%mxyz )
!!! Check (not for production)
!!  if ( any( source .ne. curgrid%op * (/ n1, n2, n3 /) ) ) then
!!print *,' Error on extract'
!!print *, n1,n2,n3,source,rho1%mxyz, &
!!curgrid%op * (/ n1, n2, n3 /), &
!!curgrid%op%rot,curgrid%op%trans, curgrid%op%grid
!!  end if

   IF ( ANY( source .LT. rho1%limits(:,1) ) .OR.	&
        ANY( source .GT. rho1%limits(:,2) ) ) THEN

    nsearch = nsearch + 1
    prevgrid => curgrid
    curgrid => startgrid
    do
      if ( .not. associated( curgrid%next ) ) then
	write(6,'(a,3i6)') &
	' Error extracting map; unable to find source point for', &
          n1,n2,n3
	stop
      end if
      curgrid => curgrid%next
      if ( associated( curgrid, prevgrid) ) nignore = nignore + 1
      if ( associated( curgrid, prevgrid) ) cycle
      source = curgrid%op * (/ n1, n2, n3 /)
      ntest = ntest + 1
   IF ( ALL( source .GE. rho1%limits(:,1) ) .AND.	&
        ALL( source .LE. rho1%limits(:,2) ) ) EXIT
    end do
    stride = MATMUL( curgrid%op%rot , (/ 1, 0, 0 /) )
  end if
!!! Check (not for production)
!!  if ( any( source .ne. curgrid%op * (/ n1, n2, n3 /) ) ) then
!!print *,' Error on extract'
!!print *, n1,n2,n3,source,curgrid%op * (/ n1, n2, n3 /)
!!  end if

  rho2%rho(n1,n2,n3) = rho1%rho( source(1), source(2), source(3) )

end do
end do
end do

write(6,'(/a/4(a40,i10/))') ' Statistics on map extraction ....', &
'Number of map points:',npoints,	&
'Number of operator searches:',nsearch,	&
'Number of tests:',ntest,		&
'Number of ignores:',nignore
return
end
