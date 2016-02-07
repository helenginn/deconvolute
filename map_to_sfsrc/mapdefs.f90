! Structure definitions for maps.
module mapdefs
use space_group
! limits	upper and lower limits to xyz
! extlim	upper and lower limits to fast med slow (as in map file)
! iuvw		index to axis order (fast,med,slow)
! mxyz		number of grid points in unit cell in x,y,z
! nsec, nw1	number of slow sections, and starting value
! nu1,nu2,nv1,nv2	lower and upper bound for fast and medium axes
! rho contains map in xyz order (and limits gives the bounds)

type map	! this mostly contains the stuff for ccp4 routines
  character	::	title*80
  integer	::	limits(3,2)
  integer	::	extlim(3,2)
  integer	::	iuvw(3),mxyz(3)
  integer	::	lspgrp,lmode
  integer	::	true_spgrp
  real		::	cell(6),rhmin,rhmax,rhmean,rhrms
  integer	::	au_box(3,2)	! required map box
  integer	::	sf_box(3,2)	! resultant structure factor box
  real		::	basis(3,3), metric(3,3), recip_metric(3,3)
  real, pointer	::	rho(:,:,:)
  real, pointer	::	sfs(:,:,:)
  type (sg_element), pointer	::	generators
  type (grid_element), pointer	::	grid_ops
end type map

end module mapdefs
