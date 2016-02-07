module space_group

! Definitions for types of coordinates

  interface assignment(=)
    module procedure op_sg_to_grid
  end interface

! Definitions for space-group stuff
! The space-group operators could be coded to say which axes interfered
! -- or for the sg as a whole. This would reduce 9 ops to 3 for orthor, etc.
! 5 possibilities :- none, x, y, z, all independent.
  type sg_operator		! this will be in fractional coords.
    integer	rot(3,3)
    integer	trans(3)	! this is scaled in 24 th/cell
  end type sg_operator

  integer, parameter		:: unit_rot(3,3) = &
		reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/))
  type (sg_operator), parameter :: unit_op = &
         sg_operator(unit_rot, (/0,0,0/))

  type grid_operator		! this will be in grid units
    integer	rot(3,3)	! 
    integer	trans(3)	! 
    integer, pointer	::	grid(:)		! Number of grid points
  end type grid_operator

  type sg_element
    type (sg_operator)		:: op
    type (sg_element), pointer  :: next
  end type sg_element
  type grid_element
    type (grid_operator)	:: op
    type (grid_element), pointer  :: next
  end type grid_element

  interface operator (*)
    module procedure sg_op_mult
  end interface

  interface operator (.eq.)
    module procedure sg_op_comp
  end interface

  interface operator (*)
    module procedure grid_op_mult
  end interface

  interface operator (.eq.)
    module procedure grid_op_comp
  end interface

  interface operator (*)
    module procedure grid_op_vector
  end interface

  interface operator (.within.)
    module procedure within_int
  end interface

  interface operator (.subgroupof.)
    module procedure subgroup_of_grid
  end interface

contains

  subroutine op_sg_to_grid(a,b)
    type (sg_operator), intent(in)  :: b
    type (grid_operator), intent(out) :: a

    a%rot = b%rot	! but this is only true if grid dims of a conform
    a%trans = modulo( (b%trans*a%grid)/24, a%grid )	! elemental
    return
  end subroutine op_sg_to_grid

subroutine cell_metric(cell, basis, metric, recip_metric)
!!!!!!!!!  so far, checked only for orthorhombic cells.

real cell(6)
real metric(3,3), basis(3,3), recip_metric(3,3), recbas(3,3)

    basis(:,1) = cell(1) * (/ 1.0,0.0,0.0 /)
    basis(:,2) = cell(2) * (/ cosd(cell(6)), sind(cell(6)), 0.0 /)
    basis(1,3) = cosd(cell(5))
    basis(2,3) = (cosd(cell(4))-cosd(cell(5))*cosd(cell(6))) &
	/sind(cell(6))
    basis(3,3) = sqrt(sind(cell(5))**2 - basis(2,3)**2)
    basis(:,3) = cell(3) * basis(:,3)

    metric = matmul(transpose(basis),basis)

! The reciprocal basis is basis^-T.
    recbas = 0.0
    recbas(1,1) = 1.0/basis(1,1)
    recbas(2,2) = 1.0/basis(2,2)
    recbas(3,3) = 1.0/basis(3,3)
    recbas(3,2) = - recbas(2,2)*basis(2,3)/basis(3,3)
    recbas(2,1) = - recbas(1,1)*basis(1,2)/basis(2,2)
    recbas(3,1) = - (recbas(1,1)*basis(1,3)+recbas(2,1)*basis(2,3))/basis(3,3)

    recip_metric = matmul(transpose(recbas),recbas)

return
end subroutine cell_metric

  function sg_op_mult(a,b)
    type (sg_operator), intent(in)	::	a,b
    type (sg_operator)		::	sg_op_mult

    sg_op_mult%rot = matmul(a%rot,b%rot)
    sg_op_mult%trans = modulo(matmul(a%rot,b%trans) + a%trans , 24)
    return
  end function sg_op_mult

  logical function sg_op_comp(a,b)
    type (sg_operator), intent(in)	::	a,b
    sg_op_comp = all(a%rot .eq. b%rot) .and. all(a%trans .eq. b%trans)
    return
  end function sg_op_comp

  function grid_op_mult(a,b)
    type (grid_operator), intent(in)	::	a,b
    type (grid_operator)		::	grid_op_mult

    if ( any(a%grid .ne. b%grid ) ) stop ' Grid_op_mult - non-equal grids.'
    grid_op_mult%grid = a%grid
    grid_op_mult%rot = matmul(a%rot,b%rot)
    grid_op_mult%trans = modulo(matmul(a%rot,b%trans) + a%trans , a%grid)
    return
  end function grid_op_mult

  function grid_op_vector(a,b)
    type (grid_operator), intent(in)	::	a
    integer, intent(in)		::	b(3)
    integer			::	grid_op_vector(3)

    grid_op_vector = modulo(matmul(a%rot,b) + a%trans , a%grid)
    return
  end function grid_op_vector

  logical function grid_op_comp(a,b)
    type (grid_operator), intent(in)	::	a, b
    grid_op_comp = all(a%rot .eq. b%rot) .and. all(a%trans .eq. b%trans) &
		.and. all(a%grid .eq. b%grid)
    return
  end function grid_op_comp

  logical function subgroup_of_grid(subg, group)
    type (grid_element), intent(in), target	::	subg, group
! this checks that each element of a is contained in b
    type (grid_element), pointer	:: 	subcur, grocur

  subgroup_of_grid = .FALSE.
  subcur => subg
  do 
    grocur => group
    do
     if ( grocur%op .eq. subcur%op ) exit		! found this element
     if ( .not. associated( grocur%next ) ) return	! not there
     grocur => grocur%next
    end do
    if ( .not. associated( subcur%next ) ) exit		! finished subg
    subcur => subcur%next
  end do
  subgroup_of_grid = .TRUE.
  return

  end function subgroup_of_grid


  logical function within_int(coord,box)
    integer, intent(in)	::	coord(:), box(:,:)
    within_int = all( coord .ge. box(:,1) ) .and. all( coord .le. box(:,2) )
    return
  end function within_int


subroutine gen_to_group( generators, group )
  type (sg_element), pointer :: group, generators, gen, current, last
  type (sg_operator)	:: newop

! this creates the group generated by the specified generators.
! put the unit element into the group
allocate (group)
group%op = unit_op
nullify(group%next)
!call print_sg_op(group%op)

current => group
last => group
!print *,' current op'
!call print_sg_op(current%op)

do		! act with all generators on current element
  gen => generators
  do
    newop = gen%op * current%op
!!call print_mult(newop,gen%op,current%op)
    if( .not. find_op( newop, group ) ) then	! search for newop in group
      allocate( last%next )
      last => last%next
      last%op = newop
      nullify(last%next)
    end if
   if ( .not. associated( gen%next ) ) exit
   gen => gen%next
  end do
!!print *,' Acting on next group element'
  if ( .not. associated( current%next ) ) exit
  current => current%next
end do
return
end subroutine gen_to_group


logical function find_op ( op, group)
  type (sg_element), pointer :: group, current
  type (sg_operator)	:: op

  current => group
  find_op = .FALSE.

  do
   if ( current%op .eq. op ) exit
   if ( .not. associated( current%next ) ) return
   current => current%next
  end do
  find_op = .TRUE.
  return
end function find_op 

integer function num_ops(group)
  type (sg_element), pointer    :: group, current

  num_ops = 0
  current => group
  do
    num_ops = num_ops + 1
    if ( .not. associated(current%next) ) exit
    current => current%next
  end do
  return
end function num_ops


subroutine print_sg(group)
  type (sg_element), pointer	:: group, current

  current => group
  do
    call print_sg_op(current%op)
    if ( .not. associated(current%next) ) exit
    current => current%next
  end do
  return
end subroutine print_sg

subroutine print_sg_op(op)
  type (sg_operator)	:: op
  integer i,j

  print '(3(1x,4i5/))', ( (op%rot(i,j),j=1,3), op%trans(i) , i = 1,3)
  return
end subroutine print_sg_op

subroutine print_mult(newop,op1,op2)
  type (sg_operator)	:: op1,op2,newop
  integer j

print '(1x,4i5,10x,4i5,5x,4i5)',(newop%rot(1,j),j=1,3), newop%trans(1), &
(op1%rot(1,j),j=1,3), op1%trans(1), &
(op2%rot(1,j),j=1,3), op2%trans(1)
print '(1x,4i5,5x,a1,4x,4i5,2x,a1,2x,4i5)',(newop%rot(2,j),j=1,3), newop%trans(2), &
'=',(op1%rot(2,j),j=1,3), op1%trans(2), &
'*',(op2%rot(2,j),j=1,3), op2%trans(2)
print '(1x,4i5,10x,4i5,5x,4i5/)',(newop%rot(3,j),j=1,3), newop%trans(3), &
(op1%rot(3,j),j=1,3), op1%trans(3), &
(op2%rot(3,j),j=1,3), op2%trans(3)

return
end subroutine print_mult

subroutine group_to_grid(ucgrid, group, grid_group)
! this creates the grid operators from the space-group ops in group.

type (sg_element), pointer	::	group
type (grid_element), pointer	::	grid_group
integer, intent(in), target	::	ucgrid(3)

type (sg_element), pointer	::	current
type (grid_element), pointer	::	curr_grid

current => group
allocate(grid_group)
nullify(grid_group%next)
curr_grid => grid_group
do
  curr_grid%op%grid => ucgrid
  curr_grid%op = current%op	! overloaded assignment
  if ( .not. associated(current%next) ) exit
  allocate(curr_grid%next)
  current => current%next
  curr_grid => curr_grid%next
  nullify(curr_grid%next)
end do

return

end subroutine group_to_grid


subroutine print_grid_ops(grid_group)
  type (grid_element), pointer	:: grid_group

  type (grid_element), pointer	:: curr_grid
  integer i,j

curr_grid => grid_group
do
  print '(3(3x,3i5,3x,i5,5x,''|'',i5/))', &
	( (curr_grid%op%rot(i,j),j=1,3), &
		curr_grid%op%trans(i), curr_grid%op%grid(i), i = 1,3)
if ( .not. associated(curr_grid%next) ) exit
curr_grid => curr_grid%next
end do

  return
end subroutine print_grid_ops

end module space_group
