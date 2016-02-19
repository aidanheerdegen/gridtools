module pathfind_functions

  use kdtree2_module
  use precision

  implicit none

  integer, parameter :: i4_huge = 2147483647
  real, parameter    :: real_huge = huge(1.0)
  real, parameter    :: earth_R = 6371000 ! units is metres
  real(kind=rd_kind),parameter :: DEG2RAD = asin(1.0_rd_kind)/90.0_rd_kind  ! PI/180

  interface make_tree
     module procedure make_tree_vectors, make_tree_latlon
  end interface make_tree

contains

  function lonlat2vector(lonlat) result(v)

    real, intent(in) :: lonlat(2)
    real :: v(3)

    real :: lonlatrad(2)

    lonlatrad = lonlat * DEG2RAD

    v(1) =  cos(lonlatrad(1))*cos(lonlatrad(2))
    v(2) =  sin(lonlatrad(1))*cos(lonlatrad(2))
    v(3) =  sin(lonlatrad(2))

    v = v * earth_R
    
  end function lonlat2vector
       
  subroutine make_tree_vectors(invectors, tree) 

    real(kdkind), intent(inout)         :: invectors(:,:)
    type(kdtree2), pointer, intent(out) :: tree
    
    ! Create a nearest-neighbour kdtree from data
    tree => kdtree2_create(invectors,sort=.false.,rearrange=.true.)

  end subroutine make_tree_vectors

  subroutine make_tree_latlon(sourcegrid, tree, vectors) 

    ! Calculate and return the vectors used in the tree. Must do this, as
    ! the tree only points to the data, so it must persist while the tree
    ! persists

    real, intent(in)                      :: sourcegrid(:,:,:)
    type(kdtree2), pointer, intent(inout) :: tree
    real(kdkind), intent(out)             :: vectors(:,:)

    integer :: i, j, ij, nlonin, nlatin
    
    nlonin = size(sourcegrid,2)
    nlatin = size(sourcegrid,3)

    ! Points on surface of spherical earth
    ij = 0
    do i = 1, nlatin
       do j = 1,  nlonin
          ij = ij + 1
          vectors(:,ij) = lonlat2vector(sourcegrid(:,j,i))
       end do
    enddo

    ! Create a nearest-neighbour kdtree from data
    ! call make_tree(pos_data, tree)
    tree => kdtree2_create(vectors,sort=.false.,rearrange=.true.)

    print *,tree%n

  end subroutine make_tree_latlon

  subroutine path_2d_tree(tree, mind, mask, origin) 

    type(kdtree2), pointer, intent(in) :: tree
    real, intent(out)             :: mind(:)
    logical, intent(in), optional :: mask(:)
    real, intent(in), optional    :: origin(2)

    real :: ohd(size(mind), size(mind))
    logical :: vectormask(size(mind))

    ! Assume a basic grid layout, so 4 possible neighbours + 1 for
    ! the point where query is made
    integer, parameter :: nneighbours = 5

    type(kdtree2_result)  :: results(nneighbours), oresult(1)
    
    real(kdkind) :: origv(3)

    integer :: npointsin, origin1d

    integer :: i, j, ij, idx

    if (present(mask)) then
       vectormask = mask
    else
       vectormask = .true.
    end if

    npointsin = size(mind)

    if (present(origin)) then
       origv = lonlat2vector(origin)
       call kdtree2_n_nearest(tree, origv, 1, oresult)
       origin1d = oresult(1)%idx
       print *,'origin1d: ',origin1d
    end if
 
    print *,"Traverse tree for all destination points"
    ij = 0

    ohd = real_huge
    do i = 1, npointsin

       ! idx = index1d(i)
       idx = i
       call kdtree2_n_nearest_around_point(tp=tree,idxin=idx,correltime=0,nn=nneighbours,results=results)
       do j = 1, nneighbours
          if (results(j)%idx == 0) then
             ! print *,"i: ",i, "idx: ",idx,"resultsidx: ",results(j)%idx
             cycle
          end if
          ! Ignore the point itself
          if (results(j)%idx == idx) then
             ohd(i,results(j)%idx) = 0.
          else if (vectormask(results(j)%idx)) then
             ohd(i,results(j)%idx) = results(j)%dis
          else
             ohd(i,results(j)%idx) = real_huge
          end if
          ! print *,'j:',j,results(j)%idx, results(j)%dis, vectormask(results(j)%idx), ohd(i,results(j)%idx)
       end do
       ! print *,idx,j,ij,results(1)%idx
    end do

    if (present(origin)) then
       ohd = cshift(cshift(ohd, origin1d, dim=1), origin1d, dim=2)
    end if

    print *,'Run dijkstras algorithm'

    call dijkstra_distance(npointsin, ohd, mind)

    print *,minval(mind)
    print *,maxval(mind)

    if (present(origin)) then
       print *,'shift mind back'
       mind = cshift(mind, -1*origin1d)
    end if

  end subroutine path_2d_tree

  subroutine path_2d_vectors(invectors, mind, mask, origin) 

    real(kdkind), intent(inout)   :: invectors(:,:)
    real, intent(out)             :: mind(size(invectors,2))
    logical, intent(in), optional :: mask(:)
    real, intent(in), optional    :: origin(2)

    real :: ohd(size(mind), size(mind))
    logical :: vectormask(size(mind))

    ! Assume a basic grid layout, so 4 possible neighbours + 1 for
    ! the point where query is made
    integer, parameter :: nneighbours = 5

    type(kdtree2_result)  :: results(nneighbours), oresult(1)
    type(kdtree2),pointer :: tree
    
    real(kdkind) :: origv(3)

    integer :: npointsin, origin1d

    integer :: i, j, ij, idx

    if (present(mask)) then
       vectormask = mask
    else
       vectormask = .true.
    end if

    npointsin = size(mind)

    print *,'Create kd tree'

    ! Create a nearest-neighbour kdtree from data
    tree => kdtree2_create(invectors,sort=.false.,rearrange=.true.)

    if (present(origin)) then
       origv = lonlat2vector(origin)
       call kdtree2_n_nearest(tree, origv, 1, oresult)
       origin1d = oresult(1)%idx
       print *,'origin1d: ',origin1d
    end if
 
    print *,"Traverse tree for all destination points"
    ij = 0

    ohd = real_huge
    do i = 1, npointsin

       ! idx = index1d(i)
       idx = i
       ! print *,"i: ",i, "idx: ",idx
       call kdtree2_n_nearest_around_point(tp=tree,idxin=idx,correltime=0,nn=nneighbours,results=results)
       do j = 1, nneighbours
          ! Ignore the point itself
          if (results(j)%idx == idx) then
             ohd(i,results(j)%idx) = 0.
          else if (vectormask(results(j)%idx)) then
             ohd(i,results(j)%idx) = results(j)%dis
          else
             ohd(i,results(j)%idx) = real_huge
          end if
          ! print *,'j:',j,results(j)%idx, results(j)%dis, vectormask(results(j)%idx), ohd(i,results(j)%idx)
       end do
       ! print *,idx,j,ij,results(1)%idx
    end do

    if (present(origin)) then
       ohd = cshift(cshift(ohd, origin1d, dim=1), origin1d, dim=2)
    end if

    print *,'Run dijkstras algorithm'

    call dijkstra_distance(npointsin, ohd, mind)

    print *,minval(mind)
    print *,maxval(mind)

    if (present(origin)) then
       print *,'shift mind back'
       mind = cshift(mind, -1*origin1d)
    end if

  end subroutine path_2d_vectors

  subroutine path_2d_distance(datain, sourcegrid, mind, origin, normalise)
  
    logical, target, intent(in)   :: datain(:,:)
    real, intent(in)              :: sourcegrid(:,:,:)
    real, intent(out)             :: mind(size(datain))

    ! Optionally specify the origin from which to calculate distance
    integer, optional, intent(in) :: origin(2)

    ! Optionally specify normalised distance to be returned
    logical, optional, intent(in) :: normalise

    real :: mindtmp(size(datain))
    real :: gridin(size(sourcegrid,1),size(sourcegrid,2),size(sourcegrid,3))
    logical :: data1d(size(datain)), normalise_mind

    real(kind=rd_kind),parameter :: DEG2RAD = asin(1.0_rd_kind)/90.0_rd_kind  ! PI/180

    ! Assume a basic grid layout, so 4 possible neighbours + 1 for
    ! the point where query is made
    integer, parameter :: nneighbours = 5

    integer :: nlonin, nlatin, npointsin, npointstrue, origin1d

    real(kdkind), allocatable :: pos_data(:,:), pos_data_true(:,:)

    integer :: i, j, ij

    if (present(normalise)) then
       normalise_mind = normalise
    else
       normalise_mind = .false.
    end if

    nlonin = size(datain,1)
    nlatin = size(datain,2)

    npointsin = nlonin*nlatin
    npointstrue = count(datain)

    data1d = reshape(datain,shape(data1d))

    ! Make all longitudes run between 0 - 360 degrees
    where(sourcegrid < 0)
       gridin = sourcegrid + 360.
    elsewhere
       gridin = sourcegrid
    end where

    ! Convert to radians
    gridin = gridin*DEG2RAD

    ! Check for origin, otherwise set to first cell in input data
    if (present(origin)) then
       ! Cycle the data and grid to have origin first. Subtract 1 from
       ! each of the shifts to convert between location and offset
       gridin = cshift(cshift(gridin,origin(1)-1,2),origin(2)-1,3)

       ! This is a 1-D offset, so only subtract 1 once (make sense?)
       origin1d = origin(1) + (origin(2)-1)*nlonin - 1

       print *,'origin1d: ',origin1d
       data1d = cshift(data1d,origin1d)
    end if

    ! Points on surface of sphere radius 1
    allocate(pos_data(3,npointsin))
    ij = 0
    do i = 1, nlatin
       do j = 1,  nlonin
          ij = ij + 1
          pos_data(1,ij) =  cos(gridin(1,j,i))*cos(gridin(2,j,i))
          pos_data(2,ij) =  sin(gridin(1,j,i))*cos(gridin(2,j,i))
          pos_data(3,ij) =  sin(gridin(2,j,i))
       end do
    enddo

    allocate(pos_data_true(3,npointstrue))

    ! Make distances relative to the origin
    forall (i=1:npointsin) mind(i) = sum((pos_data(:,i)-pos_data(:,1))**2,dim=1)

    ! Find all minimum distances without data marked as false in datain
    call path_2d_vectors(pos_data, mind, data1d)
    
    if (normalise_mind) then

       ! Re-run the distance calculations, including masked nodes
       call path_2d_vectors(pos_data, mindtmp)

       where (mind == real_huge .or. mindtmp == real_huge)
          mind = 0
          mindtmp = 0
       end where

       where (mindtmp > 0. .and. data1d)
          mind = mind / mindtmp
       elsewhere
          mind = -1
       end where

    end if

    print *,minval(mind)
    print *,maxval(mind)

    ! Check for origin and shift mind back to original position
    if (present(origin)) then
       mind = cshift(mind,-1*origin1d)
    end if

  end subroutine path_2d_distance

  subroutine path_2d(datamask, sourcegrid, mind, origin, normalise)

    ! Calculate the path from the origin to all other points in sourcegrid.
    ! Cells where datamask is false are omitted from the distance calculation.
    ! For the example of calculating paths in the ocean, the masked cells
    ! correspond to land
  
    logical, target, intent(in)   :: datamask(:,:)
    real, intent(in)              :: sourcegrid(:,:,:)
    real, intent(out)             :: mind(size(datamask))

    ! Optionally specify the origin from which to calculate distance
    real, intent(inout) :: origin(2)

    ! Optionally specify normalised distance to be returned
    logical, optional, intent(in) :: normalise

    logical :: data1d(size(datamask)), normalise_mind

    integer :: npointsin

    type(kdtree2), pointer :: tree
    real(kdkind), allocatable :: pos_data(:,:)

    real, allocatable :: mindtmp(:)

    if (present(normalise)) then
       normalise_mind = normalise
    else
       normalise_mind = .false.
    end if

    npointsin = size(datamask)

    data1d = reshape(datamask,shape(data1d))

    allocate(pos_data(3,npointsin))
    call make_tree(sourcegrid, tree, pos_data)

    ! Find all minimum distances without data marked as false in datamask
    call path_2d_tree(tree, mind, data1d, origin)
    
    if (normalise_mind) then

       allocate(mindtmp(npointsin))

       ! Re-run the distance calculations, including masked nodes
       ! call path_2d_vectors(pos_data, mindtmp, origin=origin)
       call path_2d_tree(tree, mindtmp, origin=origin)

       where (mind == real_huge .or. mindtmp == real_huge)
          mind = 0
          mindtmp = 0
       end where

       where (mind > 0. .and. data1d)
          mind = mindtmp / mind
       elsewhere
          mind = -1
       end where

       deallocate(mindtmp)

    end if

    deallocate(pos_data)

  end subroutine path_2d

  subroutine path_2d_old(datain, sourcegrid, mind, origin, normalise)
  
    logical, target, intent(in)   :: datain(:,:)
    real, intent(in)              :: sourcegrid(:,:,:)
    real, intent(out)             :: mind(size(datain))

    ! Optionally specify the origin from which to calculate distance
    real, intent(inout) :: origin(2)

    ! Optionally specify normalised distance to be returned
    logical, optional, intent(in) :: normalise

    real :: ohd(size(datain),size(datain))
    real :: mindtmp(size(datain))
    real :: gridin(size(sourcegrid,1),size(sourcegrid,2),size(sourcegrid,3))
    logical :: data1d(size(datain)), normalise_mind

    real(kind=rd_kind),parameter :: DEG2RAD = asin(1.0_rd_kind)/90.0_rd_kind  ! PI/180

    ! Assume a basic grid layout, so 4 possible neighbours + 1 for
    ! the point where query is made
    integer, parameter :: nneighbours = 5

    type(kdtree2_result)  :: results(nneighbours), oresult(1)
    type(kdtree2),pointer :: tree
    
    integer :: nlonin, nlatin, npointsin, origin1d

    real(kdkind), allocatable :: pos_data(:,:)
    real(kdkind) :: origv(3)

    integer :: i, j, ij

    if (present(normalise)) then
       normalise_mind = normalise
    else
       normalise_mind = .false.
    end if

    nlonin = size(datain,1)
    nlatin = size(datain,2)

    npointsin = nlonin*nlatin

    data1d = reshape(datain,shape(data1d))

    ! Make all longitudes run between 0 - 360 degrees
    where(sourcegrid < 0)
       gridin = sourcegrid + 360.
    elsewhere
       gridin = sourcegrid
    end where

    ! Convert to radians
    gridin = gridin*DEG2RAD
    origin = origin*DEG2RAD

    ! Points on surface of sphere radius 1
    allocate(pos_data(3,npointsin))
    ij = 0
    do i = 1, nlatin
       do j = 1,  nlonin
          ij = ij + 1
          pos_data(1,ij) =  cos(gridin(1,j,i))*cos(gridin(2,j,i))
          pos_data(2,ij) =  sin(gridin(1,j,i))*cos(gridin(2,j,i))
          pos_data(3,ij) =  sin(gridin(2,j,i))
       end do
    enddo

    ! Make distances relative to the origin
    ! forall (i=1:npointsin) mind(i) = sum((pos_data(:,i)-pos_data(:,1))**2,dim=1)
    
    ! print *,i,j,'mind = ',mind(ij)

    ! print *,pos_data(1,:)
    ! print *,pos_data(2,:)
    ! print *,pos_data(3,:)

    print *,'Create kd tree'

    ! Create a nearest-neighbour kdtree from data
    tree => kdtree2_create(pos_data,sort=.false.,rearrange=.true.)

    origv(1) =  cos(origin(1))*cos(origin(2))
    origv(2) =  sin(origin(1))*cos(origin(2))
    origv(3) =  sin(origin(2))
    
    call kdtree2_n_nearest(tree, origv, 1, oresult)

    origin1d = oresult(1)%idx

    ! Cycle the data and grid to have origin first. Subtract 1 from
    ! each of the shifts to convert between location and offset
    ! gridin = cshift(cshift(gridin,origin(1)-1,2),origin(2)-1,3)
    
    ! This is a 1-D offset, so only subtract 1 once (make sense?)
    ! origin1d = origin(1) + (origin(2)-1)*nlonin - 1

    print *,'origin1d: ',origin1d
    data1d = cshift(data1d,origin1d)
 
    print *,"Traverse tree for all destination points"
    ij = 0

    ohd = real_huge
    do i = 1, npointsin
       ! print *,"i: ",i
       if (.not. data1d(i)) cycle
       call kdtree2_n_nearest_around_point(tp=tree,idxin=i,correltime=0,nn=nneighbours,results=results)
       do j = 1, nneighbours
          ! Ignore the point itself
          if (results(j)%idx == i) then
             ohd(i,results(j)%idx) = 0.
          else if (data1d(results(j)%idx)) then
             ohd(i,results(j)%idx) = results(j)%dis
          else
             ohd(i,results(j)%idx) = real_huge
          end if
          ! print *,'j:',j,results(j)%idx, results(j)%dis, data1d(results(j)%idx)
       end do
       ! print *,i,j,ij,results(1)%idx
    end do

    print *,'Run dijkstras algorithm'

    print *,minval(mind)
    print *,maxval(mind)

    call dijkstra_distance(npointsin, ohd, mind)

    if (normalise_mind) then

       ! Re-run the distance calculations for each node, but this time don't
       ! exclude masked nodes
       do i = 1, npointsin
          call kdtree2_n_nearest_around_point(tp=tree,idxin=i,correltime=0,nn=nneighbours,results=results)
          do j = 1, nneighbours
             ! Ignore the point itself
             if (results(j)%idx == i) then
                ohd(i,results(j)%idx) = 0.
             else
                ohd(i,results(j)%idx) = results(j)%dis
             end if
          end do
       end do

       mindtmp = mind
       
       call dijkstra_distance(npointsin, ohd, mind)

       where (mind == real_huge .or. mindtmp == real_huge)
          mind = 0
          mindtmp = 0
       end where

       where (mind > 0. .and. data1d)
          mind = mindtmp / mind
       elsewhere
          mind = -1
       end where

    end if

    print *,minval(mind)
    print *,maxval(mind)

    ! Shift mind back to original position
    mind = cshift(mind,-1*origin1d)

  end subroutine path_2d_old

  subroutine dijkstra_distance(nv, ohd, mind)

    ! DIJKSTRA_DISTANCE uses Dijkstra's minimum distance algorithm.
    !
    !  Discussion:
    !
    !    We essentially build a tree.  We start with only node 0 connected
    !    to the tree, and this is indicated by setting CONNECTED(0) = TRUE.
    !
    !    We initialize MIND(I) to the one step distance from node 0 to node I.
    !    
    !    Now we search among the unconnected nodes for the node MV whose minimum
    !    distance is smallest, and connect it to the tree.  For each remaining
    !    unconnected node I, we check to see whether the distance from 0 to MV
    !    to I is less than that recorded in MIND(I), and if so, we can reduce
    !    the distance.
    !
    !    After NV-1 steps, we have connected all the nodes to 0, and computed
    !    the correct minimum distances.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Author:
    !
    !    Original C version by Norm Matloff, CS Dept, UC Davis.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) NV, the number of nodes.
    !
    !    Input, integer ( kind = 4 ) OHD(NV,NV), the distance of the direct
    !    link between nodes I and J.
    !
    !    Output, integer ( kind = 4 ) MIND(NV), the minimum 
    !    distance from node 1 to each node.
    !

    integer :: nv

    logical :: connected(nv)
    real    :: md
    real    :: mind(nv)
    integer :: mv
    real    :: ohd(nv,nv)
    integer :: step 
    !
    !  Start out with only node 1 connected to the tree.
    !
    connected(1) = .true.
    connected(2:nv) = .false.
    !
    !  Initialize the minimum distance to the one-step distance.
    !
    mind(1:nv) = ohd(1,1:nv)
    !
    !  Attach one more node on each iteration.
    !
    do step = 2, nv
       !
       !  Find the nearest unconnected node.
       !
       call find_nearest ( nv, mind, connected, md, mv )

       if ( mv == - 1 ) then
          write ( *, '(a)' )  ' '
          write ( *, '(a)' )  'DIJKSTRA_DISTANCE - Warning!'
          write ( *, '(a)' )  '  Search terminated early.'
          write ( *, '(a)' )  '  Graph might not be connected.'
          return
       end if
       !
       !  Mark this node as connected.
       !
       connected(mv) = .true.
       !
       !  Having determined the minimum distance to node MV, see if
       !  that reduces the minimum distance to other nodes.
       !
       call update_mind ( nv, connected, ohd, mv, mind )

    end do

    ! print *,mind

    return
  end subroutine dijkstra_distance

  subroutine find_nearest(nv, mind, connected, d, v)

    ! FIND_NEAREST finds the nearest unconnected node.
    ! MIND the currently computed minimum distance from node 1 to each node.
    ! CONNECTED is true for each connected node, whose minimum distance to
    !           node 1 has been determined.
    ! d the distance from node 1 to the nearest unconnected node.
    ! v the index of the nearest unconnected node.

    integer                 :: nv
    real, intent(inout)     :: mind(:)
    logical, intent(inout)  :: connected(:)
    real, intent(out)       :: d
    integer, intent(out)    :: v

    integer :: i

    d = real_huge
    v = -1

    do i = 1, nv
       if ( .not. connected(i) .and. mind(i) < d ) then
          d = mind(i)
          v = i
       end if
    end do

    return
  end subroutine find_nearest

  subroutine update_mind ( nv, connected, ohd, mv, mind )

    !*****************************************************************************80
    !
    !! UPDATE_MIND updates the minimum distance vector.
    !
    !  Discussion:
    !
    !    We've just determined the minimum distance to node MV.
    !
    !    For each node I which is not connected yet,
    !    check whether the route from node 0 to MV to I is shorter
    !    than the currently known minimum distance.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Author:
    !
    !    Original C version by Norm Matloff, CS Dept, UC Davis.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) NV, the number of nodes.
    !
    !    Input, logical ( kind = 4 ) CONNECTED(NV), is true for each connected 
    !    node, whose minimum distance to node 0 has been determined.
    !
    !    Input, integer ( kind = 4 ) OHD(NV,NV), the distance of the direct link 
    !    between nodes I and J.
    !
    !    Input, integer ( kind = 4 ) MV, the node whose minimum distance to node 20
    !    has just been determined.
    !
    !    Input/output, integer ( kind = 4 ) MIND(NV), the currently computed
    !    minimum distances from node 1 to each node.
    !

    integer :: nv

    logical :: connected(nv)
    integer :: i
    real    :: mind(nv)
    integer :: mv
    real    :: ohd(nv,nv)

    do i = 1, nv
       if ( .not. connected(i) ) then
          !
          !  If we really use the maximum integer (or something close) to indicate
          !  no link, then we'll get burned if we add it to another value
          !  Integer arithmetic can 'wrap around', so that 17 + real_huge becomes
          !  a very negative number!  So first we eliminate the possiblity that
          !  the link is infinite.
          !
          if ( ohd(mv,i) < real_huge ) then
             if ( mind(mv) + ohd(mv,i) < mind(i) ) then
                ! print *,'Updating ',i,mind(i)
                mind(i) = mind(mv) + ohd(mv,i)
                ! print *,'---->    ',mind(i)
             end if
          end if
       end if
    end do

    return
  end subroutine update_mind

end module pathfind_functions
