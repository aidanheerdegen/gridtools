module pathfind_functions

  implicit none

contains

  subroutine dijkstra_distance ( nv, ohd, mind )

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
    implicit none

    integer ( kind = 4 ) nv

    logical ( kind = 4 ) connected(nv)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) md
    integer ( kind = 4 ) mind(nv)
    integer ( kind = 4 ) mv
    integer ( kind = 4 ) ohd(nv,nv)
    integer ( kind = 4 ) step 
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

    return
  end subroutine dijkstra_distance

  subroutine find_nearest ( nv, mind, connected, d, v )

    !*****************************************************************************80
    !
    !! FIND_NEAREST finds the nearest unconnected node.
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
    !    Input, integer ( kind = 4 ) MIND(NV), the currently computed minimum 
    !    distance from node 1 to each node.
    !
    !    Input, logical ( kind = 4 ) CONNECTED(NV), is true for each connected 
    !    node, whose minimum distance to node 1 has been determined.
    !
    !    Output, integer ( kind = 4 ) D, the distance from node 1 to the nearest 
    !    unconnected node.
    !
    !    Output, integer ( kind = 4 ) V, the index of the nearest unconnected node.
    !
    implicit none

    integer ( kind = 4 ) nv

    logical ( kind = 4 ) connected(nv)
    integer ( kind = 4 ) d
    integer ( kind = 4 ) i
    integer ( kind = 4 ), parameter :: i4_huge = 2147483647
    integer ( kind = 4 ) mind(nv)
    integer ( kind = 4 ) v

    d = i4_huge
    v = -1

    do i = 1, nv
       if ( .not. connected(i) .and. mind(i) < d ) then
          d = mind(i)
          v = i
       end if
    end do

    return
  end subroutine find_nearest

  subroutine init ( nv, ohd )

    !! INIT initializes the problem data.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    01 July 2010
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
    !    Output, integer ( kind = 4 ) OHD(NV,NV), the distance of the direct
    !    link between nodes I and J.
    !
    implicit none

    integer ( kind = 4 ) nv

    integer ( kind = 4 ) i
    integer ( kind = 4 ), parameter :: i4_huge = 2147483647
    integer ( kind = 4 ) ohd(nv,nv)

    ohd(1:nv,1:nv) = i4_huge

    do i = 1, nv
       ohd(i,i) = 0
    end do

    ohd(1,2) = 40
    ohd(1,3) = 15
    ohd(2,3) = 20
    ohd(2,4) = 10
    ohd(2,5) = 25
    ohd(3,4) = 100
    ohd(2,6) = 6
    ohd(5,6) = 8

    ohd(2,1) = 40
    ohd(3,1) = 15
    ohd(3,2) = 20
    ohd(4,2) = 10
    ohd(5,2) = 25
    ohd(4,3) = 100
    ohd(6,2) = 6
    ohd(6,5) = 8

    return
  end subroutine init
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
    implicit none

    integer ( kind = 4 ) nv

    logical ( kind = 4 ) connected(nv)
    integer ( kind = 4 ) i
    integer ( kind = 4 ), parameter :: i4_huge = 2147483647
    integer ( kind = 4 ) mind(nv)
    integer ( kind = 4 ) mv
    integer ( kind = 4 ) ohd(nv,nv)

    do i = 1, nv
       if ( .not. connected(i) ) then
          !
          !  If we really use the maximum integer (or something close) to indicate
          !  no link, then we'll get burned if we add it to another value
          !  Integer arithmetic can 'wrap around', so that 17 + i4_huge becomes
          !  a very negative number!  So first we eliminate the possiblity that
          !  the link is infinite.
          !
          if ( ohd(mv,i) < i4_huge ) then
             if ( mind(mv) + ohd(mv,i) < mind(i) ) then
                mind(i) = mind(mv) + ohd(mv,i)
             end if
          end if
       end if
    end do

    return
  end subroutine update_mind

end module pathfind_functions
