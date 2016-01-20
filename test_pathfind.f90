program main

  !*****************************************************************************80
  !
  !! MAIN runs an example of Dijkstra's minimum distance algorithm.
  !
  !  Discussion:
  !
  !    Given the distance matrix that defines a graph, we seek a list
  !    of the minimum distances between node 0 and all other nodes.
  !
  !    This program sets up a small example problem and solves it.
  !
  !    The correct minimum distances are:
  !
  !      0   35   15   45   49   41
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    30 June 2010
  !
  !  Author:
  !
  !    Original C version by Norm Matloff, CS Dept, UC Davis.
  !    FORTRAN90 version by John Burkardt.
  !

  use pathfind_functions, only: dijkstra_distance, real_huge 

  implicit none

  integer, parameter :: nv = 6

  real :: mind(nv)
  real :: ohd(nv,nv)

  !  Initialize the problem data.
  call init ( nv, ohd )

  ! Call program
  call dijkstra_distance ( nv, ohd, mind )

  if (any(mind /= (/0., 35., 15., 45., 49., 41./))) then
     print *,"FAILED!"
     print *,mind
  else
     print *,"SUCCESS!"
  end if

contains

  subroutine init ( nv, ohd )

    !*****************************************************************************80
    !
    !! INIT initializes the problem data.
    !
    !  Discussion:
    !
    !    The graph uses 6 nodes, and has the following diagram and
    !    distance matrix:
    !
    !    N0--15--N2-100--N3           0   40   15  Inf  Inf  Inf
    !      \      |     /            40    0   20   10   25    6
    !       \     |    /             15   20    0  100  Inf  Inf
    !        40  20  10             Inf   10  100    0  Inf  Inf
    !          \  |  /              Inf   25  Inf  Inf    0    8
    !           \ | /               Inf    6  Inf  Inf    8    0
    !            N1
    !            / \
    !           /   \
    !          6    25
    !         /       \
    !        /         \
    !      N5----8-----N4
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

    integer :: nv

    integer :: i
    real :: ohd(nv,nv)

    ohd(1:nv,1:nv) = real_huge

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
end program main

