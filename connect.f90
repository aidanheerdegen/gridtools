module connect

  ! Regrid data using kdtree
  use kdtree2_module
  use precision

  implicit none

contains

  subroutine make_connect_2d(datain, sourcegrid, dataout)
  
    logical, target, intent(in) :: datain(:,:)
    real, intent(in)            :: sourcegrid(:,:,:)
    real, intent(out)           :: dataout(size(datain),size(datain))

    real :: gridin(size(sourcegrid,1),size(sourcegrid,2),size(sourcegrid,3))
    logical :: data1d(size(datain))

    real(kind=rd_kind),parameter :: DEG2RAD = asin(1.0_rd_kind)/90.0_rd_kind  ! PI/180

    ! Assume a basic grid layout, so 8 possible neighbours
    integer, parameter :: nneighbours = 8

    type(kdtree2_result)  :: results(nneighbours)
    type(kdtree2),pointer :: tree
    
    integer :: nlonin, nlatin, npointsin

    real(kdkind), allocatable :: pos_data(:,:)

    integer :: i, j, ij

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

    ! Create a nearest-neighbour kdtree from data
    tree => kdtree2_create(pos_data,sort=.false.,rearrange=.true.)

    print *,"Traverse tree for all destination points"
    ij = 0
    do i = 1, npointsin
       print *,i
       call kdtree2_n_nearest(tp=tree,qv=pos_data(:,i),nn=nneighbours,results=results)
       do j = 1, nneighbours
          print *,'j:',j,results(j)%idx, results(j)%dis, data1d(results(i)%idx)
          if (data1d(results(i)%idx)) then
             dataout(j,results(j)%idx) = results(j)%dis
          else
             dataout(j,results(j)%idx) = -1
          end if
       end do
       ! print *,i,j,ij,results(1)%idx
    end do

  end subroutine make_connect_2d
    
end module connect
  
