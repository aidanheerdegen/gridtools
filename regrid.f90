module regrid

  ! Regrid data using kdtree
  use kdtree2_module
  use precision

  implicit none

contains

  subroutine regrid_real_2d(datain, sourcegrid, destgrid, dataout)
  
    real, target, intent(in) :: datain(:,:)
    real, intent(in) :: sourcegrid(:,:,:)
    real, intent(in) :: destgrid(:,:,:)
    real, intent(out) :: dataout(size(destgrid,2),size(destgrid,3))

    real :: gridin(size(sourcegrid,1),size(sourcegrid,2),size(sourcegrid,3))
    real :: gridout(size(destgrid,1),size(destgrid,2),size(destgrid,3))

    real :: data1d(size(datain))

    real(kind=rd_kind),parameter :: DEG2RAD = asin(1.0_rd_kind)/90.0_rd_kind  ! PI/180

    logical :: supergrid
  
    type(kdtree2_result) :: results(1)
    type(kdtree2),pointer    :: tree
    
    integer :: nlonin, nlatin, npointsin
    integer :: nlonout, nlatout, npointsout

    real(kdkind), allocatable :: pos_data(:,:), newdist(:,:)

    integer :: i, j, ij

    nlonin = size(datain,1)
    nlatin = size(datain,2)
    npointsin = nlonin*nlatin

    nlonout = size(destgrid,2)
    nlatout = size(destgrid,3)
    npointsout = nlonout*nlatout

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

    ! Make sure destination grid is between 0 and 360
    where (gridout < 0) gridout = gridout + 360.

    ! Make all longitudes run between 0 - 360 degrees
    where(destgrid < 0)
       gridout = destgrid + 360.
    elsewhere
       gridout = destgrid
    end where

    ! Convert to radians
    gridout = gridout*DEG2RAD

    ! Cycle through the points of the new grid, find n nearest neighbours, and
    ! apply function to these neighbours, and save result at the new grid point

    deallocate(pos_data)
    allocate(pos_data(3,npointsout))
    ij = 0
    do i = 1, nlatout
       do j = 1,  nlonout
          ij = ij + 1
          print '(I4,5(X,F0.2))',ij,pos_data(:,ij),gridin(:,j,i)
          pos_data(1,ij) =  cos(gridout(1,j,i))*cos(gridout(2,j,i))
          pos_data(2,ij) =  sin(gridout(1,j,i))*cos(gridout(2,j,i))
          pos_data(3,ij) =  sin(gridout(2,j,i))
          print '(I4,5(X,F0.2))',ij,pos_data(:,ij),gridout(:,j,i)
       end do
    enddo

    print *,"Traverse tree for all destination points"
    ij = 0
    do i = 1, nlatout
       do j = 1,  nlonout
          ij = ij + 1
          call kdtree2_n_nearest(tp=tree,qv=pos_data(:,ij),nn=1,results=results)
          dataout(j,i) = data1d(results(1)%idx) 
          print *,i,j,ij,results(1)%idx
       end do
    end do

  end subroutine regrid_real_2d
    
end module regrid
  
