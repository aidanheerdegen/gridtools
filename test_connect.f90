program test_connect

  ! use connect
  use pathfind_functions
  use pnm_class, only: pnm_object, write, assignment(=)
  use ncio, only: nc_create, nc_write_dim, nc_write

  implicit none

  character(len=2000) :: outfile

  integer, parameter :: nx = 4, ny = 10

  integer :: data(nx,ny)
  real    :: grid(2,nx,ny), dist(nx*ny), knowndist(nx*ny)

  type(pnm_object) :: pnm

  integer :: i, j, offset(2)

  real :: tolerance = 0.01

  do i = 1, ny
     do j = 1, nx
        grid(1,j,i) = 5 * j
        grid(2,j,i) = 5 * i
        data(j,i) = i
     end do
     print *,grid(1,:,i)
     print *,grid(2,:,i)
  end do

  ! Hollow out an island in the middle of the data
  data(1:3,2:4) = 0

  pnm = data

  call write(pnm,'test_connect_topo.pgm') 

  ! call make_connect_2d(data>0, grid, dist)
  call path_2d(data>0, grid, dist)

  knowndist =  (/ 0., 0.007552796, 0.01510559, 0.02265838, &
       3.402823e+38, 3.402823e+38, 3.402823e+38, 0.03018265, &
       3.402823e+38, 3.402823e+38, 3.402823e+38, 0.03779325, &
       3.402823e+38, 3.402823e+38, 3.402823e+38, 0.04540385, &
       0.07176837, 0.06551707, 0.05926576, 0.05301446, &
       0.07774892, 0.07204097, 0.06633301, 0.06062507, &
       0.08355603, 0.07844924, 0.07334245, 0.06823567, &
       0.08924451, 0.08477843, 0.08031234, 0.07584626, &
       0.09487278, 0.09106748, 0.08726218, 0.08345687, &
       0.100501, 0.09735651, 0.09421199, 0.09106747 /)

  if (sum(abs(knowndist-dist),mask=pack(data,.TRUE.)>0) > tolerance) then
     print *,'Error in distance calc 2: ', sum(abs(knowndist-dist),mask=pack(data,.TRUE.)>0)
  end if

  offset=(/1,6/)

  ! call make_connect_2d(data>0, grid, dist)
  ! call path_2d(cshift(cshift(data>0,offset(1),1),offset(2),2), cshift(cshift(grid,offset(1),2),offset(2),3), dist)
  call path_2d(data>0, grid, dist, offset)

  knowndist = (/ 0.07783526, 0.07028247, 0.06272967, 0.05517688, &
       3.402823e+38, 3.402823e+38, 3.402823e+38, 0.04756627, &
       3.402823e+38, 3.402823e+38, 3.402823e+38, 0.03995567, &
       3.402823e+38, 3.402823e+38, 3.402823e+38, 0.03234507, &
       0.007610605, 0.01331856, 0.01902651, 0.02473446, &
       0.0, 0.005707953, 0.01141591, 0.01712386, &
       0.0076106, 0.01271739, 0.01782419, 0.02293097, &
       0.0152212, 0.01968729, 0.02415338, 0.02861946, &
       0.02283181, 0.02663711, 0.03044241, 0.03424772, &
       0.03044241, 0.03358693, 0.03673145, 0.03987597 /)

  if (sum(abs(knowndist-dist),mask=pack(data,.TRUE.)>0) > tolerance) then
     print *,'Error in distance calc 2: ', sum(abs(knowndist-dist),mask=pack(data,.TRUE.)>0)
  end if

  ! dist = dist - minval(dist)

  ! pnm = int(dist * (65535./maxval(dist)))

  ! call write(pnm,'test_connect_dist.pgm') 

  print *,minval(dist)
  print *,maxval(dist)

  ! Save result
  print *,"Save result"
  outfile = 'dist.nc'
  call nc_create(outfile,overwrite=.TRUE.,netcdf4=.TRUE.)
  call nc_write_dim(outfile,"i",x=(/(i,i=1,nx)/))
  call nc_write_dim(outfile,"j",x=(/(i,i=1,ny)/))
  call nc_write(outfile,"distance",reshape(dist(:),(/nx,ny/)),dim1="i",dim2="j")

end program test_connect
