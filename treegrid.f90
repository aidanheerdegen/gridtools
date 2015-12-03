program treegrid

  ! Regrid netcdf data file using kdtree

  use ncio 
  use cmdline_arguments, only: get_options, bad_options, have_args, option_exists, &
      has_value, get_value, assignment(=), next_arg, num_args
  use iso_varying_string
  use string_functions, only: join
  use file_functions, only: exists, freeunit, stderr, stdout, open

  use kdtree2_precision_module
  use kdtree2_module
  use precision

  implicit none
  
  character(len=2000) :: datafile, gridfile, fname, outfile
  character(len=200) :: latname, lonname, varname

  type (varying_string) :: myoptions(1)

  real, allocatable, target :: data(:,:)
  real, pointer :: data1d(:)

  integer, allocatable :: newdata(:,:)

  real, allocatable :: newgridx(:,:), newgridy(:,:), dst_gridx(:,:), dst_gridy(:,:)
  real, allocatable :: lon(:), lat(:)

  real(kdkind), allocatable :: pos_data(:,:), newdist(:,:)
  real(kdkind) :: oldlon, newlon

  integer :: error, id_data, id_grid, nlat, nlon, nlatnew, nlonnew, ij, nx, ny, i, j

  real(kind=rd_kind),parameter :: DEG2RAD = asin(1.0_rd_kind)/90.0_rd_kind  ! PI/180

  logical :: supergrid
  
  type(kdtree2_result) :: results(1)
  type(kdtree2),pointer    :: tree

  ! These are our accepted command line options (see subroutine usage for
  ! an explanation)

  myoptions(1) = 'help'

  ! This call parses the command line arguments for command line options
  call get_options(myoptions, error)

  ! Check we weren't passed duff options -- spit the dummy if we were
  if (error > 0) then
     write(stderr,*) 'ERROR! Unknown options: ',join(bad_options()," ")
     call usage
     STOP
  end if

  ! Check if we just want to print the usage
  if (option_exists('help')) then
     call usage
     STOP
  end if

  if (num_args() < 2) then
     write(stderr,*) 'ERROR! Must supply input data file and new grid file as command-line arguments'
     call usage
     STOP
  end if

  ! Read in data to be re-gridded
  datafile = next_arg()

  print *,'Reading in data file: ',trim(datafile)

  call nc_open(trim(datafile), id_data, writable=.false.)

  latname = 'lat'; lonname = 'lon'

  nlat = nc_size(trim(fname), trim(latname), id_data)
  nlon = nc_size(trim(fname), trim(lonname), id_data)

  print *,nlon,' x ',nlat

  allocate(lon(nlon),lat(nlat),data(nlon,nlat))

  call nc_read(trim(fname),trim(lonname),lon,ncid=id_data)
  call nc_read(trim(fname),trim(latname),lat,ncid=id_data)
  
  varname = 'elevation'

  call nc_read(trim(fname),trim(varname),data,ncid=id_data)

  ! Pointer bounds remapping allows 1d pointer to 2d array
  data1d(1:nlon*nlat) => data

  print *,"minval lon: ",minval(lon)
  print *,"maxval lon: ",maxval(lon)
 
  ! Make all longitudes run between 0 - 360 degrees
  where(lon < 0) lon = lon + 360.

  print *,"minval lon: ",minval(lon)
  print *,"maxval lon: ",maxval(lon)
 
  ! Convert to radians
  lon = lon*DEG2RAD
  lat = lat*DEG2RAD

  ! Points on surface of sphere radius 1
  allocate(pos_data(3,nlon*nlat))
  ! allocate(pos_data(3,100*100))
  ij = 0
  do i = 1, nlat
     ! print *,'y,ij ',i,ij
     do j = 1,  nlon
        ij = ij + 1
        pos_data(1,ij) =  cos(lon(j))*cos(lat(i))
        pos_data(2,ij) =  sin(lon(j))*cos(lat(i))
        pos_data(3,ij) =  sin(lat(i))
        ! if ( (ij > 1700) .and. (ij < 1900)) print *,i,j,lon(j),lat(i),pos_data(:,ij)
     end do
  enddo

  ! Create a nearest-neighbour kdtree from data
  print *,'Creating Tree'
  ! Create tree
  tree => kdtree2_create(pos_data,sort=.false.,rearrange=.true.)
  print *,'Finished Creating Tree'

  ! Read in the new grid
  
  gridfile = next_arg()

  print *,'Reading in new grid from ',trim(gridfile)

  call nc_open(trim(gridfile), id_grid, writable=.false.)

  latname = 'nyp'
  lonname = 'nxp'

  nlatnew = nc_size(trim(gridfile), trim(latname), id_grid)
  nlonnew = nc_size(trim(gridfile), trim(lonname), id_grid)
  
  print *,nlonnew,' x ',nlatnew

  allocate(newgridx(nlonnew,nlatnew))
  allocate(newgridy(nlonnew,nlatnew))

  latname = 'y'
  lonname = 'x'

  call nc_read(trim(gridfile),trim(lonname),newgridx,ncid=id_grid)
  call nc_read(trim(gridfile),trim(latname),newgridy,ncid=id_grid)

  print *,"minval newgridx: ",minval(newgridx)
  print *,"maxval newgridx: ",maxval(newgridx)
 
  ! Need to make this a command line option
  supergrid = .true.
  if (supergrid) then
     ! Supergrid, pull out every second cell
     nx = (nlonnew - 1)/2
     ny = (nlatnew - 1)/2
     allocate(dst_gridx(nx,ny), dst_gridy(nx,ny))
     dst_gridx = newgridx(2::2,2::2)
     dst_gridy = newgridy(2::2,2::2)
  else
     nx = nlonnew
     ny = nlatnew
     dst_gridx = newgridx
     dst_gridy = newgridy
  end if

  print *,"minval dst_gridx: ",minval(dst_gridx)
  print *,"maxval dst_gridx: ",maxval(dst_gridx)
 
  ! Make sure destination grid is between 0 and 360
  where (dst_gridx < 0) dst_gridx = dst_gridx + 360.
 
  print *,"minval dst_gridx: ",minval(dst_gridx)
  print *,"maxval dst_gridx: ",maxval(dst_gridx)
 
  dst_gridx = dst_gridx * DEG2RAD
  dst_gridy = dst_gridy * DEG2RAD

  ! Cycle through the points of the new grid, find n nearest neighbours, and
  ! apply function to these neighbours, and save result at the new grid point

  print *,"Prepare destination grid"

  ! Points on surface of sphere radius 1
  deallocate(pos_data)
  allocate(pos_data(3,nx*ny))
  ! allocate(pos_data(3,100*100))
  ij = 0
  do i = 1, ny
     do j = 1,  nx
        ij = ij + 1
        pos_data(1,ij) =  cos(dst_gridx(j,i))*cos(dst_gridx(j,i))
        pos_data(2,ij) =  sin(dst_gridx(j,i))*cos(dst_gridy(j,i))
        pos_data(3,ij) =  sin(dst_gridy(j,i))
        ! if ( (ij > 1700) .and. (ij < 1900)) print *,i,j,dst_gridx(j,i),dst_gridy(j,i)
     end do
  enddo

  print *,"Traverse tree for all destination points"
  allocate(newdata(nx,ny))
  allocate(newdist(nx,ny))
  ij = 0
  do i = 1, 200 !ny
     print *,i
     do j = 1,  nx
        ij = ij + 1
        ! if ( (ij > 1700) .and. (ij < 1900)) print *,i,j,lon(j),lat(i),pos_data(:,ij)
        call kdtree2_n_nearest(tp=tree,qv=pos_data(:,ij),nn=1,results=results)
        ! if (results(1)%dis > 1.) print *,i,j,dst_gridx(j,i)/DEG2RAD,dst_gridy(j,i)/DEG2RAD,pos_data(:,ij)
        newdata(j,i) = data1d(results(1)%idx) 
        newdist(j,i) = results(1)%dis
        oldlon = lon(mod(results(1)%idx,nlat))/DEG2RAD
        newlon = dst_gridx(j,i)/DEG2RAD
        if (abs(newlon-oldlon) > 1.) print *,i,j,results(1)%idx,mod(results(1)%idx,nlat),oldlon,newlon
     end do
     ! print *,newdist(:,i)
  end do
  
  ! Save result
  print *,"Save result"
  outfile = 'regridded.nc'
  call nc_create(outfile,overwrite=.TRUE.,netcdf4=.TRUE.)
  call nc_write_dim(outfile,"x",x=(/(i,i=1,nx)/))
  call nc_write_dim(outfile,"y",x=(/(i,i=1,ny)/))
  call nc_write(outfile,"newdata",newdata(:,:),dim1="x",dim2="y")

contains

  subroutine usage

    write(stderr,*)
    write(stderr,*) 'Regrid data on new grid using kd-tree'
    write(stderr,*)
    write(stderr,*) 'Usage: treegrid [--help] datafile newgridfile'
    write(stderr,*)
    write(stderr,*) '  --help    - print this message'
    write(stderr,*)

  end subroutine usage

end program treegrid
