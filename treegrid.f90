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
  use regrid, only: regrid_real_2d

  implicit none
  
  character(len=2000) :: datafile, gridfile, fname, outfile
  character(len=200) :: latname, lonname, varname

  type (varying_string) :: myoptions(1)

  real, allocatable, target :: data(:,:)

  real, allocatable :: newdata(:,:)

  real, allocatable :: newgridx(:,:), newgridy(:,:)
  real, allocatable :: src_grid(:,:,:), dst_grid(:,:,:)
  real, allocatable :: lon(:), lat(:)

  integer :: error, id_data, id_grid, nlat, nlon, nlatnew, nlonnew, ij, nx, ny, i, j

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

  nlat = nc_size(trim(datafile), trim(latname), id_data)
  nlon = nc_size(trim(datafile), trim(lonname), id_data)

  print *,nlon,' x ',nlat

  allocate(lon(nlon),lat(nlat),src_grid(2,nlon,nlat),data(nlon,nlat))

  call nc_read(trim(fname),trim(lonname),lon,ncid=id_data)
  call nc_read(trim(fname),trim(latname),lat,ncid=id_data)

  src_grid(1,:,:) = spread(lon,2,nlat)
  src_grid(2,:,:) = reshape(spread(lat,1,nlon),(/nlon,nlat/))
  
  varname = 'elevation'

  call nc_read(trim(fname),trim(varname),data,ncid=id_data)

  ! Read in data to be re-gridded
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
     allocate(dst_grid(2,nx,ny))
     dst_grid(1,:,:) = newgridx(2::2,2::2)
     dst_grid(2,:,:) = newgridy(2::2,2::2)
  else
     nx = nlonnew
     ny = nlatnew
     dst_grid(1,:,:) = newgridx
     dst_grid(2,:,:) = newgridy
  end if
  allocate(newdata(nx,ny))

  call regrid_real_2d(data, src_grid, dst_grid, newdata)

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
    write(stderr,*) 'Usage: treegrid [--help] data newgrid'
    write(stderr,*)
    write(stderr,*) '  --help    - print this message'
    write(stderr,*)

  end subroutine usage

end program treegrid
