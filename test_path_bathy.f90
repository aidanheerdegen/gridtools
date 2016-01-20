program test_path_bathy

  ! use connect
  use pathfind_functions
  use pnm_class, only: pnm_object, write, assignment(=)
  use ncio, only: ncvar, nc_read, nc_open, nc_create, nc_write_dim, nc_write, nc_get_att, &
       nc_size, nc_print_attr, nc_v_init
  use cmdline_arguments, only: get_options, bad_options, have_args, option_exists, &
      has_value, get_value, assignment(=), next_arg, num_args
  use iso_varying_string
  use string_functions, only: join
  use file_functions, only: exists, freeunit, stderr, stdout, open
  use precision

  implicit none

  character(len=2000) :: datafile
  character(len=200) :: latname, lonname, varname, fname

  real, allocatable, target :: data(:,:)

  real, allocatable :: newdata(:,:), dist(:)

  real, allocatable :: src_grid(:,:,:)
  real, allocatable :: lon(:), lat(:)

  type (varying_string) :: myoptions(6)

  integer :: i, j, offset(2)
  integer :: error, id_data, id_grid, nlat, nlon, nlatnew, nlonnew, nx, ny

  type(ncvar) :: v, xgridv, ygridv

  real :: tolerance = 0.01, depth

  logical :: normalise

  ! These are our accepted command line options (see subroutine usage for
  ! an explanation)
  myoptions(1) = 'help'
  myoptions(2) = 'var'
  myoptions(3) = 'lon'
  myoptions(4) = 'lat'
  myoptions(5) = 'depth'
  myoptions(6) = 'norm'

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

  normalise = option_exists('norm')

  if (option_exists('var')) then
     ! We have specified variable name
     if (.NOT. has_value('var')) then
        write(stderr,*) 'Option var must specify a value!'
        call usage
        stop
     end if
     varname = ""
     varname = get_value('var')
  else
     varname = "elevation"
  end if

  if (option_exists('lon')) then
     ! We have specified variable name for the longitude in dest grid
     if (.NOT. has_value('lon')) then
        write(stderr,*) 'Option lon must specify a value!'
        call usage
        stop
     end if
     lonname = ""
     lonname = get_value('lon')
  else
     lonname = "geolon_uv"
  end if

  if (option_exists('lat')) then
     ! We have specified variable name for the longitude in dest grid
     if (.NOT. has_value('lat')) then
        write(stderr,*) 'Option lat must specify a value!'
        call usage
        stop
     end if
     latname = ""
     latname = get_value('lat')
  else
     latname = "geolat_uv"
  end if

  if (option_exists('depth')) then
     ! Specified a depth
     if (.NOT. has_value('depth')) then
        write(stderr,*) 'Option depth must specify a value!'
        call usage
        stop
     end if
     depth = get_value('depth')
  else
     depth = 0.
  end if

  ! Read in data to be re-gridded
  datafile = next_arg()

  print *,'Reading in data file: ',trim(datafile)

  call nc_open(trim(datafile), id_data, writable=.false.)

  ! Initialize the netcdf variable info and load attributes
  call nc_v_init(v,trim(varname))
  call nc_get_att(id_data,v,readmeta=.TRUE.)
  call nc_print_attr(v)

  nlon = v%dlen(1)
  nlat = v%dlen(2)

  print *,nlon,'x',nlat

  allocate(src_grid(2,nlon,nlat),data(nlon,nlat),dist(nlon*nlat))

  ! Read longitude and latitude variables
  call nc_read(trim(fname),trim(lonname),src_grid(1,:,:),ncid=id_data)
  call nc_read(trim(fname),trim(latname),src_grid(2,:,:),ncid=id_data)

  call nc_read(trim(fname),trim(varname),data,ncid=id_data)

  print *,src_grid(1,1,1)
  print *,src_grid(2,1,1)

  dist = -1.

  call path_2d(data<depth, src_grid, dist, normalise=normalise)

  call save('bathy_dist.nc', dist, nlon, nlat)

contains

  subroutine save(outfile, dist, nx, ny)
    character(len=*) :: outfile
    integer :: nx, ny
    real :: dist(:)

    ! Save result
    call nc_create(outfile,overwrite=.TRUE.,netcdf4=.TRUE.)
    call nc_write_dim(outfile,"i",x=(/(i,i=1,nx)/))
    call nc_write_dim(outfile,"j",x=(/(i,i=1,ny)/))
    call nc_write(outfile,"distance",reshape(dist(:),(/nx,ny/)),dim1="i",dim2="j",missing_value=-1.)
  end subroutine save

  subroutine usage

    write(stderr,*)
    write(stderr,*) 'Test path finding on data'
    write(stderr,*)
    write(stderr,*) 'Usage: test_path_bathy [--help] data'
    write(stderr,*)
    write(stderr,*) '  --help    - print this message'
    write(stderr,*)

  end subroutine usage

end program test_path_bathy
