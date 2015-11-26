program waterdiviner

  use pnm_class, only: pnm_object, new, write
  use iso_varying_string
  use cmdline_arguments, only: get_options, bad_options, have_args, option_exists, &
      has_value, get_value, assignment(=), next_arg, num_args
  use string_functions, only: join
  use file_functions, only: exists, freeunit, stderr, stdout, open
  use ncio 
  use cluster_functions, only: cluster_image


  implicit none

  type (pnm_object) :: pnm

  integer, dimension(:,:), allocatable :: data, clusters

  integer :: i, error, ncid, nlat, nlon, nclusters

  character(len=2000) :: fname, outname, latname, lonname, varname, pnmname

  logical :: initialised = .FALSE., fixed_name = .FALSE., have_missing = .FALSE.

  type (varying_string) :: myoptions(3)

  integer, allocatable :: bathydat(:,:)

  ! These are our accepted command line options (see subroutine usage for
  ! an explanation)

  myoptions(1) = 'help'
  myoptions(2) = 'outfile'
  myoptions(3) = 'missing'

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

  if (num_args() < 1) then
     write(stderr,*) 'ERROR! Must supply GEBCO bathymetry file as a command-line argument'
     call usage
     STOP
  end if

  fname = next_arg()

  call nc_open(fname,ncid,writable=.false.)

  varname = 'elevation'
  latname = 'lat'
  lonname = 'lon'

  nlat = nc_size(trim(fname), trim(latname), ncid)
  nlon = nc_size(trim(fname), trim(lonname), ncid)
  
  print *,nlon,' x ',nlat

  allocate(bathydat(nlon,nlat))

  call nc_read(fname,varname,bathydat,ncid=ncid)

  print *,"minval: ",minval(bathydat)
  print *,"maxval: ",maxval(bathydat)
 
  allocate(clusters(nlon,nlat))

  ! Return integer array where all the pixels which are 'true' in the input
  ! are labelled by cluster number. In this case all the true pixels
  ! will be those which exceed our mask value
  clusters = cluster_image( bathydat < 0, sort=.true.) 

  ! The clusters are numbered from 1 .. n, so the maximum value will
  ! be the number of different clusters
  nclusters = maxval(clusters)

  print *,'Num clusters found: ',nclusters

  do i = 1, 5
     print *,i,count(clusters==i)
  end do

  call new(pnm, clusters, origin='bl')
  pnmname = 'clusters.pgm'
  call write(pnm, trim(pnmname))

contains

  subroutine usage

    write(stderr,*)
    write(stderr,*) 'Find contiguous water bodies in GEBCO bathymetry file'
    write(stderr,*)
    write(stderr,*) 'Usage: waterdiviner [--help] bathymetryfile'
    write(stderr,*)
    write(stderr,*) '  --help    - print this message'
    write(stderr,*)

  end subroutine usage


end program waterdiviner
