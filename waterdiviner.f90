program waterdiviner

  use cmdline_arguments, only: get_options, bad_options, have_args, option_exists, &
      has_value, get_value, assignment(=), next_arg, num_args
  use iso_varying_string
  use string_functions, only: join
  use file_functions, only: exists, freeunit, stderr, stdout, open
  use ncio 
  use cluster_functions, only: cluster_image


  implicit none

  integer, dimension(:,:), allocatable :: clusters

  integer :: i, error, ncid, nlat, nlon, nclusters

  character(len=2000) :: fname, latname, lonname, varname

  type (varying_string) :: myoptions(4)

  integer, allocatable :: bathydat(:,:)

  logical :: maskwater, changesign
  integer :: maskid

  ! These are our accepted command line options (see subroutine usage for
  ! an explanation)

  myoptions(1) = 'help'
  myoptions(2) = 'outfile'
  myoptions(3) = 'maskid'
  myoptions(4) = 'changesign'


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

  if (option_exists('maskid')) then
     ! Mask out (make into land) all clusterids > value
     if (.NOT. has_value('maskid')) then
        write(stderr,*) 'Option maskid must specify a value!'
        call usage
        stop
     end if
     maskwater = .TRUE.
     maskid = get_value('maskid')
  else
     maskwater = .FALSE.
  end if

  fname = next_arg()

  call nc_open(fname,ncid,writable=.true.)

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


  ! 43200  x        21600   = 933120000 points
  ! minval:       -10977
  ! maxval:         8685
  ! Num clusters found:        17717
  !           1   613311378 = 65.73 %  -- World Oceans
  !           2     1143179 = 0.123 %  -- Caspian
  !           3      757185 = 0.081 %  -- Black Sea
  !           4       41545 = 0.004 %
  !           5       26679 = 0.003 %

  do i = 1, 5
     print *,i,count(clusters==i)
  end do

  if (maskwater) then
     where (clusters > maskid) clusters = 0
     where (clusters < 1 ) bathydat = 0
  end if
  if (changesign) then
     where (clusters > 0) bathydat = -1 * bathydat
  end if
  if (maskwater .or. changesign) then
     call nc_write(trim(fname),trim(varname),bathydat,dim1="lon",dim2="lat", ncid=ncid)
  end if

  call nc_write(fname,"clusterid",clusters,dim1="lon",dim2="lat", ncid=ncid)
  call nc_write_attr(trim(fname), "clusterid", "standard_name", "ocean_cluster_id")
  call nc_write_attr(trim(fname), "clusterid", "long_name", "Ocean cluster id, 0=land, 1=ocean, 2+=inland water bodies")

contains

  subroutine usage

    write(stderr,*)
    write(stderr,*) 'Find contiguous water bodies in GEBCO bathymetry file'
    write(stderr,*)
    write(stderr,*) 'Usage: waterdiviner [--help] bathymetryfile'
    write(stderr,*)
    write(stderr,*) '  --help          - print this message'
    write(stderr,*) '  --maskid=value  - mask (make into land) all clusters with id greater than value'
    write(stderr,*)

  end subroutine usage


end program waterdiviner
