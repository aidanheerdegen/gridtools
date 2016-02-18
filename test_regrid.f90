program test_regrid

  use regrid

  implicit none

  real, allocatable :: sdata(:,:), sgrid(:,:,:), ogrid(:,:,:), odata(:,:)
  integer :: nlon, nlat, i, j
  real :: lat, lon, datavalue

  nlon = 10
  nlat = 6

  allocate(sdata(nlon,nlat),sgrid(2,nlon,nlat),ogrid(2,nlon,nlat),odata(nlon,nlat))

  lat = -80.
  datavalue = 0.
  do i = 1, nlat
     lon = -120.
     lat = lat + 20.
     do j = 1, nlon
        lon = lon + 10.
        datavalue = datavalue + 1.
        sgrid(1,j,i) = lon
        sgrid(2,j,i) = lat
        sdata(j,i) = datavalue
     end do
  end do

  ogrid = sgrid + 2.

  do i = 1, nlat
     print *,'nlat=',i
     print '(10(F7.2,X))',sgrid(1,:,i)
     print '(10(F7.2,X))',ogrid(1,:,i)
     print '(10(F7.2,X))',sgrid(2,:,i)
     print '(10(F7.2,X))',ogrid(2,:,i)
  end do
  call regrid_real_2d(sdata, sgrid, ogrid, odata)

  do i = 1, nlat
     print *,i
     print '(10(F0.2,X))',sdata(:,i)
     print '(10(F0.2,X))',odata(:,i)
  end do
  
  ogrid = sgrid + 8.

  do i = 1, nlat
     print *,'nlat=',i
     print '(10(F7.2,X))',sgrid(1,:,i)
     print '(10(F7.2,X))',ogrid(1,:,i)
     print '(10(F7.2,X))',sgrid(2,:,i)
     print '(10(F7.2,X))',ogrid(2,:,i)
  end do
  call regrid_real_2d(sdata, sgrid, ogrid, odata)

  do i = 1, nlat
     print *,i
     print '(10(F0.2,X))',sdata(:,i)
     print '(10(F0.2,X))',odata(:,i)
  end do
  
end program test_regrid
