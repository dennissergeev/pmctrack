program read_netcdf
  use netcdf
  implicit none
  
  ! This is the name of the data file we will read.
  character (len = *), parameter :: FILE_NAME = "../../../phd/reanalysis/ERA5/era5.an.sfc.2011.01.vo.nc"
  integer :: ncid

  character (len = nf90_max_name), parameter :: LVL_NAME = "level"
  character (len = nf90_max_name), parameter :: LAT_NAME = "latitude"
  character (len = nf90_max_name), parameter :: LON_NAME = "longitude"
  character (len = nf90_max_name), parameter :: REC_NAME = "time"
  character (len = nf90_max_name), parameter :: VAR_NAME = "vo"
  integer :: lvl_dimid, lon_dimid, lat_dimid, rec_dimid
  integer :: ntime, nlevs, nlats, nlons

  real(4), allocatable :: vort(:, :, :, :) 
  real(4), allocatable :: lats(:) 
  real(4), allocatable :: lons(:) 

  integer :: lon_varid, lat_varid
  ! Loop indices
  integer :: lvl, lat, lon, t, i

  ! Open the file. 
  call check( nf90_open(FILE_NAME, nf90_nowrite, ncid) )

  ! Get the varids of the latitude and longitude coordinate variables.
  call check( nf90_inq_dimid(ncid, REC_NAME, rec_dimid) )
  call check( nf90_inq_dimid(ncid, LVL_NAME, lvl_dimid) )
  call check( nf90_inq_dimid(ncid, LAT_NAME, lat_dimid) )
  call check( nf90_inq_dimid(ncid, LON_NAME, lon_dimid) )

  call check( nf90_inquire_dimension(ncid, rec_dimid, len=ntime) )
  call check( nf90_inquire_dimension(ncid, lvl_dimid, len=nlevs) )
  call check( nf90_inquire_dimension(ncid, lat_dimid, len=nlats) )
  call check( nf90_inquire_dimension(ncid, lon_dimid, len=nlons) )

  print*, ntime, nlevs, nlats, nlons

  allocate(lats(nlats))
  allocate(lons(nlons))
  allocate(vort(ntime, nlevs, nlats, nlons))
  ! Read the latitude and longitude data.
  call check( nf90_inq_varid(ncid, LAT_NAME, lat_varid) )
  call check( nf90_inq_varid(ncid, LON_NAME, lon_varid) )
  call check( nf90_get_var(ncid, lat_varid, lats) )
  call check( nf90_get_var(ncid, lon_varid, lons) )
  
  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file.
  call check( nf90_close(ncid) )

  ! If we got this far, everything worked as expected. Yipee! 
  print *,"end"

contains
  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print*, 'status:', status
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  
end program read_netcdf
