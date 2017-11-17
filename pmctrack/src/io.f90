program read_netcdf
  use netcdf
  implicit none
  
  integer, parameter :: wp = kind(0.0)

  ! This is the name of the data file we will read.
  character (len = *), parameter :: FILE_NAME = "../../../phd/reanalysis/ERA5/era5.an.sfc.2011.01.vo.nc"
  integer :: ncid

  character (len = nf90_max_name), parameter :: LVL_NAME = "level"
  character (len = nf90_max_name), parameter :: LAT_NAME = "latitude"
  character (len = nf90_max_name), parameter :: LON_NAME = "longitude"
  character (len = nf90_max_name), parameter :: REC_NAME = "time"
  character (len = nf90_max_name), parameter :: VORT_NAME = "vo"
  integer :: lvl_dimid, lon_dimid, lat_dimid, rec_dimid
  integer :: ntime, nlevs, nlats, nlons
  integer :: xtype, ndims

  real(wp), allocatable :: vort(:, :, :, :) 
  integer, allocatable :: time(:) 
  integer, allocatable :: levs(:) 
  real(wp), allocatable :: lats(:) 
  real(wp), allocatable :: lons(:) 

  integer :: rec_varid, lvl_varid, lat_varid, lon_varid
  integer :: vort_varid
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

  allocate(time(ntime))
  allocate(levs(nlevs))
  allocate(lats(nlats))
  allocate(lons(nlons))
  allocate(vort(nlons, nlats, nlevs, ntime))
  ! Read the latitude and longitude data.
  call check( nf90_inq_varid(ncid, REC_NAME, rec_varid) )
  call check( nf90_inq_varid(ncid, LVL_NAME, lvl_varid) )
  call check( nf90_inq_varid(ncid, LAT_NAME, lat_varid) )
  call check( nf90_inq_varid(ncid, LON_NAME, lon_varid) )
  call check( nf90_get_var(ncid, lvl_varid, levs) )
  call check( nf90_get_var(ncid, lat_varid, lats) )
  call check( nf90_get_var(ncid, lon_varid, lons) )
  print*,levs
  print*,maxval(lats), minval(lats)
  print*,maxval(lons), minval(lons)


  call read_data(ncid, VORT_NAME, 4, vort)
  print*, minval(vort), maxval(vort)
  
  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file.
  call check( nf90_close(ncid) )

  ! If we got this far, everything worked as expected. Yipee! 
  print *,"end"

contains
  subroutine read_data(ncid, var_name, ndim, var_data)
    implicit none

    real(wp), intent(out) :: var_data
    character (len=nf90_max_name), intent(in) :: var_name
    integer, intent(in) :: ncid
    integer, intent(in) :: ndim

    integer :: var_id
    real(wp) :: scale_factor, add_offset

    call check( nf90_inq_varid(ncid, var_name, var_id) )
    !call check( nf90_inquire_variable(ncid, var_id, xtype=xtype, ndims=ndims) )
    call check( nf90_get_var(ncid, var_id, var_data) )
    call check( nf90_get_att(ncid, var_id, "scale_factor", scale_factor) )
    call check( nf90_get_att(ncid, var_id, "add_offset", add_offset) )
  
    var_data = scale_factor * var_data + add_offset
    
  end subroutine read_data


  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print*, 'status:', status
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  
end program read_netcdf
