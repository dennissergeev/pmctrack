program read_netcdf
  use netcdf
  !use iso_fortran_env
  implicit none
  
  integer, parameter :: wp = kind(0.0)

  ! This is the name of the data file we will read.
  character (len=*), parameter :: CONFIG_FILE = "config.txt"
  character (len=256) :: nc_file_name
  character (len=256) :: datadir


  character (len=*), parameter :: LVL_NAME = "level"
  character (len=*), parameter :: LAT_NAME = "latitude"
  character (len=*), parameter :: LON_NAME = "longitude"
  character (len=*), parameter :: REC_NAME = "time"
  character (len=*), parameter :: VORT_NAME = "vo"
  integer :: ntime, nlvls, nlats, nlons

  real(wp), allocatable :: vort(:, :, :, :) 
  real(wp), allocatable :: vort3d(:, :, :) 
  integer, allocatable :: time(:) 
  integer, allocatable :: lvls(:) 
  real(wp), allocatable :: lats(:) 
  real(wp), allocatable :: lons(:) 

  integer :: vort_varid
  ! Loop indices
  integer :: lvl, lat, lon, t, i

  integer :: proj
  integer :: year
  integer :: month
  integer :: thelevel
  integer :: lvl_idx

  character(len=nf90_max_name), dimension(4) :: DIM_NAMES 
  integer, dimension(4) :: DIMS

!call get_command_argument(1, FILE_NAME)

  DIM_NAMES(1) = trim(REC_NAME)
  DIM_NAMES(2) = trim(LVL_NAME) 
  DIM_NAMES(3) = trim(LAT_NAME)
  DIM_NAMES(4) = trim(LON_NAME)

  call get_config_params(CONFIG_FILE, datadir, year, month, thelevel, proj)

  write(nc_file_name, '(A,A,A,I4.4,A,I2.2,A,A,A)') trim(datadir), '/', &
                                                 & 'era5.an.sfc.', year, '.', &
                                                 & month, '.', VORT_NAME, '.nc'
  print*, nc_file_name

  call get_dims(nc_file_name, DIM_NAMES, ntime, nlvls, nlats, nlons)

  allocate(time(ntime))
  allocate(lvls(nlvls))
  allocate(lats(nlats))
  allocate(lons(nlons))
  allocate(vort(nlons, nlats, nlvls, ntime))
  allocate(vort3d(nlons, nlats, ntime))

  call get_coords(nc_file_name, DIM_NAMES, time, lvls, lats, lons)

  call get_data_4d(nc_file_name, VORT_NAME, vort)
  
  lvl_idx = minloc(abs(lvls-thelevel), 1)
  call get_one_level(nc_file_name, VORT_NAME, vort3d, lvl_idx)
  print*, minval(vort3d), maxval(vort3d)

  ! If we got this far, everything worked as expected. Yipee! 
  print *,"end"
  deallocate(time)
  deallocate(lvls)
  deallocate(lats)
  deallocate(lons)
  deallocate(vort)
  deallocate(vort3d)

contains
  
  subroutine get_dims(nc_file_name, dim_names, ntime, nlvls, nlats, nlons)
    implicit none

    character(len=*)              , intent(in)    :: nc_file_name
    character(len=*), dimension(4), intent(in)    :: dim_names
    integer                       , intent(inout) :: ntime 
    integer                       , intent(inout) :: nlvls 
    integer                       , intent(inout) :: nlats 
    integer                       , intent(inout) :: nlons 
    ! Local variables
    integer                                       :: ncid
    integer                                       :: dimid
    integer                       , dimension(4)  :: dims

    call check( nf90_open(nc_file_name, nf90_nowrite, ncid) )

    do i = 1, size(dim_names)
      call check( nf90_inq_dimid(ncid, dim_names(i), dimid) )
      call check( nf90_inquire_dimension(ncid, dimid, len=dims(i)) )
    end do
    ntime = dims(1)
    nlvls = dims(2)
    nlats = dims(3)
    nlons = dims(4)

    call check( nf90_close(ncid) )
  end subroutine get_dims

  subroutine get_coords(nc_file_name, dim_names, time, lvls, lats, lons)
    implicit none

    character(len=*)              , intent(in)    :: nc_file_name
    character(len=*), dimension(4), intent(in)    :: dim_names
    integer                       , intent(inout) :: time(:) 
    integer                       , intent(inout) :: lvls(:) 
    real(wp)                      , intent(inout) :: lats(:) 
    real(wp)                      , intent(inout) :: lons(:) 
    ! Local variables
    integer                                       :: ncid
    integer                                       :: var_id

    call check( nf90_open(nc_file_name, nf90_nowrite, ncid) )

    call check( nf90_inq_varid(ncid, dim_names(1), var_id) )
    call check( nf90_get_var(ncid, var_id, time) )
    call check( nf90_inq_varid(ncid, dim_names(2), var_id) )
    call check( nf90_get_var(ncid, var_id, lvls) )
    call check( nf90_inq_varid(ncid, dim_names(3), var_id) )
    call check( nf90_get_var(ncid, var_id, lats) )
    call check( nf90_inq_varid(ncid, dim_names(4), var_id) )
    call check( nf90_get_var(ncid, var_id, lons) )

    call check( nf90_close(ncid) )
  end subroutine get_coords

  subroutine get_data_4d(nc_file_name, var_name, var_data)
    implicit none

    character(len=*)              , intent(in)    :: nc_file_name
    character(len=*)              , intent(in)    :: var_name
    real(wp)                      , intent(inout) :: var_data(:, :, :, :)
    ! Local variables
    integer                                       :: ncid
    integer                                       :: var_id
    real(wp)                                      :: scale_factor, add_offset

    call check( nf90_open(nc_file_name, nf90_nowrite, ncid) )

    call check( nf90_inq_varid(ncid, var_name, var_id) )
    !call check( nf90_inquire_variable(ncid, var_id, xtype=xtype, ndims=ndims) )
    call check( nf90_get_var(ncid, var_id, var_data) )
    call check( nf90_get_att(ncid, var_id, "scale_factor", scale_factor) )
    call check( nf90_get_att(ncid, var_id, "add_offset", add_offset) )
  
    var_data = scale_factor * var_data + add_offset

    call check( nf90_close(ncid) )
  end subroutine get_data_4d

  subroutine get_data_3d(nc_file_name, var_name, var_data)
    implicit none

    character(len=*)              , intent(in)    :: nc_file_name
    character(len=*)              , intent(in)    :: var_name
    real(wp)                      , intent(inout) :: var_data(:, :, :)
    ! Local variables
    integer                                       :: ncid
    integer                                       :: var_id
    real(wp)                                      :: scale_factor, add_offset

    call check( nf90_open(nc_file_name, nf90_nowrite, ncid) )

    call check( nf90_inq_varid(ncid, var_name, var_id) )
    call check( nf90_get_var(ncid, var_id, var_data) )
    call check( nf90_get_att(ncid, var_id, "scale_factor", scale_factor) )
    call check( nf90_get_att(ncid, var_id, "add_offset", add_offset) )
  
    var_data = scale_factor * var_data + add_offset

    call check( nf90_close(ncid) )
  end subroutine get_data_3d

  subroutine get_one_level(nc_file_name, var_name, var_data, lvl_idx)
    implicit none

    character(len=*)              , intent(in)    :: nc_file_name
    character(len=*)              , intent(in)    :: var_name
    real(wp)                      , intent(inout) :: var_data(:, :, :)
    integer                       , intent(in)    :: lvl_idx
    ! Local variables
    integer                                       :: ncid
    integer                                       :: var_id
    real(wp)                                      :: scale_factor, add_offset
    integer                                       :: shp(3)

    shp = shape(var_data)

    call check( nf90_open(nc_file_name, nf90_nowrite, ncid) )

    call check( nf90_inq_varid(ncid, var_name, var_id) )
    call check( nf90_get_var(ncid, var_id, var_data,  &
                             start=(/1, 1, lvl_idx, 1/), &
                             count=(/shp(1), shp(2), 1, shp(3)/)) )
    call check( nf90_get_att(ncid, var_id, "scale_factor", scale_factor) )
    call check( nf90_get_att(ncid, var_id, "add_offset", add_offset) )
  
    var_data = scale_factor * var_data + add_offset

    call check( nf90_close(ncid) )
  end subroutine get_one_level


  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print*, 'status:', status
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  


  subroutine get_config_params(config_file, datadir, year, month, &
    & thelevel, proj)

    implicit none
  
    character(len=*), intent(in) :: config_file
    character(len=*), intent(inout) :: datadir
    integer         , intent(inout) :: year
    integer         , intent(inout) :: month
    integer         , intent(inout) :: thelevel
    integer         , intent(inout) :: proj
    ! Local variables
    integer, parameter :: fh = 999
    character (len=256) :: buffer
    character (len=256) :: label
    integer :: line
    integer :: ios
    integer :: sep
    
    ios = 0
    line = 0
  
    open (fh, file=CONFIG_FILE, form='formatted', status='old', action='read')
    do while (ios == 0)
      read(fh, '(A)', iostat=ios) buffer
      if (ios == 0) then
        line = line + 1
  
        ! Find the first instance of whitespace.  Split label and data.
        sep = scan(buffer, '=')
        if (sep == 0) then
          label = buffer
        else
          label = buffer(1:sep-1)
          buffer = buffer(sep+1:len_trim(buffer))
        end if
  
        select case (trim(label))
        case ('datadir'); read(buffer, *, iostat=ios) datadir
        case ('year'); read(buffer, *, iostat=ios) year
        case ('month'); read(buffer, *, iostat=ios) month
        case ('thelevel'); read(buffer, *, iostat=ios) thelevel
        case ('proj'); read(buffer, *, iostat=ios) proj
        case default
          if (index (trim(label), "#") == 1) then
            print*, buffer
          else
            write(*, *) 'ConfigParseWarning: Skipping invalid line', line
          end if
        end select
      end if
    end do
    close(fh)
  end subroutine get_config_params
end program read_netcdf
