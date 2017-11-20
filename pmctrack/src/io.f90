program read_netcdf
  use netcdf
  use types, only : wp
  !use iso_fortran_env
  implicit none
  
  character(len=*), parameter :: CONFIG_FILE = "config.ini"
  character(len=256) :: datadir
  character(len=256) :: vort_name
  character(len=256) :: u_name
  character(len=256) :: v_name
  character(len=256) :: psea_name
  character(len=256) :: land_name
  integer            :: year
  integer            :: month
  integer            :: thelevel
  integer            :: proj
  integer            :: vert_grid
  integer            :: nx1, nx2, ny1, ny2
  integer            :: nt
  real(wp)           :: del_t
  integer            :: smth_type
  integer            :: nsmth_x, nsmth_y
  real(wp)           :: r_smth
  real(wp)           :: zeta_max0, zeta_min0
  real(wp)           :: int_zeta_min0, gamma
  real(wp)           :: d_cf_min, size_synop
  real(wp)           :: del_psea_min, distance_ec
  integer            :: steering_type
  integer            :: n_steering_x, n_steering_y
  real(wp)           :: r_steering
  integer            :: track_type
  real(wp)           :: del_lon, del_lat,del_r
  integer            :: period_min

  character (len=256)            :: nc_file_name
  character (len=*), parameter   :: LVL_NAME = "level"
  character (len=*), parameter   :: LAT_NAME = "latitude"
  character (len=*), parameter   :: LON_NAME = "longitude"
  character (len=*), parameter   :: REC_NAME = "time"
  integer                        :: ntime, nlvls, nlats, nlons
  integer          , allocatable :: time(:) 
  integer          , allocatable :: lvls(:) 
  real(wp)         , allocatable :: lats(:) 
  real(wp)         , allocatable :: lons(:) 

  real(wp)         , allocatable :: vor(:, :, :) 
  real(wp)         , allocatable :: psea(:, :, :) 
  real(wp)         , allocatable :: u(:, :, :, :) 
  real(wp)         , allocatable :: v(:, :, :, :) 
  real(wp)         , allocatable :: land_mask(:, :) 

  ! Loop indices
  integer :: lvl, lat, lon, t, i

  integer :: lvl_idx

  character(len=nf90_max_name), dimension(4) :: DIM_NAMES 
  integer, dimension(4) :: DIMS

  !call get_command_argument(1, FILE_NAME)

  DIM_NAMES(1) = trim(REC_NAME)
  DIM_NAMES(2) = trim(LVL_NAME) 
  DIM_NAMES(3) = trim(LAT_NAME)
  DIM_NAMES(4) = trim(LON_NAME)

  call get_config_params(CONFIG_FILE, datadir, year, month,       &
    & vort_name, u_name, v_name, psea_name, land_name,            &
    & thelevel, proj, vert_grid, nx1, nx2, ny1, ny2, nt, del_t,   &
    & smth_type, nsmth_x, nsmth_y, r_smth,                        &
    & zeta_max0, zeta_min0, int_zeta_min0, gamma,                 &
    & d_cf_min, size_synop, del_psea_min, distance_ec,            &
    & steering_type, n_steering_x, n_steering_y, r_steering,      &
    & track_type, del_lon, del_lat, del_r, period_min)

  write(nc_file_name, '(A,A,A,I4.4,A,I2.2,A,A,A)') trim(datadir), '/', &
                                                 & 'era5.an.sfc.', year, '.', &
                                                 & month, '.', &
                                                 & trim(vort_name), '.nc'

  ! Get dimensions from the vorticity file
  call get_dims(nc_file_name, DIM_NAMES, ntime, nlvls, nlats, nlons)

  allocate(time(ntime))
  allocate(lvls(nlvls))
  allocate(lats(nlats))
  allocate(lons(nlons))
  allocate(vor      (nlons, nlats, ntime))
  allocate(u        (nlons, nlats, nlvls, ntime))
  allocate(v        (nlons, nlats, nlvls, ntime))
  allocate(psea     (nlons, nlats, ntime))
  allocate(land_mask(nlons, nlats))

  call get_coords(nc_file_name, DIM_NAMES, time, lvls, lats, lons)

  ! Read vorticity at the specified level
  lvl_idx = minloc(abs(lvls - thelevel), 1)
  call get_one_level(nc_file_name, vort_name, vor, lvl_idx)
  print*, minval(vor), maxval(vor)

  write(nc_file_name, '(A,A,A,I4.4,A,I2.2,A,A,A)') trim(datadir), '/', &
                                                 & 'era5.an.sfc.', year, '.', &
                                                 & month, '.', &
                                                 & trim(u_name), '.nc'
  call get_data_4d(nc_file_name, u_name, u)
  print*, minval(u), maxval(u)

  write(nc_file_name, '(A,A,A,I4.4,A,I2.2,A,A,A)') trim(datadir), '/', &
                                                 & 'era5.an.sfc.', year, '.', &
                                                 & month, '.', &
                                                 & trim(v_name), '.nc'
  call get_data_4d(nc_file_name, v_name, v)
  print*, minval(v), maxval(v)

  write(nc_file_name, '(A,A,A,I4.4,A,I2.2,A,A,A)') trim(datadir), '/', &
                                                 & 'era5.an.sfc.', year, '.', &
                                                 & month, '.', &
                                                 & trim(psea_name), '.nc'
  call get_data_3d(nc_file_name, psea_name, psea)
  print*, minval(psea), maxval(psea)

  ! If we got this far, everything worked as expected. Yipee! 
  print *,"end"
  deallocate(time)
  deallocate(lvls)
  deallocate(lats)
  deallocate(lons)
  deallocate(vor)
  deallocate(u)
  deallocate(v)
  deallocate(psea)
  deallocate(land_mask)

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
    & vort_name, u_name, v_name, psea_name, land_name,            &
    & thelevel, proj, vert_grid, nx1, nx2, ny1, ny2, nt, del_t,   &
    & smth_type, nsmth_x, nsmth_y, r_smth,                        &
    & zeta_max0, zeta_min0, int_zeta_min0, gamma,                 &
    & d_cf_min, size_synop, del_psea_min, distance_ec,            &
    & steering_type, n_steering_x, n_steering_y, r_steering,      &
    & track_type, del_lon, del_lat, del_r, period_min)

    implicit none
  
    character(len=*), intent(in)    :: config_file
    character(len=*), intent(inout) :: datadir
    character(len=*), intent(inout) :: vort_name
    character(len=*), intent(inout) :: u_name
    character(len=*), intent(inout) :: v_name
    character(len=*), intent(inout) :: psea_name
    character(len=*), intent(inout) :: land_name
    integer         , intent(inout) :: year
    integer         , intent(inout) :: month
    integer         , intent(inout) :: thelevel
    integer         , intent(inout) :: proj
    integer         , intent(inout) :: vert_grid
    integer         , intent(inout) :: nx1, nx2, ny1, ny2
    integer         , intent(inout) :: nt
    real(wp)        , intent(inout) :: del_t
    ! parameter for smoothing of vorticity
    integer         , intent(inout) :: smth_type
    integer         , intent(inout) :: nsmth_x, nsmth_y
    real(wp)        , intent(inout) :: r_smth
    ! parameter for detecting vortex
    real(wp)        , intent(inout) :: zeta_max0, zeta_min0
    real(wp)        , intent(inout) :: int_zeta_min0, gamma
    ! parameter for excluding the synoptic scale disturbances
    real(wp)        , intent(inout) :: d_cf_min, size_synop
    real(wp)        , intent(inout) :: del_psea_min, distance_ec
    ! parameter for calculating steering winds
    integer         , intent(inout) :: steering_type
    integer         , intent(inout) :: n_steering_x, n_steering_y
    real(wp)        , intent(inout) :: r_steering
    ! parameter for linking vortex
    integer         , intent(inout) :: track_type
    real(wp)        , intent(inout) :: del_lon, del_lat,del_r
    ! parameter for checking the track
    integer         , intent(inout) :: period_min

    ! Local variables
    integer, parameter              :: fh = 999
    character(len=256)              :: buffer
    character(len=256)              :: label
    integer                         :: line
    integer                         :: ios
    integer                         :: sep
    
    ios = 0
    line = 0
  
    open (fh, file=CONFIG_FILE, form='formatted', status='old', &
      &       iostat=ios, action='read')
    if (ios == 0) then
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
          case ('vort_name'); read(buffer, *, iostat=ios) vort_name
          case ('u_name'); read(buffer, *, iostat=ios) u_name
          case ('v_name'); read(buffer, *, iostat=ios) v_name
          case ('psea_name'); read(buffer, *, iostat=ios) psea_name
          case ('land_name'); read(buffer, *, iostat=ios) land_name
          case ('proj'); read(buffer, *, iostat=ios) proj
          case ('vert_grid'); read(buffer, *, iostat=ios) vert_grid
          case ('nx1'); read(buffer, *, iostat=ios) nx1
          case ('nx2'); read(buffer, *, iostat=ios) nx2
          case ('ny1'); read(buffer, *, iostat=ios) ny1
          case ('ny2'); read(buffer, *, iostat=ios) ny2
          !case ('lons'); read(buffer, *, iostat=ios) lons
          !case ('lats'); read(buffer, *, iostat=ios) lats
          !case ('lonin'); read(buffer, *, iostat=ios) lonin
          !case ('latin'); read(buffer, *, iostat=ios) latin
          case ('del_t'); read(buffer, *, iostat=ios) del_t
          case ('smth_type'); read(buffer, *, iostat=ios) smth_type
          case ('nsmth_x'); read(buffer, *, iostat=ios) nsmth_x
          case ('nsmth_y'); read(buffer, *, iostat=ios) nsmth_y
          case ('r_smth'); read(buffer, *, iostat=ios) r_smth
          case ('zeta_max0'); read(buffer, *, iostat=ios) zeta_max0
          case ('zeta_min0'); read(buffer, *, iostat=ios) zeta_min0
          case ('int_zeta_min0'); read(buffer, *, iostat=ios) int_zeta_min0
          case ('gamma'); read(buffer, *, iostat=ios) gamma
          case ('d_cf_min'); read(buffer, *, iostat=ios) d_cf_min
          case ('size_synop'); read(buffer, *, iostat=ios) size_synop
          case ('distance_ec'); read(buffer, *, iostat=ios) distance_ec
          case ('del_psea_min'); read(buffer, *, iostat=ios) del_psea_min
          case ('steering_type'); read(buffer, *, iostat=ios) steering_type
          case ('n_steering_x'); read(buffer, *, iostat=ios) n_steering_x
          case ('n_steering_y'); read(buffer, *, iostat=ios) n_steering_y
          case ('r_steering'); read(buffer, *, iostat=ios) r_steering
          case ('track_type'); read(buffer, *, iostat=ios) track_type
          case ('del_lon'); read(buffer, *, iostat=ios) del_lon
          case ('del_lat'); read(buffer, *, iostat=ios) del_lat
          case ('del_r'); read(buffer, *, iostat=ios) del_r
          case ('period_min'); read(buffer, *, iostat=ios) period_min
          case default
            if (index (trim(label), "#") == 1) then
              print*, buffer
            else
              write(*, *) 'ConfigParseWarning: Skipping invalid line', line
            end if
          end select
        end if
      end do
    end if
    close(fh)
  end subroutine get_config_params
end program read_netcdf
