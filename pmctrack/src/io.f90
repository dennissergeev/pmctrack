program read_netcdf
  use netcdf
  use types, only : wp
  use nc_io, only : get_dims, get_coords, get_one_level, &
    & get_data_4d, get_data_3d, get_data_2d
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
  real(wp)           :: lon0
  real(wp)           :: lat0
  real(wp)           :: lonin
  real(wp)           :: latin
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
    & thelevel, proj, vert_grid, nx1, nx2, ny1, ny2,              &
    & smth_type, nsmth_x, nsmth_y, r_smth,                        &
    & zeta_max0, zeta_min0, int_zeta_min0, gamma,                 &
    & d_cf_min, size_synop, del_psea_min, distance_ec,            &
    & steering_type, n_steering_x, n_steering_y, r_steering,      &
    & track_type, del_lon, del_lat, del_r, period_min)
#ifdef debug
  print*, datadir, year, month,       &
    & vort_name, u_name, v_name, psea_name, land_name,            &
    & thelevel, proj, vert_grid, nx1, nx2, ny1, ny2,              &
    & smth_type, nsmth_x, nsmth_y, r_smth,                        &
    & zeta_max0, zeta_min0, int_zeta_min0, gamma,                 &
    & d_cf_min, size_synop, del_psea_min, distance_ec,            &
    & steering_type, n_steering_x, n_steering_y, r_steering,      &
    & track_type, del_lon, del_lat, del_r, period_min
#endif

  write(nc_file_name, '(A,A,A,I4.4,A,I2.2,A,A,A)') trim(datadir), '/', &
                                                 & 'era5.an.sfc.', year, '.', &
                                                 & month, '.', &
                                                 & trim(vort_name), '.nc'

  ! Get dimensions from the vorticity file
  call get_dims(nc_file_name, DIM_NAMES, ntime, nlvls, nlats, nlons)

  allocate(time(0:ntime-1))
  allocate(lvls(0:nlvls-1))
  allocate(lats(0:nlats-1))
  allocate(lons(0:nlons-1))
  allocate(vor      (0:nlons-1, 0:nlats-1, 0:ntime-1))
  allocate(u        (0:nlons-1, 0:nlats-1, 0:nlvls-1, 0:ntime-1))
  allocate(v        (0:nlons-1, 0:nlats-1, 0:nlvls-1, 0:ntime-1))
  allocate(psea     (0:nlons-1, 0:nlats-1, 0:ntime-1))
  allocate(land_mask(0:nlons-1, 0:nlats-1))

  call get_coords(nc_file_name, DIM_NAMES, time, lvls, lats, lons)
  lvls = lvls(nlvls-1:0:-1)
  lats = lats(nlats-1:0:-1)
  lon0 = lons(0)
  lat0 = lats(0)
  ! Calculate grid spacing and time increment assuming the grid is uniform
  lonin = lons(1) - lons(0)
  latin = lats(1) - lats(0)
  del_t = (time(1) - time(0)) * 3600

  ! Read vorticity at the specified level
  lvl_idx = minloc(abs(lvls - thelevel), 1)
  call get_one_level(nc_file_name, vort_name, vor, lvl_idx)
  vor = vor(:, nlats-1:0:-1, :)

  write(nc_file_name, '(A,A,A,I4.4,A,I2.2,A,A,A)') trim(datadir), '/', &
                                                 & 'era5.an.sfc.', year, '.', &
                                                 & month, '.', &
                                                 & trim(u_name), '.nc'
  call get_data_4d(nc_file_name, u_name, u)
  u = u(:, nlats-1:0:-1, :, :)

  write(nc_file_name, '(A,A,A,I4.4,A,I2.2,A,A,A)') trim(datadir), '/', &
                                                 & 'era5.an.sfc.', year, '.', &
                                                 & month, '.', &
                                                 & trim(v_name), '.nc'
  call get_data_4d(nc_file_name, v_name, v)
  v = v(:, nlats-1:0:-1, :, :)

  write(nc_file_name, '(A,A,A,I4.4,A,I2.2,A,A,A)') trim(datadir), '/', &
                                                 & 'era5.an.sfc.', year, '.', &
                                                 & month, '.', &
                                                 & trim(psea_name), '.nc'
  call get_data_3d(nc_file_name, psea_name, psea)
  psea = psea(:, nlats-1:0:-1, :)

  write(nc_file_name, '(A,A,A,I4.4,A,I2.2,A,A,A)') trim(datadir), '/', &
                                                 & 'era5.an.sfc.', year, '.', &
                                                 & month, '.', &
                                                 & trim(land_name), '.nc'
  call get_data_2d(nc_file_name, land_name, land_mask)
  land_mask = land_mask(:, nlats-1:0:-1)

  deallocate(time)
  deallocate(lvls)
  deallocate(lats)
  deallocate(lons)
  deallocate(vor)
  deallocate(u)
  deallocate(v)
  deallocate(psea)
  deallocate(land_mask)

end program read_netcdf
