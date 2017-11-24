program main
  use netcdf
  use datetime_module

  use types, only : wp
  use params, only :  get_config_params, datadir,                 &
    & year_start, month_start, day_start,                         &
    & year_end, month_end, day_end,                               &
    & vort_name, u_name, v_name, psea_name, land_name,            &
    & thelevel, proj, vert_grid, nx1, nx2, ny1, ny2,              &
    & smth_type, nsmth_x, nsmth_y, r_smth,                        &
    & zeta_max0, zeta_min0, int_zeta_min0, gamma,                 &
    & d_cf_min, size_synop, del_psea_min, distance_ec,            &
    & steering_type, n_steering_x, n_steering_y, r_steering,      &
    & track_type, del_lon, del_lat, del_r, period_min
  use nc_io, only : get_dims, get_coords, get_one_level, &
    & get_data_4d, get_data_3d, get_data_2d
  use utils, only : apply_mask_2d

  implicit none
  
  real(wp)           :: lon0
  real(wp)           :: lat0
  real(wp)           :: lonin
  real(wp)           :: latin
  real(wp)           :: del_t
  integer            :: ntimes

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
  integer :: kt

  character(len=nf90_max_name), dimension(4) :: DIM_NAMES 
  ! integer, dimension(4) :: DIMS

  type(datetime) :: a, b, calendar_start  
  type(timedelta) :: td

  DIM_NAMES(1) = trim(REC_NAME)
  DIM_NAMES(2) = trim(LVL_NAME) 
  DIM_NAMES(3) = trim(LAT_NAME)
  DIM_NAMES(4) = trim(LON_NAME)

  call get_config_params()
#ifdef debug
  print*, datadir, &
    & year_start, month_start, day_start,                         &
    & year_end, month_end, day_end,                               &
    & vort_name, u_name, v_name, psea_name, land_name,            &
    & thelevel, proj, vert_grid, nx1, nx2, ny1, ny2,              &
    & smth_type, nsmth_x, nsmth_y, r_smth,                        &
    & zeta_max0, zeta_min0, int_zeta_min0, gamma,                 &
    & d_cf_min, size_synop, del_psea_min, distance_ec,            &
    & steering_type, n_steering_x, n_steering_y, r_steering,      &
    & track_type, del_lon, del_lat, del_r, period_min
#endif

  write(nc_file_name, '(A,A,A,I4.4,A,I2.2,A,A,A)') trim(datadir), '/', &
                                                 & 'era5.an.pl.', year_start, '.', &
                                                 & month_start, '.', &
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

  a = datetime(year_start, month_start, day_start)
  b = datetime(year_end, month_end, day_end)
  td = b - a
  ntimes = int(td%total_seconds() / del_t )
  print*, 'ntimes=', ntimes

  calendar_start = datetime(1900, 1, 1)
  td = timedelta(hours=time(0))
  a = calendar_start + td
  print*, 'Start date: ', a

  if (month_end > month_start) then
   print*,'' 
  end if


  ! Read vorticity at the specified level
  lvl_idx = minloc(abs(lvls - thelevel), 1)
  call get_one_level(nc_file_name, vort_name, vor, lvl_idx)
  vor = vor(:, nlats-1:0:-1, :)

  write(nc_file_name, '(A,A,A,I4.4,A,I2.2,A,A,A)') trim(datadir), '/', &
                                                 & 'era5.an.pl.', year_start, '.', &
                                                 & month_start, '.', &
                                                 & trim(u_name), '.nc'
  call get_data_4d(nc_file_name, u_name, u)
  u = u(:, nlats-1:0:-1, :, :)
  write(nc_file_name, '(A,A,A,I4.4,A,I2.2,A,A,A)') trim(datadir), '/', &
                                                 & 'era5.an.pl.', year_start, '.', &
                                                 & month_start, '.', &
                                                 & trim(v_name), '.nc'
  call get_data_4d(nc_file_name, v_name, v)
  v = v(:, nlats-1:0:-1, :, :)
  write(nc_file_name, '(A,A,A,I4.4,A,I2.2,A,A,A)') trim(datadir), '/', &
                                                 & 'era5.an.sfc.', year_start, '.', &
                                                 & month_start, '.', &
                                                 & trim(psea_name), '.nc'
  call get_data_3d(nc_file_name, psea_name, psea)
  psea = psea(:, nlats-1:0:-1, :)

  write(nc_file_name, '(A,A,A,A)') trim(datadir), '/', trim(land_name), '.nc'
  call get_data_2d(nc_file_name, land_name, land_mask)
  land_mask = land_mask(:, nlats-1:0:-1)

  do kt = 1, ntimes
    call apply_mask_2d(vor(0:nlons-1, 0:nlats-1, kt), nlons-1, nlats-1, land_mask)
  end do

  call tracking_main(vor, u, v, psea, &
    & nlons-1, nlats-1, nlvls, lvls, ntimes, &
    & lons, lats, lonin, latin, del_t)

  deallocate(time)
  deallocate(lvls)
  deallocate(lats)
  deallocate(lons)
  deallocate(vor)
  deallocate(u)
  deallocate(v)
  deallocate(psea)
  deallocate(land_mask)

end program main
