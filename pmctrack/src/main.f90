program main
  use datetime_module

  use types, only : wp
  use params, only : get_config_params, dbg, datadir,              &
    & year_start, month_start, day_start, hour_start,              &
    & year_end, month_end, day_end, hour_end,                      &
    & vort_name, u_name, v_name, psea_name, land_name,             &
    & vor_lvl, steer_lvl_btm, steer_lvl_top
  use nc_io, only : get_dims, get_time, get_coords, get_one_level, &
    & get_data_4d, get_data_3d, get_data_2d
  use utils, only : apply_mask_2d

  implicit none

  character(len=*)                , parameter :: LVL_NAME = "level"
  character(len=*)                , parameter :: LAT_NAME = "latitude"
  character(len=*)                , parameter :: LON_NAME = "longitude"
  character(len=*)                , parameter :: REC_NAME = "time"
  character(len=256), dimension(4)            :: DIM_NAMES
  character(len=256)                          :: nc_file_name
  real(wp)                                    :: lonin
  real(wp)                                    :: latin
  real(wp)                                    :: del_t
  integer                                     :: nt_per_file
  ! Coordinate arrays
  integer                                     :: ntime, nlvls, nlats, nlons
  integer          , allocatable              :: time_temp(:)
  integer          , allocatable              :: time(:)
  integer          , allocatable              :: lvls(:)
  real(wp)         , allocatable              :: lats(:)
  real(wp)         , allocatable              :: lons(:)
  ! Data arrays
  real(wp)         , allocatable              :: vor(:, :, :)
  real(wp)         , allocatable              :: psea(:, :, :)
  real(wp)         , allocatable              :: u(:, :, :, :)
  real(wp)         , allocatable              :: v(:, :, :, :)
  real(wp)         , allocatable              :: land_mask(:, :)

  ! Local variables
  ! Indices
  integer                                     :: lvl_idx
  integer                                     :: steer_idx_btm, steer_idx_top
  integer                                     :: nsteer_lvl
  integer                                     :: kt

  ! Time and date variables
  type(datetime)                              :: calendar_start
  type(datetime)                              :: dt_start, dt_end
  type(datetime)                              :: dt_min ! within a file
  type(timedelta)                             :: td
  integer                                     :: time_idx


  ! Store dimension names in one array
  DIM_NAMES(1) = trim(REC_NAME)
  DIM_NAMES(2) = trim(LVL_NAME)
  DIM_NAMES(3) = trim(LAT_NAME)
  DIM_NAMES(4) = trim(LON_NAME)

  call get_config_params()

  ! Get dimensions from the vorticity file using first year and first month
  ! Assume all the other files are organised in the same way
  write(nc_file_name, '(A,A,A,I4.4,A,I2.2,A,A,A)') trim(datadir), '/', &
                                                 & 'era5.an.pl.', year_start, '.', &
                                                 & month_start, '.', &
                                                 & trim(vort_name), '.nc'
  call get_dims(nc_file_name, DIM_NAMES, nt_per_file, nlvls, nlats, nlons)

  ! Time & calendar
  allocate(time_temp(0:nt_per_file-1))
  call get_time(nc_file_name, REC_NAME, time_temp)

  calendar_start = datetime(1900, 1, 1)  ! TODO: define automatically
  del_t = (time_temp(1) - time_temp(0)) * 3600
  print*, 'time_temp', time_temp(0)
  print*, 'del_t', del_t
  dt_start = datetime(year_start, month_start, day_start, hour_start)
  dt_end = datetime(year_end, month_end, day_end, hour_end)
  print*,'dt_start', dt_start
  print*,'dt_end', dt_end
  td = dt_end - dt_start
  ntime = int(td%total_seconds() / del_t)
  print*, 'ntime=', ntime

  td = timedelta(hours=time_temp(0))
  dt_min = calendar_start + td
  td = dt_start - dt_min
  time_idx = int(td%total_seconds() / del_t) + 1 ! time_temp(0) +
  print*, 'Start date: ', time_idx

  if (month_end > month_start) then
   print*,''
  end if
  deallocate(time_temp)

  allocate(time(ntime))
  allocate(lvls(nlvls))
  allocate(lats(0:nlats-1))
  allocate(lons(0:nlons-1))

  call get_coords(nc_file_name, DIM_NAMES, lons, lats, lvls, &
    & time, time_idx, ntime)
  lvl_idx = minloc(abs(lvls - vor_lvl), 1)
  steer_idx_btm = minloc(abs(lvls - steer_lvl_btm), 1)
  steer_idx_top = minloc(abs(lvls - steer_lvl_top), 1)
  nsteer_lvl = steer_idx_btm - steer_idx_top + 1
  !lvls = lvls(nlvls:1:-1)
  lats = lats(nlats-1:0:-1)
  ! Calculate grid spacing assuming the grid is uniform
  lonin = lons(1) - lons(0)
  latin = lats(1) - lats(0)

  ! Define array sizes
  allocate(vor      (0:nlons-1, 0:nlats-1, ntime))
  allocate(u        (0:nlons-1, 0:nlats-1, nsteer_lvl, ntime))
  allocate(v        (0:nlons-1, 0:nlats-1, nsteer_lvl, ntime))
  allocate(psea     (0:nlons-1, 0:nlats-1, ntime))
  allocate(land_mask(0:nlons-1, 0:nlats-1))
  ! Read vorticity at the specified level
  call get_one_level(nc_file_name, vort_name, lvl_idx, time_idx, ntime, vor)
  vor = vor(:, nlats-1:0:-1, :)
  ! print*, shape(vor), minval(vor(:, :, 0)), maxval(vor(:, :, 0))
  ! print*, shape(vor), minval(vor), maxval(vor)

  write(nc_file_name, '(A,A,A,I4.4,A,I2.2,A,A,A)') trim(datadir), '/', &
                                                 & 'era5.an.pl.', year_start, '.', &
                                                 & month_start, '.', &
                                                 & trim(u_name), '.nc'
  call get_data_4d(nc_file_name, u_name, time_idx, ntime, &
    & steer_idx_top, nsteer_lvl, u)
  u = u(:, nlats-1:0:-1, :, :)
  ! print*, shape(u), minval(u(:, :, :, 0)), maxval(u(:, :, :, 0))

  write(nc_file_name, '(A,A,A,I4.4,A,I2.2,A,A,A)') trim(datadir), '/', &
                                                 & 'era5.an.pl.', year_start, '.', &
                                                 & month_start, '.', &
                                                 & trim(v_name), '.nc'
  call get_data_4d(nc_file_name, v_name, time_idx, ntime, &
    & steer_idx_top, nsteer_lvl, v)
  v = v(:, nlats-1:0:-1, :, :)

  write(nc_file_name, '(A,A,A,I4.4,A,I2.2,A,A,A)') trim(datadir), '/', &
                                                 & 'era5.an.sfc.', year_start, '.', &
                                                 & month_start, '.', &
                                                 & trim(psea_name), '.nc'
  call get_data_3d(nc_file_name, psea_name, time_idx, ntime, psea)
  psea = psea(:, nlats-1:0:-1, :)
  psea = 1e-2 * psea

  write(nc_file_name, '(A,A,A,A)') trim(datadir), '/', trim(land_name), '.nc'
  call get_data_2d(nc_file_name, land_name, land_mask)
  land_mask = land_mask(:, nlats-1:0:-1)

  do kt = 1, ntime
    call apply_mask_2d(vor(:, :, kt), nlons-1, nlats-1, land_mask)
  end do
  ! print*, shape(v), minval(v), maxval(v)
  ! print*, shape(psea), minval(psea), maxval(psea)
  ! print*, lvls

  !write (*,*)'nx=', nx, 'ny=', ny,' nt=', nt, 'nz=', nz
  !write (*,*)'nx1=', nx1, 'nx2=', nx2, 'ny1=', ny1, 'ny2=', ny2
  call tracking_main(vor, u, v, psea, &
    & nlons-1, nlats-1, nsteer_lvl, lvls(1:nsteer_lvl), ntime, &
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
