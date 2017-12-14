program main
  use datetime_module

  use types, only : wp
  use constants, only : fillval, steer_nt
  use params, only : get_config_params, dbg,                                  &
    & datadir, prefix_sfc, prefix_lvl,                                        &
    & year_start, month_start, day_start, hour_start,                         &
    & year_end, month_end, day_end, hour_end,                                 &
    & vort_name, u_name, v_name, psea_name, land_name,                        &
    & vor_lvl, steer_lvl_btm, steer_lvl_top,                                  &
    & nx1, nx2, ny1, ny2
  use nc_io, only : get_dims, get_time, get_coords,                           &
    & get_xy_from_xyzt, get_xy_from_xyt, get_xyz_from_xyzt,                   &
    & get_data_2d
  use utils, only : apply_mask_2d, make_nc_file_name

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
  real(wp)         , allocatable              :: vor(:, :)
  real(wp)         , allocatable              :: psea(:, :)
  real(wp)         , allocatable              :: u(:, :, :, :)
  real(wp)         , allocatable              :: v(:, :, :, :)
  real(wp)         , allocatable              :: vor_ft(:, :, :)
  real(wp)         , allocatable              :: psea_ft(:, :, :)
  real(wp)         , allocatable              :: u_ft(:, :, :, :)
  real(wp)         , allocatable              :: v_ft(:, :, :, :)
  real(wp)         , allocatable              :: land_mask(:, :)

  ! Local variables
  ! Indices
  integer                                     :: lvl_idx
  integer                                     :: steer_idx_btm, steer_idx_top
  integer                                     :: nsteer_lvl
  integer                                     :: kt, kt2

  ! Time and date variables
  type(datetime)                              :: calendar_start
  type(datetime)                              :: dt_start, dt_end
  type(datetime)                              :: dt_min ! within a file
  type(timedelta)                             :: td
  integer                                     :: time_idx
  type(datetime)                              :: idt
  type(datetime)                              :: idt_pair(steer_nt)


  ! Store dimension names in one array
  DIM_NAMES(1) = trim(REC_NAME)
  DIM_NAMES(2) = trim(LVL_NAME)
  DIM_NAMES(3) = trim(LAT_NAME)
  DIM_NAMES(4) = trim(LON_NAME)

  call get_config_params()

  calendar_start = datetime(1900, 1, 1)  ! TODO: define automatically
  dt_start = datetime(year_start, month_start, day_start, hour_start)
  if (.not. dt_start%isValid()) then
    write(*, *) 'Start date ', dt_start, ' is not valid'
    stop
  endif
  dt_end = datetime(year_end, month_end, day_end, hour_end)
  if (.not. dt_end%isValid()) then
    write(*, *) 'End date ', dt_end, ' is not valid'
    stop
  endif

  ! Get dimensions from the vorticity file using first year and first month
  ! Assume all the other files are organised in the same way
  idt = dt_start
  call make_nc_file_name(nc_file_name, datadir, prefix_lvl, &
                       & idt%year, idt%month, vort_name)
  call get_dims(nc_file_name, DIM_NAMES, nt_per_file, nlvls, nlats, nlons)
  if (nx1 == -1) nx1 = 0
  if (nx2 == -1) nx2 = nlons-1
  if (ny1 == -1) ny1 = 0
  if (ny2 == -1) ny2 = nlats-1

  ! Time & calendar
  allocate(time_temp(0:nt_per_file-1))
  call get_time(nc_file_name, REC_NAME, time_temp)

  del_t = (time_temp(1) - time_temp(0)) * 3600
  print*, 'time_temp', time_temp(0)
  print*, 'del_t', del_t
  td = dt_end - dt_start
  ntime = int(td%total_seconds() / del_t) + 1
  print*, 'ntime=', ntime

  td = timedelta(hours=time_temp(0))
  dt_min = calendar_start + td
  td = dt_start - dt_min
  time_idx = int(td%total_seconds() / del_t) + 1 ! time_temp(0) +
  deallocate(time_temp)

  ! Assume space coordinates are the same for all files
  allocate(time(1))
  allocate(lvls(nlvls))
  allocate(lats(0:nlats-1))
  allocate(lons(0:nlons-1))

  call get_coords(nc_file_name, DIM_NAMES, lons, lats, lvls, &
    & time, 1, 1)
  print*,lvls
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
  allocate(vor_ft      (0:nlons-1, 0:nlats-1, ntime))
  allocate(u_ft        (0:nlons-1, 0:nlats-1, nsteer_lvl, ntime))
  allocate(v_ft        (0:nlons-1, 0:nlats-1, nsteer_lvl, ntime))
  allocate(psea_ft     (0:nlons-1, 0:nlats-1, ntime))
  allocate(vor      (0:nlons-1, 0:nlats-1))
  allocate(u        (0:nlons-1, 0:nlats-1, nsteer_lvl, steer_nt))
  allocate(v        (0:nlons-1, 0:nlats-1, nsteer_lvl, steer_nt))
  allocate(psea     (0:nlons-1, 0:nlats-1))
  allocate(land_mask(0:nlons-1, 0:nlats-1))

  vor_ft = fillval
  u_ft = fillval
  v_ft = fillval
  psea_ft = fillval

  write(nc_file_name, '(A,A,A,A)') trim(datadir), '/', trim(land_name), '.nc'
  call get_data_2d(nc_file_name, land_name, land_mask)
  land_mask = land_mask(:, nlats-1:0:-1)

  do kt = 1, ntime ! including both start and end dates
    call make_nc_file_name(nc_file_name, datadir, prefix_lvl, &
                         & idt%year, idt%month, vort_name)
    call get_dims(nc_file_name, DIM_NAMES, nt_per_file, nlvls, nlats, nlons)
    allocate(time_temp(0:nt_per_file-1))
    call get_time(nc_file_name, REC_NAME, time_temp)
    del_t = (time_temp(1) - time_temp(0)) * 3600
    td = timedelta(hours=time_temp(0))
    deallocate(time_temp)
    dt_min = calendar_start + td ! Start date time of each file
    td = idt - dt_min
    time_idx = int(td%total_seconds() / del_t) + 1
    print*, kt, 'idt=', idt, 'time_idx=', time_idx
    ! Read vorticity at the specified level
    call get_xy_from_xyzt(nc_file_name, vort_name, lvl_idx, time_idx, vor)
    vor_ft(:, :, kt) = vor(:, nlats-1:0:-1)
    call apply_mask_2d(vor_ft(:, :, kt), nlons-1, nlats-1, land_mask)

    ! Read sea level pressure
    call make_nc_file_name(nc_file_name, datadir, prefix_sfc, &
                         & idt%year, idt%month, psea_name)
    call get_xy_from_xyt(nc_file_name, psea_name, time_idx, psea)
    psea_ft(:, :, kt) = psea(:, nlats-1:0:-1)

    if (kt > 1 .and. mod(kt, steer_nt) == 0) then
      ! Read u- and v-winds
      do kt2 = 1, steer_nt
        call make_nc_file_name(nc_file_name, datadir, prefix_lvl, &
                             & idt_pair(kt2)%year, idt_pair(kt2)%month, u_name)
        call get_xyz_from_xyzt(nc_file_name, u_name, time_idx, &
                             & steer_idx_top, nsteer_lvl, u(:, :, :, kt2))
        u_ft(:, :, :, kt+kt2-steer_nt) = u(:, nlats-1:0:-1, :, kt2)

        call make_nc_file_name(nc_file_name, datadir, prefix_lvl, &
                             & idt_pair(kt2)%year, idt_pair(kt2)%month, v_name)
        call get_xyz_from_xyzt(nc_file_name, v_name, time_idx, &
                             & steer_idx_top, nsteer_lvl, v(:, :, :, kt2))
        v_ft(:, :, :, kt+kt2-steer_nt) = v(:, nlats-1:0:-1, :, kt2)
      enddo
    end if

    idt_pair(1) = idt
    idt = idt + timedelta(hours=del_t / 3600)
    idt_pair(2) = idt
  enddo
  print*, minval(vor_ft), maxval(vor_ft)
  print*, minval(u_ft), maxval(u_ft)
  print*, minval(v_ft), maxval(v_ft)
  print*, minval(psea_ft), maxval(psea_ft)

  ! print*, shape(vor), minval(vor(:, :, 0)), maxval(vor(:, :, 0))
  ! print*, shape(vor), minval(vor), maxval(vor)

  !! print*, shape(u), minval(u(:, :, :, 0)), maxval(u(:, :, :, 0))

  !write(nc_file_name, '(A,A,A,I4.4,A,I2.2,A,A,A)') trim(datadir), '/', &
  !                                               & 'era5.an.pl.', year_start, '.', &
  !                                               & month_start, '.', &
  !                                               & trim(v_name), '.nc'
  !call get_data_4d(nc_file_name, v_name, time_idx, ntime, &
  !  & steer_idx_top, nsteer_lvl, v)
  !v = v(:, nlats-1:0:-1, :, :)

  !write(nc_file_name, '(A,A,A,I4.4,A,I2.2,A,A,A)') trim(datadir), '/', &
  !                                               & 'era5.an.sfc.', year_start, '.', &
  !                                               & month_start, '.', &
  !                                               & trim(psea_name), '.nc'
  !call get_data_3d(nc_file_name, psea_name, time_idx, ntime, psea)
  !psea = psea(:, nlats-1:0:-1, :)
  !psea = 1e-2 * psea

  ! print*, shape(v), minval(v), maxval(v)
  ! print*, shape(psea), minval(psea), maxval(psea)
  ! print*, lvls

  !write (*,*)'nx=', nx, 'ny=', ny,' nt=', nt, 'nz=', nz
  !write (*,*)'nx1=', nx1, 'nx2=', nx2, 'ny1=', ny1, 'ny2=', ny2
  call tracking_main(vor_ft, u_ft, v_ft, psea_ft, &
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
  deallocate(vor_ft)
  deallocate(u_ft)
  deallocate(v_ft)
  deallocate(psea_ft)
  deallocate(land_mask)

end program main
