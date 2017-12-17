program main
  use datetime_module

  use types, only : wp
  use constants, only : fillval, steer_nt, nmax, rkilo
  use params, only : get_config_params, dbg,                                  &
    & datadir, outdir, prefix_sfc, prefix_lvl,                                &
    & year_start, month_start, day_start, hour_start,                         &
    & year_end, month_end, day_end, hour_end,                                 &
    & vort_name, u_name, v_name, psea_name, land_name,                        &
    & vor_lvl, steer_lvl_btm, steer_lvl_top,                                  &
    & nx1, nx2, ny1, ny2,                                                     &
    & smth_type, proj, steering_type
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
  character(len=256)                          :: fname_bin
  character(len=256)                          :: fname_vormaxloc
  real(wp)                                    :: lonin
  real(wp)                                    :: latin
  real(wp)                                    :: del_t
  real(wp)                                    :: time_step_s
  integer                                     :: nt_per_file
  real(wp)                                    :: lon0
  real(wp)                                    :: lat0
  ! Coordinate arrays
  integer                                     :: ntime, nlvls, nlats, nlons
  integer                                     :: ny, nx
  integer                                     :: ny12, nx12
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
  !real(wp)         , allocatable              :: vor_ft(:, :, :)
  !real(wp)         , allocatable              :: psea_ft(:, :, :)
  !real(wp)         , allocatable              :: u_ft(:, :, :, :)
  !real(wp)         , allocatable              :: v_ft(:, :, :, :)
  real(wp)         , allocatable              :: land_mask(:, :)

  ! Local arrays
  real(wp)         , allocatable              :: vor_smth(:, :)
  integer(4)       , allocatable              :: vor_part(:, :)
  real(wp)         , allocatable              :: dummy(:, :)
  real(wp)         , allocatable              :: mlat(:), mlon(:)
  integer(4)       , allocatable              :: mi(:), mj(:)
  real(wp)         , allocatable              :: max_vor(:)
  real(wp)         , allocatable              :: minlat(:), minlon(:)
  real(wp)         , allocatable              :: z_min(:)
  real(wp)         , allocatable              :: z_min_size(:)
  real(wp)         , allocatable              :: s_part(:)
  real(wp)         , allocatable              :: mtype(:)
  real(wp)         , allocatable              :: u_vor_f(:)
  real(wp)         , allocatable              :: v_vor_f(:)

  ! Local variables
  ! work
  integer                                     :: n_min, n_max
  ! Indices
  integer                                     :: lvl_idx
  integer                                     :: steer_idx_btm, steer_idx_top
  integer                                     :: nsteer_lvl
  integer                                     :: kt, kt2
  integer                                     :: i, j
  integer                                     :: i_max

  ! Time and date variables
  type(datetime)                              :: cal_start
  type(datetime)                              :: dt_start, dt_end
  type(datetime)                              :: dt_min ! within a file
  type(timedelta)                             :: td
  integer                                     :: time_idx
  type(datetime)                              :: idt
  type(datetime)                              :: idt_pair(steer_nt)

  ! IO units
  integer                         , parameter :: fh_bin = 998
  integer                         , parameter :: fh_maxloc = 997


  ! Store dimension names in one array
  DIM_NAMES(1) = trim(REC_NAME)
  DIM_NAMES(2) = trim(LVL_NAME)
  DIM_NAMES(3) = trim(LAT_NAME)
  DIM_NAMES(4) = trim(LON_NAME)

  call get_config_params()

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
  nx = nlons - 1
  ny = nlats - 1
  if (nx1 == -1) nx1 = 0
  if (nx2 == -1) nx2 = nx
  if (ny1 == -1) ny1 = 0
  if (ny2 == -1) ny2 = ny
  nx12 = nx2 - nx1
  ny12 = ny2 - ny1

  ! Time & calendar
  allocate(time_temp(0:nt_per_file-1))
  call get_time(nc_file_name, REC_NAME, time_temp, time_step_s, cal_start)

  del_t = (time_temp(1) - time_temp(0)) * time_step_s
  print*, 'time_temp', time_temp(0)
  print*, 'del_t', del_t
  td = dt_end - dt_start
  ntime = int(td%total_seconds() / del_t) + 1
  print*, 'ntime=', ntime

  td = timedelta(hours=time_temp(0))
  dt_min = cal_start + td
  td = dt_start - dt_min
  time_idx = int(td%total_seconds() / del_t) + 1 ! time_temp(0) +
  deallocate(time_temp)

  ! Assume space coordinates are the same for all files
  allocate(time(1))
  allocate(lvls(nlvls))
  allocate(lats(0:ny))
  allocate(lons(0:nx))

  call get_coords(nc_file_name, DIM_NAMES, lons, lats, lvls, &
    & time, 1, 1)
  print*,lvls
  lvl_idx = minloc(abs(lvls - vor_lvl), 1)
  steer_idx_btm = minloc(abs(lvls - steer_lvl_btm), 1)
  steer_idx_top = minloc(abs(lvls - steer_lvl_top), 1)
  nsteer_lvl = steer_idx_btm - steer_idx_top + 1
  !lvls = lvls(nlvls:1:-1)
  lats = lats(ny:0:-1)
  ! Calculate grid spacing assuming the grid is uniform
  lon0 = lons(0)
  lat0 = lats(0)
  lonin = lons(1) - lons(0)
  latin = lats(1) - lats(0)

  ! Define input array sizes
  allocate(vor      (0:nx, 0:ny))
  allocate(u        (0:nx, 0:ny, nsteer_lvl, steer_nt))
  allocate(v        (0:nx, 0:ny, nsteer_lvl, steer_nt))
  allocate(psea     (0:nx, 0:ny))
  allocate(land_mask(0:nx, 0:ny))

  vor = fillval
  u = fillval
  v = fillval
  psea = fillval

  print*,'size(vor)', shape(vor)
  print*,'nx', nx, 'ny', ny

  ! Allocate work arrays
  allocate(vor_smth  (nx1:nx2, ny1:ny2      ))
  allocate(vor_part  (nx1:nx2, ny1:ny2      ))
  allocate(dummy     (nx1:nx2, ny1:ny2      ))
  allocate(mlat      (                  nmax))
  allocate(mlon      (                  nmax))
  allocate(max_vor   (                  nmax))
  allocate(s_part    (                  nmax))
  allocate(mtype     (                  nmax))
  allocate(minlat    (                  nmax))
  allocate(minlon    (                  nmax))
  allocate(z_min     (                  nmax))
  allocate(z_min_size(                  nmax))
  allocate(mi        (                  nmax));      mi = 0
  allocate(mj        (                  nmax));      mj = 0
  allocate(u_vor_f   (                  nmax)); u_vor_f = fillval
  allocate(v_vor_f   (                  nmax)); v_vor_f = fillval

  write(nc_file_name, '(A,A,A,A)') trim(datadir), '/', trim(land_name), '.nc'
  call get_data_2d(nc_file_name, land_name, land_mask)
  land_mask = land_mask(:, ny:0:-1)

  do kt = 1, ntime ! including both start and end dates
    call make_nc_file_name(nc_file_name, datadir, prefix_lvl, &
                         & idt%year, idt%month, vort_name)
    call get_dims(nc_file_name, DIM_NAMES, nt_per_file, nlvls, nlats, nlons)
    allocate(time_temp(0:nt_per_file-1))
    call get_time(nc_file_name, REC_NAME, time_temp, time_step_s, cal_start)
    del_t = (time_temp(1) - time_temp(0)) * time_step_s
    td = timedelta(hours=time_temp(0))
    deallocate(time_temp)
    dt_min = cal_start + td ! Start date time of each file
    td = idt - dt_min
    time_idx = int(td%total_seconds() / del_t) + 1
    print*, kt, 'idt=', idt, 'time_idx=', time_idx
    ! Read vorticity at the specified level
    call get_xy_from_xyzt(nc_file_name, vort_name, lvl_idx, time_idx, vor)
    vor(:, :) = vor(:, ny:0:-1)
    call apply_mask_2d(vor, nx, ny, land_mask)

    ! Read sea level pressure
    call make_nc_file_name(nc_file_name, datadir, prefix_sfc, &
                         & idt%year, idt%month, psea_name)
    call get_xy_from_xyt(nc_file_name, psea_name, time_idx, psea)
    psea(:, :) = 1e-2 * psea(:, ny:0:-1)

    if (kt > 1 .and. mod(kt, steer_nt) == 0) then
      ! TODO: ensure all times are read in
      ! Read u- and v-winds
      do kt2 = 1, steer_nt
        call make_nc_file_name(nc_file_name, datadir, prefix_lvl, &
                             & idt_pair(kt2)%year, idt_pair(kt2)%month, u_name)
        call get_xyz_from_xyzt(nc_file_name, u_name, time_idx, &
                             & steer_idx_top, nsteer_lvl, u(:, :, :, kt2))

        call make_nc_file_name(nc_file_name, datadir, prefix_lvl, &
                             & idt_pair(kt2)%year, idt_pair(kt2)%month, v_name)
        call get_xyz_from_xyzt(nc_file_name, v_name, time_idx, &
                             & steer_idx_top, nsteer_lvl, v(:, :, :, kt2))
      enddo
      u(:, :, :, :) = u(:, ny:0:-1, :, :)
      v(:, :, :, :) = v(:, ny:0:-1, :, :)
    endif

    if (smth_type == 1) then
      call smth(vor(0:nx, 0:ny), nx, ny, vor_smth(nx1:nx2, ny1:ny2))
    elseif (smth_type == 2) then
      call smth_r(vor(0:nx, 0:ny), nx, ny, lons(0:nx), lats(0:ny),            &
                & vor_smth(nx1:nx2, ny1:ny2))
    else
      ! No smoothing
      vor_smth(nx1:nx2, ny1:ny2) = vor(nx1:nx2, ny1:ny2)
    endif


    write(fname_bin, '(A,A,A,I4.4,A)') trim(outdir), '/',                     &
                                     & 'vor_out_', kt, '.dat'
    open(unit=fh_bin, file=fname_bin, form='unformatted', access='sequential')
    write(unit=fh_bin) vor(nx1:nx2,ny1:ny2)
    write(unit=fh_bin) vor_smth(nx1:nx2,ny1:ny2)

    write (*,'(A,I4.4)') 'Detecting vortex at kt = ', kt

    call vor_partition(vor_smth(nx1:nx2, ny1:ny2),                            &
                     & nx12, ny12,                                            &
                     & mlat(:), mlon(:),                                      &
                     & max_vor(:), mtype(:), n_max,                           &
                     & lats(ny1:ny2), lons(nx1:nx2),                          &
                     & vor_part(nx1:nx2, ny1:ny2), s_part(:))
    if (n_max >= 1) then
      call min_z(psea(nx1:nx2, ny1:ny2),                                      &
               & nx12, ny12,                                                  &
               & minlat(:), minlon(:),                                        &
               & z_min(:), n_min, lats(ny1:ny2), lons(nx1:nx2),               &
               & z_min_size(:))
    else
      n_min = 0
    endif

    if (maxval(mtype(:)) >= 1) then
      call synop_check(mlon(:), mlat(:), n_max,                               &
                     & minlon(:), minlat(:), n_min, mtype(:))
    endif

    ! Save vor_part to a dummy array and write it to unformatted output
    dummy = fillval
    ! dummy(nx1, ny2)=-1.
    do j = ny1, ny2
      do i = nx1, nx2
        if (vor_part(i, j) /= 0.) then
          dummy(i, j) = vor_part(i, j)
        endif
      enddo
    enddo
    write(unit=fh_bin) dummy(nx1:nx2, ny1:ny2)


    write(fname_vormaxloc, '(A,A,A,I4.4,A)') trim(outdir), '/',               &
                                           & 'vormax_loc_', kt, '.txt'
    open(unit=fh_maxloc, file=fname_vormaxloc, form='formatted')
    do i_max = 1, n_max
      mi(i_max) = nint((mlon(i_max) - lon0) / lonin)
      mj(i_max) = nint((mlat(i_max) - lat0) / latin)
      if (proj == 1) then
        write(unit=fh_maxloc, fmt=*) mlon(i_max),                             &
                                  & mlat(i_max),                              &
                                  & max_vor(i_max)*rkilo,                     &
                                  & nint(s_part(i_max)),                      &
                                  & mtype(i_max)
        if (dbg) then
          write(*, *) mlon(i_max),                                            &
                    & mlat(i_max),                                            &
                    & max_vor(i_max)*rkilo,                                   &
                    & nint(s_part(i_max)),                                    &
                    & mtype(i_max)
        endif
      elseif (proj == 2) then
        write(unit=fh_maxloc, fmt=*) mlon(i_max)/rkilo,                       &
                                  & mlat(i_max)/rkilo,                        &
                                  & max_vor(i_max)*rkilo,                     &
                                  & nint(s_part(i_max)),                      &
                                  & mtype(i_max)
        if (dbg) then
          write(*, *) mlon(i_max)/rkilo,                                      &
                    & mlat(i_max)/rkilo,                                      &
                    & max_vor(i_max)*rkilo,                                   &
                    & nint(s_part(i_max)),                                    &
                    & mtype(i_max)
        endif
      endif
      !
      ! ----- calculate steering wind -----!
      !
      if (steering_type == 1) then
        write(*, *) "NotImplementedError"; stop
        !call steering_wind_f(u,v,&
        !     &levs,nx,ny,nz,nt,kt,&
        !     &mi(i_max,kt),mj(i_max,kt),&
        !     &u_vor_f(i_max,kt),v_vor_f(i_max,kt))

        !call steering_wind_b(u,v,&
        !     &levs,nx,ny,nz,nt,kt,&
        !     &mi(i_max,kt),mj(i_max,kt),&
        !     &u_vor_b(i_max,kt),v_vor_b(i_max,kt))
      elseif (steering_type == 2) then
        if (kt < ntime) then
          ! Forward
          call steering_wind_r(u, v,                                          &
                             & lvls(1:nsteer_lvl), lons, lats,                &
                             & nx, ny, nsteer_lvl,                            &
                             & steer_nt, 1, 2,                                &
                             & mi(i_max), mj(i_max),                          &
                             & u_vor_f(i_max), v_vor_f(i_max))
        endif
        !if (kt > 1) then
        !  ! Backward
        !  call steering_wind_r(u, v,     &
        !       &lvls, lons,lats,nx,ny,nz,nt,kt,kt-1,&
        !       &mi(i_max,kt),mj(i_max,kt),&
        !       &u_vor_b(i_max,kt),v_vor_b(i_max,kt))
        !endif
      endif

    enddo ! i_max loop
    close(unit=fh_maxloc)

    !---- output steeering wind ----!
    dummy = fillval
    do i_max=1, n_max
      dummy(mi(i_max), mj(i_max)) = u_vor_f(i_max)
    end do
    write(unit=fh_bin) dummy(nx1:nx2, ny1:ny2)
    dummy = fillval
    do i_max=1, n_max
      dummy(mi(i_max), mj(i_max)) = v_vor_f(i_max)
    end do
    write(unit=fh_bin) dummy(nx1:nx2, ny1:ny2)

    !     dummy=fillval
    
    !     do i_max=1,n_max(kt)
    !       dummy(mi(i_max,kt),mj(i_max,kt))=-u_vor_b(i_max,kt)
    !     end do
    
    !     write (12)dummy(nx1:nx2,ny1:ny2)
    
    !     dummy=fillval
    
    !     do i_max=1,n_max(kt)
    !       dummy(mi(i_max,kt),mj(i_max,kt))=-v_vor_b(i_max,kt)
    !     end do
    
    !     write (12)dummy(nx1:nx2,ny1:ny2)
    ! SLP output
    write(unit=fh_bin) psea(nx1:nx2, ny1:ny2)

    close(unit=fh_bin)
    idt_pair(1) = idt
    idt = idt + timedelta(hours=del_t / time_step_s)
    idt_pair(2) = idt
  enddo

  !write (*,*)'nx=', nx, 'ny=', ny,' nt=', nt, 'nz=', nz
  !write (*,*)'nx1=', nx1, 'nx2=', nx2, 'ny1=', ny1, 'ny2=', ny2
  !call tracking_main(vor_ft, u_ft, v_ft, psea_ft, &
  !  & nx-1, ny-1, nsteer_lvl, lvls(1:nsteer_lvl), ntime, &
  !  & lons, lats, lonin, latin, del_t)

  deallocate(time)
  deallocate(lvls)
  deallocate(lats)
  deallocate(lons)
  deallocate(vor)
  deallocate(u)
  deallocate(v)
  deallocate(psea)
  deallocate(land_mask)

  deallocate(vor_smth  )
  deallocate(vor_part  )
  deallocate(dummy     )
  deallocate(mlat      )
  deallocate(mlon      )
  deallocate(max_vor   )
  deallocate(s_part    )
  deallocate(mtype     )
  deallocate(minlat    )
  deallocate(minlon    )
  deallocate(z_min     )
  deallocate(z_min_size)
  deallocate(mi        )
  deallocate(mj        )
  deallocate(u_vor_f   )
  deallocate(v_vor_f   )

end program main
