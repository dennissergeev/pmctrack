program main
  use datetime_module

  use types, only: wp
  use constants, only: ifillval, fillval, missval,                             &
    & steer_nt, nmax, pmax, rkilo, fh_bin, fh_maxloc
  use params, only: get_config_params, copy_config_file,                       &
    & set_lonlat_bounds_auto, dbg,                                             &
    & datadir, outdir, fname_sfc, fname_lvl, dt_start, dt_end,                 &
    & t_dim,                                                                   &
    & vort_name, u_name, v_name, psea_name, land_name,                         &
    & vor_lvl, steer_lvl_btm, steer_lvl_top, tfreq,                            &
    & nx1, nx2, ny1, ny2,                                                      &
    & land_mask_type, halo_r, smth_type, proj, steering_type, track_type,      &
    & vor_out_on
  use nc_io, only: get_dims, get_time, get_coords,                             &
    & get_xy_from_xyzt, get_xy_from_xyt, get_xyz_from_xyzt,                    &
    & get_data_2d, get_units
  use utils, only: extend_mask_2d, make_nc_file_name, write_vortrack,          &
    & makedirs_p, lower

  implicit none

  character(len=256)                          :: nc_file_name
  character(len=256)                          :: fname_bin
  character(len=256)                          :: fname_vormaxloc
  character(len=256)                          :: psea_units
  real     (wp)                               :: lonin
  real     (wp)                               :: latin
  real     (wp)                               :: del_t
  real     (wp)                               :: data_del_t
  real     (wp)                               :: time_step_s
  integer                                     :: nt_per_file
  real     (wp)                               :: lon0
  real     (wp)                               :: lat0
  ! Coordinate arrays
  integer                                     :: ntime, nlvls, nlats, nlons
  integer                                     :: ny, nx
  integer                                     :: ny12, nx12
  integer          , allocatable              :: time_temp(:)
  integer          , allocatable              :: time(:)
  real     (wp)    , allocatable              :: lvls(:)
  real     (wp)    , allocatable              :: lats(:)
  real     (wp)    , allocatable              :: lons(:)
  ! Data arrays
  real     (wp)    , allocatable              :: vor(:, :)
  real     (wp)    , allocatable              :: psea(:, :)
  real     (wp)    , allocatable              :: u(:, :, :, :)
  real     (wp)    , allocatable              :: v(:, :, :, :)
  real     (wp)    , allocatable              :: land_mask(:, :)

  ! Local arrays
  real     (wp)    , allocatable              :: vor_smth(:, :)
  integer          , allocatable              :: vor_part(:, :)
  real     (wp)    , allocatable              :: dummy(:, :)
  real     (wp)    , allocatable              :: mlat(:), mlon(:)
  real     (wp)    , allocatable              :: mlat_prev(:), mlon_prev(:)
  integer          , allocatable              :: mi(:), mj(:)
  real     (wp)    , allocatable              :: max_vor(:)
  real     (wp)    , allocatable              :: minlat(:), minlon(:)
  real     (wp)    , allocatable              :: z_min(:)
  integer          , allocatable              :: z_min_size(:)
  real     (wp)    , allocatable              :: s_part(:)
  integer          , allocatable              :: mtype(:)
  real     (wp)    , allocatable              :: u_vor_f(:)
  real     (wp)    , allocatable              :: v_vor_f(:)
  real     (wp)    , allocatable              :: u_vor_f_prev(:)
  real     (wp)    , allocatable              :: v_vor_f_prev(:)
  integer          , allocatable              :: vor_merge(:)
  integer          , allocatable              :: vor_index(:)
  integer          , allocatable              :: merged_count(:)
  real     (wp)    , allocatable              :: vortex(:, :)
  integer          , allocatable              :: vor_merge_num(:)

  ! Local scalars
  real                                        :: factor
  ! work
  integer                                     :: n_min
  integer                                     :: n_max
  integer                                     :: n_max_prev
  integer                                     :: vor_num
  ! Indices
  integer                                     :: lvl_idx
  integer                                     :: steer_idx_btm, steer_idx_top
  integer                                     :: nsteer_lvl
  integer                                     :: kt, kt2
  integer                                     :: i, j
  integer                                     :: i_max
  integer                                     :: i_vor_num
  integer                                     :: ix, iy
  ! Time and date variables
  type   (datetime)                           :: cal_start
  type   (datetime)                           :: dt_min ! within a file
  type   (timedelta)                          :: td
  integer                                     :: time_idx
  type   (datetime)                           :: idt
  type   (datetime)                           :: idt_pair(steer_nt)


  ! Read configs from settings.conf file
  call get_config_params()

  ! Create output directory if it does not exist
  call makedirs_p(outdir, overwrite=.true.)

  ! Copy config file to the output directory
  call copy_config_file()

  ! Get dimensions from the vorticity file using first year and first month
  ! Assume all the other files are organised in the same way
  idt = dt_start
  call make_nc_file_name(nc_file_name, datadir, fname_lvl, &
                       & idt%year, idt%month, idt%day, vort_name)

  call get_dims(nc_file_name, nt_per_file, nlvls, nlats, nlons)
  nx = nlons - 1
  ny = nlats - 1

  ! Time & calendar
  allocate(time_temp(0:nt_per_file-1))
  call get_time(nc_file_name, t_dim, time_temp, time_step_s, cal_start)

  ! Time resolution of the input data
  data_del_t = (time_temp(1) - time_temp(0)) * time_step_s
  ! Time step of tracking
  del_t = data_del_t * tfreq
  td = dt_end - dt_start
  ntime = int(td%total_seconds() / del_t) + 1

  td = timedelta(hours=time_temp(0))
  dt_min = cal_start + td
  td = dt_start - dt_min
  time_idx = int(td%total_seconds() / data_del_t) + 1
  deallocate(time_temp)

  ! Assume space coordinates are the same for all files
  allocate(time(1))
  allocate(lvls(nlvls))
  allocate(lats(0:ny))
  allocate(lons(0:nx))

  call get_coords(nc_file_name, time, lvls, lats, lons, 1, 1)
  deallocate(time)  ! time array is not needed and is handled by get_time() instead

  lvl_idx = minloc(abs(lvls - vor_lvl), 1)
  steer_idx_btm = minloc(abs(lvls - steer_lvl_btm), 1)
  steer_idx_top = minloc(abs(lvls - steer_lvl_top), 1)
  nsteer_lvl = abs(steer_idx_top - steer_idx_btm) + 1
  steer_idx_btm = min(steer_idx_btm, steer_idx_top) ! steer_idx_top is not used below

  lats = lats(ny:0:-1)
  ! Calculate grid spacing assuming the grid is uniform
  lon0 = lons(0)
  lat0 = lats(0)
  lonin = lons(1) - lons(0)
  latin = lats(1) - lats(0)
  ! Calculate sub domain boundaries
  call set_lonlat_bounds_auto(nx, ny, lons, lats)
  nx12 = nx2 - nx1
  ny12 = ny2 - ny1

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
  land_mask = fillval

  ! Allocate work arrays
  allocate(vor_smth     (nx1:nx2, ny1:ny2      ))
  allocate(vor_part     (nx1:nx2, ny1:ny2      )); vor_part = ifillval
  allocate(dummy        (nx1:nx2, ny1:ny2      )); dummy = fillval
  allocate(mlat         (                  nmax)); mlat = fillval
  allocate(mlon         (                  nmax)); mlon = fillval
  allocate(mlat_prev    (                  nmax)); mlat_prev = fillval
  allocate(mlon_prev    (                  nmax)); mlon_prev = fillval
  allocate(max_vor      (                  nmax))
  allocate(s_part       (                  nmax))
  allocate(mtype        (                  nmax)); mtype = ifillval
  allocate(minlat       (                  nmax))
  allocate(minlon       (                  nmax))
  allocate(z_min        (                  nmax))
  allocate(z_min_size   (                  nmax))
  allocate(mi           (                  nmax));      mi = ifillval
  allocate(mj           (                  nmax));      mj = ifillval
  allocate(u_vor_f      (                  nmax)); u_vor_f = fillval
  allocate(v_vor_f      (                  nmax)); v_vor_f = fillval
  allocate(u_vor_f_prev (                  nmax)); u_vor_f_prev = fillval
  allocate(v_vor_f_prev (                  nmax)); v_vor_f_prev = fillval
  allocate(vor_index    (                  pmax)); vor_index = ifillval
  allocate(merged_count (                  pmax)); merged_count = ifillval
  allocate(vor_merge    (                  pmax)); vor_merge = ifillval
  allocate(vor_merge_num(                  pmax));

  vor_num = 0
  vor_merge_num(:) = 1

 if (land_mask_type == 1) then
  nc_file_name = trim(datadir) // '/' // trim(land_name) // '.nc'
  call get_data_2d(nc_file_name, land_name, land_mask)
  land_mask = land_mask(:, ny:0:-1)
 !the land mask is set to 1 if it has a value larger than 1, so for lakes, islands, etc
  where ( land_mask(:, :) > 0 ) land_mask(:, :) = 1

  if (halo_r > 0) then
    call extend_mask_2d(nx, ny, land_mask, lons(0:nx), lats(0:ny),            &
      &                 proj, halo_r)
  endif
 endif

  ! MAIN TIME LOOP ------------------------------------------------------------
  do kt = 1, ntime ! including both start and end dates
    call make_nc_file_name(nc_file_name, datadir, fname_lvl, &
                         & idt%year, idt%month, idt%day, vort_name)
    call get_dims(nc_file_name, nt_per_file, nlvls, nlats, nlons)
    allocate(time_temp(0:nt_per_file-1))
    call get_time(nc_file_name, t_dim, time_temp, time_step_s, cal_start)
    ! Time resolution of the input data
    data_del_t = (time_temp(1) - time_temp(0)) * time_step_s
    ! Time step of tracking
    del_t = data_del_t * tfreq
    td = timedelta(hours=time_temp(0))
    deallocate(time_temp)
    dt_min = cal_start + td ! Start date time of each file
    td = idt - dt_min
    time_idx = int(td%total_seconds() / data_del_t) + 1
    write(*, *) ''
    write(*, *) '============================================================='
    write(*, *) 'kt=', kt, 'idt=', trim(idt%strftime('%Y-%m-%d %H:%M'))
    write(*, *) '============================================================='
    ! Read vorticity at the specified level

    call get_xy_from_xyzt(nc_file_name, vort_name, lvl_idx, time_idx, vor)
    vor(:, :) = vor(:, ny:0:-1)

    if (land_mask_type == 1) then
      ! Apply land mask
      where ( land_mask(:, :) == 1 ) vor(:, :) = missval
    endif

    ! Read sea level pressure
    call make_nc_file_name(nc_file_name, datadir, fname_sfc, &
                         & idt%year, idt%month, idt%day, psea_name)
    call get_xy_from_xyt(nc_file_name, psea_name, time_idx, psea)
    psea_units = get_units(nc_file_name, psea_name)
    if (lower(psea_units) == 'pa') then
      ! Convert to hPa
      factor = 1e-2
    else
      factor = 1.0
    endif
    psea(:, :) = factor * psea(:, ny:0:-1)

    ! Read 2 time steps (forward)
    idt_pair(1) = idt
    idt_pair(2) = idt + timedelta(hours=del_t / time_step_s)
    ! if (kt > 1 .and. mod(kt, steer_nt) == 0) then
    if (kt < ntime) then
      ! Read u- and v-winds
      do kt2 = 1, steer_nt
!        if (idt_pair(2)%month /= idt_pair(1)%month .and. kt2 == 2) then
        if (idt_pair(2)%day /= idt_pair(1)%day .and. kt2 == 2) then
          time_idx = 0
        endif
        call make_nc_file_name(nc_file_name, datadir, fname_lvl, &
                             & idt_pair(kt2)%year, idt_pair(kt2)%month, idt_pair(kt2)%day, u_name)
        call get_xyz_from_xyzt(nc_file_name, u_name, time_idx+(kt2-1)*tfreq, &
                             & steer_idx_btm, nsteer_lvl, u(:, :, :, kt2))

        call make_nc_file_name(nc_file_name, datadir, fname_lvl, &
                             & idt_pair(kt2)%year, idt_pair(kt2)%month, idt_pair(kt2)%day, v_name)
        call get_xyz_from_xyzt(nc_file_name, v_name, time_idx+(kt2-1)*tfreq, &
                             & steer_idx_btm, nsteer_lvl, v(:, :, :, kt2))

      enddo
      u(:, :, :, :) = u(:, ny:0:-1, :, :)
      v(:, :, :, :) = v(:, ny:0:-1, :, :)
    endif
    ! END OF INPUT

    if (smth_type == 1) then
      call smth(vor(0:nx, 0:ny), nx, ny, vor_smth(nx1:nx2, ny1:ny2))
    elseif (smth_type == 2) then
      call smth_r(vor(0:nx, 0:ny), nx, ny, lons(0:nx), lats(0:ny),            &
                & vor_smth(nx1:nx2, ny1:ny2))
    else
      ! No smoothing
      vor_smth(nx1:nx2, ny1:ny2) = vor(nx1:nx2, ny1:ny2)
    endif


    if (vor_out_on==1) then
      write(fname_bin, '(A,A,A,A,A)') trim(outdir), '/',                     &
                        & 'vor_out_', trim(idt%strftime('%Y%m%d%H%M')), '.dat'
      open(unit=fh_bin, file=fname_bin, form='unformatted', access='sequential')
      write(unit=fh_bin) vor(nx1:nx2,ny1:ny2)
      write(unit=fh_bin) vor_smth(nx1:nx2,ny1:ny2)
    endif

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

    write(fname_vormaxloc, '(A,A,A,A,A)') trim(outdir), '/',               &
                      & 'vormax_loc_', trim(idt%strftime('%Y%m%d%H%M')), '.txt'
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
      if (kt < ntime) then
        if (steering_type == 1) then
          ! Forward
          call steering_wind_f(u, v,                                          &
                             & lvls(1:nsteer_lvl),                            &
                             & nx, ny, nsteer_lvl, steer_nt,                  &
                             & mi(i_max), mj(i_max),                          &
                             & u_vor_f(i_max), v_vor_f(i_max))
        elseif (steering_type == 2) then
          ! Forward
          call steering_wind_r(u, v,                                          &
                             & lvls(1:nsteer_lvl), lons, lats,                &
                             & nx, ny, nsteer_lvl, steer_nt,                  &
                             & mi(i_max), mj(i_max),                          &
                             & u_vor_f(i_max), v_vor_f(i_max))
        endif
      endif
    enddo ! i_max loop
    close(unit=fh_maxloc)

    if (vor_out_on==1) then
      ! Save vor_part to a dummy array and write it to unformatted output
      ! dummy(nx1, ny2)=-1.
      do j = ny1, ny2
        do i = nx1, nx2
          if (vor_part(i, j) /= 0.) then
            dummy(i, j) = vor_part(i, j)
          endif
        enddo
      enddo
      write(unit=fh_bin) dummy(nx1:nx2, ny1:ny2)
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
    endif

    ! Link vortices
    if (kt > 1) then
      if (track_type == 1) then
        write(*, *) "NotImplementedError"; stop
      elseif (track_type == 2) then
       call link_vort_rad(nx12, ny12, lons(nx1:nx2), lats(ny1:ny2),           &
                        & del_t, mtype, mlon_prev, mlat_prev, mlon, mlat,     &
                        & u_vor_f_prev, v_vor_f_prev, s_part,                 &
                        & vor_part(nx1:nx2, ny1:ny2),                         &
                        & n_max_prev, n_max,                                  &
                        & vor_num, vor_index, vor_merge)
      endif
    else
      vor_num = n_max
      do i_max = 1, n_max
        vor_index(i_max) = i_max
      enddo
    endif

    allocate(vortex     (6, vor_num)); vortex = fillval

    do i_vor_num = 1, vor_num
      if (vor_index(i_vor_num) > 0) then
        if (proj == 1) then
          vortex(1, i_vor_num) = mlon(vor_index(i_vor_num))
          vortex(2, i_vor_num) = mlat(vor_index(i_vor_num))
        elseif (proj == 2) then
          vortex(1, i_vor_num) = mlon(vor_index(i_vor_num)) / rkilo
          vortex(2, i_vor_num) = mlat(vor_index(i_vor_num)) / rkilo
        endif
        vortex(3, i_vor_num) = max_vor(vor_index(i_vor_num)) * rkilo
        vortex(4, i_vor_num) = s_part(vor_index(i_vor_num))
        vortex(5, i_vor_num) = mtype(vor_index(i_vor_num))

        ! Output SLP at the vortex centre
        ! TODO: do it within a radius? or use z_min (but some values are 0)
        ix = minloc(abs(lons-mlon(vor_index(i_vor_num))), 1) - 1
        iy = minloc(abs(lats-mlat(vor_index(i_vor_num))), 1) - 1
        ! vortex(6, i_vor_num) = z_min(vor_index(i_vor_num))
        vortex(6, i_vor_num) = psea(ix, iy)

      end if

      ! --- check the track ---
      ! call check_track(vortex(i_vor_num, 1:3), vortex_flag(i_vor_num))
      ! disabled, because requires full track in time
      ! So it's easier to filter them out when analysing the full output
    enddo

    !------------ vortrack out put ----------------------
    do i_vor_num = 1, vor_num
      if (vor_index(i_vor_num) > 0) then! .and. merged_count(i_vor_num) /= 1) then
        if (vor_merge(i_vor_num) > 0) then
          if (merged_count(i_vor_num) /= 1) then
            vor_merge_num(vor_merge(i_vor_num)) = &
              & vor_merge_num(vor_merge(i_vor_num)) + 1
            merged_count(i_vor_num) = 1
            ! Write current time step to the merged vortex file too
            call write_vortrack(outdir, &
                              & vortex(:, i_vor_num), i_vor_num, 1, idt)
          endif
          call write_vortrack(outdir, &
                            & vortex(:, i_vor_num), vor_merge(i_vor_num), &
                            & vor_merge_num(vor_merge(i_vor_num)), idt)
        else
          call write_vortrack(outdir, vortex(:, i_vor_num), i_vor_num, 1, idt)
        endif
      endif
    enddo
    deallocate(vortex)

    ! Save for the next time step
    n_max_prev = n_max
    mlon_prev = mlon
    mlat_prev = mlat
    u_vor_f_prev = u_vor_f
    v_vor_f_prev = v_vor_f
    ! idt_pair(1) = idt
    idt = idt + timedelta(hours=del_t / time_step_s)
    ! idt_pair(2) = idt
  enddo ! MAIN TIME LOOP ------------------------------------------------------

  deallocate(lvls)
  deallocate(lats)
  deallocate(lons)
  deallocate(vor)
  deallocate(u)
  deallocate(v)
  deallocate(psea)
  deallocate(land_mask)

  deallocate(vor_smth     )
  deallocate(vor_part     )
  deallocate(dummy        )
  deallocate(mlat         )
  deallocate(mlon         )
  deallocate(mlat_prev    )
  deallocate(mlon_prev    )
  deallocate(max_vor      )
  deallocate(s_part       )
  deallocate(mtype        )
  deallocate(minlat       )
  deallocate(minlon       )
  deallocate(z_min        )
  deallocate(z_min_size   )
  deallocate(mi           )
  deallocate(mj           )
  deallocate(u_vor_f      )
  deallocate(v_vor_f      )
  deallocate(u_vor_f_prev )
  deallocate(v_vor_f_prev )
  deallocate(vor_merge    )
  deallocate(vor_index    )
  deallocate(merged_count )

  write(*, *) 'Owari'
end program main
