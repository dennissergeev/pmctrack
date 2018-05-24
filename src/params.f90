module params
  use datetime_module

  use types, only: wp
  use constants, only: fillval
  use utils, only: check_exist

  implicit none

  character(len=256), protected :: config_file
  type(datetime)    , protected :: dt_start
  type(datetime)    , protected :: dt_end
  character(len=256), protected :: datadir
  character(len=256), protected :: outdir
  character(len=256), protected :: vort_name
  character(len=256), protected :: u_name
  character(len=256), protected :: v_name
  character(len=256), protected :: psea_name
  character(len=256), protected :: land_name
  character(len=256), protected :: prefix_lvl
  character(len=256), protected :: prefix_sfc
  integer           , protected :: vor_lvl
  integer           , protected :: steer_lvl_btm, steer_lvl_top
  integer           , protected :: tfreq
  integer           , protected :: proj
  integer           , protected :: vert_grid
  real   (wp)       , protected :: lon1, lon2, lat1, lat2
  integer           , protected :: nx1, nx2, ny1, ny2
  !integer          , protected:: nt
  !real(wp)         , protected:: del_t
  ! parameters for masking land
  real   (wp)       , protected :: halo_r
  ! parameter for smoothing of vorticity
  integer           , protected :: smth_type
  integer           , protected :: nsmth_x, nsmth_y
  real   (wp)       , protected :: r_smth
  ! parameter for detecting vortex
  real   (wp)       , protected :: zeta_max0, zeta_min0
  real   (wp)       , protected :: int_zeta_min0, gamma
  ! parameter for excluding the synoptic scale disturbances
  real   (wp)       , protected :: d_cf_min, size_synop
  real   (wp)       , protected :: del_psea_min, distance_ec
  ! parameter for calculating steering winds
  integer           , protected :: steering_type
  integer           , protected :: n_steering_x, n_steering_y
  real   (wp)       , protected :: r_steering
  ! parameter for linking vortex
  integer           , protected :: track_type
  real   (wp)       , protected :: del_lon, del_lat, del_r
  integer           , protected :: merge_opt
  ! Output flags
  integer           , protected :: vor_out_on
  ! Debug flag
  logical           , protected :: dbg

contains
  subroutine get_config_file_name()
    implicit none
    integer :: narg

    narg = iargc()

    if (narg == 0) then
      config_file = "settings.conf"
    else
      call getarg(1, config_file)
    endif
    call check_exist( trim( config_file ) )
  end subroutine get_config_file_name


  subroutine copy_config_file()
      call system('cp '//trim(config_file)//' '//trim(outdir))
  end subroutine copy_config_file


  subroutine get_config_params()

    use constants, only : fh_conf

    implicit none

    ! Local variables
    character(len=256)            :: buffer
    character(len=256)            :: label
    integer                       :: line
    integer                       :: ios
    integer                       :: sep
    character(len=14)             :: dt_start_str, dt_end_str

#ifdef debug
    dbg = .true.
#else
    dbg = .false.
#endif

    ios = 0
    line = 0

    lon1 = fillval
    lon2 = fillval
    lat1 = fillval
    lat2 = fillval

    ! Default time frequency is one
    tfreq = 1

    call get_config_file_name()
    open (fh_conf, file=trim(config_file), form='formatted', status='old', &
      &   iostat=ios, action='read')
    if (ios == 0) then
      if (dbg) write(*, *) '--------------------------------'
      if (dbg) write(*, *) 'Reading configuration parameters'
      do while (ios == 0)
        read(fh_conf, '(A)', iostat=ios) buffer
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
          case('dt_start'); read(buffer, *, iostat=ios) dt_start_str; if (dbg) write(*, *) dt_start_str
          case('dt_end'); read(buffer, *, iostat=ios) dt_end_str; if (dbg) write(*, *) dt_end_str
          case('vor_lvl'); read(buffer, *, iostat=ios) vor_lvl; if (dbg) write(*, *) vor_lvl
          case('steer_lvl_btm'); read(buffer, *, iostat=ios) steer_lvl_btm; if (dbg) write(*, *) steer_lvl_btm
          case('steer_lvl_top'); read(buffer, *, iostat=ios) steer_lvl_top; if (dbg) write(*, *) steer_lvl_top
          case('datadir'); read(buffer, *, iostat=ios) datadir; if (dbg) write(*, *) datadir
          case('outdir'); read(buffer, *, iostat=ios) outdir; if (dbg) write(*, *) outdir
          case('vort_name'); read(buffer, *, iostat=ios) vort_name; if (dbg) write(*, *) vort_name
          case('u_name'); read(buffer, *, iostat=ios) u_name; if (dbg) write(*, *) u_name
          case('v_name'); read(buffer, *, iostat=ios) v_name; if (dbg) write(*, *) v_name
          case('psea_name'); read(buffer, *, iostat=ios) psea_name; if (dbg) write(*, *) psea_name
          case('land_name'); read(buffer, *, iostat=ios) land_name; if (dbg) write(*, *) land_name
          case('prefix_lvl'); read(buffer, *, iostat=ios) prefix_lvl; if (dbg) write(*, *) prefix_lvl
          case('prefix_sfc'); read(buffer, *, iostat=ios) prefix_sfc; if (dbg) write(*, *) prefix_sfc
          case('tfreq'); read(buffer, *, iostat=ios) tfreq; if (dbg) write(*, *) tfreq
          case('proj'); read(buffer, *, iostat=ios) proj; if (dbg) write(*, *) proj
          case('vert_grid'); read(buffer, *, iostat=ios) vert_grid; if (dbg) write(*, *) vert_grid
          case('lon1'); read(buffer, *, iostat=ios) lon1; if (dbg) write(*, *) lon1
          case('lon2'); read(buffer, *, iostat=ios) lon2; if (dbg) write(*, *) lon2
          case('lat1'); read(buffer, *, iostat=ios) lat1; if (dbg) write(*, *) lat1
          case('lat2'); read(buffer, *, iostat=ios) lat2; if (dbg) write(*, *) lat2
          case('smth_type'); read(buffer, *, iostat=ios) smth_type; if (dbg) write(*, *) smth_type
          case('nsmth_x'); read(buffer, *, iostat=ios) nsmth_x; if (dbg) write(*, *) nsmth_x
          case('nsmth_y'); read(buffer, *, iostat=ios) nsmth_y; if (dbg) write(*, *) nsmth_y
          case('r_smth'); read(buffer, *, iostat=ios) r_smth; if (dbg) write(*, *) r_smth
          case('halo_r'); read(buffer, *, iostat=ios) halo_r; if (dbg) write(*, *) halo_r
          case('zeta_max0'); read(buffer, *, iostat=ios) zeta_max0; if (dbg) write(*, *) zeta_max0
          case('zeta_min0'); read(buffer, *, iostat=ios) zeta_min0; if (dbg) write(*, *) zeta_min0
          case('int_zeta_min0'); read(buffer, *, iostat=ios) int_zeta_min0; if (dbg) write(*, *) int_zeta_min0
          case('gamma'); read(buffer, *, iostat=ios) gamma; if (dbg) write(*, *) gamma
          case('d_cf_min'); read(buffer, *, iostat=ios) d_cf_min; if (dbg) write(*, *) d_cf_min
          case('size_synop'); read(buffer, *, iostat=ios) size_synop; if (dbg) write(*, *) size_synop
          case('distance_ec'); read(buffer, *, iostat=ios) distance_ec; if (dbg) write(*, *) distance_ec
          case('del_psea_min'); read(buffer, *, iostat=ios) del_psea_min; if (dbg) write(*, *) del_psea_min
          case('steering_type'); read(buffer, *, iostat=ios) steering_type; if (dbg) write(*, *) steering_type
          case('n_steering_x'); read(buffer, *, iostat=ios) n_steering_x; if (dbg) write(*, *) n_steering_x
          case('n_steering_y'); read(buffer, *, iostat=ios) n_steering_y; if (dbg) write(*, *) n_steering_y
          case('r_steering'); read(buffer, *, iostat=ios) r_steering; if (dbg) write(*, *) r_steering
          case('track_type'); read(buffer, *, iostat=ios) track_type; if (dbg) write(*, *) track_type
          case('del_lon'); read(buffer, *, iostat=ios) del_lon; if (dbg) write(*, *) del_lon
          case('del_lat'); read(buffer, *, iostat=ios) del_lat; if (dbg) write(*, *) del_lat
          case('del_r'); read(buffer, *, iostat=ios) del_r; if (dbg) write(*, *) del_r
          case('merge_opt'); read(buffer, *, iostat=ios) merge_opt; if (dbg) write(*, *) merge_opt
          case('vor_out_on'); read(buffer, *, iostat=ios) vor_out_on; if (dbg) write(*, *) vor_out_on
          case default
            if (index (trim(label), "#") /= 1) then
              write(*, *) 'ConfigParseWarning: Skipping invalid line', line
            end if
          end select
        end if
      end do
      if (dbg) write(*, *) '--------------------------------'
    end if
    close(fh_conf)

    dt_start = strptime(trim(dt_start_str), '%Y_%m_%d_%H%M')
    !dt_start%minute = 0
    dt_start%second = 0
    dt_start%millisecond = 0
    if (.not. dt_start%isValid()) then
      write(*, *) 'Start date ', dt_start, ' is not valid'
      stop
    endif
    dt_end = strptime(trim(dt_end_str), '%Y_%m_%d_%H%M')
    !dt_end%minute = 0
    dt_end%second = 0
    dt_end%millisecond = 0
    if (.not. dt_end%isValid()) then
      write(*, *) 'End date ', dt_end, ' is not valid'
      stop
    endif

    if (tfreq < 1) then
      write(*, *) 'ConfigError: tfreq must be a positive integer; given:', tfreq
    endif

    if(proj==2)then
      write (*,*)'Cartesian coordinate'
    elseif(proj==1)then
      write (*,*)'Geographical coordinate'
    else
      write (*,*)'Coordinate system not supported'
    end if

    if(vert_grid==1)then
      write (*,*)'pressure levels'
    elseif(vert_grid==2)then
      write (*,*)'height levels'
    else
      write (*,*)'vertical coordinate system not supported'
    end if

    if(smth_type==1)then
      write (*,*)'Smoothing in lat-lon rectangle'
    elseif(smth_type==2)then
      write (*,*)'Smooting in radius'
    end if

    if(steering_type==1)then
      write (*,*)'calculate steering wind in latlon coordinate'
    elseif(steering_type==2)then
      write (*,*)'calculate steering wind in radius'
    end if

    if(track_type==1)then
      write (*,*)'Use del_lon and del_lat to connect vortices'

    elseif(track_type==2)then
      write (*,*)'Use radius to connect vortices'
    end if
  end subroutine get_config_params


  subroutine set_lonlat_bounds_auto(nx, ny, lons, lats)

    integer    , intent(in)    :: nx, ny
    real   (wp), intent(in)    :: lons(0:nx      )
    real   (wp), intent(in)    :: lats(      0:ny)

    if (lon1 == fillval) lon1 = lons(0)
    if (lon2 == fillval) lon2 = lons(nx)
    if (lat1 == fillval) lat1 = lats(0)
    if (lat2 == fillval) lat2 = lats(ny)

    nx1 = minloc(abs(lons-lon1), 1) - 1 ! -1 because arrays start at 0
    nx2 = minloc(abs(lons-lon2), 1) - 1
    ny1 = minloc(abs(lats-lat1), 1) - 1
    ny2 = minloc(abs(lats-lat2), 1) - 1

    ! print*, lon1, lon2, lat1, lat2
    ! print*, nx1, nx2, ny1, ny2
    ! stop

  end subroutine set_lonlat_bounds_auto
end module params
