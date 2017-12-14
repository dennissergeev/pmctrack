module params
  use types, only : wp

  implicit none

  character(len=*), parameter :: CONFIG_FILE = "settings.conf"
  character(len=256)          :: datadir
  character(len=256)          :: outdir
  character(len=256)          :: vort_name
  character(len=256)          :: u_name
  character(len=256)          :: v_name
  character(len=256)          :: psea_name
  character(len=256)          :: land_name
  character(len=256)          :: prefix_lvl
  character(len=256)          :: prefix_sfc
  integer                     :: year_start, month_start, day_start, hour_start
  integer                     :: year_end, month_end, day_end, hour_end
  integer                     :: vor_lvl
  integer                     :: steer_lvl_btm, steer_lvl_top
  integer                     :: proj
  integer                     :: vert_grid
  integer                     :: nx1, nx2, ny1, ny2
  !integer                    :: nt
  !real(wp)                   :: del_t
  ! parameter for smoothing of vorticity
  integer                     :: smth_type
  integer                     :: nsmth_x, nsmth_y
  real(wp)                    :: r_smth
  ! parameter for detecting vortex
  real(wp)                    :: zeta_max0, zeta_min0
  real(wp)                    :: int_zeta_min0, gamma
  ! parameter for excluding the synoptic scale disturbances
  real(wp)                    :: d_cf_min, size_synop
  real(wp)                    :: del_psea_min, distance_ec
  ! parameter for calculating steering winds
  integer                     :: steering_type
  integer                     :: n_steering_x, n_steering_y
  real(wp)                    :: r_steering
  ! parameter for linking vortex
  integer                     :: track_type
  real(wp)                    :: del_lon, del_lat, del_r
  ! parameter for checking the track
  integer                     :: period_min

  ! Debug flag
  logical                     :: dbg

contains
  subroutine get_config_params()

    implicit none

    ! Local variables
    integer, parameter            :: fh = 999
    character(len=256)            :: buffer
    character(len=256)            :: label
    integer                       :: line
    integer                       :: ios
    integer                       :: sep

#ifdef debug
    dbg = .true.
#else
    dbg = .false.
#endif

    ios = 0
    line = 0

    open (fh, file=trim(config_file), form='formatted', status='old', &
      &   iostat=ios, action='read')
    if (ios == 0) then
      if (dbg) write(*, *) '--------------------------------'
      if (dbg) write(*, *) 'Reading configuration parameters'
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
          case('year_start'); read(buffer, *, iostat=ios) year_start; if (dbg) write(*, *) year_start
          case('month_start'); read(buffer, *, iostat=ios) month_start; if (dbg) write(*, *) month_start
          case('day_start'); read(buffer, *, iostat=ios) day_start; if (dbg) write(*, *) day_start
          case('hour_start'); read(buffer, *, iostat=ios) hour_start; if (dbg) write(*, *) hour_start
          case('year_end'); read(buffer, *, iostat=ios) year_end; if (dbg) write(*, *) year_end
          case('month_end'); read(buffer, *, iostat=ios) month_end; if (dbg) write(*, *) month_end
          case('day_end'); read(buffer, *, iostat=ios) day_end; if (dbg) write(*, *) day_end
          case('hour_end'); read(buffer, *, iostat=ios) hour_end; if (dbg) write(*, *) hour_end
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
          case('proj'); read(buffer, *, iostat=ios) proj; if (dbg) write(*, *) proj
          case('vert_grid'); read(buffer, *, iostat=ios) vert_grid; if (dbg) write(*, *) vert_grid
          case('nx1'); read(buffer, *, iostat=ios) nx1; if (dbg) write(*, *) nx1
          case('nx2'); read(buffer, *, iostat=ios) nx2; if (dbg) write(*, *) nx2
          case('ny1'); read(buffer, *, iostat=ios) ny1; if (dbg) write(*, *) ny1
          case('ny2'); read(buffer, *, iostat=ios) ny2; if (dbg) write(*, *) ny2
          case('smth_type'); read(buffer, *, iostat=ios) smth_type; if (dbg) write(*, *) smth_type
          case('nsmth_x'); read(buffer, *, iostat=ios) nsmth_x; if (dbg) write(*, *) nsmth_x
          case('nsmth_y'); read(buffer, *, iostat=ios) nsmth_y; if (dbg) write(*, *) nsmth_y
          case('r_smth'); read(buffer, *, iostat=ios) r_smth; if (dbg) write(*, *) r_smth
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
          case('period_min'); read(buffer, *, iostat=ios) period_min; if (dbg) write(*, *) period_min
          case default
            if (index (trim(label), "#") /= 1) then
              write(*, *) 'ConfigParseWarning: Skipping invalid line', line
            end if
          end select
        end if
      end do
      if (dbg) write(*, *) '--------------------------------'
    end if
    close(fh)

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
end module params
