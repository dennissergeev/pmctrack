module config_params
  use types, only : wp

  implicit none

  character(len=*), parameter :: CONFIG_FILE = "config.ini"
  character(len=256) :: datadir
  character(len=256) :: vort_name
  character(len=256) :: u_name
  character(len=256) :: v_name
  character(len=256) :: psea_name
  character(len=256) :: land_name
  integer            :: year_start
  integer            :: month_start
  integer            :: day_start
  integer            :: year_end
  integer            :: month_end
  integer            :: day_end
  integer            :: thelevel
  integer            :: proj
  integer            :: vert_grid
  integer            :: nx1, nx2, ny1, ny2
  !integer            :: nt
  !real(wp)           :: del_t
  ! parameter for smoothing of vorticity
  integer            :: smth_type
  integer            :: nsmth_x, nsmth_y
  real(wp)           :: r_smth
  ! parameter for detecting vortex
  real(wp)           :: zeta_max0, zeta_min0
  real(wp)           :: int_zeta_min0, gamma
  ! parameter for excluding the synoptic scale disturbances
  real(wp)           :: d_cf_min, size_synop
  real(wp)           :: del_psea_min, distance_ec
  ! parameter for calculating steering winds
  integer            :: steering_type
  integer            :: n_steering_x, n_steering_y
  real(wp)           :: r_steering
  ! parameter for linking vortex
  integer            :: track_type
  real(wp)           :: del_lon, del_lat,del_r
  ! parameter for checking the track
  integer            :: period_min
  
contains
  subroutine get_config_params()
  
    implicit none
  
    ! Local variables
    integer, parameter              :: fh = 999
    character(len=256)              :: buffer
    character(len=256)              :: label
    integer                         :: line
    integer                         :: ios
    integer                         :: sep
    
    ios = 0
    line = 0
  
    open (fh, file=trim(config_file), form='formatted', status='old', &
      &   iostat=ios, action='read')
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
          case ('year_start'); read(buffer, *, iostat=ios) year_start
          case ('month_start'); read(buffer, *, iostat=ios) month_start
          case ('day_start'); read(buffer, *, iostat=ios) day_start
          case ('year_end'); read(buffer, *, iostat=ios) year_end
          case ('month_end'); read(buffer, *, iostat=ios) month_end
          case ('day_end'); read(buffer, *, iostat=ios) day_end
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
          !case ('del_t'); read(buffer, *, iostat=ios) del_t
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
            if (index (trim(label), "#") /= 1) then
              write(*, *) 'ConfigParseWarning: Skipping invalid line', line
            end if
          end select
        end if
      end do
    end if
    close(fh)
  end subroutine get_config_params
end module config_params
