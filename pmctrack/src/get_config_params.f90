subroutine get_config_params(config_file, datadir, year, month, &
  & vort_name, u_name, v_name, psea_name, land_name,            &
  & thelevel, proj, vert_grid, nx1, nx2, ny1, ny2,              &
  & smth_type, nsmth_x, nsmth_y, r_smth,                        &
  & zeta_max0, zeta_min0, int_zeta_min0, gamma,                 &
  & d_cf_min, size_synop, del_psea_min, distance_ec,            &
  & steering_type, n_steering_x, n_steering_y, r_steering,      &
  & track_type, del_lon, del_lat, del_r, period_min)

  use types, only : wp

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
  !integer         , intent(inout) :: nt
  !real(wp)        , intent(inout) :: del_t
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
