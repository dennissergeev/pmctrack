module nc_io

  use netcdf
  use datetime_module
  use types, only : wp

contains
  subroutine get_dims(nc_file_name, dim_names, ntime, nlvls, nlats, nlons)
    implicit none

    character(len=*)              , intent(in)    :: nc_file_name
    character(len=*), dimension(4), intent(in)    :: dim_names
    integer                       , intent(inout) :: ntime 
    integer                       , intent(inout) :: nlvls 
    integer                       , intent(inout) :: nlats 
    integer                       , intent(inout) :: nlons 
    ! Local variables
    integer                                       :: ncid
    integer                                       :: dimid
    integer                       , dimension(4)  :: dims
    integer                                       :: i

    call check( nf90_open(nc_file_name, nf90_nowrite, ncid) )

    do i = 1, size(dim_names)
      call check( nf90_inq_dimid(ncid, dim_names(i), dimid) )
      call check( nf90_inquire_dimension(ncid, dimid, len=dims(i)) )
    end do
    ntime = dims(1)
    nlvls = dims(2)
    nlats = dims(3)
    nlons = dims(4)

    call check( nf90_close(ncid) )
  end subroutine get_dims


  subroutine get_time(nc_file_name, time_name, time, time_step_s, cal_start)
    implicit none

    character(len=*)              , intent(in)    :: nc_file_name
    character(len=*)              , intent(in)    :: time_name
    integer                       , intent(out)   :: time(:) 
    real(wp)                      , intent(out)   :: time_step_s
    type(datetime)                , intent(out)   :: cal_start
    ! Local variables
    integer                                       :: ncid
    integer                                       :: var_id
    character(len=256)                            :: units
    character(len=256)                            :: units_date
    character(len=*), parameter                   :: cal_fmt = "since %Y-%m-%d %H:%M:%S"
    integer                                       :: i

    call check( nf90_open(nc_file_name, nf90_nowrite, ncid) )

    call check( nf90_inq_varid(ncid, time_name, var_id) )
    call check( nf90_get_var(ncid, var_id, time) )

    call check( nf90_get_att(ncid, var_id, "units", units) )

    ! "hours since 1900-01-01 00:00:0.0"
    units_date = trim(units)
    i = scan(units_date, ' ')
    units = trim(units_date(1:i-1))
    if (units == "days") then
      time_step_s = 24 * 3600.
    elseif (units == "hours") then
      time_step_s = 3600.
    elseif (units == "seconds") then
      time_step_s = 1.
    else
      write(*, *) 'CalendarParseError: unrecognised units'
      stop
    endif
   
    i = scan(trim(units_date(i+1:)), ' ')
    units_date = trim(units_date(i+1:))
    cal_start = strptime(trim(units_date), cal_fmt)
    ! Calendars usually start from 00:00
    cal_start%hour = 0
    cal_start%minute = 0
    cal_start%second = 0
    cal_start%millisecond = 0
    if (.not. cal_start%isValid()) then
      write(*, *) 'CalendarParseError: calendar start', cal_start, ' is not valid'
      stop
    endif

    call check( nf90_close(ncid) )
  end subroutine get_time


  subroutine get_coords(nc_file_name, dim_names, lons, lats, lvls, &
    & time, time_idx, nt)
    implicit none

    character(len=*)              , intent(in)    :: nc_file_name
    character(len=*), dimension(4), intent(in)    :: dim_names
    integer                       , intent(inout) :: time(:) 
    real(wp)                      , intent(inout) :: lvls(:) 
    real(wp)                      , intent(inout) :: lats(:) 
    real(wp)                      , intent(inout) :: lons(:) 
    integer                       , intent(in)    :: time_idx
    integer                       , intent(in)    :: nt
    ! Local variables
    integer                                       :: ncid
    integer                                       :: var_id

    call check( nf90_open(nc_file_name, nf90_nowrite, ncid) )

    call check( nf90_inq_varid(ncid, dim_names(1), var_id) )
    call check( nf90_get_var(ncid, var_id, time,  &
                             start=(/time_idx/), &
                             count=(/      nt/)) )
    call check( nf90_inq_varid(ncid, dim_names(2), var_id) )
    call check( nf90_get_var(ncid, var_id, lvls) )
    call check( nf90_inq_varid(ncid, dim_names(3), var_id) )
    call check( nf90_get_var(ncid, var_id, lats) )
    call check( nf90_inq_varid(ncid, dim_names(4), var_id) )
    call check( nf90_get_var(ncid, var_id, lons) )

    call check( nf90_close(ncid) )
  end subroutine get_coords


  subroutine get_xyzt_from_xyzt(nc_file_name, var_name, time_idx, nt, &
                              & lvl_idx, nz, var_data)
    implicit none

    character(len=*)              , intent(in)    :: nc_file_name
    character(len=*)              , intent(in)    :: var_name
    integer                       , intent(in)    :: time_idx
    integer                       , intent(in)    :: nt
    integer                       , intent(in)    :: lvl_idx
    integer                       , intent(in)    :: nz
    real(wp)                      , intent(inout) :: var_data(:, :, :, :)
    ! Local variables
    integer                                       :: ncid
    integer                                       :: var_id
    real(wp)                                      :: scale_factor, add_offset
    integer                                       :: shp(4)

    shp = shape(var_data)

    call check( nf90_open(nc_file_name, nf90_nowrite, ncid) )
    call check( nf90_inq_varid(ncid, var_name, var_id) )
    call check( nf90_get_var(ncid, var_id, var_data,  &
                             start=(/     1,      1, lvl_idx, time_idx/), &
                             count=(/shp(1), shp(2),      nz,       nt/)) )

    call check( nf90_get_att(ncid, var_id, "scale_factor", scale_factor) )
    call check( nf90_get_att(ncid, var_id, "add_offset", add_offset) )
  
    var_data = scale_factor * var_data + add_offset

    call check( nf90_close(ncid) )
  end subroutine get_xyzt_from_xyzt


  subroutine get_xyz_from_xyzt(nc_file_name, var_name, time_idx, &
                             & lvl_idx, nz, var_data)
    implicit none

    character(len=*)              , intent(in)    :: nc_file_name
    character(len=*)              , intent(in)    :: var_name
    integer                       , intent(in)    :: time_idx
    integer                       , intent(in)    :: lvl_idx
    integer                       , intent(in)    :: nz
    real(wp)                      , intent(inout) :: var_data(:, :, :)
    ! Local variables
    integer                                       :: ncid
    integer                                       :: var_id
    real(wp)                                      :: scale_factor, add_offset
    integer                                       :: shp(3)

    shp = shape(var_data)

    call check( nf90_open(nc_file_name, nf90_nowrite, ncid) )
    call check( nf90_inq_varid(ncid, var_name, var_id) )
    call check( nf90_get_var(ncid, var_id, var_data,  &
                             start=(/     1,      1, lvl_idx, time_idx/), &
                             count=(/shp(1), shp(2),      nz,        1/)) )

    call check( nf90_get_att(ncid, var_id, "scale_factor", scale_factor) )
    call check( nf90_get_att(ncid, var_id, "add_offset", add_offset) )
  
    var_data = scale_factor * var_data + add_offset

    call check( nf90_close(ncid) )
  end subroutine get_xyz_from_xyzt


  subroutine get_data_3d(nc_file_name, var_name, time_idx, nt, var_data)
    implicit none

    character(len=*)              , intent(in)    :: nc_file_name
    character(len=*)              , intent(in)    :: var_name
    integer                       , intent(in)    :: time_idx
    integer                       , intent(in)    :: nt
    real(wp)                      , intent(inout) :: var_data(:, :, :)
    ! Local variables
    integer                                       :: ncid
    integer                                       :: var_id
    real(wp)                                      :: scale_factor, add_offset
    integer                                       :: shp(3)

    shp = shape(var_data)

    call check( nf90_open(nc_file_name, nf90_nowrite, ncid) )
    call check( nf90_inq_varid(ncid, var_name, var_id) )
    call check( nf90_get_var(ncid, var_id, var_data,  &
                             start=(/     1,      1, time_idx/), &
                             count=(/shp(1), shp(2),       nt/)) )

    call check( nf90_get_att(ncid, var_id, "scale_factor", scale_factor) )
    call check( nf90_get_att(ncid, var_id, "add_offset", add_offset) )
  
    var_data = scale_factor * var_data + add_offset

    call check( nf90_close(ncid) )
  end subroutine get_data_3d


  subroutine get_data_2d(nc_file_name, var_name, var_data)
    implicit none

    character(len=*)              , intent(in)    :: nc_file_name
    character(len=*)              , intent(in)    :: var_name
    real(wp)                      , intent(inout) :: var_data(:, :)
    ! Local variables
    integer                                       :: ncid
    integer                                       :: var_id
    real(wp)                                      :: scale_factor, add_offset

    call check( nf90_open(nc_file_name, nf90_nowrite, ncid) )

    call check( nf90_inq_varid(ncid, var_name, var_id) )
    call check( nf90_get_var(ncid, var_id, var_data) )
    call check( nf90_get_att(ncid, var_id, "scale_factor", scale_factor) )
    call check( nf90_get_att(ncid, var_id, "add_offset", add_offset) )
  
    var_data = scale_factor * var_data + add_offset

    call check( nf90_close(ncid) )
  end subroutine get_data_2d


  subroutine get_xy_from_xyzt(nc_file_name, var_name, lvl_idx, time_idx, &
                            & var_data)
    implicit none

    character(len=*)              , intent(in)    :: nc_file_name
    character(len=*)              , intent(in)    :: var_name
    integer                       , intent(in)    :: lvl_idx
    integer                       , intent(in)    :: time_idx
    real(wp)                      , intent(inout) :: var_data(:, :)
    ! Local variables
    integer                                       :: ncid
    integer                                       :: var_id
    real(wp)                                      :: scale_factor, add_offset
    integer                                       :: shp(2)

    shp = shape(var_data)

    call check( nf90_open(nc_file_name, nf90_nowrite, ncid) )

    call check( nf90_inq_varid(ncid, var_name, var_id) )
    call check( nf90_get_var(ncid, var_id, var_data,  &
                             start=(/     1,      1, lvl_idx, time_idx/), &
                             count=(/shp(1), shp(2),       1,       1/)) )
    call check( nf90_get_att(ncid, var_id, "scale_factor", scale_factor) )
    call check( nf90_get_att(ncid, var_id, "add_offset", add_offset) )
  
    var_data = scale_factor * var_data + add_offset

    call check( nf90_close(ncid) )
  end subroutine get_xy_from_xyzt


  subroutine get_xy_from_xyt(nc_file_name, var_name, time_idx, var_data)
    implicit none

    character(len=*)              , intent(in)    :: nc_file_name
    character(len=*)              , intent(in)    :: var_name
    integer                       , intent(in)    :: time_idx
    real(wp)                      , intent(inout) :: var_data(:, :)
    ! Local variables
    integer                                       :: ncid
    integer                                       :: var_id
    real(wp)                                      :: scale_factor, add_offset
    integer                                       :: shp(2)

    shp = shape(var_data)

    call check( nf90_open(nc_file_name, nf90_nowrite, ncid) )

    call check( nf90_inq_varid(ncid, var_name, var_id) )
    call check( nf90_get_var(ncid, var_id, var_data,  &
                             start=(/     1,      1, time_idx/), &
                             count=(/shp(1), shp(2),       1/)) )
    call check( nf90_get_att(ncid, var_id, "scale_factor", scale_factor) )
    call check( nf90_get_att(ncid, var_id, "add_offset", add_offset) )
  
    var_data = scale_factor * var_data + add_offset

    call check( nf90_close(ncid) )
  end subroutine get_xy_from_xyt


  subroutine get_xyt_from_xyzt(nc_file_name, var_name, lvl_idx, &
                             & time_idx, nt, var_data)
    implicit none

    character(len=*)              , intent(in)    :: nc_file_name
    character(len=*)              , intent(in)    :: var_name
    integer                       , intent(in)    :: lvl_idx
    integer                       , intent(in)    :: time_idx
    integer                       , intent(in)    :: nt
    real(wp)                      , intent(inout) :: var_data(:, :, :)
    ! Local variables
    integer                                       :: ncid
    integer                                       :: var_id
    real(wp)                                      :: scale_factor, add_offset
    integer                                       :: shp(3)

    shp = shape(var_data)

    call check( nf90_open(nc_file_name, nf90_nowrite, ncid) )

    call check( nf90_inq_varid(ncid, var_name, var_id) )
    call check( nf90_get_var(ncid, var_id, var_data,  &
                             start=(/     1,      1, lvl_idx, time_idx/), &
                             count=(/shp(1), shp(2),       1,       nt/)) )
    call check( nf90_get_att(ncid, var_id, "scale_factor", scale_factor) )
    call check( nf90_get_att(ncid, var_id, "add_offset", add_offset) )
  
    var_data = scale_factor * var_data + add_offset

    call check( nf90_close(ncid) )
  end subroutine get_xyt_from_xyzt


  subroutine check(status)
    implicit none
    integer, intent (in) :: status

    if (status /= nf90_noerr) then 
      write(*, *) "status:", status
      write(*, *) trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  

end module nc_io
