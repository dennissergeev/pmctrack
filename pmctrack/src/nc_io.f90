module nc_io

  use netcdf
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

  subroutine get_coords(nc_file_name, dim_names, time, lvls, lats, lons)
    implicit none

    character(len=*)              , intent(in)    :: nc_file_name
    character(len=*), dimension(4), intent(in)    :: dim_names
    integer                       , intent(inout) :: time(:) 
    integer                       , intent(inout) :: lvls(:) 
    real(wp)                      , intent(inout) :: lats(:) 
    real(wp)                      , intent(inout) :: lons(:) 
    ! Local variables
    integer                                       :: ncid
    integer                                       :: var_id

    call check( nf90_open(nc_file_name, nf90_nowrite, ncid) )

    call check( nf90_inq_varid(ncid, dim_names(1), var_id) )
    call check( nf90_get_var(ncid, var_id, time) )
    call check( nf90_inq_varid(ncid, dim_names(2), var_id) )
    call check( nf90_get_var(ncid, var_id, lvls) )
    call check( nf90_inq_varid(ncid, dim_names(3), var_id) )
    call check( nf90_get_var(ncid, var_id, lats) )
    call check( nf90_inq_varid(ncid, dim_names(4), var_id) )
    call check( nf90_get_var(ncid, var_id, lons) )

    call check( nf90_close(ncid) )
  end subroutine get_coords


  subroutine get_data_4d(nc_file_name, var_name, var_data)
    implicit none

    character(len=*)              , intent(in)    :: nc_file_name
    character(len=*)              , intent(in)    :: var_name
    real(wp)                      , intent(inout) :: var_data(:, :, :, :)
    ! Local variables
    integer                                       :: ncid
    integer                                       :: var_id
    real(wp)                                      :: scale_factor, add_offset

    call check( nf90_open(nc_file_name, nf90_nowrite, ncid) )

    call check( nf90_inq_varid(ncid, var_name, var_id) )
    !call check( nf90_inquire_variable(ncid, var_id, xtype=xtype, ndims=ndims) )
    call check( nf90_get_var(ncid, var_id, var_data) )
    call check( nf90_get_att(ncid, var_id, "scale_factor", scale_factor) )
    call check( nf90_get_att(ncid, var_id, "add_offset", add_offset) )
  
    var_data = scale_factor * var_data + add_offset

    call check( nf90_close(ncid) )
  end subroutine get_data_4d


  subroutine get_data_3d(nc_file_name, var_name, var_data)
    implicit none

    character(len=*)              , intent(in)    :: nc_file_name
    character(len=*)              , intent(in)    :: var_name
    real(wp)                      , intent(inout) :: var_data(:, :, :)
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


  subroutine get_one_level(nc_file_name, var_name, var_data, lvl_idx)
    implicit none

    character(len=*)              , intent(in)    :: nc_file_name
    character(len=*)              , intent(in)    :: var_name
    real(wp)                      , intent(inout) :: var_data(:, :, :)
    integer                       , intent(in)    :: lvl_idx
    ! Local variables
    integer                                       :: ncid
    integer                                       :: var_id
    real(wp)                                      :: scale_factor, add_offset
    integer                                       :: shp(3)

    shp = shape(var_data)

    call check( nf90_open(nc_file_name, nf90_nowrite, ncid) )

    call check( nf90_inq_varid(ncid, var_name, var_id) )
    call check( nf90_get_var(ncid, var_id, var_data,  &
                             start=(/1, 1, lvl_idx, 1/), &
                             count=(/shp(1), shp(2), 1, shp(3)/)) )
    call check( nf90_get_att(ncid, var_id, "scale_factor", scale_factor) )
    call check( nf90_get_att(ncid, var_id, "add_offset", add_offset) )
  
    var_data = scale_factor * var_data + add_offset

    call check( nf90_close(ncid) )
  end subroutine get_one_level


  subroutine get_time_slice_3d(nc_file_name, var_name, var_data, time_idx)
    implicit none

    character(len=*)              , intent(in)    :: nc_file_name
    character(len=*)              , intent(in)    :: var_name
    real(wp)                      , intent(inout) :: var_data(:, :, :)
    integer                       , intent(in)    :: time_idx
    ! Local variables
    integer                                       :: ncid
    integer                                       :: var_id
    real(wp)                                      :: scale_factor, add_offset
    integer                                       :: shp(3)

    shp = shape(var_data)

    call check( nf90_open(nc_file_name, nf90_nowrite, ncid) )

    call check( nf90_inq_varid(ncid, var_name, var_id) )
    call check( nf90_get_var(ncid, var_id, var_data,  &
                             start=(/1, 1, 1, time_idx/), &
                             count=(/shp(1), shp(2), shp(3), 1/)) )
    call check( nf90_get_att(ncid, var_id, "scale_factor", scale_factor) )
    call check( nf90_get_att(ncid, var_id, "add_offset", add_offset) )
  
    var_data = scale_factor * var_data + add_offset

    call check( nf90_close(ncid) )
  end subroutine get_time_slice_3d


  subroutine get_time_slice_2d(nc_file_name, var_name, var_data, time_idx)
    implicit none

    character(len=*)              , intent(in)    :: nc_file_name
    character(len=*)              , intent(in)    :: var_name
    real(wp)                      , intent(inout) :: var_data(:, :)
    integer                       , intent(in)    :: time_idx
    ! Local variables
    integer                                       :: ncid
    integer                                       :: var_id
    real(wp)                                      :: scale_factor, add_offset
    integer                                       :: shp(2)

    shp = shape(var_data)

    call check( nf90_open(nc_file_name, nf90_nowrite, ncid) )

    call check( nf90_inq_varid(ncid, var_name, var_id) )
    call check( nf90_get_var(ncid, var_id, var_data,  &
                             start=(/1, 1, time_idx/), &
                             count=(/shp(1), shp(2), 1/)) )
    call check( nf90_get_att(ncid, var_id, "scale_factor", scale_factor) )
    call check( nf90_get_att(ncid, var_id, "add_offset", add_offset) )
  
    var_data = scale_factor * var_data + add_offset

    call check( nf90_close(ncid) )
  end subroutine get_time_slice_2d


  subroutine check(status)
    implicit none
    integer, intent (in) :: status

    if (status /= nf90_noerr) then 
      write(*, *), "status:", status
      write(*, *), trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  

end module nc_io