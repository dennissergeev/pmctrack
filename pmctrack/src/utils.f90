module utils

use types, only : wp
use constants, only : missval, pi, deg2rad

implicit none

contains

  subroutine make_nc_file_name(nc_file_name, datadir, prefix, y, m, var_name)

    character(len=*)              , intent(inout) :: nc_file_name
    character(len=*)              , intent(in)    :: datadir
    character(len=*)              , intent(in)    :: prefix
    character(len=*)              , intent(in)    :: var_name
    integer                       , intent(in)    :: y
    integer                       , intent(in)    :: m
    
    write(nc_file_name, '(A,A,A,I4.4,A,I2.2,A,A,A)') trim(datadir), '/',      &
                                                   & trim(prefix),            &
                                                   & y, '.',                  &
                                                   & m, '.',                  &
                                                   & trim(var_name), '.nc'
  end subroutine make_nc_file_name


  subroutine apply_mask_2d(var, nx, ny, flag)
  
    integer(4), intent (in)    :: nx, ny
    real   (4), intent (inout) :: var (0:nx, 0:ny)
    real   (4), intent (in)    :: flag(0:nx, 0:ny)
    
    integer(4) :: i, j
  
    do j = 0, ny
      do i = 0, nx
        !  print*, i, j, flag(i, j)
        if (flag(i, j) == 1.) then
          var(i, j) = missval
        endif
      enddo
    enddo
  
  end subroutine apply_mask_2d
  
  
  function integral_p(var, p, nz)
    integer , intent (in)  :: nz
    real(wp), intent (in)  :: var(nz)
    real(wp), intent (in)  :: p  (nz)
    real(wp)               :: integral_p
    integer                :: k
    
    integral_p = 0.
  
    do k = 1, nz-1
      integral_p = integral_p + 0.5 * (var(k) + var(k+1)) * (p(k) - p(k+1))
    end do
    integral_p = integral_p / (p(1) - p(nz))
  
  end function integral_p


  function sind(x)
    real(wp) :: x
    real(wp) :: sind
    sind = sin(deg2rad * x)
  end function sind 


  function cosd(x)
    real(wp) :: x
    real(wp) :: cosd
    cosd = cos(deg2rad * x)
  end function cosd 


  function great_circle(lon1, lat1, lon2, lat2, ra)
    real(wp) :: lon1, lat1, lon2, lat2 ! in degrees
    real(wp) :: ra
    real(wp) :: ang_cos
    real(wp) :: great_circle

    ang_cos = cosd(lat1) * cosd(lat2) * cosd(lon2 - lon1)             &
                      & + sind(lat1) * sind(lat2)

    if (abs(ang_cos) < 1.0) then
      great_circle = ra * acos(ang_cos)
    else
      great_circle = 0.
      write(*, *) 'ValueError: ang_cos is', ang_cos, '; Setting great_circle to 0'
    endif
  end function

end module utils
