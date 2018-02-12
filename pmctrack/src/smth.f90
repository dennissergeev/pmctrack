subroutine smth(var, nx, ny, var_smth)

  use types, only : wp
  use params, only : nsmth_x, nsmth_y, nx1, nx2, ny1, ny2

  implicit none

  integer    , intent(in) :: nx
  integer    , intent(in) :: ny
  real   (wp), intent(in) :: var     (0:nx,       0:ny)
  real   (wp), intent(out):: var_smth(nx1:nx2, ny1:ny2)
  ! Local variables
  real   (wp)             :: tmp     (nx1:nx2, ny1:ny2)
  integer                 :: nc
  integer                 :: i, j, ii, jj
  real   (wp), parameter  :: var_thresh = -10.0


  var_smth(nx1:nx2, ny1:ny2) = var(nx1:nx2, ny1:ny2)

  do j = ny1, ny2
    if (j - nsmth_y >=0 .and. j + nsmth_y <= ny) then
      do i = nx1, nx2
        if (i - nsmth_x >= 0 .and. i + nsmth_x <= nx) then
          nc = 0
          tmp(i, j) = 0.
          if (var(i, j) > var_thresh) then
            do jj = -nsmth_y, nsmth_y
              do ii = -nsmth_x, nsmth_x
                if (var(i+ii, j+jj) > var_thresh) then
                  tmp(i, j) = tmp(i, j) + var(i+ii, j+jj)
                  nc = nc + 1
                endif
              enddo
            enddo
            tmp(i, j) = tmp(i, j) / ((2 * nsmth_x + 1) * (2 * nsmth_y + 1))
            var_smth(i, j) = tmp(i, j)
          endif
        endif
      enddo
    endif
  enddo
end subroutine smth


subroutine smth_r(var, nx, ny, lon, lat, var_smth)

  use types, only : wp
  use constants, only : pi, ra, rkilo, deg2rad
  use params, only : proj, r_smth, nx1, nx2, ny1, ny2
  use utils, only : great_circle, cosd

  implicit none

  integer    , intent(in) :: nx
  integer    , intent(in) :: ny
  real   (wp), intent(in) :: var     (0:nx,       0:ny)
  real   (wp), intent(in) :: lon     (0:nx            )
  real   (wp), intent(in) :: lat     (            0:ny)
  real   (wp), intent(out):: var_smth(nx1:nx2, ny1:ny2)
  ! Local variables
  real   (wp)             :: tmp     (nx1:nx2, ny1:ny2)
  integer                 :: nc
  integer                 :: i, j, ii, jj
  real   (wp)             :: lonin, latin
  integer                 :: x_smth, y_smth
  real   (wp)             :: dist
  real   (wp), parameter  :: var_thresh = -10.0

  x_smth = 0
  y_smth = 0

  dist = 0.

  lonin = lon(2) - lon(1)
  latin = lat(2) - lat(1)

  if (proj==1) then
    x_smth = nint(r_smth / (ra * lonin * deg2rad * cosd(lat(ny1+ny2/2)) / rkilo)) + 1
    y_smth = nint(r_smth / (ra * latin * deg2rad / rkilo)) + 1
  elseif (proj==2) then
    x_smth = nint(r_smth / lonin / rkilo) + 2
    y_smth = nint(r_smth / latin / rkilo) + 2 ! TODO: check why +2
  endif

  var_smth(nx1:nx2, ny1:ny2) = var(nx1:nx2, ny1:ny2)

  do j = ny1, ny2
    if (j - y_smth >= 0 .and. j + y_smth <= ny) then
      do i = nx1, nx2
        if (i - x_smth >= 0 .and. i + x_smth <= nx) then
          nc = 0
          tmp(i, j) = 0.
          if (var(i, j) > var_thresh) then
            do jj = -y_smth, y_smth
              do ii = -x_smth, x_smth
                if (proj == 1) then
                  dist = great_circle(lon(i), lon(i+ii), lat(j), lat(j+jj), ra)

                elseif(proj==2)then
                  dist = sqrt((ii * lonin)**2 + (jj * latin)**2)
                endif

                if (dist <= r_smth * rkilo) then
                  nc = nc + 1
                  if (var(i+ii, j+jj) > var_thresh) then
                    tmp(i, j) = tmp(i, j) + var(i+ii, j+jj)
                  endif
                endif
              enddo
            enddo
            var_smth(i, j) = tmp(i, j) / nc
          endif
        endif
      enddo
    endif
  enddo
end subroutine smth_r
