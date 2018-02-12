subroutine cf_synop_check(vor_in, vor_part, &
                        & n_part, nx, ny,   &
                        & lon, lat, mtype_part)

  use types, only: wp
  use constants, only: pi, ra, rkilo, pmax, deg2rad, pi_thirds, pi_ninths
  use params, only: proj, d_cf_min, size_synop
  use utils, only: cosd, sind, great_circle

  implicit none

  real   (wp), intent(in)  :: vor_in  (0:nx, 0:ny)
  integer    , intent(in)  :: vor_part(0:nx, 0:ny)
  integer    , intent(in)  :: n_part
  integer    , intent(in)  :: nx
  integer    , intent(in)  :: ny
  real   (wp), intent(in)  :: lon     (0:nx      )
  real   (wp), intent(in)  :: lat     (      0:ny)
  integer    , intent(out) :: mtype_part
  ! Local variables
  real   (wp)              :: lonin, latin
  integer                  :: pnum !----cold front
  integer    , parameter   :: pnum_thresh = 6
  integer                  :: j_n, j_s, i_n, i_s
  logical                  :: flag_one
  real   (wp)              :: d_cf, theta_d_cf
  real   (wp)              :: one(pmax, 1:2)
  integer                  :: one_num
  real   (wp)              :: a, b, c, k, r2
  real   (wp)              :: size_vor ! synoptic
  integer                  :: i, j


  lonin = lon(1)-lon(0)
  latin = lat(1)-lat(0)

  theta_d_cf = 0. ! to avoid warnings
  d_cf = 0. ! just in case

  i_n = 0
  i_s = nx
  j_n = 0
  j_s = ny

  size_vor = 0.

  one_num = 0
  one(:, :) = 0

  do j = 1, ny-1
    do i = 1, nx-1
      if (vor_part(i, j) == n_part) then
        if (proj == 1) then
          size_vor = size_vor &
                   & + (lonin * deg2rad * ra / rkilo * cosd(lat(j)) &
                   &    * latin * deg2rad * ra / rkilo)
        elseif (proj == 2) then
          size_vor = size_vor + lonin * latin * rkilo**(-2)
        endif

        call median_check(vor_in(i-1:i+1, j-1:j+1), flag_one, pnum)

        if (flag_one) then
          if (pnum >= pnum_thresh) then
            one_num = one_num + 1
            one(one_num, 1) = lon(i)
            one(one_num, 2) = lat(j)
          endif
        endif

        if (j < j_s) then
          j_s = j
          i_s = i
        endif
        if (j >= j_n) then
          j_n = j
          i_n = i
        endif
      endif
    enddo
  enddo

  !if (dbg) then
  !    print*, 'cf_synop_check(99):', 'i_s', i_s, 'i_n', i_n
  !endif

  if (proj == 1) then
    d_cf = great_circle(lon(i_n), lon(i_s), lat(j_n), lat(j_s), ra)
    theta_d_cf = acos((ra * (lat(j_n) - lat(j_s)) * deg2rad) / d_cf)
  elseif(proj==2)then
    d_cf = sqrt((lon(i_s) - lon(i_n))**2 + (lat(j_s) - lat(j_n))**2) * 1.0e-3
    theta_d_cf = (lat(j_s) - lat(j_n)) * 1.0e-3 / d_cf
  endif

  if (d_cf >= d_cf_min * rkilo) then
    if (lon(i_n) > lon(i_s) .and. theta_d_cf <= pi_thirds)then
      call quadric_fit(one(1:one_num, 1:2), one_num, a, b, c, k, r2)

      if(r2 >= 0.8 .and. k <= 0.1) then
        mtype_part=1
      endif

    elseif (theta_d_cf <= pi_ninths) then
      call quadric_fit(one(1:one_num, 1:2), one_num, a, b, c, k, r2)

      if(r2 >= 0.8 .and. k<= 0.1) then
        mtype_part=1
      endif
    endif
  endif

!--------------- Check synoptic low ----------------
  if (size_vor >= size_synop) then
    mtype_part = mtype_part + 2
  endif
end subroutine cf_synop_check


subroutine median_check(val, flag, pnum)

  use types, only: wp
  use constants, only: mx, my

  implicit none

  real   (wp), intent (in)  :: val(-1:1, -1:1)
  logical    , intent (out) :: flag
  integer    , intent (out) :: pnum
  ! Local variables
  integer                   :: m
  integer                   :: mmax

  flag = .false.
  pnum = 0
  mmax = size(mx)

  do m = 1, mmax
    if (val (0, 0) > val(mx(m), my(m))) then
      pnum = pnum + 1
    endif
  enddo

  if (pnum >= 5) flag = .true.
end subroutine median_check


subroutine quadric_fit(one, one_num, a, b, c, k, r2)

  use types, only: wp

  implicit none

  real   (wp), intent (in) :: one(1:one_num, 1:2)
  integer    , intent (in) :: one_num
  real   (wp), intent (out):: a, b, c, k, r2
  ! Local variables
  integer                  :: i
  real   (wp)              :: sx, sx2, sx3, sx4, sx2y, sxy, sy, s1
  real   (wp)              :: xbar, ybar
  real   (wp)              :: x  (1:one_num)
  real   (wp)              :: y  (1:one_num)
  real   (wp)              :: fit(1:one_num)

  sx = 0.
  sx2 = 0.
  sx3 = 0.
  sx4 = 0.
  sx2y = 0.
  sxy = 0.
  sy = 0.
  s1 = 0.

  xbar = 0.
  ybar = 0.

  do i = 1, one_num
    xbar = xbar + one(i, 2)
    ybar = ybar + one(i, 1)
  enddo

  xbar = xbar / one_num
  ybar = ybar / one_num

  do i = 1, one_num
    x(i) = one(i, 2) - xbar
    y(i) = one(i, 1) - ybar
  enddo

  do i = 1, one_num
    sx4 = sx4 + x(i)**4
    sx3 = sx3 + x(i)**3
    sx2 = sx2 + x(i)**2
    sx = sx + x(i)**1
    sx2y = sx2y + x(i)**2 * y(i)
    sxy = sxy + x(i) * y(i)
    sy = sy + y(i)
    s1 = s1 + 1.0
  enddo

  a = (s1*sx2*sx2y-sx*sx*sx2y+sx*sx2*sxy-s1*sx3*sxy+sx*sx3*sy-sx2*sx2*sy)&
       &/(2*sx*sx2*sx3+s1*sx2*sx4-sx*sx*sx4-s1*sx3*sx3-sx2*sx2*sx2)

  b = (sx*sx2*sx2y-s1*sx3*sx2y+s1*sx4*sxy-sx2*sx2*sxy+sx2*sx3*sy-sx*sx4*sy)&
       &/(2*sx*sx2*sx3+s1*sx2*sx4-sx*sx*sx4-s1*sx3*sx3-sx2*sx2*sx2)

  c = (-sx2*sx2*sx2y+sx*sx3*sx2y-sx*sx4*sxy+sx2*sx3*sxy-sx3*sx3*sy+sx2*sx4*sy)&
       &/(2*sx*sx2*sx3+s1*sx2*sx4-sx*sx*sx4-s1*sx3*sx3-sx2*sx2*sx2)

  do i = 1, one_num
    fit(i) = a * x(i)**2 + b * x(i) + c
  enddo

  call correlation(y(1:one_num), fit(1:one_num), one_num, r2)

  k = 0.
  do i = 1, one_num
    k = k + a / (1 + (2 * a * x(i) + b)**2.) ** 1.5
  enddo

  k = k / one_num
end subroutine quadric_fit


subroutine correlation(x, y, n, r2)

  use types, only: wp

  implicit none

  integer    , intent (in) :: n
  real   (wp), intent (in) :: x(n)
  real   (wp), intent (in) :: y(n)
  real   (wp), intent (out):: r2
  ! Local variables
  real   (wp)              :: xm, ym, sigx, sigy, cor
  real   (wp)              :: sx, sy, sxy, sx2, sy2
  integer                  :: i


  sx = 0; sy = 0; sxy = 0; sx2 = 0; sy2 = 0

  if (n > 2) then
    do i = 1, n
      sx = sx + x(i)
      sy = sy + y(i)
    enddo
    xm = sx / n
    ym = sy / n
    do i = 1, n
      sx2 = sx2 + (x(i) - xm) * (x(i) - xm)
      sy2 = sy2 + (y(i) - ym) * (y(i) - ym)
      sxy = sxy + (x(i) - xm) * (y(i) - ym)
    enddo

    r2 = sxy * sxy / (sx2 * sy2)

    sigx = sqrt(sx2/n)
    sigy = sqrt(sy2/n)
    if (sigx > 0 .and. sigy > 0) then
      cor = (sxy / n) / (sigx * sigy)
    endif
  endif
end subroutine correlation
