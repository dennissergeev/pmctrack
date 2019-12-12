subroutine synop_check(mlon, mlat, n_max, minlon, minlat, n_min, mtype)

  use types, only: wp
  use constants, only: pi, ra, rkilo
  use params, only: proj, r_ec
  use utils, only: great_circle

  implicit none

  integer    , intent(in)    :: n_max
  integer    , intent(in)    :: n_min
  real   (wp), intent(in)    :: mlon(n_max)
  real   (wp), intent(in)    :: mlat(n_max)
  real   (wp), intent(in)    :: minlon(n_min)
  real   (wp), intent(in)    :: minlat(n_min)
  integer    , intent(inout) :: mtype(n_max)
  ! Local variables
  integer                    :: i_max, i_min
  real   (wp)                :: dist
  logical                    :: flag_synop(n_max)

  dist = 0.

  flag_synop = .false.

  do i_max = 1, n_max
    if (mtype(i_max) >= 1) then

      do i_min = 1, n_min
        if (proj==1) then
          dist = great_circle(mlon(i_max), minlon(i_min), &
                            & mlat(i_max), minlat(i_min), ra)
        elseif(proj==2)then
          dist = sqrt((mlon(i_max) - minlon(i_min))**2 &
                   & +(mlat(i_max) - minlat(i_min))**2)
        endif

        if (dist < r_ec * rkilo) flag_synop(i_max) = .true.
      enddo

      if (mtype(i_max) >= 2 .and. .not. flag_synop(i_max)) then
        mtype(i_max) = mtype(i_max) - 2
      endif
    endif
  enddo

end subroutine synop_check
