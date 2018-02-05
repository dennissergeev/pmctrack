subroutine link_vort_rad(nx, ny, lon, lat, del_t, mtype,                      &
                       & mlon_prev, mlat_prev, mlon, mlat,                    &
                       & u_vor_f_prev, v_vor_f_prev,                          &
                       & vor_part, n_max_prev, n_max,                         &
                       & vor_num, vor_index, vor_merge)

  use types, only: wp
  use constants, only: ra, fillval, nmax, pmax, rad2deg, deg2rad
  use params, only: proj, del_r !, dbg
  use utils, only: sind, cosd, great_circle

  implicit none

  integer    , intent(in)    :: nx, ny
  real   (wp), intent(in)    :: lon          (     0:nx        )
  real   (wp), intent(in)    :: lat          (           0:ny  )
  real   (wp), intent(in)    :: del_t
  integer    , intent(in)    :: mtype        (nmax             )
  integer    , intent(in)    :: n_max_prev, n_max
  real   (wp), intent(in)    :: mlon_prev    (nmax             )
  real   (wp), intent(in)    :: mlat_prev    (nmax             )
  real   (wp), intent(in)    :: mlon         (nmax             )
  real   (wp), intent(in)    :: mlat         (nmax             )
  real   (wp), intent(in)    :: u_vor_f_prev (nmax             )
  real   (wp), intent(in)    :: v_vor_f_prev (nmax             )
  integer    , intent(in)    :: vor_part     (     0:nx, 0:ny  )
  integer    , intent(inout) :: vor_index    (pmax             ) ! 2? ! TODO: change to nmax?
  integer    , intent(out)   :: vor_merge    (pmax             )
  integer    , intent(inout) :: vor_num
  ! Local variables
  integer                    :: i_max, i_max1
  integer                    :: i_next       (nmax             )
  integer                    :: vor_index_old(pmax             )
  real   (wp)                :: r_next       (nmax             ) ! TODO: check!
  real   (wp)                :: r_next_tmp
  integer                    :: vor_part_s   (nmax             )
  integer                    :: i, j
  integer                    :: i_vor_num, i_vor_num2
  logical                    :: vor_prev_flag(nmax             )
  real   (wp)                :: e_mv_lon, e_mv_lat
  real   (wp)                :: e_mlon, e_mlat
  real   (wp)                :: r_c_min
  real   (wp)                :: r_tmp
  ! NEW VARIABLES
  real   (wp)                :: max_dist ! Search radius for a vortex at the next time step


  max_dist = del_r * 1.e3

  r_tmp = 0.
  e_mv_lon = 0.
  e_mv_lat = 0.

  r_next = fillval

  vor_merge = -999
  vor_index_old = vor_index
  vor_index = -999

  i_next = -999

  vor_prev_flag(1:nmax) = .false.

  print*, '###################################################################'
  !print*, '>>> u_vor_f_prev', u_vor_f_prev(1:5)
  !print*, '>>> mlon_prev', mlon_prev(1:5)
  !print*, '>>> mlon', mlon(1:5)
  print*, '>>> n_max_prev, n_max',  n_max_prev, n_max
  print*, '>>> vor_index_old',  vor_index_old(1:5)
  print*, '>>> vor_num',  vor_num

  ! Skipped on the first iteration: vor_num=0
  do i_max = 1, n_max_prev
    do i_vor_num = 1, vor_num
      if (i_max == vor_index_old(i_vor_num)) then
        ! The vortex labeled as i_max at kt  existed at kt-1
        print*, '>>>vor_prev_flag(i_max) = .true.', i_max, i_vor_num
        vor_prev_flag(i_max) = .true.
      endif
    enddo
  enddo
  ! vor_num = 0

  ! t=kt -> t=kt+1
  do i_max = 1, n_max_prev ! Loop over vortices at the previous time step

    ! Estimate the distance travelled by i_max'th vortex from the steering winds
    if (proj==1) then
      e_mv_lon = (u_vor_f_prev(i_max) * del_t / (ra * cosd(mlat_prev(i_max)))) * rad2deg
      e_mv_lat = (v_vor_f_prev(i_max) * del_t /   ra                         ) * rad2deg
    elseif (proj==2) then
      e_mv_lon = u_vor_f_prev(i_max) * del_t
      e_mv_lat = v_vor_f_prev(i_max) * del_t
    endif

    ! Estimated lon and lat
    e_mlon = mlon_prev(i_max) + e_mv_lon
    e_mlat = mlat_prev(i_max) + e_mv_lat

    ! -------- Tracking by estimation of movement ------

    r_c_min = max_dist
    do i_max1 = 1, n_max ! Loop over vortices at the current time step
      ! Compare the distance travelled by vortex
      ! with the one estimated by steering winds
      if (proj == 1) then
        r_tmp = great_circle(mlon(i_max1), e_mlon,                  &
                           & mlat(i_max1), e_mlat, ra)
      elseif (proj == 2) then
        r_tmp=sqrt( (mlon(i_max1) - e_mlon)**2                         &
                   +(mlat(i_max1) - e_mlat)**2)
      endif 

      if (r_tmp <= r_c_min) then 
        i_next(i_max) = i_max1  ! i_max1'th vortex is i_max'th in i_next array
        r_c_min = r_tmp !TODO: check if this is valid
        r_next(i_max) = r_tmp
      endif
      print*, '>>>', i_max, i_max1, r_c_min, r_tmp, mlon_prev(i_max), mlat_prev(i_max), &
                   & e_mlon, e_mlat, mlon(i_max1), mlat(i_max1)
    enddo

    ! ------ Tracking by part of vortex --------

    !print*, '*****************************************************************'
    !print*, e_mlon, e_mlat
    !print*, i_next(1:5)
    !print*, r_next(1:5)
    !print*, '*****************************************************************'

    vor_part_s = 0

    if (i_next(i_max) < 1) then
    ! i_max'th vortex on the previous time step does not have a location
    ! within the r_c_min radius on the current time step
    ! Then the second requirement is considered, which is that a part of 
    ! an isolated vortex area at the next time step overlaps with the
    ! estimated area
      do j =0, ny
        do i = 0, nx
          if (proj == 1) then
            r_tmp = great_circle(lon(i), e_mlon,                            &
                               & lat(j), e_mlat, ra)
          elseif (proj == 2) then
            r_tmp = sqrt((lon(i)-e_mlon)**2 + (lat(j)-e_mlat)**2)
          endif

          if (r_tmp <= max_dist) then
            if (vor_part(i, j) > 0) then
              ! The vortex (i_max) is partly within the radius
              print*, i, j, lon(i), lat(j)
              vor_part_s(vor_part(i, j)) = vor_part_s(vor_part(i, j)) + 1
            endif
          endif
        enddo
      enddo

      print*, 'vor_part_s', vor_part_s(1:5)

      do ! "while loop"
        if (maxval(vor_part_s) == 0) exit

        ! If there are more than two isolated vortex areas for which part of 
        ! the area overlaps with the estimated area, the vortex at the previous
        ! time step is linked to the vortex with an isolated vortex area 
        ! having the largest amount of overlap with the estimated area
        i_next(i_max) = maxloc(vor_part_s, dim=1)

        if (proj == 1) then
          r_next(i_max) = great_circle(mlon(i_next(i_max)), e_mlon, &
                                     & mlat(i_next(i_max)), e_mlat, ra)
        elseif (proj == 2) then
          r_next(i_max) = sqrt((e_mlon - mlon(i_next(i_max)))**2      &
                            & +(e_mlat - mlat(i_next(i_max)))**2)
        endif


        if (r_next(i_max) <= 2.0 * max_dist) exit
        if (mtype(i_next(i_max)) >= 1) exit
        if (r_next(i_max) > 2.0 * max_dist) then
          vor_part_s(i_next(i_max)) = -999
          i_next(i_max) = -999
        endif
      enddo ! "while loop"
    endif  ! i_next(i_max) == 0
    !print*, '*****************************************************************'
    !print*, r_next(1:5)
    !print*, '*****************************************************************'
  enddo ! i_max loop

  print*, n_max_prev
  print*, i_next(1:5)
  print*, vor_prev_flag(1:5)

  ! ------------- Connecting vortex(kt)to vortex(kt+1) ------------------
  do i_max = 1, n_max ! ?????????????????
    !if (i_next(i_max) > 0) then
      write (*,'(A,I4,A,I4,A,I4,A,I4)') "Vortex labeled as ", i_max, &
            & ' is connected to vortex ', i_next(i_max), ' at the next time step'

      if (vor_prev_flag(i_max)) then  ! The vortex existed at kt-1
        write (*,'(A,I4,A,I4,A)') "Vortex labeled as ", i_max, ' existed at previous time step'
        do i_vor_num = 1, vor_num
          if (vor_index_old(i_vor_num) == i_max) then
            vor_index(i_vor_num) = i_next(i_max)
          endif
        enddo
      else ! The vortex appear at kt
        print*, '->', vor_num, i_max, i_next(i_max)
        vor_num = vor_num + 1
        vor_index_old(vor_num) = i_max
        vor_index(vor_num) = i_next(i_max)
        ! print*, '--->', vor_num, i_max, i_next(i_max)
      endif
    !endif ! i_next > 0
  enddo ! i_max loop

  print*, '+++vor_num', vor_num
  print*,vor_index_old
  !stop


  !------- check the merger of the vortices -------
  do i_vor_num = 1, vor_num
    if (vor_index_old(i_vor_num) > 0) then ! TODO: check this!!!
      r_next_tmp = r_next(vor_index_old(i_vor_num))

      do i_vor_num2 = 1, vor_num
        if (vor_index(i_vor_num) == vor_index(i_vor_num2)         &
          & .and. vor_index(i_vor_num) > 0                              &
          & .and. i_vor_num /= i_vor_num2) then
          if (vor_merge(i_vor_num) /= i_vor_num2 &
            & .and. vor_merge(i_vor_num2) == 0) then

            if (r_next_tmp > r_next(vor_index_old(i_vor_num2))) then
              vor_merge(i_vor_num) = i_vor_num2
              r_next_tmp = r_next(vor_index_old(i_vor_num2))
            endif

            if (vor_merge(i_vor_num) > 0) then
              write (*, *) 'vortex', i_vor_num, &
                   &' merged with vortex', vor_merge(i_vor_num)
            endif
          endif
        endif
      enddo ! i_vor_num2
    endif
  enddo ! i_vor_num

  print*, '>>> vor_index',  vor_index(1:5)
  print*, '###################################################################'

end subroutine link_vort_rad
