! temporary reference:
! call linkin_vort2(mlon,mlat,mtype,u_vor_f,v_vor_f,&
!      &nt,n_max,vor_index,vor_num,vor_merge,&
!      &vor_part(nx1:nx2,ny1:ny2,1:nt),nx12,ny12,lon(nx1:nx2),lat(ny1:ny2),&
!      del_t)
subroutine linkin_vort2(mlon, mlat, mtype, u_vor_f, v_vor_f, nt,              &
  &                     n_max, vor_index, vor_num, vor_merge, vor_part,       &
  &                     nx, ny, lon, lat, del_t)

  use types, only: wp
  use constants, only: ra, fillval, nmax, pmax, rad2deg, deg2rad
  use params, only: proj, del_r, dbg 
  use utils, only: sind, cosd, great_circle

  implicit none 

  integer    , intent(in)  :: nt, nx, ny
  real   (wp), intent(in)  :: del_t
  integer    , intent(in)  :: n_max        (                 nt)
  real   (wp), intent(in)  :: mlon         (nmax,            nt)
  real   (wp), intent(in)  :: mlat         (nmax,            nt)
  integer    , intent(in)  :: mtype        (nmax,            nt)
  real   (wp), intent(in)  :: u_vor_f      (nmax,            nt)
  real   (wp), intent(in)  :: v_vor_f      (nmax,            nt)
  real   (wp), intent(in)  :: lon          (     0:nx          )
  real   (wp), intent(in)  :: lat          (           0:ny    )
  integer    , intent(in)  :: vor_part     (     0:nx, 0:ny, nt)
  integer    , intent(out) :: vor_index    (pmax,            nt) ! TODO: change to nmax?
  integer    , intent(out) :: vor_merge    (pmax               )
  integer    , intent(out) :: vor_num
  ! Local variables
  integer                  :: kt
  integer                  :: i_max, i_max1
  integer                  :: i_next       (nmax               )
  real   (wp)              :: r_next     (0:nmax               ) ! TODO: check!
  real   (wp)              :: r_next_tmp
  integer                  :: vor_part_s   (nmax               )
  integer                  :: i, j
  integer                  :: i_vor_num, i_vor_num2
  logical                  :: vor_prev_flag(nmax)
  ! integer                  :: vor_prev_idx (nmax,            nt)
  real   (wp)              :: e_mv_lon, e_mv_lat
  real   (wp)              :: e_mlon, e_mlat
  real   (wp)              :: r_c_min
  real   (wp)              :: r_tmp, theta_tmp
  ! NEW VARIABLES
  real   (wp)              :: max_dist


  max_dist = del_r * 1.e3    

  r_tmp = 0.
  theta_tmp = 0.
  e_mv_lon = 0.
  e_mv_lat = 0.

  r_next = fillval

  vor_merge = 0
  vor_index = 0
  vor_num = 0

  i_next = 0

  vor_prev_flag(1:nmax) = .false.
  ! vor_prev_idx = 0

  do kt = 1, nt-1 ! Time loop

    ! Skipped on the first iteration: vor_num=0
    do i_max = 1, n_max(kt)
      do i_vor_num = 1, vor_num           
        if (i_max == vor_index(i_vor_num, kt)) then 
          ! The vortex labeled as i_max at kt  existed at kt-1
          vor_prev_flag(i_max) = .true.
          ! vor_prev_idx(i_max, kt) = vor_index(i_vor_num, kt-1)
        endif
      enddo
    enddo

    ! t=kt -> t=kt+1
    do i_max = 1, n_max(kt)

      ! Calculate distance travelled by the vortex
      if (proj==1) then
        e_mv_lon = (u_vor_f(i_max, kt) * del_t / (ra * cosd(mlat(i_max,kt)))) * rad2deg        
        e_mv_lat = (v_vor_f(i_max, kt) * del_t / ra                         ) * rad2deg
      elseif (proj==2) then
        e_mv_lon = u_vor_f(i_max, kt) * del_t
        e_mv_lat = v_vor_f(i_max, kt) * del_t
      endif

      e_mlon = mlon(i_max, kt) + e_mv_lon        
      e_mlat = mlat(i_max, kt) + e_mv_lat
        
      ! -------- Tracking by estimation of movement ------

      r_c_min = max_dist
      do i_max1 = 1, n_max(kt+1)

        ! Compare the distance travelled by vortex
        ! with the one estimated by steering winds
        if (proj == 1) then
          r_tmp = great_circle(mlon(i_max1, kt+1), e_mlon,                  &
                             & mlat(i_max1, kt+1), e_mlat, ra)
            ! dist = cosd(e_mlat)*cosd(mlat(i_max1, kt+1))                    &
            !      & * cosd((e_mlon-mlon(i_max1,kt+1)))                       &
            !      & + sind(e_mlat)*sind(mlat(i_max1, kt+1))
            ! if(abs(dist)<1.0)then
            !   theta_tmp = acos(dist)
            ! else
            !   theta_tmp = 0.0
            ! endif
            ! r_tmp = ra * theta_tmp
        elseif (proj == 2) then
          r_tmp=sqrt( (mlon(i_max1,kt+1) - e_mlon)**2                         &
                     +(mlat(i_max1,kt+1) - e_mlat)**2)
        endif

        if (r_tmp <= r_c_min) then
          i_next(i_max) = i_max1
          r_c_min = r_tmp
          r_next(i_max) = r_tmp
        endif
      enddo

      ! ------ Tracking by part of vortex --------

      vor_part_s = 0

      if (i_next(i_max) == 0) then
        
        do j =0, ny
          do i = 0, nx
            if (proj == 1) then
              r_tmp = great_circle(lon(i), e_mlon,                            &
                                 & lat(j), e_mlat, ra)
              ! if(abs(cos(pi/180*e_mlat)*cos(pi/180*lat(j))&
              !        &*cos(pi/180*(e_mlon-lon(i)))&
              !        &+sin(pi/180*e_mlat)*sin(pi/180*lat(j)))<1.0)then

              !   theta_tmp=acos(cos(pi/180*e_mlat)*cos(pi/180*lat(j))&
              !        &*cos(pi/180*(e_mlon-lon(i)))&
              !        &+sin(pi/180*e_mlat)*sin(pi/180*lat(j)))

              ! else
              !   theta_tmp=0.0
              ! endif
              ! r_tmp=ra*theta_tmp
 
            elseif (proj == 2) then
              r_tmp = sqrt((lon(i)-e_mlon)**2 + (lat(j)-e_mlat)**2)
            endif
          
            if (r_tmp <= max_dist) then
              if (vor_part(i, j, kt+1) > 0) then
                vor_part_s(vor_part(i, j, kt+1)) = vor_part_s(vor_part(i, j, kt+1)) + 1
              endif
            endif
          enddo
        enddo

        
        do ! "while loop"
          if (maxval(vor_part_s) == 0) exit
          
          i_next(i_max) = maxloc(vor_part_s, 1)          
          
          if (proj == 1) then
            if (e_mlon /= mlon(i_next(i_max), kt+1))then
              r_next(i_max) = great_circle(mlon(i_next(i_max), kt+1), e_mlon, &
                                         & mlat(i_next(i_max), kt+1), e_mlat, ra)
              ! r_next(i_max, kt) =ra*&
              !      &acos(cos(pi/180*e_mlat)*cos(pi/180*mlat(i_next(i_max),kt+1))&
              !      &*cos(pi/180*(e_mlon-mlon(i_next(i_max),kt+1)))&
              !      &+sin(pi/180*e_mlat)*sin(pi/180*mlat(i_next(i_max),kt+1)))
            else
              r_next(i_max) = ra * deg2rad * (e_mlat - mlat(i_next(i_max), kt+1))
            endif
          elseif (proj == 2) then
            r_next(i_max) = sqrt((e_mlon - mlon(i_next(i_max), kt+1))**2      &
                              & +(e_mlat - mlat(i_next(i_max), kt+1))**2)
          endif


          if (r_next(i_max) <= 2.0 * max_dist) exit
          if (mtype(i_next(i_max), kt+1) >= 1) exit
          if (r_next(i_max) > 2.0 * max_dist) then
            vor_part_s(i_next(i_max)) = 0
            i_next(i_max) = 0
          endif
        enddo ! "while loop"
      endif  ! i_next(i_max) == 0
    enddo ! i_max loop

    ! ------------- Connecting vortex(kt)to vortex(kt+1) ------------------
    do i_max=1, n_max(kt)
      if (i_next(i_max) > 0) then
        write (*,'(A,I4,A,I4,A,I4,A,I4)') "Vortex labeled as ", i_max, ' at ' &
             & , kt, ' is connected to vortex ', i_next(i_max), ' at ', kt + 1

        if (vor_prev_flag(i_max)) then  ! The vortex existed at kt-1
          write (*,'(A,I4,A,I4,A)') "Vortex labeled as ", i_max, ' at ', kt, ' existed at kt-1'
          do i_vor_num = 1, vor_num
            if (vor_index(i_vor_num, kt) == i_max) then 
              vor_index(i_vor_num, kt+1) = i_next(i_max)
            endif
          enddo
        else ! The vortex appear at kt
          vor_num = vor_num + 1
          vor_index(vor_num, kt) = i_max
          vor_index(vor_num, kt+1) = i_next(i_max)
        endif
      endif ! i_next > 0
    enddo ! i_max loop
    

    !------- check the merger of the vortices -------  
    do i_vor_num = 1, vor_num
      if (dbg) then
        print*, 'linkin_vort2(229):', i_vor_num, kt, vor_index(i_vor_num, kt)
      endif
      r_next_tmp = r_next(vor_index(i_vor_num, kt))
          
      do i_vor_num2 = 1, vor_num
        if (vor_index(i_vor_num, kt+1) == vor_index(i_vor_num2, kt+1)         &
          & .and. vor_index(i_vor_num, kt+1) > 0                              &
          & .and. i_vor_num /= i_vor_num2) then
          if (vor_merge(i_vor_num) /= i_vor_num2 &
            & .and. vor_merge(i_vor_num2) == 0) then
            
            if (r_next_tmp > r_next(vor_index(i_vor_num2, kt))) then 
              vor_merge(i_vor_num) = i_vor_num2
              r_next_tmp = r_next(vor_index(i_vor_num2, kt))
            endif
            
            if (vor_merge(i_vor_num) > 0) then
              write (*, *) 'vortex', i_vor_num, &
                   &' merged with vortex', vor_merge(i_vor_num)
            endif
          endif
        endif
      enddo ! i_vor_num2
    enddo ! i_vor_num

  enddo  ! Time loop 

end subroutine linkin_vort2
