subroutine tracking_main(vor, u, v, psea, &
    & proj, vert_grid, &
    & nx1, nx2, ny1, ny2, levs, &
    & lons, lats, lonin, latin, del_t, &
    & nsmth_x, nsmth_y, r_smth, smth_type, &
    & zeta_max0, zeta_min0, int_zeta_min0, gamma, &
    & n_steering_x, n_steering_y, r_steering, steering_type, &
    & del_lon, del_lat, del_r, track_type, &
    & period_min, d_cf_min, size_synop, del_psea_min, distance_ec, &
    & nx, ny, nz, nt)
   
  use kind 

  implicit none

  integer :: nx, ny, nx1, nx2, ny1, ny2, nz, nt

  real(4) :: vor (0:nx, 0:ny,       1:nt)
  real(4) :: u   (0:nx, 0:ny, 1:nz, 1:nt)
  real(4) :: v   (0:nx, 0:ny, 1:nz, 1:nt)
  real(4) :: psea(0:nx, 0:ny,       1:nt)
  real(4) :: lons, lats, lonin, latin, levs, del_t

  integer (4) :: proj, vert_grid  ! 20

! parameter for smoothing of vorticity
  integer :: smth_type
  integer :: nsmth_x, nsmth_y
  real(4) :: r_smth

! parameter for detecting vortex
  real(4) :: zeta_max0, zeta_min0, int_zeta_min0, gamma

! parameter for excluding the synoptic scale disturbances
  real(4) :: d_cf_min, size_synop, del_psea_min, distance_ec

! parameter for calculating steering winds
  integer :: steering_type
  integer :: n_steering_x, n_steering_y
  real(4) :: r_steering

! parameter for linking vortex
  integer :: track_type
  real(4) :: del_lon, del_lat, del_r

! parameter for checking the track
  integer :: period_min ! 41 inputs

  !f2py intent(in) nx, ny, nx1, nx2, ny1, ny2, nz, nt
  !f2py intent(in) vor, u, v, psea
  !f2py intent(in) lons, lats, lonin, latin, levs, del_t
  !f2py intent(in) proj, vert_grid
  !f2py intent(in) smth_type, nsmth_x, nsmth_y, r_smth
  !f2py intent(in) zeta_max0, zeta_min0, int_zeta_min0, gamma
  !f2py intent(in) d_cf_min, size_synop, del_psea_min, distance_ec
  !f2py intent(in) steering_type, n_steering_x, n_steering_y, r_steering
  !f2py intent(in) track_type, del_lon, del_lat, del_r
  !f2py intent(in) period_min


  ! LOCAL VARIABLES
  real    (4), parameter  :: undef=9.99e20

  real    (4), allocatable :: lon(:), lat(:)
  real    (4), allocatable :: vor_smth  (:, :, :)
  real    (4), allocatable :: dummy     (:, :)
  integer (4), allocatable :: vor_part  (:, :, :)
  real    (4), allocatable :: vor_part_r(:, :)

  integer (4) :: i_rec

  integer (4) :: nx_out, ny_out
  integer (4) :: nx12, ny12
  integer (4) :: i, j, k
  integer (4) :: kt
  integer (4) :: ios
  real    (4), allocatable :: mlat  (:, :), mlon  (:, :), max_vor(:, :)
  real    (4), allocatable :: minlat(:, :), minlon(:, :), z_min  (:, :)
  integer  :: min_i, min_j
  real    (4), allocatable :: s_part(:, :)
  integer (4), allocatable :: mi(:, :), mj(:, :)
  integer (4), allocatable :: mtype(:, :), z_min_size(:, :)
  
  real    (4), allocatable :: onelat(:, :), onelon(:, :), vor_one(:, :)

  real    (4), allocatable :: u_vor_f(:, :), v_vor_f(:, :)
  real    (4), allocatable :: u_vor_b(:, :), v_vor_b(:, :)
  integer (4), allocatable :: n_max(:), n_min(:)
  integer (4) :: i_max, i_min

  integer (4) :: vor_num, i_vor_num, vor_num_out=0
  integer (4) :: vor_merge(10000), vor_merge_num(10000)
  integer (4), allocatable :: vor_index(:, :)
  integer (4) :: kt_appear(1000)
  integer (4) :: vor_period(10000)
  real    (4) :: vor_max_track(10000)

  real    (4) :: vor_movement

  real    (4), allocatable :: vortex(:, :, :)
  integer (4), allocatable :: vortex_type(:, :)
  integer (4), allocatable :: vortex_flag(:)

  character (50) :: fname_out, fname_loc, fname_track, fname_in


  write (*, *) 'nx=', nx, 'ny=', ny, 'nt=', nt, 'nz=', nz
  write (*, *) 'nx1=', nx1, 'nx2=', nx2, 'ny1=', ny1, 'ny2=', ny2

  allocate(lon(0:nx))
  allocate(lat(0:ny))

  do i=0, nx
    lon(i) = lons + lonin * i
  end do
  do j=0, ny
    lat(j) = lats + latin * j
  end do


  nx12 = nx2 - nx1
  ny12 = ny2 - ny1


  !Read data

!  write (*, *)"allocate time" 
  ! Time

  if(proj==2)then
    write (*, *)'Cartesian coordinate'
  elseif(proj==1)then
    write (*, *)'Geographical coordinate'
  else
    write (*, *)'Coordinate system not supported'
  end if

  if(vert_grid==1)then
    write (*, *)'pressure levels'
  elseif(vert_grid==2)then
    write (*, *)'height levels'
  else
    write (*, *)'vertical coordinate system not supported'
  end if


  if(smth_type==1)then 
    write (*, *)'smth latlon'
  elseif(smth_type==2)then
    write (*, *)'smth radius', r_smth, 'km'
  end if

  if(steering_type==1)then 
    write (*, *)'calculate steering wind in latlon coordinate'
  elseif(steering_type==2)then
    write (*, *)'calculate steering wind in radius', r_steering, 'km'
  end if

  if(track_type==1)then
    write (*, *)'Use del_lon, del_lat', del_lon, del_lat

  elseif(track_type==2)then
    write (*, *)'Use radious ', del_r, 'km'
  end if

  allocate(vor_smth(nx1:nx2, ny1:ny2, nt))
  allocate(vor_part(nx1:nx2, ny1:ny2, nt))
  allocate(vor_part_r(nx1:nx2, ny1:ny2))
  allocate(n_max(nt), n_min(nt))
  allocate(mlat(100, nt), mlon(100, nt), max_vor(100, nt), s_part(100, nt), mtype(100, nt))
  allocate(minlat(100, nt), minlon(100, nt), z_min(100, nt), z_min_size(100, nt))

  allocate (onelat(1000, nt), onelon(1000, nt), vor_one(1000, nt))


  allocate (mi(100, nt), mj(100, nt))
  allocate (u_vor_f(100, nt), v_vor_f(100, nt))
  allocate (u_vor_b(100, nt), v_vor_b(100, nt))

  allocate (vor_index(10000, nt))


  allocate(dummy(nx1:nx2, ny1:ny2))


  do kt=1, nt


    if(smth_type==1)then 
!      write (*, *)'smth latlon'
      call smth(vor(0:nx, 0:ny, kt), nx, ny, vor_smth(nx1:nx2, ny1:ny2, kt), &
           &nx1, nx2, ny1, ny2, nsmth_x, nsmth_y)
    elseif(smth_type==2)then
      call smth_r(vor(0:nx, 0:ny, kt), nx, ny, lon(0:nx), lat(0:ny), vor_smth(nx1:nx2, ny1:ny2, kt), &
             &nx1, nx2, ny1, ny2, r_smth, proj)
    end if


    write (fname_out, '(A, I4.4, A)')'vor_out_', kt, '.dat'
    open(12, file=fname_out, form='unformatted', access='sequential')    
    write (12)vor(nx1:nx2, ny1:ny2, kt)
    write (12)vor_smth(nx1:nx2, ny1:ny2, kt)
  
  
  
    write (*, '(A, I4.4, A, I2.2, A, I2.2, A, I2.2, A, I2.2)')'Detecting vortex at kt = ', kt

    call vor_partition(vor_smth(nx1:nx2, ny1:ny2, kt), &
         &nx12, ny12, proj, &
         &mlat(:, kt), mlon(:, kt), &
         &max_vor(:, kt), mtype(:, kt), n_max(kt), lat(ny1:ny2), lon(nx1:nx2), &
         &vor_part(nx1:nx2, ny1:ny2, kt), s_part(:, kt), &
         &zeta_max0, zeta_min0, int_zeta_min0, gamma, &
         &d_cf_min, size_synop)


    if(n_max(kt)>=1)then
      call min_z(psea(nx1:nx2, ny1:ny2, kt), &
           &nx12, ny12, proj, &
           &minlat(:, kt), minlon(:, kt), &
           &z_min(:, kt), n_min(kt), lat(ny1:ny2), lon(nx1:nx2), &
           &z_min_size(:, kt), del_psea_min)
    else
      n_min(kt)=0
    endif


    if(maxval(mtype(:, kt))>=1)then

     call synop_check(mlon(:, kt), mlat(:, kt), n_max(kt), minlon(:, kt), minlat(:, kt), n_min(kt), mtype(:, kt), proj, distance_ec)
    end if

 

    vor_part_r=9.99e20
    vor_part_r(nx1, ny2)=-1.

    do j=ny1, ny2
      do i=nx1, nx2
        if(vor_part(i, j, kt)/=0)vor_part_r(i, j)=vor_part(i, j, kt)
      end do
    end do

    write (12)vor_part_r(nx1:nx2,ny1:ny2)

    write (fname_loc,'(A,I4.4,A)')'vormax_loc_',kt,'.txt'
    open(82,file=fname_loc,form='formatted')    


    do i_max=1,n_max(kt)

      mi(i_max,kt)=nint((mlon(i_max,kt)-lons)/lonin)
      mj(i_max,kt)=nint((mlat(i_max,kt)-lats)/latin)
      if(proj==1)then
        write (82,*)mlon(i_max,kt),mlat(i_max,kt),&
             &max_vor(i_max,kt)*1000,nint(s_part(i_max,kt)),mtype(i_max,kt)
      elseif(proj==2)then
        write (82,*)mlon(i_max,kt)/1000,mlat(i_max,kt)/1000,&
             &max_vor(i_max,kt)*1000,nint(s_part(i_max,kt)),mtype(i_max,kt)
      end if
 
     if(proj==1)then
        write (*,*)mlon(i_max,kt),mlat(i_max,kt),&
             &max_vor(i_max,kt)*1000,nint(s_part(i_max,kt)),mtype(i_max,kt)
      elseif(proj==2)then
        write (*,*)mlon(i_max,kt)/1000,mlat(i_max,kt)/1000,&
             &max_vor(i_max,kt)*1000,nint(s_part(i_max,kt)),mtype(i_max,kt)
      end if
 


   end do

    close(82)


! ----- calculate steering wind -----!


    do i_max=1,n_max(kt)


    if(steering_type==1)then 
 
      call steering_wind_f(u,v,&
           &levs,nx,ny,nz,nt,kt,&
           &mi(i_max,kt),mj(i_max,kt),&
           &u_vor_f(i_max,kt),v_vor_f(i_max,kt),&
           &n_steering_x,n_steering_y)

      call steering_wind_b(u,v,&
           &levs,nx,ny,nz,nt,kt,&
           &mi(i_max,kt),mj(i_max,kt),&
           &u_vor_b(i_max,kt),v_vor_b(i_max,kt),&
           &n_steering_x,n_steering_y)


    elseif(steering_type==2)then
      if(kt<nt)then
        call steering_wind_r(u,v,&
             &levs,lon,lat,proj,nx,ny,nz,nt,kt,kt+1,&
             &mi(i_max,kt),mj(i_max,kt),&
             &u_vor_f(i_max,kt),v_vor_f(i_max,kt),&
             &r_steering)
      endif
      if(kt>1)then
        call steering_wind_r(u,v,&
             &levs,lon,lat,proj,nx,ny,nz,nt,kt,kt-1,&
             &mi(i_max,kt),mj(i_max,kt),&
             &u_vor_b(i_max,kt),v_vor_b(i_max,kt),&
             &r_steering)
      endif

    end if


  end do


!---- output steeering wind ----!
    dummy=undef

    do i_max=1,n_max(kt)
      dummy(mi(i_max,kt),mj(i_max,kt))=u_vor_f(i_max,kt)
    end do

    write (12)dummy(nx1:nx2,ny1:ny2)

    dummy=undef

    do i_max=1,n_max(kt)
      dummy(mi(i_max,kt),mj(i_max,kt))=v_vor_f(i_max,kt)
    end do

    write (12)dummy(nx1:nx2,ny1:ny2)


!     dummy=undef

!     do i_max=1,n_max(kt)
!       dummy(mi(i_max,kt),mj(i_max,kt))=-u_vor_b(i_max,kt)
!     end do

!     write (12)dummy(nx1:nx2,ny1:ny2)

!     dummy=undef

!     do i_max=1,n_max(kt)
!       dummy(mi(i_max,kt),mj(i_max,kt))=-v_vor_b(i_max,kt)
!     end do

!     write (12)dummy(nx1:nx2,ny1:ny2)

    write (12)psea(nx1:nx2,ny1:ny2,kt)



  end do

  write (*,*)'Tracking vortex'
  if(track_type==1)then

    call tracking(mlon,mlat,max_vor,mtype,u_vor_f,v_vor_f,u_vor_f,v_vor_f,&
         &nt,n_max,vor_index,vor_num,vor_merge,&
         &vor_part(nx1:nx2,ny1:ny2,1:nt),nx12,ny12,proj,lon(nx1:nx2),lat(ny1:ny2),&
         &del_lon,del_lat,del_t)
  elseif(track_type==2)then
    call tracking2(mlon,mlat,max_vor,mtype,u_vor_f,v_vor_f,u_vor_f,v_vor_f,&
         &nt,n_max,vor_index,vor_num,vor_merge,&
         &vor_part(nx1:nx2,ny1:ny2,1:nt),nx12,ny12,proj,lon(nx1:nx2),lat(ny1:ny2),&
         &del_r,del_t)
  end if
    

  allocate (vortex(vor_num,1:nt,4))
  allocate (vortex_type(vor_num,1:nt))
  allocate (vortex_flag(vor_num))
  vortex=0.

  do i_vor_num=1,vor_num
    do kt=1,nt
      if(vor_index(i_vor_num,kt)>0)then
        if(proj==1)then
          vortex(i_vor_num,kt,1)=mlon(vor_index(i_vor_num,kt),kt)
          vortex(i_vor_num,kt,2)=mlat(vor_index(i_vor_num,kt),kt)
        elseif(proj==2)then
          vortex(i_vor_num,kt,1)=mlon(vor_index(i_vor_num,kt),kt)/1000
          vortex(i_vor_num,kt,2)=mlat(vor_index(i_vor_num,kt),kt)/1000
        endif
        vortex(i_vor_num,kt,3)=max_vor(vor_index(i_vor_num,kt),kt)
        vortex(i_vor_num,kt,4)=s_part(vor_index(i_vor_num,kt),kt)
        vortex_type(i_vor_num,kt)=mtype(vor_index(i_vor_num,kt),kt)
      end if
    end do
  end do

  ! --- check the track ---

  write (*,*)'Check the track'
  
  do i_vor_num=1,vor_num
    call track_check2(vortex(i_vor_num,:,:),vortex_flag(i_vor_num),nt,period_min)
 end do

  
!------------ out put ----------------------


    vor_merge_num=1
    
    do i_vor_num=1,vor_num
      if(vortex_flag(i_vor_num)==1)then
        
        if(vor_merge(i_vor_num)==0)then
          !        write (fname_track,'(A,I4.4,A)')'vortrack_',vor_num_out,'.txt'
          write (fname_track,'(A,I4.4,A,I4.4,A)')'vortrack_',i_vor_num,&
             &'_',1,'.txt'
          !          if(vor_merge_num(i_vor_num)==1)then
          !            vor_merge_num(i_vor_num)=vor_merge_num(i_vor_num)+1          
          !          end if
          
        else
          vor_merge_num(vor_merge(i_vor_num))=&
               &vor_merge_num(vor_merge(i_vor_num))+1          
          
          write (fname_track,'(A,I4.4,A,I4.4,A)')'vortrack_',&
               &vor_merge(i_vor_num),'_'&
               &,vor_merge_num(vor_merge(i_vor_num)),'.txt'
          !          write (*,*)vor_merge(i_vor_num),vor_merge_num(vor_merge(i_vor_num))                    
        end if
      
    
        open(67,file=fname_track,form='formatted')    
        
        do kt=1,nt
          if(vortex(i_vor_num,kt,3)>0.0000001)then
            write (67,'(3f12.5,I6,f15.5,I3)')vortex(i_vor_num,kt,1),&
                 &vortex(i_vor_num,kt,2),vortex(i_vor_num,kt,3)*1000.,&
                 &kt,vortex(i_vor_num,kt,4),vortex_type(i_vor_num,kt)
            
          end if
          
        end do
    
        close (67)

      else
        write (*,*)'vortex ',i_vor_num,'not satisfy the criterion'

      end if
    end do

  write (*, *) 'Owari'

  return
end subroutine tracking_main