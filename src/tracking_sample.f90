program main
  implicit none   
  real(4)::p0(100)
  real(4),allocatable::levs(:),lon(:),lat(:)
  real(4),allocatable::vor(:,:,:),psea(:,:,:),vor_in(:,:,:)
  real(4),allocatable::u(:,:,:,:),v(:,:,:,:)

  integer (4)::i_rec

  real(4)::lats,lons
  real(4)::latin,lonin
  real(4)::del_t
  integer(4)::nx,ny,nz
  integer ::proj,vert_grid
  integer (4)::nt
  integer (4)::nx1,nx2,ny1,ny2
  integer(4)::i,j,k
  integer (4)::ios

  integer,allocatable ::yyyy(:),mm(:),dd(:),hh(:),mn(:)
  integer ::kt

  character (50)::fname
  character (50)::fname_in="msm_pv_slp_"
  
 
! parameter for smoothing of vorticity
  integer::smth_type
  integer::nsmth_x,nsmth_y
  real(4)::r_smth

! parameter for detecting vortex
  real(4)::zeta_max0,zeta_min0,int_zeta_min0,gamma

! parameter for excluding the synoptic scale disturbances
  real(4)::d_cf_min,size_synop,del_psea_min,distance_ec
  
! parameter for calculating steering winds
  integer::steering_type
  integer::n_steering_x,n_steering_y
  real(4)::r_steering

! parameter for linking vortex
  integer::track_type
  real(4)::del_lon,del_lat,del_r

! parameter for checking the track
  integer::period_min



  proj=1
  vert_grid=1

  nx=240
  ny=252  
  nt=248
  nz=8

  nx1=60
  nx2=176
  ny1=106
  ny2=246

  lons=120.0
  lats=22.4
  lonin=0.125
  latin=0.1  

  del_t=10800.0
  
  ! parameter for smoothing of vorticity
  smth_type=2
  nsmth_x=10
  nsmth_y=10
  r_smth=30.0
  
  ! parameter for detecting vortex
  zeta_max0=2.0e-4
  zeta_min0=1.5e-4
  int_zeta_min0=0.02e-4
  gamma=0.25
  
! parameter for excluding the synoptic scale disturbances
  d_cf_min=400.0
  size_synop=40000.0
  distance_ec=300.0
  del_psea_min=0.5
  
! parameter for calculating steering winds
  steering_type=2
  n_steering_x=20
  n_steering_y=20
  r_steering=200.0

! parameter for linking vortex
  track_type=2
  del_lon=1.0
  del_lat=0.8
  del_r=120.0

! parameter for checking the track
  period_min=3


  ! Grid 
  allocate(lon(0:nx))
  allocate(lat(0:ny))
  allocate(levs(1:nz))

  do i=0,nx
    lon(i)=lons+lonin*i
  end do
  do j=0,ny
    lat(j)=lats+latin*j
  end do

  levs=(/ 1000, 975, 950, 925, 900, 850, 800, 700 /)


  ! Time
  allocate (yyyy(nt),mm(nt),dd(nt),hh(nt),mn(nt))

  write (*,*)"Read timecard"
  open (40,file='timecard',form='formatted')
  do kt=1,nt
    read(40,*)yyyy(kt),mm(kt),dd(kt),hh(kt),mn(kt)
  end do
  
  close (40)
  

  ! Read data

  allocate(vor(0:nx,0:ny,nt))
  allocate(u(0:nx,0:ny,nz,nt))
  allocate(v(0:nx,0:ny,nz,nt))
  allocate(psea(0:nx,0:ny,nt))



  write (*,*)'Read MSM data',nt

  do kt=1,nt
    i_rec=1

    write (fname,'(A,I4.4,I2.2,I2.2,I2.2,I2.2,A)')trim(fname_in),yyyy(kt),mm(kt),dd(kt),hh(kt),mn(kt),'.dat'


    open(11,file=fname,form='unformatted',access='direct',status='old',iostat=ios,recl=4*(nx+1)*(ny+1))    

    if(ios==0)then

      write (0,'(A,I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2)')&
           &'Reading data at ',yyyy(kt),' ',mm(kt),' ',dd(kt),' ',hh(kt),' ',mn(kt)
      
      read(11,rec=i_rec) psea(0:nx,0:ny,kt)

      i_rec=i_rec+1
      do k=1,16
        i_rec=i_rec+1

      end do

      do k=1,16
        if(k<=nz)read(11,rec=i_rec) u(0:nx,0:ny,k,kt)
        i_rec=i_rec+1
      end do
      do k=1,16
        if(k<=nz)read(11,rec=i_rec) v(0:nx,0:ny,k,kt)
        i_rec=i_rec+1
      end do
      !        write (0,*)'Read v  '

      do k=1,16    
        i_rec=i_rec+1   !vvel
      end do
      do k=1,16
        i_rec=i_rec+1   !t
      end do
      do k=1,12
        i_rec=i_rec+1   !rh
      end do
      do k=1,16
        i_rec=i_rec+1   !pt
      end do
      do k=1,16
        i_rec=i_rec+1   !sigmai
      end do

      do k=1,16
        if(k==3)read(11,rec=i_rec) vor(0:nx,0:ny,kt)
        i_rec=i_rec+1
      end do
      !        write (0,*)'Read vor  '


      close(11)    

    end if
  end do




  call tracking_main(vor,u,v,psea,&
     &proj,vert_grid,&
     &nx,ny,nx1,nx2,ny1,ny2,nz,levs,nt,&
     &lons,lats,lonin,latin,del_t,&
     &nsmth_x,nsmth_y,r_smth,smth_type,&
     &zeta_max0,zeta_min0,int_zeta_min0,gamma,&
     &n_steering_x,n_steering_y,r_steering,steering_type,&
     &del_lon,del_lat,del_r,track_type,&
     &period_min,d_cf_min,size_synop,del_psea_min,distance_ec)



  stop
end program main

subroutine land(var,nx,ny,land_flag,sl,sl_land)
  implicit none 
  integer ,intent (in)::nx,ny
  real (4),intent (in)::var(0:nx,0:ny)
  real (4),intent (in)::sl(0:nx,0:ny)  
  real(4),intent (in)::sl_land
  logical,intent (out)::land_flag(0:nx,0:ny)
  real(4)::sl_tmp
  

  integer ::i,j,ii,jj
  real (4),parameter ::undef=9.99e20

  land_flag=.false.

  do j=0,ny
    do i=0,nx
      !        sl_tmp=0.
      !        do ii=-2,2
      !           do  jj=-2,2
      !              sl_tmp=sl_tmp+sl(i+ii,j+jj)
      !           end do
      !        end do
      !        if(sl_tmp>0.0001)then

      if(sl(i,j)>sl_land)then
        land_flag(i,j)=.true.
!        var(i,j)=undef
      end if

    end do
  end do

  return
end subroutine land

subroutine japan_sea(var,nx,ny,land_flag,lon,lat)
  implicit none 
  integer ,intent (in)::nx,ny
  real (4),intent (in)::var(0:nx,0:ny)
  real (4),intent (in)::lon(0:nx),lat(0:ny)  
  logical,intent (out)::land_flag(0:nx,0:ny)
  real(4)::sl_tmp

  integer ::i,j,ii,jj
  real (4),parameter ::undef=9.99e20

  !  land_flag=.false.

  do j=0,ny
    do i=0,nx

      if(lon(i)>=132..and.lat(j)<=34)then
        land_flag(i,j)=.true.
!        var(i,j)=undef
      elseif(lon(i)>=136..and.lat(j)<=35)then
        land_flag(i,j)=.true.
!        var(i,j)=undef
      elseif(lon(i)>=141..and.lat(j)<=42)then
        land_flag(i,j)=.true.
!        var(i,j)=undef
      elseif(lon(i)>=138..and.lat(j)<=36.9)then
        land_flag(i,j)=.true.
!        var(i,j)=undef
      elseif(lon(i)>=140..and.lat(j)>=42.and.lat(j)<=42.5)then
        land_flag(i,j)=.true.
!        var(i,j)=undef
      elseif(lon(i)>=131..and.lat(j)<=34)then
        land_flag(i,j)=.true.
!        var(i,j)=undef
      elseif(lon(i)>=133..and.lat(j)<=35)then
        land_flag(i,j)=.true.

      elseif(lon(i)>=140.2.and.lat(j)<=42)then
        land_flag(i,j)=.true.

      elseif(lon(i)>=130..and.lon(i)<=135.and.lat(j)>=44)then
        land_flag(i,j)=.true.

      elseif(lon(i)<129..and.lat(j)<=35)then
        land_flag(i,j)=.true.

      elseif(lat(j)<33.5)then
        land_flag(i,j)=.true.

        



!        var(i,j)=undef



      end if

      if(lon(i)>=130.0.and.lon(i)<=132.0.and.lat(j)>=37.0.and.lat(j)<=38.)then
        land_flag(i,j)=.false.
      elseif(lon(i)>=132.5.and.lon(i)<=134.0.and.lat(j)>=35.8.and.lat(j)<=37.)then
        land_flag(i,j)=.false.
      elseif(lon(i)>=129.6.and.lon(i)<=129.9.and.lat(j)>=33.58.and.lat(j)<=34.)then
        land_flag(i,j)=.false.
      elseif(lon(i)>=128.95.and.lon(i)<=130.0.and.lat(j)>=33.9.and.lat(j)<=34.8)then
        land_flag(i,j)=.false.
      elseif(lon(i)>=136.6.and.lon(i)<=137.4.and.lat(j)>=36.8.and.lat(j)<=38.0)then
        land_flag(i,j)=.false.
      elseif(lon(i)>=138.0.and.lon(i)<=138.65.and.lat(j)>=37.75.and.lat(j)<=38.4)then
        land_flag(i,j)=.false.
      elseif(lon(i)>=139.and.lon(i)<=139.8.and.lat(j)>=41.95.and.lat(j)<=42.5)then
        land_flag(i,j)=.false.
      elseif(lon(i)>=140.5.and.lon(i)<=141.5.and.lat(j)>=44.9.and.lat(j)<=45.6)then
        land_flag(i,j)=.false.
      end if

    end do
  end do

  return
end subroutine japan_sea



subroutine undef_in(var,nx,ny,land_flag)
  implicit none 
  integer ,intent (in)::nx,ny
  real (4),intent (inout)::var(0:nx,0:ny)
  logical,intent (in)::land_flag(0:nx,0:ny)
  
  integer (4)::i,j

  real (4),parameter ::undef=-100.0


  do j=0,ny
    do i=0,nx
      if(land_flag(i,j))var(i,j)=undef
    end do
  end do

  return
end subroutine undef_in


