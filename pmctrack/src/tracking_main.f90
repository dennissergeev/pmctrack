subroutine tracking_main(vor,u,v,psea,&
     &nx,ny,nz,levs,nt,&
     &lons,lats,lonin,latin,del_t)

  use types, only : wp
  use constants, only : ikilo, rkilo, fillval, nmax, kmax, pmax
  use params, only : proj, vert_grid, nx1, nx2, ny1, ny2, &
    & smth_type, steering_type, track_type, outdir

  implicit none

  integer,intent (in) ::nx,ny,nz,nt
  real(4),intent (in)::vor(0:nx,0:ny,1:nt)
  real(4),intent (in)::u(0:nx,0:ny,1:nz,1:nt),v(0:nx,0:ny,1:nz,1:nt)
  real(4),intent (in)::psea(0:nx,0:ny,1:nt)
  real(4),intent (in)::lons,lats,lonin,latin,levs,del_t

! LOCAL VARIABLES
  real(4),allocatable::lon(:),lat(:)

  real(4),allocatable::vor_smth(:,:,:)
  real(4),allocatable::dummy(:,:)
  integer (4),allocatable::vor_part(:,:,:)
  real (4),allocatable::vor_part_r(:,:)

  integer (4)::nx12,ny12
  integer(4)::i,j
  integer (4)::kt
  real(4),allocatable::mlat(:,:),mlon(:,:),max_vor(:,:)
  real(4),allocatable::minlat(:,:),minlon(:,:),z_min(:,:)
  real(4),allocatable::s_part(:,:)
  integer (4),allocatable::mi(:,:),mj(:,:)
  integer (4),allocatable::mtype(:,:),z_min_size(:,:)

  real(4),allocatable::onelat(:,:),onelon(:,:),vor_one(:,:)

  real(4),allocatable::u_vor_f(:,:),v_vor_f(:,:)
  real(4),allocatable::u_vor_b(:,:),v_vor_b(:,:)
  integer (4),allocatable::n_max(:),n_min(:)
  integer (4)::i_max

  integer (4)::vor_num,i_vor_num
  integer (4)::vor_merge(pmax),vor_merge_num(pmax)
  integer (4),allocatable::vor_index(:,:)

  real(4),allocatable::vortex(:,:,:)
  integer (4),allocatable::vortex_type(:,:)
  integer (4),allocatable::vortex_flag(:)

  character (100) :: fname_out, fname_loc, fname_track


  allocate(lon(0:nx))
  allocate(lat(0:ny))

  do i=0,nx
    lon(i)=lons+lonin*i
  end do
  do j=0,ny
    lat(j)=lats+latin*j
  end do


  nx12=nx2-nx1
  ny12=ny2-ny1



  allocate(vor_smth(nx1:nx2,ny1:ny2,nt))
  allocate(vor_part(nx1:nx2,ny1:ny2,nt))
  allocate(vor_part_r(nx1:nx2,ny1:ny2))
  allocate (n_max(nt),n_min(nt))
  allocate (mlat(nmax,nt),mlon(nmax,nt),max_vor(nmax,nt),s_part(nmax,nt),mtype(nmax,nt))
  allocate (minlat(nmax,nt),minlon(nmax,nt),z_min(nmax,nt),z_min_size(nmax,nt))

  allocate (onelat(kmax,nt),onelon(kmax,nt),vor_one(kmax,nt))


  allocate (mi(nmax,nt),mj(nmax,nt))
  allocate (u_vor_f(nmax,nt),v_vor_f(nmax,nt))
  allocate (u_vor_b(nmax,nt),v_vor_b(nmax,nt))

  allocate (vor_index(pmax,nt))


  allocate(dummy(nx1:nx2,ny1:ny2))


  do kt=1,nt


    if(smth_type==1)then
!      write (*,*)'smth latlon'
      call smth(vor(0:nx,0:ny,kt),nx,ny,vor_smth(nx1:nx2,ny1:ny2,kt))
    elseif(smth_type==2)then
      call smth_r(vor(0:nx,0:ny,kt),nx,ny,lon(0:nx),lat(0:ny),vor_smth(nx1:nx2,ny1:ny2,kt))
    else
      vor_smth(nx1:nx2,ny1:ny2,kt) = vor(nx1:nx2,ny1:ny2,kt)
    end if


    write (fname_out,'(A,A,A,I4.4,A)')trim(outdir),'/','vor_out_',kt,'.dat'
    open(12,file=fname_out,form='unformatted',access='sequential')
    write (12)vor(nx1:nx2,ny1:ny2,kt)
    write (12)vor_smth(nx1:nx2,ny1:ny2,kt)



    write (*,'(A,I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2)')'Detecting vortex at kt = ',kt

    call vor_partition(vor_smth(nx1:nx2,ny1:ny2,kt),&
         &nx12,ny12,&
         &mlat(:,kt),mlon(:,kt),&
         &max_vor(:,kt),mtype(:,kt),n_max(kt),lat(ny1:ny2),lon(nx1:nx2),&
         &vor_part(nx1:nx2,ny1:ny2,kt),s_part(:,kt))
    if(n_max(kt)>=1)then
      call min_z(psea(nx1:nx2,ny1:ny2,kt),&
           &nx12,ny12,&
           &minlat(:,kt),minlon(:,kt),&
           &z_min(:,kt),n_min(kt),lat(ny1:ny2),lon(nx1:nx2),&
           &z_min_size(:,kt))
    else
      n_min(kt)=0
    endif


    if(maxval(mtype(:,kt))>=1)then

     call synop_check(mlon(:,kt),mlat(:,kt),n_max(kt),minlon(:,kt),minlat(:,kt),n_min(kt),mtype(:,kt))
    end if



    vor_part_r = fillval

    vor_part_r(nx1,ny2)=-1.

    do j=ny1,ny2
      do i=nx1,nx2
        if(vor_part(i,j,kt)/=0)vor_part_r(i,j)=vor_part(i,j,kt)
      end do
    end do

    write (12)vor_part_r(nx1:nx2,ny1:ny2)

    write (fname_loc,'(A,A,A,I4.4,A)')trim(outdir),'/','vormax_loc_',kt,'.txt'
    open(82,file=fname_loc,form='formatted')


    do i_max=1,n_max(kt)

      mi(i_max,kt)=nint((mlon(i_max,kt)-lons)/lonin)
      mj(i_max,kt)=nint((mlat(i_max,kt)-lats)/latin)
      if(proj==1)then
        write (82,*)mlon(i_max,kt),mlat(i_max,kt),&
             &max_vor(i_max,kt)*ikilo,nint(s_part(i_max,kt)),mtype(i_max,kt)
      elseif(proj==2)then
        write (82,*)mlon(i_max,kt)/ikilo,mlat(i_max,kt)/ikilo,&
             &max_vor(i_max,kt)*ikilo,nint(s_part(i_max,kt)),mtype(i_max,kt)
      end if

     if(proj==1)then
        write (*,*)mlon(i_max,kt),mlat(i_max,kt),&
             &max_vor(i_max,kt)*ikilo,nint(s_part(i_max,kt)),mtype(i_max,kt)
      elseif(proj==2)then
        write (*,*)mlon(i_max,kt)/ikilo,mlat(i_max,kt)/ikilo,&
             &max_vor(i_max,kt)*ikilo,nint(s_part(i_max,kt)),mtype(i_max,kt)
      end if



   end do

    close(82)


! ----- calculate steering wind -----!


    do i_max=1,n_max(kt)


    if(steering_type==1)then

      call steering_wind_f(u,v,&
           &levs,nx,ny,nz,nt,kt,&
           &mi(i_max,kt),mj(i_max,kt),&
           &u_vor_f(i_max,kt),v_vor_f(i_max,kt))

      call steering_wind_b(u,v,&
           &levs,nx,ny,nz,nt,kt,&
           &mi(i_max,kt),mj(i_max,kt),&
           &u_vor_b(i_max,kt),v_vor_b(i_max,kt))


    elseif(steering_type==2)then
      if(kt<nt)then
        call steering_wind_r(u,v,&
             &levs,lon,lat,nx,ny,nz,nt,kt,kt+1,&
             &mi(i_max,kt),mj(i_max,kt),&
             &u_vor_f(i_max,kt),v_vor_f(i_max,kt))
      endif
      if(kt>1)then
        call steering_wind_r(u,v,&
             &levs,lon,lat,nx,ny,nz,nt,kt,kt-1,&
             &mi(i_max,kt),mj(i_max,kt),&
             &u_vor_b(i_max,kt),v_vor_b(i_max,kt))
      endif

    end if


  end do


!---- output steeering wind ----!
    dummy=fillval

    do i_max=1,n_max(kt)
      dummy(mi(i_max,kt),mj(i_max,kt))=u_vor_f(i_max,kt)
    end do

    write (12)dummy(nx1:nx2,ny1:ny2)

    dummy=fillval

    do i_max=1,n_max(kt)
      dummy(mi(i_max,kt),mj(i_max,kt))=v_vor_f(i_max,kt)
    end do

    write (12)dummy(nx1:nx2,ny1:ny2)


!     dummy=fillval

!     do i_max=1,n_max(kt)
!       dummy(mi(i_max,kt),mj(i_max,kt))=-u_vor_b(i_max,kt)
!     end do

!     write (12)dummy(nx1:nx2,ny1:ny2)

!     dummy=fillval

!     do i_max=1,n_max(kt)
!       dummy(mi(i_max,kt),mj(i_max,kt))=-v_vor_b(i_max,kt)
!     end do

!     write (12)dummy(nx1:nx2,ny1:ny2)

    write (12)psea(nx1:nx2,ny1:ny2,kt)



  end do

  write (*, *) 'Linking vortices'
  if(track_type==1)then
    call linkin_vort(mlon,mlat,mtype,u_vor_f,v_vor_f,&
         &nt,n_max,vor_index,vor_num,vor_merge,&
         &vor_part(nx1:nx2,ny1:ny2,1:nt),nx12,ny12,lon(nx1:nx2),lat(ny1:ny2),&
         del_t)
  elseif(track_type==2)then
    call linkin_vort2(mlon,mlat,mtype,u_vor_f,v_vor_f,&
         &nt,n_max,vor_index,vor_num,vor_merge,&
         &vor_part(nx1:nx2,ny1:ny2,1:nt),nx12,ny12,lon(nx1:nx2),lat(ny1:ny2),&
         del_t)
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
          vortex(i_vor_num,kt,1)=mlon(vor_index(i_vor_num,kt),kt)/ikilo
          vortex(i_vor_num,kt,2)=mlat(vor_index(i_vor_num,kt),kt)/ikilo
        endif
        vortex(i_vor_num,kt,3)=max_vor(vor_index(i_vor_num,kt),kt)
        vortex(i_vor_num,kt,4)=s_part(vor_index(i_vor_num,kt),kt)
        vortex_type(i_vor_num,kt)=mtype(vor_index(i_vor_num,kt),kt)
      end if
    end do
  end do

  ! --- check the track ---

  write (*,*) 'Check the track'

  do i_vor_num=1,vor_num
    call track_check2(vortex(i_vor_num,:,:),vortex_flag(i_vor_num),nt)
 end do


!------------ out put ----------------------


    vor_merge_num=1

    do i_vor_num=1,vor_num
      if(vortex_flag(i_vor_num)==1)then

        if(vor_merge(i_vor_num)==0)then
          !        write (fname_track,'(A,I4.4,A)')'vortrack_',vor_num_out,'.txt'
          write (fname_track,'(A,A,A,I4.4,A,I4.4,A)')trim(outdir),'/','vortrack_',i_vor_num,&
             &'_',1,'.txt'
          !          if(vor_merge_num(i_vor_num)==1)then
          !            vor_merge_num(i_vor_num)=vor_merge_num(i_vor_num)+1
          !          end if

        else
          vor_merge_num(vor_merge(i_vor_num))=&
               &vor_merge_num(vor_merge(i_vor_num))+1

          write (fname_track,'(A,A,A,I4.4,A,I4.4,A)')trim(outdir),'/','vortrack_',&
               &vor_merge(i_vor_num),'_'&
               &,vor_merge_num(vor_merge(i_vor_num)),'.txt'
          !          write (*,*)vor_merge(i_vor_num),vor_merge_num(vor_merge(i_vor_num))
        end if


        open(67,file=fname_track,form='formatted')

        do kt=1,nt
          if(vortex(i_vor_num,kt,3)>0.0000001)then
            write (67,'(3f12.5,I6,f15.5,I3)')vortex(i_vor_num,kt,1),&
                 &vortex(i_vor_num,kt,2),vortex(i_vor_num,kt,3)*rkilo,&
                 &kt,vortex(i_vor_num,kt,4),vortex_type(i_vor_num,kt)

          end if

        end do

        close (67)

      else
        write (*,*)'vortex ',i_vor_num,'not satisfy the criterion'

      end if
    end do

  write (*,*)'Owari'

  return
end subroutine tracking_main
