subroutine cf_synop_check(vor_in,vor_part,n_part,nx,ny,proj,lon,lat,mtype_part,d_cf_min,size_synop)

  use constants
  use params

  implicit none 
  integer (4),intent (in)::nx,ny
  integer (4),intent (in)::proj
  integer (4),intent (in)::vor_part(0:nx,0:ny)
  integer (4),intent (in)::n_part

  integer (4),intent (out)::mtype_part

  real(4),intent (in)::lon(0:nx),lat(0:ny)
  real(4),intent (in)::vor_in(0:nx,0:ny)
  real(4),intent (in)::d_cf_min,size_synop


  real(4)::lonin,latin

  !----cold front
  integer (4)::pnum
  integer (4)::j_n,j_s,i_w,i_e,i_n,i_s
  logical ::flag_coldfront
  logical::flag_one
  character (len=80)::fname_one
  real (4)::d_cf,theta_d_cf

  real(4)::slope,width,length
  real(4)::one(pmax,1:2)
  integer (4)::one_num

  real(4)::a,b,c,k,r2

  !----synoptic
  real (4)::size_vor

  integer (4)::i,j

  lonin=lon(1)-lon(0)
  latin=lat(1)-lat(0)


!--------------- Check cold front ----------------
!    flag_coldfront=.false.


!    write (*,*)fname_one,kt,yyyy(kt),mm(kt),dd(kt),hh(kt),mn(kt)


    j_n=0
    j_s=ny

  
    size_vor=0.

    one_num=0
    one(:,:)=0

    do j=0,ny
      do i=0,nx
        

        if(vor_part(i,j)==n_part)then
          if(proj==1)then
            size_vor=size_vor&
                 &+lonin*(pi/180.)*ra/rkilo*cos(lat(j)*pi/180.)*latin*(pi/180.)*ra/rkilo
          elseif(proj==2)then
            size_vor=size_vor&
                 &+lonin*latin*1.0e-6
          end if


          call median_check(vor_in(i-1:i+1,j-1:j+1),flag_one,pnum)
          if(flag_one)then
            if(pnum>=6)then
!              write (81,*)lon(i),lat(j),pnum,n_part
              one_num=one_num+1
              one(one_num,1)=lon(i)
              one(one_num,2)=lat(j)
            end if
          end if



          if(j<j_s)then
            j_s=j
            i_s=i
          end if
          if(j>=j_n)then
            j_n=j
            i_n=i
          end if
        end if
      end do
    end do

    if(proj==1)then
      if(i_s/=i_n)then
        d_cf=ra*acos(cos(pi/180*lat(j_n))*cos(pi/180*lat(j_s))&
             &*cos(pi/180*(lon(i_s)-lon(i_n)))&
             &+sin(pi/180*lat(j_n))*sin(pi/180*lat(j_s)))
        
        theta_d_cf=acos((ra*(lat(j_n)-lat(j_s))*pi/180.)/d_cf)
      else
        d_cf=ra*pi/180.0*(lat(j_n)-lat(j_s))
        theta_d_cf=0.0
        
      end if

    elseif(proj==2)then
      d_cf=sqrt((lon(i_s)-lon(i_n))**2+ (lat(j_s)-lat(j_n))**2)*1.0e-3
      theta_d_cf=(lat(j_s)-lat(j_n))*1.0e-3/d_cf      
    end if

!      write (82,*)n_part
!      write (82,*)d_cf*1.0e-3,theta_d_cf/3.1415*180.,size_vor,size_vor/(d_cf*1.0e-3)
!      write (82,*)lon(i_n),lat(j_n),lon(i_s),lat(j_s)

    if(d_cf>=d_cf_min*1.0e3)then
      if(lon(i_n)>lon(i_s).and.theta_d_cf<=1.*pi/3.)then
        call quadric_fit(one(1:one_num,1:2),one_num,a,b,c,k,r2)
!        write (82,*)'quad',a,b,c,k,r2

        if(r2>=0.8.and.k<=0.1)then

!          write (20,*)'cold front' ,0,d_cf,theta_d_cf*180./pi,lon(i_n),lat(j_n),lon(i_s),lat(j_s)
          mtype_part=1
        end if

      elseif(theta_d_cf<=pi/9.)then
        call quadric_fit(one(1:one_num,1:2),one_num,a,b,c,k,r2)
!        write (82,*)'quad',a,b,c,k,r2

        if(r2>=0.8.and.k<=0.1)then
!          write (20,*)'cold front' ,1,d_cf,theta_d_cf*180./pi,lon(i_n),lat(j_n),lon(i_s),lat(j_s)
          mtype_part=1
        end if

      end if

    end if


!--------------- Check synoptic low ----------------


!  write (*,*)size_vor


  if(size_vor>=size_synop)then
!    mtype(n_part,2)=2
    mtype_part=mtype_part+2
!    write (20,*)'synoptic low'
  endif
!    if(s_part>=s_part_min)then

  return
end subroutine cf_synop_check

subroutine median_check(val,flag,pnum)
  implicit none 
  real (4),intent (in)::val(-1:1,-1:1)
  logical (4),intent (out)::flag
  integer (4),intent (out)::pnum
  integer (4)::m
  integer (4)::mx(8),my(8)

  mx(1:8)=(/ 1, 1, 0,-1,-1,-1, 0, 1/)
  my(1:8)=(/ 0, 1, 1, 1, 0,-1,-1,-1/)


  flag=.false.
  pnum=0
!  mnum=0
  
  do m=1,8
    if(val(0,0)>val(mx(m),my(m)))then
      pnum=pnum+1
!    elseif(val(0,0)<val(mx(m),my(m)))then
!      mnum=mnum+1
    end if
  end do

!  if(pnum/=4.or.mnum/=4)flag=.true.
  if(pnum>=5)flag=.true.

  return
end subroutine median_check


subroutine quadric_fit(one,one_num,a,b,c,k,r2)
  implicit none 
  integer (4),intent (in)::one_num
  real(4),intent (in)::one(1:one_num,1:2)
  real (4),intent (out)::a,b,c,k,r2
  real(4)::r
  integer (4)::i
  
  real(8)::sx,sx2,sx3,sx4,sx2y,sxy,sy,s1

  real(4)::xbar,ybar

  real(4)::x(1:one_num),y(1:one_num)

  real(4)::fit(1:one_num)

  sx=0.
  sx2=0.
  sx3=0.
  sx4=0.
  sx2y=0.
  sxy=0.
  sy=0.
  s1=0.


  xbar=0.
  ybar=0.
  
  do i=1,one_num
    xbar=xbar+one(i,2)
    ybar=ybar+one(i,1)
  end do
  
  xbar=xbar/one_num
  ybar=ybar/one_num

  do i=1,one_num
    x(i)=one(i,2)-xbar
    y(i)=one(i,1)-ybar
  end do

  do i = 1,one_num
    sx4=sx4+x(i)**4
    sx3=sx3+x(i)**3
    sx2=sx2+x(i)**2
    sx=sx+x(i)**1
    sx2y=sx2y+x(i)**2*y(i)
    sxy=sxy+x(i)*y(i)
    sy=sy+y(i)
    s1=s1+1.0
  end do

!  write (71,*)kt,yyyy(kt),mm(kt),dd(kt),hh(kt),mn(kt)
!  write (71,*)xbar,ybar

!  write (71,*)sx,sx2,sx3,sx4,sx2y,sxy,sy,s1



  a=(s1*sx2*sx2y-sx*sx*sx2y+sx*sx2*sxy-s1*sx3*sxy+sx*sx3*sy-sx2*sx2*sy)&
       &/(2*sx*sx2*sx3+s1*sx2*sx4-sx*sx*sx4-s1*sx3*sx3-sx2*sx2*sx2)

  b=(sx*sx2*sx2y-s1*sx3*sx2y+s1*sx4*sxy-sx2*sx2*sxy+sx2*sx3*sy-sx*sx4*sy)&
       &/(2*sx*sx2*sx3+s1*sx2*sx4-sx*sx*sx4-s1*sx3*sx3-sx2*sx2*sx2)

  c=(-sx2*sx2*sx2y+sx*sx3*sx2y-sx*sx4*sxy+sx2*sx3*sxy-sx3*sx3*sy+sx2*sx4*sy)&
       &/(2*sx*sx2*sx3+s1*sx2*sx4-sx*sx*sx4-s1*sx3*sx3-sx2*sx2*sx2)

!  write (71,*)a,b,c
  
  do i=1,one_num
    fit(i)=a*x(i)**2+b*x(i)+c
  end do


  call correlation(y(1:one_num), fit(1:one_num), one_num, r2)

!  write (71,*)fit(1)+ybar,y(1)+ybar,r2



!  r2=xy*xy/(x2*y2)
  
!  r2=r*r
  
  k=0.
  do i=1,one_num
    k=k+a/(1+(2*a*x(i)+b)**2.)**1.5
  end do

  k=k/one_num
  
  return
end subroutine quadric_fit

subroutine correlation(x, y, n, r2)

  implicit none
  integer,intent (in) :: n
  real(4),intent (in)    :: x(n), y(n)
  real(4),intent (out)::r2

  real :: xm, ym, sigx, sigy, cor

  real(kind=8) :: sx, sy, sxy, sx2, sy2
  integer (4)::i


  sx=0; sy=0; sxy=0; sx2=0; sy2=0
  
  
  if( n > 2 ) then
    do i=1, n
      sx = sx +x(i)
      sy = sy +y(i)
    end do
    xm = sx/n; ym = sy/n
    do i=1, n
      sx2 = sx2 + (x(i)-xm)*(x(i)-xm)
      sy2 = sy2 + (y(i)-ym)*(y(i)-ym)
      sxy = sxy + (x(i)-xm)*(y(i)-ym)
    end do
    
    r2=sxy*sxy/(sx2*sy2)

!    write (71,*)'cor',xm,ym,sx2,sy2,sxy

    sigx = sqrt(sx2/n); sigy = sqrt(sy2/n)
    if( sigx > 0 .AND. sigy > 0 ) then
      cor = (sxy/n)/(sigx*sigy)
    end if
  end if
end subroutine correlation


