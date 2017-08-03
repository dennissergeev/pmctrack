subroutine steering_wind_f(u,v,p,nx,ny,nz,nt,kt,mi,mj,u_vor,v_vor,n_steering_x,n_steering_y)

  implicit none

  integer ,intent (in)::nx,ny,nz,nt,kt
  integer ,intent (in)::n_steering_x,n_steering_y
  real(4) ,intent (in)::u(0:nx,0:ny,nz,nt),v(0:nx,0:ny,nz,nt)
  real(4) ,intent (in)::p(nz)
  
  integer (4),intent (in)::mi,mj
  real (4),intent (out)::u_vor,v_vor
  real(4) ::u_t0t1(0:nx,0:ny,nz),v_t0t1(0:nx,0:ny,nz)
  real(4) ::u_int(0:nx,0:ny),v_int(0:nx,0:ny)

  integer ::n_steering_s
  integer ::ii,jj


  u_vor=0.
  v_vor=0.
  n_steering_s=0


    if(kt>=1.and.kt<nt)then
       u_t0t1(0:nx,0:ny,1:nz)=0.5*(u(0:nx,0:ny,1:nz,kt)+u(0:nx,0:ny,1:nz,kt+1))
       v_t0t1(0:nx,0:ny,1:nz)=0.5*(v(0:nx,0:ny,1:nz,kt)+v(0:nx,0:ny,1:nz,kt+1))
    else
       u_t0t1(0:nx,0:ny,1:nz)=u(0:nx,0:ny,1:nz,kt)
       v_t0t1(0:nx,0:ny,1:nz)=v(0:nx,0:ny,1:nz,kt)
    end if


  do jj=-n_steering_y,n_steering_y
     do ii=-n_steering_x,n_steering_x
        if(mi+ii>=0.and.mi+ii<=nx.and.mj+jj>=0.and.mj+jj<=ny)then
           call integral_p(u_t0t1(mi+ii,mj+jj,1:nz),u_int(mi+ii,mj+jj),nz,p)
           call integral_p(v_t0t1(mi+ii,mj+jj,1:nz),v_int(mi+ii,mj+jj),nz,p)
           u_vor=u_vor+u_int(mi+ii,mj+jj)
           v_vor=v_vor+v_int(mi+ii,mj+jj)
           n_steering_s=n_steering_s+1
           
        end if
     end do
  end do

  u_vor=u_vor/n_steering_s
  v_vor=v_vor/n_steering_s

!  write (*,*)u_vor,v_vor

  return
end subroutine steering_wind_f


subroutine steering_wind_b(u,v,p,nx,ny,nz,nt,kt,mi,mj,u_vor,v_vor,n_steering_x,n_steering_y)

  implicit none

  integer ,intent (in)::nx,ny,nz,nt,kt
  real(4) ,intent (in)::u(0:nx,0:ny,nz,nt),v(0:nx,0:ny,nz,nt)
  real(4) ,intent (in)::p(nz)
  integer (4),intent (in)::mi,mj
  real (4),intent (out)::u_vor,v_vor
  real(4) ::u_t0t1(0:nx,0:ny,nz),v_t0t1(0:nx,0:ny,nz)
  real(4) ::u_int(0:nx,0:ny),v_int(0:nx,0:ny)

  integer ,intent (in)::n_steering_x,n_steering_y
  integer ::n_steering_s
  integer ::ii,jj

  
  u_vor=0.
  v_vor=0.
  n_steering_s=0


    if(kt>1.and.kt<=nt)then
       u_t0t1(0:nx,0:ny,1:nz)=0.5*(u(0:nx,0:ny,1:nz,kt)+u(0:nx,0:ny,1:nz,kt-1))
       v_t0t1(0:nx,0:ny,1:nz)=0.5*(v(0:nx,0:ny,1:nz,kt)+v(0:nx,0:ny,1:nz,kt-1))
    else
       u_t0t1(0:nx,0:ny,1:nz)=u(0:nx,0:ny,1:nz,kt)
       v_t0t1(0:nx,0:ny,1:nz)=v(0:nx,0:ny,1:nz,kt)
    end if

  do jj=-n_steering_y,n_steering_y
     do ii=-n_steering_x,n_steering_x
        if(mi+ii>=0.and.mi+ii<=nx.and.mj+jj>=0.and.mj+jj<=ny)then
           call integral_p(u_t0t1(mi+ii,mj+jj,1:nz),u_int(mi+ii,mj+jj),nz,p)
           call integral_p(v_t0t1(mi+ii,mj+jj,1:nz),v_int(mi+ii,mj+jj),nz,p)
           u_vor=u_vor+u_int(mi+ii,mj+jj)
           v_vor=v_vor+v_int(mi+ii,mj+jj)
           n_steering_s=n_steering_s+1
           
        end if
     end do
  end do

  u_vor=u_vor/n_steering_s
  v_vor=v_vor/n_steering_s

!  write (*,*)u_vor,v_vor

  return
end subroutine steering_wind_b

subroutine steering_wind_r(u,v,p,lon,lat,proj,nx,ny,nz,nt,kt1,kt2,mi,mj,&
     &u_vor,v_vor,r_steering)

  use constants, only: pi, ra

  implicit none 
  integer ,intent (in)::nx,ny,nz,nt,kt1,kt2
  integer ,intent (in)::proj
  real(4),intent (in)::r_steering
  real(4) ,intent (in)::u(0:nx,0:ny,nz,nt),v(0:nx,0:ny,nz,nt)
  real(4) ,intent (in)::p(nz),lon(0:nx),lat(0:ny)
  integer (4),intent (in)::mi,mj
  real (4),intent (out)::u_vor,v_vor
  real(4) ::u_t0t1(0:nx,0:ny,nz),v_t0t1(0:nx,0:ny,nz)
  real(4) ::u_int(0:nx,0:ny),v_int(0:nx,0:ny)

  real(4) ::s_tot
  integer ::ii,jj

  real(4)::d,theta_d

  real(4)::lonin,latin
  integer::n_steering_x,n_steering_y

  lonin=lon(2)-lon(1)
  latin=lat(2)-lat(1)

  n_steering_x=nint(r_steering/(ra*lonin*pi/180.0*cos(lat(mj)*pi/180.0)*1.0e-3)+5)
  n_steering_y=nint(r_steering/(ra*latin*pi/180.0*1.0e-3)+5)

 
  u_vor=0.
  v_vor=0.
  s_tot=0.0


  u_t0t1(0:nx,0:ny,1:nz)=0.5*(u(0:nx,0:ny,1:nz,kt1)+u(0:nx,0:ny,1:nz,kt2))
  v_t0t1(0:nx,0:ny,1:nz)=0.5*(v(0:nx,0:ny,1:nz,kt1)+v(0:nx,0:ny,1:nz,kt2))


  do jj=-n_steering_y,n_steering_y
     do ii=-n_steering_x,n_steering_x
       if(mi+ii>=0.and.mi+ii<=nx.and.mj+jj>=0.and.mj+jj<=ny)then
       if(proj==1)then
         if(ii/=0.and.jj/=0)then
           theta_d=acos(cos(pi/180*lat(mj))*cos(pi/180*lat(mj+jj))&
                &*cos(pi/180*(lon(mi)-lon(mi+ii)))&
                &+sin(pi/180*lat(mj))*sin(pi/180*lat(mj+jj)))
           d=ra*theta_d
         else
           d=0.0
         endif

       elseif(proj==2)then
         d=sqrt((ii*lonin)**2+(jj*latin)**2)
       end if
       if(d<=r_steering*1.0e3)then
         call integral_p(u_t0t1(mi+ii,mj+jj,1:nz),u_int(mi+ii,mj+jj),nz,p)
         call integral_p(v_t0t1(mi+ii,mj+jj,1:nz),v_int(mi+ii,mj+jj),nz,p)
         if(proj==1)then
           u_vor=u_vor+u_int(mi+ii,mj+jj)*cos(pi/180*lat(mj+jj))
           v_vor=v_vor+v_int(mi+ii,mj+jj)*cos(pi/180*lat(mj+jj))
           s_tot=s_tot+cos(pi/180*lat(mj+jj))
         elseif(proj==2)then
           u_vor=u_vor+u_int(mi+ii,mj+jj)
           v_vor=v_vor+v_int(mi+ii,mj+jj)
           s_tot=s_tot+1
         end if
       end if
     end if
     end do
  end do

  u_vor=u_vor/s_tot
  v_vor=v_vor/s_tot

!  write (*,*)u_vor,v_vor

  return
end subroutine steering_wind_r


subroutine integral_p(var,int,nz,p)
  implicit none 
  integer ,intent (in)::nz
  real (4),intent (in)::var(nz)
  real (4),intent (in)::p(nz)
  real (4),intent (out)::int
  integer ::k
  
  int=0.

  do k=1,nz-1
     int=int+0.5*(var(k)+var(k+1))*(p(k)-p(k+1))
  end do
  int=int/(p(1)-p(nz))

  return
end subroutine integral_p


