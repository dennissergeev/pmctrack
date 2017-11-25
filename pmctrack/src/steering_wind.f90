subroutine steering_wind_f(u,v,p,nx,ny,nz,nt,kt,mi,mj,u_vor,v_vor)

  use types, only : wp
  use params, only : n_steering_x, n_steering_y
  use utils, only : integral_p

  implicit none

  integer ,intent (in)::nx,ny,nz,nt,kt
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


subroutine steering_wind_b(u,v,p,nx,ny,nz,nt,kt,mi,mj,u_vor,v_vor)

  use types, only : wp
  use params, only : n_steering_x, n_steering_y
  use utils, only : integral_p

  implicit none

  integer ,intent (in)::nx,ny,nz,nt,kt
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

subroutine steering_wind_r(u,v,p,lon,lat,nx,ny,nz,nt,kt1,kt2,mi,mj,&
     &u_vor,v_vor)

  use types, only : wp
  use constants, only: pi, ra
  use params, only : r_steering, proj
  use utils, only : integral_p

  implicit none 
  integer ,intent (in)::nx,ny,nz,nt,kt1,kt2
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
  integer::x_steer,y_steer

  lonin=lon(2)-lon(1)
  latin=lat(2)-lat(1)

  x_steer=nint(r_steering/(ra*lonin*pi/180.0*cos(lat(mj)*pi/180.0)*1.0e-3)+5)
  y_steer=nint(r_steering/(ra*latin*pi/180.0*1.0e-3)+5)

 
  u_vor=0.
  v_vor=0.
  s_tot=0.0
  d = 0.
  theta_d = 0.


  do jj=-y_steer,y_steer
     do ii=-x_steer,x_steer
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
