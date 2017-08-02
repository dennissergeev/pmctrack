subroutine smth(var,nx,ny,var_smth,nx1,nx2,ny1,ny2,nsmth_x,nsmth_y)
  implicit none 
  integer ,intent (in)::nx,ny,nx1,nx2,ny1,ny2
  integer ,intent (in)::nsmth_x,nsmth_y
  real (4),intent (in)::var(0:nx,0:ny)
  real (4),intent (out)::var_smth(nx1:nx2,ny1:ny2)
  real (4)::tmp(nx1:nx2,ny1:ny2)
  integer ::nc
  integer ::i,j,ii,jj


  var_smth(nx1:nx2,ny1:ny2)=var(nx1:nx2,ny1:ny2)

!  if(nsmth_x>0.and.nsmth_y>0)then

  do j=ny1,ny2
    if(j-nsmth_y>=0.and.j+nsmth_y<=ny)then
      do i=nx1,nx2
        if(i-nsmth_x>=0.and.i+nsmth_x<=nx)then


          nc=0
          tmp(i,j)=0.
          if(var(i,j)>-10.0)then
            !           write (*,*)'smth at ',i,j
            do jj=-nsmth_y,nsmth_y
              do ii=-nsmth_x,nsmth_x
                if(var(i+ii,j+jj)>-10.0)then
                  tmp(i,j)=tmp(i,j)+var(i+ii,j+jj)
                  nc=nc+1
                end if
              end do
            end do
            !          tmp(i,j)=tmp(i,j)/nc
            tmp(i,j)=tmp(i,j)/((2*nsmth_x+1)*(2*nsmth_y+1))
            var_smth(i,j)=tmp(i,j)
          end if
        end if

      end do
    end if
  end do
!end if

  return 
end subroutine smth




subroutine smth_r(var,nx,ny,lon,lat,var_smth,nx1,nx2,ny1,ny2,r_smth,proj)

  use constants, only: pi, ra

  implicit none

  integer ,intent (in)::nx,ny,nx1,nx2,ny1,ny2
  integer ,intent (in)::proj
  real (4),intent (in)::var(0:nx,0:ny)
  real (4),intent (in)::lon(0:nx),lat(0:ny)
  real (4),intent (in)::r_smth
  real (4),intent (out)::var_smth(nx1:nx2,ny1:ny2)
  real (4)::tmp(nx1:nx2,ny1:ny2)
  integer ::nc
  integer ::i,j,ii,jj
  real(4)::lonin,latin
  integer ::nsmth_x,nsmth_y


  real(4)::d,d_ave,theta_d

  lonin=lon(2)-lon(1)
  latin=lat(2)-lat(1)

  if(proj==1)then
    nsmth_x=r_smth/nint(ra*lonin*pi/180.0*cos(lat(ny1+ny2/2)*pi/180.0)*1.0e-3)+1
    nsmth_y=r_smth/nint(ra*latin*pi/180.0*1.0e-3)+1
  elseif(proj==2)then
    nsmth_x=nint(r_smth/lonin*1.0e-3)+2
    nsmth_y=nint(r_smth/latin*1.0e-3)+2
  end if

!  write (*,*)nsmth_x,nsmth_y

  var_smth(nx1:nx2,ny1:ny2)=var(nx1:nx2,ny1:ny2)

  do j=ny1,ny2
    if(j-nsmth_y>=0.and.j+nsmth_y<=ny)then
      do i=nx1,nx2
        if(i-nsmth_x>=0.and.i+nsmth_x<=nx)then
          
          nc=0
          tmp(i,j)=0.
          if(var(i,j)>-10.0)then
            !           write (*,*)'smth at ',i,j
            do jj=-nsmth_y,nsmth_y
              do ii=-nsmth_x,nsmth_x
                if(proj==1)then
                  if(ii/=0.and.jj/=0)then
                    theta_d=acos(cos(pi/180*lat(j))*cos(pi/180*lat(j+jj))&
                         &*cos(pi/180*(lon(i)-lon(i+ii)))&
                         &+sin(pi/180*lat(j))*sin(pi/180*lat(j+jj)))
                    d=ra*theta_d
                  else
                    d=0.0
                  end if
                  
                elseif(proj==2)then
                  d=sqrt((ii*lonin)**2+(jj*latin)**2)
                end if
              

                  
                if(d<=r_smth*1.0e3)then


                  nc=nc+1
                  if(var(i+ii,j+jj)>-10.0)then
                    tmp(i,j)=tmp(i,j)+var(i+ii,j+jj)
                  end if
                end if
              end do
            end do
  
          var_smth(i,j)=tmp(i,j)/nc
          end if
          
        end if

 !       if(i==(nx1+nx2)/2.and.j==(ny1+ny2)/2)write (*,*)lon(i),lat(j),nc

      end do
    end if
  end do

  return
end subroutine smth_r

