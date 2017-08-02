subroutine min_z(z, nx, ny, minlat, minlon, &
  &              z_min, n_min, lat, lon, type_min, del_z_min)
 
  use params, only: fillval, nmax, pmax4, mx, my

  implicit none 
  integer(4),intent (in)::nx,ny
  real(4),intent (in)::lon(0:nx),lat(0:ny)
  real(4),intent (in)::z(0:nx,0:ny)
  real(4),intent (out)::minlat(nmax),minlon(nmax),z_min(nmax)
  real(4),intent (in)::del_z_min
  integer ,intent (out)::n_min,type_min(nmax)
  integer (4)::i,j,ii,jj,m
  integer (4)::i_min


  real(4)::z_tmp(0:nx,0:ny),z_part_r(0:nx,0:ny)
  integer (4)::z_part(0:nx,0:ny),z_part_tmp(0:nx,0:ny)

  integer (4)::flag
  real(4)::r,theta
  integer (4)::mi,mj,mij(2),mi_tmp,mj_tmp

  real(4)::min,min0,max

  real(4)::zmin_tmp
  
  integer (4)::n_part,n_part_tmp
  integer (4)::s_part


  integer (4)::var_part_tmp(1:2,1:pmax4)
  integer (4)::p=0
  real(4)::surround8(8)
  integer (4)::buf_mij(2,nmax)
  logical(4)::mij_flag
 
  

  minlat=0
  minlon=0
  z_min=0
  type_min=0

  n_min=1

  n_part=0
 
  z_part=0
  z_part_tmp=0

  buf_mij=-1


  do j=0,ny
    do i=0,nx
        z_tmp(i,j)=z(i,j)
    end do
  end do

  max=maxval(z)
  min0=minval(z)

  write (*,*)'psea max min0',max,min0

  zmin_tmp=max

  do
    zmin_tmp=zmin_tmp-0.1

    do
      s_part=0

      p=0
      min=minval(z_tmp)
      mij=minloc(z_tmp)
      mi=mij(1)-1
      mj=mij(2)-1
      
      if(min>=zmin_tmp)exit

      n_part=n_part+1
      z_part_tmp(mi,mj)=n_part
      z_tmp(mi,mj)=fillval

      s_part=s_part+1


      do m=1,8
        if(mi+mx(m)>=0.and.mi+mx(m)<=nx.and.mj+my(m)>=0.and.mj+my(m)<=ny)then
          if(z_tmp(mi+mx(m),mj+my(m))<zmin_tmp)then  
            z_part_tmp(mi+mx(m),mj+my(m))=n_part
            z_tmp(mi+mx(m),mj+my(m))=fillval
            
            s_part=s_part+1

            
            if (p < pmax4) then
              p = p + 1 
              var_part_tmp(1,p) = mi+mx(m)
              var_part_tmp(2,p) = mj+my(m)
            end if
          end if
        end if
      end do
        
        
      do 
        if(p==0)exit
        mi_tmp=var_part_tmp(1,p)
        mj_tmp=var_part_tmp(2,p)
        p=p-1
        
        do m=1,8
          if(mi_tmp+mx(m)>=0.and.mi_tmp+mx(m)<=nx.and.mj_tmp+my(m)>=0.and.mj_tmp+my(m)<=ny)then
            
            if(z_tmp(mi_tmp+mx(m),mj_tmp+my(m))<zmin_tmp)then
              if(z_part_tmp(mi_tmp+mx(m),mj_tmp+my(m))==0)then
                z_part_tmp(mi_tmp+mx(m),mj_tmp+my(m))=n_part
                z_tmp(mi_tmp+mx(m),mj_tmp+my(m))=fillval

            s_part=s_part+1

                               
                if (p < pmax4) then
                  p = p + 1 
                  var_part_tmp(1,p) = mi_tmp+mx(m)
                  var_part_tmp(2,p) = mj_tmp+my(m)
                end if
              end if
            end if
          end if
        end do
      end do

      
      if(s_part>200)then
      mij_flag=.true.

      do i_min=1,n_min
      
        if(buf_mij(1,i_min)==mi.and.buf_mij(2,i_min)==mj)mij_flag=.false.

      end do
      
      if(mij_flag)then
        if(zmin_tmp-min>del_z_min)then
          

          do j=0,ny
            do i=0,nx
              if(z_part_tmp(i,j)==n_part.and.n_part/=1) z_part(i,j)=n_min
            end do
          end do



          buf_mij(1,n_min)=mi
          buf_mij(2,n_min)=mj         

          minlon(n_min)=lon(mi)
          minlat(n_min)=lat(mj)
          z_min(n_min)=min
          type_min(n_min)=s_part


          !      write (*,*)(var_part(ii,1,n_min),ii=-l,l),(var_part(ii,2,n_min),ii=-l,l)
          !     write (*,*)n_min,max,mlon(n_min),mlat(n_min)
!          write (99,*)minlon(n_min),minlat(n_min),z_min(n_min)
       !   write (*,*)minlon(n_min),minlat(n_min),z_min(n_min),zmin_tmp,min,del_z_min

!          write (*,*)n_min,minlon(n_min),minlat(n_min)
          n_min=n_min+1         
          

        end if
      end if
    end if
  !    write (*,*)'test'
    end do



   z_tmp=fillval
   do j=0,ny
     do i=0,nx
!       if(z_part_tmp(i,j)>0.and.zs(i,j)<500.)then
!        if(z_part_tmp(i,j)>0.and..not.land_flag(i,j))then
        z_tmp(i,j)=z(i,j)
         
!       end if
     end do
   end do
!   var_part_max=0
   z_part_tmp=0


   if(zmin_tmp<min0)exit

  end do

  n_min=n_min-1

!-------------------------------------------------------------!
  if(n_min>1)then

  do j=0,ny
    do i=0,nx
!      if(zs(i,j)<500)then
!      if(.not.land_flag(i,j))then
        z_tmp(i,j)=z(i,j)
!      else
 !       z_tmp(i,j)=fillval
!      end if
    end do
  end do





  zmin_tmp=z_min(1)


  z_part_tmp=0
  do

    zmin_tmp=zmin_tmp+1.

    if(zmin_tmp>max)exit

!    write (*,*)zmin_tmp


!   write (*,*)zmin_tmp

    s_part=0

    p=0
    min=z_min(1)
    mi=buf_mij(1,1)
    mj=buf_mij(2,1)
    
!    write (*,*)minlon(1),minlat(1),z_min(1)
!    write (*,*)lon(mi),lat(mj),z_tmp(mi,mj)
    !      if(min>=zmin_tmp)exit
    
    n_part=1
    z_part_tmp(mi,mj)=n_part
    z_tmp(mi,mj)=fillval
    
    s_part=s_part+1
    

    
    do m=1,8
      if(mi+mx(m)>=0.and.mi+mx(m)<=nx.and.mj+my(m)>=0.and.mj+my(m)<=ny)then
        if(z_tmp(mi+mx(m),mj+my(m))<zmin_tmp)then  
          z_part_tmp(mi+mx(m),mj+my(m))=n_part
          z_tmp(mi+mx(m),mj+my(m))=fillval
          
          s_part=s_part+1
          
          
          if (p < pmax4) then
            p = p + 1 
            var_part_tmp(1,p) = mi+mx(m)
            var_part_tmp(2,p) = mj+my(m)
          end if
        end if
      end if
    end do
    
        
    do 
      if(p==0)exit
      mi_tmp=var_part_tmp(1,p)
      mj_tmp=var_part_tmp(2,p)
      p=p-1
      
      do m=1,8
        if(mi_tmp+mx(m)>=0.and.mi_tmp+mx(m)<=nx.and.mj_tmp+my(m)>=0.and.mj_tmp+my(m)<=ny)then
          
          if(z_tmp(mi_tmp+mx(m),mj_tmp+my(m))<zmin_tmp)then
            if(z_part_tmp(mi_tmp+mx(m),mj_tmp+my(m))==0)then
              z_part_tmp(mi_tmp+mx(m),mj_tmp+my(m))=n_part
              z_tmp(mi_tmp+mx(m),mj_tmp+my(m))=fillval
              
              s_part=s_part+1
              
              
              if (p < pmax4) then
                p = p + 1 
                var_part_tmp(1,p) = mi_tmp+mx(m)
                var_part_tmp(2,p) = mj_tmp+my(m)
              end if
            end if
          end if
        end if
      end do
    end do

      
    mij_flag=.true.
      
    do j=0,ny
      do i=0,nx
        if(z_part_tmp(i,j)==1)then

!          write (*,*)lon(i),lat(j)


          do i_min=2,n_min
            if(buf_mij(1,i_min)==i.and.buf_mij(2,i_min)==j)mij_flag=.false.
          end do

        end if
      end do
    end do
    
    if(.not.mij_flag)exit



    

    type_min(1)=s_part


    do j=0,ny
      do i=0,nx
        if(z_part_tmp(i,j)==1.and.z_part(i,j)<=1)z_part(i,j)=1
      end do
    end do
    


!   z_tmp=fillval
   do j=0,ny
     do i=0,nx
!       if(z_part_tmp(i,j)>0.and.zs(i,j)<500.)then
!        if(z_part_tmp(i,j)>0.and..not.land_flag(i,j))then
        z_tmp(i,j)=z(i,j)
         
!       end if
     end do
   end do
!   var_part_max=0
   z_part_tmp=0

 end do


else
  z_part=1


end if




  z_part_r=fillval

  do j=0,ny
    do i=0,nx
      if(z_part(i,j)>0) z_part_r(i,j)=real(z_part(i,j))
    end do
  end do

!   if(type_min(1)<200)then
!     minlon(1)=0.
!     minlat(1)=0.
!     z_min(1)=0.
!     do j=0,ny
!       do i=0,nx
!         if(z_part(i,j)==1) z_part_r(i,j)=fillval
!       end do
!     end do
  

!   end if




!  write (*,*)'output z_part'

!  write (97)z_part_r
!  write (97)z_tmp
  




  return 
end subroutine min_z

