subroutine vor_partition(vor_in,nx,ny,proj,&
     &mlat_out,mlon_out,vor_max_out,mtype_out,&
     &n_max_out,lat,lon,vor_part_max,s_part_out,&
     &vor_min,vor_part_min,del_vor_min,del_vor_max_coeff,d_cf_min,size_synop)
  use const
  implicit none 
  integer(4),intent (in)::nx,ny
  integer(4),intent (in)::proj
  real(4),intent (in)::lon(0:nx),lat(0:ny)
  real(4),intent (in)::vor_in(0:nx,0:ny)
  real(4),intent (out)::mlat_out(100),mlon_out(100),vor_max_out(100)
  real (4),intent (out)::s_part_out(100)
  integer ,intent (out)::n_max_out
  integer ,intent (out)::vor_part_max(0:nx,0:ny)
  integer ,intent (out)::mtype_out(100)

  real(4),intent (in)::vor_min,vor_part_min,del_vor_min,del_vor_max_coeff
  real(4),intent (in)::d_cf_min,size_synop


  real(4)::vor_out(0:nx,0:ny)

  real(4)::mlat(100),mlon(100),vor_max(100)
  integer (4)::n_max
  real(4)::lonin,latin

  real(4)::vor(0:nx,0:ny)
  integer (4)::i,j,ii,jj,m
  integer (4)::i_max,i_max2
  real(4)::max_tmp
  real (4)::vor_min_tmp
  integer (4)::mi,mj,mij(2),mi_tmp,mj_tmp
  integer (4)::i_vor_min,n_vor_min
  logical(4)::mij_flag,topo_flag
  integer (4)::n_part,i_part
  integer (4)::vor_part(0:nx,0:ny),vor_part_tmp(0:nx,0:ny)
  real(4)::vor_part_r(0:nx,0:ny),vor_part_max_r(0:nx,0:ny)
  integer (4)::buf_mij(2,100),buf_mij_out(2,100)

  integer (4)::part_num(100),part_num_out(100)
  integer (4)::mtype_part(100),mtype(100,2)
  real (4)::partmax(100)

  real(4)::del_vor_max_tmp

  real (4),parameter ::undef=9.99e20

  real (4)::d,theta_d
  real (4)::d_cf,theta_d_cf
  logical (4)::remove_flag(100)
  integer (4)::remove_num

  real (4)::l(100),l_min(100),theta_l
  integer (4)::l_min_loc

  real(4)::cdot,theta
  real(4)::gradx,grady

  real(4)::max0_divide


  !----cold front
  integer (4)::width,width_max,pnum
  integer (4)::j_n,j_s,i_w,i_e,i_n,i_s
  logical ::flag_coldfront
  logical::flag_one


  character (len=80)::fname_one,fname_cf

  !----synoptic
  real (4)::size_vor

  integer (4),parameter ::pmax=10000
  integer (4)::surround8_buf(1:2,1:pmax)
  integer (4)::p=0
  integer (4)::mx(8),my(8)

  
  mx(1:8)=(/1,1,0,-1,-1,-1,0,1/)
  my(1:8)=(/0,1,1,1,0,-1,-1,-1/)




  lonin=lon(1)-lon(0)
  latin=lat(1)-lat(0)


  n_max=0
  n_part=0
  max_tmp=0.

  buf_mij=-1

  vor_min_tmp=vor_min

  vor_out(0:nx,0:ny)=undef

  vor_out(0,ny)=-1

  vor_part=0
  part_num=0
  vor_part_max=0
  partmax=0.

  mtype(:,:)=0

  mtype_out(:)=0
  mtype_part(:)=0
  do j=0,ny
    do i=0,nx
      vor(i,j)=vor_in(i,j)
    end do
  end do

  max_tmp=maxval(vor)
  n_vor_min=int((max_tmp-vor_min)/del_vor_min)+1

  do
    p=0
   

    max_tmp=maxval(vor)
    mij=maxloc(vor)
    mi=mij(1)-1
    mj=mij(2)-1
    if(max_tmp<=vor_min)exit
    n_part=n_part+1
    vor_part(mi,mj)=n_part
    partmax(n_part)=vor(mi,mj)


    vor(mi,mj)=0.


    do m=1,8
      if(vor(mi+mx(m),mj+my(m))>vor_part_min)then     
        vor_part(mi+mx(m),mj+my(m))=n_part
        vor(mi+mx(m),mj+my(m))=0.

        if (p < pmax) then
          p = p + 1 
          surround8_buf(1,p) = mi+mx(m)
          surround8_buf(2,p) = mj+my(m)
        end if
      end if
    end do

    do 
      if(p==0)exit
      mi_tmp=surround8_buf(1,p)
      mj_tmp=surround8_buf(2,p)
      p=p-1

      do m=1,8
        if(vor(mi_tmp+mx(m),mj_tmp+my(m))>vor_part_min)then
          if(vor_part(mi_tmp+mx(m),mj_tmp+my(m))==0)then
            vor_part(mi_tmp+mx(m),mj_tmp+my(m))=n_part

            vor(mi_tmp+mx(m),mj_tmp+my(m))=0.
            if (p < pmax) then
              p = p + 1 
              surround8_buf(1,p) = mi_tmp+mx(m)
              surround8_buf(2,p) = mj_tmp+my(m)
            end if
          end if
        end if
      end do
    end do




!--------------- Check cold front and synoptic low----------------
    call cf_synop_check(vor_in(0:nx,0:ny),vor_part(0:nx,0:ny),n_part,nx,ny,proj,&
         &lon(0:nx),lat(0:ny),mtype_part(n_part),d_cf_min,size_synop)

      n_max=n_max+1         
      buf_mij(1,n_max)=mi
      buf_mij(2,n_max)=mj         


      mlon(n_max)=lon(mi)
      mlat(n_max)=lat(mj)
      vor_max(n_max)=max_tmp
      part_num(n_max)=vor_part(mi,mj)
      mtype(n_max,1)=mtype_part(part_num(n_max))



    end do

  vor=0.
  do j=0,ny
    do i=0,nx
      if(vor_part(i,j)>0)then
        vor(i,j)=vor_in(i,j)
        !           if(vor_out(i,j)==0)vor_out(i,j)=vor_part_max(i,j)
      end if
    end do
  end do

  vor_part_tmp=0


    !---- divide several peaks in a single part -----

  !-------------------------------------------------------------

  do i_vor_min=0,n_vor_min
    vor_min_tmp=vor_part_min+i_vor_min*del_vor_min
    do
      p=0

      max_tmp=maxval(vor)
      mij=maxloc(vor)
      mi=mij(1)-1
      mj=mij(2)-1
      if(max_tmp<=vor_min)exit



      n_part=n_part+1
      vor_part_tmp(mi,mj)=n_part
      vor(mi,mj)=0.


      do m=1,8
        if(vor(mi+mx(m),mj+my(m))>vor_min_tmp)then     
          vor_part_tmp(mi+mx(m),mj+my(m))=n_part
          vor(mi+mx(m),mj+my(m))=0.

          if (p < pmax) then
            p = p + 1 
            surround8_buf(1,p) = mi+mx(m)
            surround8_buf(2,p) = mj+my(m)
          end if
        end if
      end do

      do 
        if(p==0)exit
        mi_tmp=surround8_buf(1,p)
        mj_tmp=surround8_buf(2,p)
        p=p-1

        do m=1,8
          if(vor(mi_tmp+mx(m),mj_tmp+my(m))>vor_min_tmp)then
            if(vor_part_tmp(mi_tmp+mx(m),mj_tmp+my(m))==0)then
              vor_part_tmp(mi_tmp+mx(m),mj_tmp+my(m))=n_part
              vor(mi_tmp+mx(m),mj_tmp+my(m))=0.

              if (p < pmax) then
                p = p + 1 
                surround8_buf(1,p) = mi_tmp+mx(m)
                surround8_buf(2,p) = mj_tmp+my(m)
              end if
            end if
          end if
        end do
      end do



      !----- check new vortex ----
      mij_flag=.true.
      do i_max=1,n_max
        if(buf_mij(1,i_max)==mi.and.buf_mij(2,i_max)==mj)mij_flag=.false.
      end do

      if(mij_flag)then
        del_vor_max_tmp=del_vor_max_coeff*max_tmp

        if(max_tmp-vor_min_tmp>del_vor_max_tmp)then

          n_max=n_max+1         
          buf_mij(1,n_max)=mi
          buf_mij(2,n_max)=mj         
          mlon(n_max)=lon(mi)
          mlat(n_max)=lat(mj)
          vor_max(n_max)=max_tmp
          part_num(n_max)=vor_part(mi,mj)
          mtype(n_max,1)=mtype_part(part_num(n_max))

        end if
      end if

    end do

    vor=0.
    do j=0,ny
      do i=0,nx
        if(vor_part_tmp(i,j)>0)then
          vor(i,j)=vor_in(i,j)

        end if
      end do
    end do
    vor_part_tmp=0

  end do




  do j=0,ny
    do i=0,nx
      vor(i,j)=vor_in(i,j)
    end do
  end do


!--------------------------------------------------------
!----- check the topographical effect ------
  remove_flag=.false.

!   do i_max=1,n_max
!     if(vor_max(i_max)<0.3e-3)then
     
!       if(mlon(i_max)<135.and.mlon(i_max)>125.and.&
!            &mlat(i_max)<44.and.mlat(i_max)>38)then
!         do j=0,ny
!           do i=0,nx
!             if(abs(lon(i)-mlon(i_max))<0.5.and.abs(lat(j)-mlat(i_max))<0.4)then
!               if(land_flag(i,j))remove_flag(i_max)=.true.
!             end if
!           end do
!         end do
!       elseif(mlon(i_max)>135.and.mlon(i_max)<139.and.&
!            &mlat(i_max)<47.and.mlat(i_max)>44)then
!         do j=0,ny
!           do i=0,nx
!             if(abs(lon(i)-mlon(i_max))<0.5.and.abs(lat(j)-mlat(i_max))<0.4)then
!               if(land_flag(i,j))remove_flag(i_max)=.true.
!             end if
!           end do
!         end do
!       end if
  
!     end if
!   end do
  

!   do i_part=1,maxval(part_num,1)
!     topo_flag=.true.
!     do i_max=1,n_max
!       if(part_num(i_max)==i_part.and..not.remove_flag(i_max))then
!         topo_flag=.false.
!       end if
!     end do
!     if(topo_flag)then
!       do j=0,ny
!         do i=0,nx
!           if(vor_part(i,j)==i_part)vor_part(i,j)=0
!         end do
!       end do
!     end if
!   end do





  remove_num=0
  n_max_out=0
  do i_max=1,n_max
    if(remove_flag(i_max))then
      remove_num=remove_num+1
      do j=0,ny
        do i=0,nx
          if(vor_part(i,j)==i_max)vor_part(i,j)=-10
        end do
      end do
    else
      n_max_out=n_max_out+1
      mlon_out(n_max_out)=mlon(i_max)
      mlat_out(n_max_out)=mlat(i_max)
      mtype_out(n_max_out)=mtype(i_max,1)
!      mtype_out(n_max_out)=mtype(i_max,1)+mtype(i_max,2)
      vor_max_out(n_max_out)=vor_max(i_max)
      buf_mij_out(1:2,n_max_out)=buf_mij(1:2,i_max)
      part_num_out(n_max_out)=part_num(i_max)
    end if
  end do

    
  
  

      
!---------------------------------------------------------

!----- Partation ----

  do i_max=1,n_max_out

    vor_min_tmp=vor_max_out(i_max)


    do
      vor_part_tmp=0

      vor_min_tmp=vor_min_tmp-del_vor_min

      if(vor_min_tmp<vor_part_min)exit

      

      p=0
      max_tmp=vor_max_out(i_max)
      mi=buf_mij_out(1,i_max)
      mj=buf_mij_out(2,i_max)


      n_part=i_max
      vor(mi,mj)=0.


      do m=1,8
        if(mi+mx(m)>=0.and.mi+mx(m)<=nx.and.mj+my(m)>=0.and.mj+my(m)<=ny)then
          if(vor(mi+mx(m),mj+my(m))>vor_min_tmp)then  
            vor_part_tmp(mi+mx(m),mj+my(m))=n_part
            vor(mi+mx(m),mj+my(m))=0



            if (p < pmax) then
              p = p + 1 
              surround8_buf(1,p) = mi+mx(m)
              surround8_buf(2,p) = mj+my(m)
            end if
          end if
        end if
      end do

      do 
        if(p==0)exit
        mi_tmp=surround8_buf(1,p)
        mj_tmp=surround8_buf(2,p)
        p=p-1

        do m=1,8
          if(mi_tmp+mx(m)>=0.and.mi_tmp+mx(m)<=nx.and.mj_tmp+my(m)>=0.and.mj_tmp+my(m)<=ny)then

            if(vor(mi_tmp+mx(m),mj_tmp+my(m))>vor_min_tmp)then
              if(vor_part_tmp(mi_tmp+mx(m),mj_tmp+my(m))==0)then
                vor_part_tmp(mi_tmp+mx(m),mj_tmp+my(m))=n_part
                vor(mi_tmp+mx(m),mj_tmp+my(m))=0.



                if (p < pmax) then
                  p = p + 1 
                  surround8_buf(1,p) = mi_tmp+mx(m)
                  surround8_buf(2,p) = mj_tmp+my(m)
                end if
              end if
            end if
          end if
        end do
      end do


      mij_flag=.true.

      do j=0,ny
        do i=0,nx
          if(vor_part_tmp(i,j)==i_max)then

            do i_max2=1,n_max_out
              if(buf_mij_out(1,i_max2)==i.and.buf_mij_out(2,i_max2)==j&
                   &.and.i_max2/=i_max)   mij_flag=.false.
            end do

          end if
        end do
      end do

      if(.not.mij_flag)exit

      do j=0,ny
        do i=0,nx
          if(vor_part_tmp(i,j)==i_max.and.vor_part_max(i,j)==0)vor_part_max(i,j)=i_max
        end do
      end do

      do j=0,ny
        do i=0,nx
            vor(i,j)=vor_in(i,j)
        end do
      end do

    end do

  end do


!--------------------------------------------------

!----- Nearest part ---------
  vor_part_tmp=0


  do j=0,ny
    do i=0,nx
      if(vor_part(i,j)>0.and.vor_part_max(i,j)==0)then

        l_min(:)=1.e9
        do jj=0,ny
          do ii=0,nx


            if(vor_part(ii,jj)==vor_part(i,j).and.vor_part_max(ii,jj)/=0)then
              if(proj==1)then
                if(i/=ii.and.j/=jj)then
                  theta_l=acos(cos(pi/180*lat(j))*cos(pi/180*lat(jj))&
                       &*cos(pi/180*(lon(i)-lon(ii)))&
                       &+sin(pi/180*lat(j))*sin(pi/180*lat(jj)))
                  l(vor_part_max(ii,jj))=ra*theta_l
                else
                  l(vor_part_max(ii,jj))=0.0
                endif
              elseif(proj==2)then
                l(vor_part_max(ii,jj))=sqrt((lon(i)-lon(ii))**2+(lat(j)-lat(jj))**2)
              end if
              if(l(vor_part_max(ii,jj))<l_min(vor_part_max(ii,jj)))then
                l_min(vor_part_max(ii,jj))=l(vor_part_max(ii,jj))
              end if
            end if
          end do
        end do
        


        if(minval(l_min)<0.9e9)then
          vor_part_tmp(i,j)=minloc(l_min,1)
        end if


      end if

      if(vor_part(i,j)==-10.and.vor_part_max(i,j)==0)then

        l_min(:)=1.e9
        do jj=0,ny
          do ii=0,nx
            if(proj==1)then
              if(i/=ii.and.j/=jj)then
                theta_l=acos(cos(pi/180*lat(j))*cos(pi/180*lat(jj))&
                     &*cos(pi/180*(lon(i)-lon(ii)))&
                     &+sin(pi/180*lat(j))*sin(pi/180*lat(jj)))
                l(vor_part_max(ii,jj))=ra*theta_l
              else
                l(vor_part_max(ii,jj))=0.0
              endif
            elseif(proj==2)then
              l(vor_part_max(ii,jj))=sqrt((lon(i)-lon(ii))**2+(lat(j)-lat(jj))**2)
            end if
              
            if(l(vor_part_max(ii,jj))<l_min(vor_part_max(ii,jj)))then
              l_min(vor_part_max(ii,jj))=l(vor_part_max(ii,jj))
            end if
            
          end do
        end do
        


        if(minval(l_min)<0.9e9)then
          vor_part_tmp(i,j)=minloc(l_min,1)
        end if


      end if


    end do
  end do

  do j=0,ny
    do i=0,nx
      if(vor_part_tmp(i,j)/=0.and.vor_part_max(i,j)==0)then
        vor_part_max(i,j)=vor_part_tmp(i,j)

      end if
    end do
  end do

!--------------------------------------------------



  

!--------------------------------------------

!------ calclate size of vortex -----
  

  s_part_out=0.
  do j=0,ny
    do i=0,nx
      if(vor_part_max(i,j)>0)then
        
        if(proj==1)then
          s_part_out(vor_part_max(i,j))=s_part_out(vor_part_max(i,j))&
               &+lonin*(pi/180.)*ra/1000.*cos(lat(j)*pi/180.)*latin*(pi/180.)*ra/1000.
        elseif(proj==2)then
          s_part_out(vor_part_max(i,j))=s_part_out(vor_part_max(i,j))&
               &+lonin*latin*1.0e-6
        end if
      end if
    end do
  end do


  do j=0,ny
    do i=0,nx
      if(vor_part_max(i,j)>0)then
        vor_out(i,j)=real(vor_part_max(i,j))
      end if
    end do
  end do


  vor_part_r=undef
  vor_part_max_r=undef

  vor_part_r(0,ny)=-1.
  vor_part_max_r(0,ny)=-1.
  

  do j=0,ny
    do i=0,nx
      if(vor_part_max(i,j)>0)then
        vor_part_max_r(i,j)=real(vor_part_max(i,j))
      end if
    end do
  end do

  do j=0,ny
    do i=0,nx
      if(vor_part(i,j)/=0)then
        vor_part_r(i,j)=real(vor_part(i,j))
      end if
    end do
  end do

!  write (98)vor_part_r
!  write (98)vor_part_max_r


!  write (*,*)'vor_partition finish!'



  return
end subroutine vor_partition




