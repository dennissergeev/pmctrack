subroutine vor_partition(vor_in, nx, ny, &
     & mlat_out, mlon_out, vor_max_out, mtype_out, &
     & n_max_out, lat, lon, vor_part_max, s_part_out)

use types, only: wp
use constants, only: pi, ra, rkilo, fillval, nmax, pmax, mx, my
use params, only : proj, zeta_max0, zeta_min0, int_zeta_min0, gamma

implicit none

integer    , intent (in)  :: nx
integer    , intent (in)  :: ny
real   (wp), intent (in)  :: lon         (0:nx           )
real   (wp), intent (in)  :: lat         (      0:ny     )
real   (wp), intent (in)  :: vor_in      (0:nx, 0:ny     )
real   (wp), intent (out) :: mlat_out    (           nmax)
real   (wp), intent (out) :: mlon_out    (           nmax)
real   (wp), intent (out) :: vor_max_out (           nmax)
real   (wp), intent (out) :: s_part_out  (           nmax)
integer    , intent (out) :: n_max_out
integer    , intent (out) :: vor_part_max(0:nx, 0:ny     )
integer    , intent (out) :: mtype_out   (           nmax)
! Local variables
logical                   :: mij_flag
logical                   :: remove_flag(nmax)
integer                   :: remove_num
integer                   :: n_max
integer                   :: i, j, ii, jj, m
integer                   :: mmax
integer                   :: i_max, i_max2
integer                   :: mi, mj, mij(2), mi_tmp, mj_tmp
integer                   :: i_vor_min, n_vor_min
integer                   :: n_part
integer                   :: vor_part    (0:nx, 0:ny      )
integer                   :: vor_part_tmp(0:nx, 0:ny      )
integer                   :: buf_mij     (2,          nmax)
integer                   :: buf_mij_out (2,          nmax)
integer                   :: part_num    (            nmax)
integer                   :: mtype_part  (            nmax)
integer                   :: mtype       (            nmax, 2)
integer                   :: surr8_buf   (2,        1:pmax)
integer                   :: p
real   (wp)               :: mlat        (            nmax)
real   (wp)               :: mlon        (            nmax)
real   (wp)               :: vor_max     (            nmax)
real   (wp)               :: lonin, latin
real   (wp)               :: vor         (0:nx, 0:ny      )
real   (wp)               :: max_tmp
real   (wp)               :: vor_min_tmp
real   (wp)               :: del_vor_max_tmp
real   (wp)               :: l(nmax), l_min(nmax), theta_l


lonin = lon(1) - lon(0)
latin = lat(1) - lat(0)

mmax = size(mx)
n_max = 0
n_part = 0
max_tmp = 0
buf_mij = -1
vor_min_tmp = zeta_max0
vor_part = 0
part_num = 0
vor_part_max = 0
mtype(:, :) = 0
mtype_out(:) = 0
mtype_part(:) = 0

do j = 0, ny
  do i = 0, nx
    vor(i, j) = vor_in(i, j)
  enddo
enddo

max_tmp = maxval(vor)
n_vor_min = int((max_tmp - zeta_max0) / int_zeta_min0) + 1

do
  p = 0

  max_tmp = maxval(vor)
  mij = maxloc(vor)
  mi = mij(1) - 1
  mj = mij(2) - 1
  if (max_tmp <= zeta_max0) exit
  n_part = n_part + 1
  vor_part(mi, mj) = n_part

  vor(mi, mj) = 0.

  do m = 1, mmax
    if (       mi + mx(m) >= 0  &
       & .and. mi + mx(m) <= nx &
       & .and. mj + my(m) >= 0  &
       & .and. mj + my(m) <= ny) then
      if (vor(mi+mx(m), mj+my(m)) > zeta_min0) then
        vor_part(mi+mx(m), mj+my(m)) = n_part
        vor(mi+mx(m), mj+my(m)) = 0.

        if (p < pmax) then
          p = p + 1
          surr8_buf(1, p) = mi + mx(m)
          surr8_buf(2, p) = mj + my(m)
        endif
      endif
    endif
  enddo

  do
    if (p == 0) exit
    mi_tmp = surr8_buf(1, p)
    mj_tmp = surr8_buf(2, p)
    p = p - 1

    do m = 1, mmax
      if (       mi_tmp + mx(m) >= 0  &
         & .and. mi_tmp + mx(m) <= nx &
         & .and. mj_tmp + my(m) >= 0  &
         & .and. mj_tmp + my(m) <= ny) then
        if (vor(mi_tmp+mx(m), mj_tmp+my(m)) > zeta_min0) then
          if (vor_part(mi_tmp+mx(m), mj_tmp+my(m)) == 0) then
            vor_part(mi_tmp+mx(m), mj_tmp+my(m)) = n_part
            vor(mi_tmp+mx(m), mj_tmp+my(m)) = 0.

            if (p < pmax) then
              p = p + 1
              surr8_buf(1, p) = mi_tmp + mx(m)
              surr8_buf(2, p) = mj_tmp + my(m)
            endif
          endif
        endif
      endif
    enddo
  enddo


!--------------- Check cold front and synoptic low----------------
  call cf_synop_check(vor_in(0:nx, 0:ny), vor_part(0:nx,0:ny), &
                    & n_part, nx, ny, &
                    & lon(0:nx), lat(0:ny), mtype_part(n_part))

  n_max = n_max + 1
  buf_mij(1, n_max) = mi
  buf_mij(2, n_max) = mj


      mlon(n_max)=lon(mi)
      mlat(n_max)=lat(mj)
      vor_max(n_max)=max_tmp
      part_num(n_max)=vor_part(mi,mj)
      mtype(n_max,1)=mtype_part(part_num(n_max))



    enddo

  vor=0.
  do j=0,ny
    do i=0,nx
      if(vor_part(i,j)>0)then
        vor(i,j)=vor_in(i,j)
      endif
    enddo
  enddo

  vor_part_tmp=0


    !---- divide several peaks in a single part -----

  !-------------------------------------------------------------

  do i_vor_min=0,n_vor_min
    vor_min_tmp=zeta_min0+i_vor_min*int_zeta_min0
    do
      p=0

      max_tmp=maxval(vor)
      mij=maxloc(vor)
      mi=mij(1)-1
      mj=mij(2)-1
      if(max_tmp<=zeta_max0)exit



      n_part=n_part+1
      vor_part_tmp(mi,mj)=n_part
      vor(mi,mj)=0.


      do m=1,8
        if ( mi+mx(m) >= 0  .and. &
           & mi+mx(m) <= nx .and. &
           & mj+my(m) >= 0  .and. &
           & mj+my(m) <= ny       ) then
          if(vor(mi+mx(m),mj+my(m))>vor_min_tmp)then
            vor_part_tmp(mi+mx(m),mj+my(m))=n_part
            vor(mi+mx(m),mj+my(m))=0.

            if (p < pmax) then
              p = p + 1
              surr8_buf(1,p) = mi+mx(m)
              surr8_buf(2,p) = mj+my(m)
            endif
          endif
        endif
      enddo

      do
        if(p==0)exit
        mi_tmp=surr8_buf(1,p)
        mj_tmp=surr8_buf(2,p)
        p=p-1

        do m=1,8
          if ( mi_tmp+mx(m) >= 0  .and. &
             & mi_tmp+mx(m) <= nx .and. &
             & mj_tmp+my(m) >= 0  .and. &
             & mj_tmp+my(m) <= ny       ) then
            if(vor(mi_tmp+mx(m),mj_tmp+my(m))>vor_min_tmp)then
              if(vor_part_tmp(mi_tmp+mx(m),mj_tmp+my(m))==0)then
                vor_part_tmp(mi_tmp+mx(m),mj_tmp+my(m))=n_part
                vor(mi_tmp+mx(m),mj_tmp+my(m))=0.

                if (p < pmax) then
                  p = p + 1
                  surr8_buf(1,p) = mi_tmp+mx(m)
                  surr8_buf(2,p) = mj_tmp+my(m)
                endif
              endif
            endif
          endif
        enddo
      enddo



      !----- check new vortex ----
      mij_flag=.true.
      do i_max=1,n_max
        if(buf_mij(1,i_max)==mi.and.buf_mij(2,i_max)==mj)mij_flag=.false.
      enddo

      if(mij_flag)then
        del_vor_max_tmp=gamma*max_tmp

        if(max_tmp-vor_min_tmp>del_vor_max_tmp)then

          n_max=n_max+1
          buf_mij(1,n_max)=mi
          buf_mij(2,n_max)=mj
          mlon(n_max)=lon(mi)
          mlat(n_max)=lat(mj)
          vor_max(n_max)=max_tmp
          part_num(n_max)=vor_part(mi,mj)
          mtype(n_max,1)=mtype_part(part_num(n_max))

        endif
      endif

    enddo

    vor=0.
    do j=0,ny
      do i=0,nx
        if(vor_part_tmp(i,j)>0)then
          vor(i,j)=vor_in(i,j)

        endif
      enddo
    enddo
    vor_part_tmp=0

  enddo




  do j=0,ny
    do i=0,nx
      vor(i,j)=vor_in(i,j)
    enddo
  enddo


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
!             endif
!           enddo
!         enddo
!       elseif(mlon(i_max)>135.and.mlon(i_max)<139.and.&
!            &mlat(i_max)<47.and.mlat(i_max)>44)then
!         do j=0,ny
!           do i=0,nx
!             if(abs(lon(i)-mlon(i_max))<0.5.and.abs(lat(j)-mlat(i_max))<0.4)then
!               if(land_flag(i,j))remove_flag(i_max)=.true.
!             endif
!           enddo
!         enddo
!       endif

!     endif
!   enddo


  remove_num=0
  n_max_out=0
  do i_max=1,n_max
    if(remove_flag(i_max))then
      remove_num=remove_num+1
      do j=0,ny
        do i=0,nx
          if(vor_part(i,j)==i_max)vor_part(i,j)=-10
        enddo
      enddo
    else
      n_max_out=n_max_out+1
      mlon_out(n_max_out)=mlon(i_max)
      mlat_out(n_max_out)=mlat(i_max)
      mtype_out(n_max_out)=mtype(i_max,1)
!      mtype_out(n_max_out)=mtype(i_max,1)+mtype(i_max,2)
      vor_max_out(n_max_out)=vor_max(i_max)
      buf_mij_out(1:2,n_max_out)=buf_mij(1:2,i_max)
    endif
  enddo


!----- Partation ----

  do i_max=1,n_max_out

    vor_min_tmp=vor_max_out(i_max)


    do
      vor_part_tmp=0

      vor_min_tmp=vor_min_tmp-int_zeta_min0

      if(vor_min_tmp<zeta_min0)exit



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
              surr8_buf(1,p) = mi+mx(m)
              surr8_buf(2,p) = mj+my(m)
            endif
          endif
        endif
      enddo

      do
        if(p==0)exit
        mi_tmp=surr8_buf(1,p)
        mj_tmp=surr8_buf(2,p)
        p=p-1

        do m=1,8
          if(mi_tmp+mx(m)>=0.and.mi_tmp+mx(m)<=nx.and.mj_tmp+my(m)>=0.and.mj_tmp+my(m)<=ny)then

            if(vor(mi_tmp+mx(m),mj_tmp+my(m))>vor_min_tmp)then
              if(vor_part_tmp(mi_tmp+mx(m),mj_tmp+my(m))==0)then
                vor_part_tmp(mi_tmp+mx(m),mj_tmp+my(m))=n_part
                vor(mi_tmp+mx(m),mj_tmp+my(m))=0.



                if (p < pmax) then
                  p = p + 1
                  surr8_buf(1,p) = mi_tmp+mx(m)
                  surr8_buf(2,p) = mj_tmp+my(m)
                endif
              endif
            endif
          endif
        enddo
      enddo


      mij_flag=.true.

      do j=0,ny
        do i=0,nx
          if(vor_part_tmp(i,j)==i_max)then

            do i_max2=1,n_max_out
              if(buf_mij_out(1,i_max2)==i.and.buf_mij_out(2,i_max2)==j&
                   &.and.i_max2/=i_max)   mij_flag=.false.
            enddo

          endif
        enddo
      enddo

      if(.not.mij_flag)exit

      do j=0,ny
        do i=0,nx
          if(vor_part_tmp(i,j)==i_max.and.vor_part_max(i,j)==0)vor_part_max(i,j)=i_max
        enddo
      enddo

      do j=0,ny
        do i=0,nx
            vor(i,j)=vor_in(i,j)
        enddo
      enddo

    enddo

  enddo


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
              endif
              if(l(vor_part_max(ii,jj))<l_min(vor_part_max(ii,jj)))then
                l_min(vor_part_max(ii,jj))=l(vor_part_max(ii,jj))
              endif
            endif
          enddo
        enddo



        if(minval(l_min)<0.9e9)then
          vor_part_tmp(i,j)=minloc(l_min,1)
        endif


      endif

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
            endif

            if(l(vor_part_max(ii,jj))<l_min(vor_part_max(ii,jj)))then
              l_min(vor_part_max(ii,jj))=l(vor_part_max(ii,jj))
            endif

          enddo
        enddo



        if(minval(l_min)<0.9e9)then
          vor_part_tmp(i,j)=minloc(l_min,1)
        endif


      endif


    enddo
  enddo

  do j=0,ny
    do i=0,nx
      if(vor_part_tmp(i,j)/=0.and.vor_part_max(i,j)==0)then
        vor_part_max(i,j)=vor_part_tmp(i,j)

      endif
    enddo
  enddo

!--------------------------------------------------





!--------------------------------------------

!------ calclate size of vortex -----


  s_part_out=0.
  do j=0,ny
    do i=0,nx
      if(vor_part_max(i,j)>0)then

        if(proj==1)then
          s_part_out(vor_part_max(i,j))=s_part_out(vor_part_max(i,j))&
               &+lonin*(pi/180.)*ra/rkilo*cos(lat(j)*pi/180.)*latin*(pi/180.)*ra/rkilo
        elseif(proj==2)then
          s_part_out(vor_part_max(i,j))=s_part_out(vor_part_max(i,j))&
               &+lonin*latin*1.0e-6
        endif
      endif
    enddo
  enddo

end subroutine vor_partition
