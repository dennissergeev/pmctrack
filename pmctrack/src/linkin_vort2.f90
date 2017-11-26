subroutine linkin_vort2(mlon,mlat,mtype,u_vor_f,v_vor_f,nt,&
  &                     n_max,vor_index,vor_num,vor_merge,vor_part,&
  &                     nx,ny,lon,lat,del_t)

  use types, only: wp
  use constants, only: pi, ra, fillval, nmax, pmax 
  use params, only: proj, del_r 

  implicit none 

  integer    , intent(in)  :: nt,nx,ny
  integer    , intent(in)  :: n_max(nt)
  real   (wp), intent(in)  :: mlon(nmax, nt), mlat(nmax, nt)
  integer    , intent(in)  :: mtype(nmax, nt)
  real   (wp), intent(in)  :: u_vor_f(nmax, nt), v_vor_f(nmax, nt)
  ! real (4),intent (in)::u_vor_b(nmax,nt),v_vor_b(nmax,nt)
  real   (wp), intent(in)  :: lon(0:nx),lat(0:ny)
  real   (wp), intent(in)  :: del_t
  integer    , intent(in)  ::   vor_part(0:nx,0:ny,nt)
  integer    , intent(out) :: vor_index(pmax,nt)
  integer    , intent(out) :: vor_num
  integer    , intent(out) :: vor_merge(pmax)
  ! Local variables
  integer                  :: kt
  integer                  :: i_max,i_max1
  integer                  :: i_next(nmax, nt)
  real   (wp)              :: r_next(0:nmax, nt), r_next_tmp
  integer                  :: vor_part_s(nmax)
  integer                  :: i, j

  integer                  :: i_vor_num, i_vor_num2

  logical                  :: vor_previous_flag(nmax, nt)
  integer                  :: vor_previous_index(nmax, nt)

  integer                  :: vor_c_f_index(nmax, nmax), n_c_f(nmax)

  real   (wp)              :: e_mv_lon, e_mv_lat
  real   (wp)              :: e_mlon, e_mlat

  integer                  :: land_num, all_num

  real   (wp)              :: r_c(nmax, nmax), r_c_min

  real   (wp)              :: r_tmp, theta_tmp


  r_tmp = 0.
  theta_tmp = 0.
  e_mv_lon = 0.
  e_mv_lat = 0.

  r_next = fillval

  vor_merge=0

  vor_index=0
  vor_num=0

  vor_c_f_index=0
  i_next=0

  n_c_f=0

  vor_previous_flag(1:nmax,1:nt)=.false.
  vor_previous_index=0

  do kt = 1, nt-1 ! Time loop


    do i_max=1,n_max(kt)
      do i_vor_num=1,vor_num           
        if(i_max==vor_index(i_vor_num,kt))then 
          !The vortex labeled as i_max at kt  existed at kt-1
          vor_previous_flag(i_max,kt) = .true.
          vor_previous_index(i_max,kt) = vor_index(i_vor_num,kt-1)
        end if
      end do
    end do


    ! t=kt -> t=kt+1

    n_c_f=0
    r_c=del_r*1.0e3
    do i_max=1,n_max(kt)

      if(proj==1)then
        e_mv_lon=(u_vor_f(i_max,kt)*del_t/(ra*cos(mlat(i_max,kt)*pi/180.)))*(180/pi)        
        e_mv_lat=(v_vor_f(i_max,kt)*del_t/ra)*(180/pi)
      elseif(proj==2)then
        e_mv_lon=u_vor_f(i_max,kt)*del_t
        e_mv_lat=v_vor_f(i_max,kt)*del_t
      end if

      e_mlon=mlon(i_max,kt)+e_mv_lon        
      e_mlat=mlat(i_max,kt)+e_mv_lat
        
!        if(kt>=541.and.kt<=542)then
!          write (77,*)kt,mlon(i_max,kt),mlat(i_max,kt),e_mlon,e_mlat
!        end if


        ! -------- Tracking by estimation of movement ------

      r_c_min=del_r*1.0e3    
      do i_max1=1,n_max(kt+1)

          if(proj==1)then
            if(abs(cos(pi/180*e_mlat)*cos(pi/180*mlat(i_max1,kt+1))&
                 &*cos(pi/180*(e_mlon-mlon(i_max1,kt+1)))&
                 &+sin(pi/180*e_mlat)*sin(pi/180*mlat(i_max1,kt+1)))<1.0)then
              theta_tmp=acos(cos(pi/180*e_mlat)*cos(pi/180*mlat(i_max1,kt+1))&
                   &*cos(pi/180*(e_mlon-mlon(i_max1,kt+1)))&
                   &+sin(pi/180*e_mlat)*sin(pi/180*mlat(i_max1,kt+1)))
            else
              theta_tmp=0.0
            endif
            
          
            r_tmp=ra*theta_tmp
          elseif(proj==2)then
            r_tmp=sqrt((mlon(i_max1,kt+1)-e_mlon)**2+(mlat(i_max1,kt+1)-e_mlat)**2)
          end if

        if(r_tmp<=r_c_min)then
!          write (78,*)kt,mlon(i_max,kt),mlat(i_max,kt),mlon(i_max1,kt+1),mlat(i_max1,kt+1),r_tmp*1.0e-3

          i_next(i_max,kt)=i_max1
          r_c_min=r_tmp
          r_next(i_max,kt)=r_tmp
          
          
        end if
      end do





      ! ------ Tracking by part of vortex --------

      land_num=0
      all_num=0
      vor_part_s=0

      if(i_next(i_max,kt)==0)then
        
        do j=0,ny
          do i=0,nx
            !if(abs(e_mlon-lon(i))<=del_lon.and.abs(e_mlat-lat(j))<=del_lat)then
          if(proj==1)then
            if(abs(cos(pi/180*e_mlat)*cos(pi/180*lat(j))&
                   &*cos(pi/180*(e_mlon-lon(i)))&
                   &+sin(pi/180*e_mlat)*sin(pi/180*lat(j)))<1.0)then

              theta_tmp=acos(cos(pi/180*e_mlat)*cos(pi/180*lat(j))&
                   &*cos(pi/180*(e_mlon-lon(i)))&
                   &+sin(pi/180*e_mlat)*sin(pi/180*lat(j)))

            else
              theta_tmp=0.0
            endif
            r_tmp=ra*theta_tmp
 
          elseif(proj==2)then
            r_tmp=sqrt((lon(i)-e_mlon)**2+(lat(j)-e_mlat)**2)
          end if
          
          
          if(r_tmp<=del_r*1.0e3)then
            
            
            if(vor_part(i,j,kt+1)>0)then
              vor_part_s(vor_part(i,j,kt+1))=vor_part_s(vor_part(i,j,kt+1))+1
            end if
          end if
          
          !  end if
          
        end do
      end do



        
        do 
          if(maxval(vor_part_s)==0)exit
          
          i_next(i_max,kt)=maxloc(vor_part_s,1)          
          
          if(proj==1)then
            if(e_mlon/=mlon(i_next(i_max,kt),kt+1))then
              r_next(i_max,kt)=ra*&
                   &acos(cos(pi/180*e_mlat)*cos(pi/180*mlat(i_next(i_max,kt),kt+1))&
                   &*cos(pi/180*(e_mlon-mlon(i_next(i_max,kt),kt+1)))&
                   &+sin(pi/180*e_mlat)*sin(pi/180*mlat(i_next(i_max,kt),kt+1)))
            else
              r_next(i_max,kt)=ra*pi/180*(e_mlat-mlat(i_next(i_max,kt),kt+1))
            end if
          elseif(proj==2)then
            r_next(i_max,kt)=sqrt((e_mlon-mlon(i_next(i_max,kt),kt+1))**2+&
                                  &(e_mlat-mlat(i_next(i_max,kt),kt+1))**2)
          end if

!          write (76,*)kt,r_next(i_max,kt)*1.0e-3,e_mlon,mlon(i_next(i_max,kt),kt+1),e_mlat,mlat(i_next(i_max,kt),kt+1)


          if(r_next(i_max,kt)<=2.0*del_r*1.0e3)exit
          if(mtype(i_next(i_max,kt),kt+1)>=1)exit
          if(r_next(i_max,kt)>2.0*del_r*1.0e3)then
            vor_part_s(i_next(i_max,kt))=0
            i_next(i_max,kt)=0
          end if
           
        end do
         
      end if
    end do
    ! ------------- Connecting vortex(kt)to vortex(kt+1) ------------------



    do i_max=1,n_max(kt)
      if(i_next(i_max,kt)>0)then
        
        write (*,'(A,I4,A,I4,A,I4,A,I4)')"Vortex labeled as ",i_max,' at '&
             &,kt,' is connected to vortex ', i_next(i_max,kt),' at ',kt+1

        if(vor_previous_flag(i_max,kt))then !The vortex existed at kt-1
          write (*,'(A,I4,A,I4,A)')"Vortex labeled as ",i_max,' at ',kt,' existed at kt-1'
       


          do i_vor_num=1,vor_num           
            if(vor_index(i_vor_num,kt)==i_max)then 
              

              vor_index(i_vor_num,kt+1)=i_next(i_max,kt)
            end if
          end do

        end if



        if(.not.vor_previous_flag(i_max,kt))then ! The vortex appear at kt
          vor_num=vor_num+1
          vor_index(vor_num,kt)=i_max
          vor_index(vor_num,kt+1)=i_next(i_max,kt)
        end if



      end if



    end do
    

      !------- check the merger of the vortices -------  
    do i_vor_num=1,vor_num
#ifdef debug
      print*, 'linkin_vort2(259):', i_vor_num, kt, vor_index(i_vor_num, kt)
#endif
      r_next_tmp=r_next(vor_index(i_vor_num,kt),kt)
          
      
      do i_vor_num2=1,vor_num
        if(vor_index(i_vor_num,kt+1)==vor_index(i_vor_num2,kt+1).and.&
             &vor_index(i_vor_num,kt+1)>0.and.&
             &i_vor_num/=i_vor_num2)then
          if(vor_merge(i_vor_num)/=i_vor_num2.and.vor_merge(i_vor_num2)==0)then
            
            if(r_next_tmp>r_next(vor_index(i_vor_num2,kt),kt))then 
              vor_merge(i_vor_num)=i_vor_num2
              r_next_tmp=r_next(vor_index(i_vor_num2,kt),kt)
            end if
            
            if(vor_merge(i_vor_num)>0)then
              write (*,*)'vortex', i_vor_num,&
                   &' merged with vortex',vor_merge(i_vor_num)
            end if
            
          end if
                    
        end if
        
      end do
            
    end do

  end do  ! Time loop 

  return
end subroutine linkin_vort2
