subroutine synop_check(mlon,mlat,n_max,minlon,minlat,n_min,mtype,proj,d_min)
  use const
  implicit none 
  integer (4),intent (in)::n_max,n_min
  real(4),intent (in)::mlon(n_max),mlat(n_max)
  real(4),intent (in)::minlon(n_min),minlat(n_min)
  real(4),intent (in) ::d_min
  integer (4),intent (inout)::mtype(n_max)
  integer(4),intent (in)::proj

  integer (4)::i_max,i_min
  real (4)::d,theta_d
  logical (4)::flag_synop(n_max)


  
  flag_synop=.false.

!  write (*,*)'synop_check'

  do i_max=1,n_max
    
    if(mtype(i_max)>=1)then

      do i_min=1,n_min
        if(proj==1)then
          theta_d=acos(cos(pi/180*mlat(i_max))*cos(pi/180*minlat(i_min))&
               &*cos(pi/180*(mlon(i_max)-minlon(i_min)))&
               &+sin(pi/180*mlat(i_max))*sin(pi/180*minlat(i_min)))
          d=ra*theta_d
        elseif(proj==2)then
          d=sqrt((mlon(i_max)-minlon(i_min))**2+(mlat(i_max)-minlat(i_min))**2)
        end if
 !       write (*,*)'d=',d
        
        if(d<d_min*1.0e3)flag_synop(i_max)=.true.
      end do

      if(mtype(i_max)>=2.and..not.flag_synop(i_max))then
        mtype(i_max)=mtype(i_max)-2


!      elseif(mtype(i_max)==1.and.flag_synop(i_max))then
!        mtype(i_max)=mtype(i_max)+2
      end if


    end if
  end do

  return
end subroutine synop_check
