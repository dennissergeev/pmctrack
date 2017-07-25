subroutine track_check2(vortex,vortex_flag,nt,period_min)
  use const
  implicit none 
  integer(4),intent (in)::nt
  integer(4),intent (in)::period_min
  real(4),intent (in)::vortex(1:nt,3)
  integer(4),intent (out)::vortex_flag
  
  integer ::kt
  integer ::vor_period,vor_period_st
  real(4)::vor_max_track,vor_movement

  

  vortex_flag=0

  ! --- check the lifetime of vortex and vor. max---

  vor_period=0
  vor_period_st=0
  vor_max_track=0.
  vor_movement=0.

  do kt=1,nt
    if(vortex(kt,3)>0.00001)then
      vor_period=vor_period+1

!      if(vortex(kt,3)>=vor_max_track0)vor_period_st=vor_period_st+1
!      if(vortex(kt,3)>vor_max_track) vor_max_track=vortex(kt,3)

      
    end if
  end do

  

  if(vor_period>=period_min)vortex_flag=1
    
  return
end subroutine track_check2
