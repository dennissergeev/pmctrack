module utils

use types, only : wp
use constants, only: missval

implicit none

contains

  subroutine apply_mask_2d(var, nx, ny, flag)
  
    implicit none
  
    integer(4), intent (in)    :: nx, ny
    real   (4), intent (inout) :: var (0:nx, 0:ny)
    real   (4), intent (in)    :: flag(0:nx, 0:ny)
    
    integer(4) :: i, j
  
    do j = 0, ny
      do i = 0, nx
        !  print*, i, j, flag(i, j)
        if (flag(i, j) == 1.) then
          var(i, j) = missval
        endif
      enddo
    enddo
  
  end subroutine apply_mask_2d
  
  
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

end module utils
