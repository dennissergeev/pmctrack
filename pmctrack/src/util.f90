module util

use params, only: missval

implicit none

contains

subroutine apply_mask_2d(var, nx, ny, flag)

  implicit none

  integer(4), intent (in)    :: nx, ny
  real   (4), intent (inout) :: var      (0:nx, 0:ny)
  real(4), intent (in)    :: flag(0:nx, 0:ny)
  
  integer(4) :: i, j

  do j = 0, ny
    do i = 0, nx
        print*, i, j, flag(i, j)
      if (flag(i, j) == 1.) then
        var(i, j) = missval
      endif
    enddo
  enddo

end subroutine apply_mask_2d

end module util
