module params
  
  use types, only : wp 

  implicit none 

  real   (kind=wp), parameter :: fillval = 9.99e20_wp
  real   (kind=wp), parameter :: missval = -999.9_wp

  integer, parameter :: nmax = 300
  integer, parameter :: kmax = 1000
  integer, parameter :: pmax = 10000 
  integer, parameter :: pmax4 = pmax * 4

  integer, parameter :: mx(8) = (/1, 1, 0, -1, -1, -1,  0,  1/)
  integer, parameter :: my(8) = (/0, 1, 1,  1,  0, -1, -1, -1/)
end module params

