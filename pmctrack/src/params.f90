module params
  implicit none 

  real   (4), parameter :: fillval = 9.99e20

  integer(4), parameter :: nmax = 100
  integer(4), parameter :: kmax = 1000
  integer(4), parameter :: pmax = 10000 
  integer(4), parameter :: pmax4 = pmax * 4

  integer(4), parameter :: mx(8) = (/1, 1, 0, -1, -1, -1, 0, 1/)
  integer(4), parameter :: my(8) = (/1, 1, 0, -1, -1, -1, 0, 1/)
end module params

