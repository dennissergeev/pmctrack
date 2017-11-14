module params
  implicit none 

  real   (4), parameter :: fillval = 9.99e20
  real   (4), parameter :: missval = -999.9

  integer(4), parameter :: nmax = 300
  integer(4), parameter :: kmax = 1000
  integer(4), parameter :: pmax = 10000 
  integer(4), parameter :: pmax4 = pmax * 4

  integer(4), parameter :: mx(8) = (/1, 1, 0, -1, -1, -1,  0,  1/)
  integer(4), parameter :: my(8) = (/0, 1, 1,  1,  0, -1, -1, -1/)
end module params

