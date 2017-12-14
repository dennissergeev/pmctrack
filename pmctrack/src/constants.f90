module constants

use types, only : wp

implicit none

! Numerical parameters
real   (kind=wp), parameter :: fillval = 9.99e20_wp
real   (kind=wp), parameter :: missval = -999.9_wp

integer, parameter :: steer_nt = 2
integer, parameter :: nmax = 300
integer, parameter :: kmax = 1000
integer, parameter :: pmax = 10000 
integer, parameter :: pmax4 = pmax * 4

integer, parameter :: mx(8) = (/1, 1, 0, -1, -1, -1,  0,  1/)
integer, parameter :: my(8) = (/0, 1, 1,  1,  0, -1, -1, -1/)

! Physical constants
real(wp), parameter :: pi = 4.0 * atan(1.0_wp)
real(wp), parameter :: ra = 6378.0e3_wp

! Metric parameters
integer , parameter :: ikilo = 1000
real(wp), parameter :: rkilo = 1000.0_wp

end module constants

