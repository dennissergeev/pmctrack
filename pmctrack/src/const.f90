module constants

use types, only : wp

implicit none

! Physical constants
real(wp), parameter :: pi = 4.0 * atan(1.0_wp)
real(wp), parameter :: ra = 6378.0e3_wp

! Metric parameters
integer , parameter :: ikilo = 1000
real(wp), parameter :: rkilo = 1000.0_wp

end module constants

