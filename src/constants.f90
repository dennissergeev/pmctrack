module constants

use types, only : wp

implicit none

! Numerical parameters
real   (kind=wp), parameter :: fillval = 9.99e20_wp
integer         , parameter :: ifillval = -999
real   (kind=wp), parameter :: missval = -999.9_wp

integer, parameter :: steer_nt = 2
integer, parameter :: nmax = 2000
integer, parameter :: pmax = 20000 
integer, parameter :: pmax4 = pmax * 4

integer, parameter :: mx(8) = (/1, 1, 0, -1, -1, -1,  0,  1/)
integer, parameter :: my(8) = (/0, 1, 1,  1,  0, -1, -1, -1/)

! Physical constants
real(wp), parameter :: pi = 4.0 * atan(1.0_wp)
real(wp), parameter :: deg2rad = pi / 180.
real(wp), parameter :: rad2deg = 1. / deg2rad
real(wp), parameter :: ra = 6378.0e3_wp

! Metric parameters
integer , parameter :: ikilo = 1000
real(wp), parameter :: rkilo = 1000.0_wp

! IO units
integer, parameter :: fh_conf = 999
integer, parameter :: fh_bin = 998
integer, parameter :: fh_maxloc = 997
integer, parameter :: fh_track = 996

end module constants

