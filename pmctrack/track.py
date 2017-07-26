#!/usr/bin/env python
# coding: utf-8
"""
Track PMCs
"""
# import warnings
# warnings.filterwarnings('ignore')
import iris
import numpy as np
from path import Path
from os import getenv

from core import tracking_main

iris.FUTURE.netcdf_promote = True
iris.FUTURE.cell_datetime_objects = True

dirname = Path(getenv('HOME') / 'phd' / 'reanalysis' / 'ERA5'

cl = iris.load(dirname.listdir('*.nc'))

cl

slp = cl.extract_strict('air_pressure_at_sea_level')[..., ::-1, :]
vo = cl.extract_strict('atmosphere_relative_vorticity')[..., ::-1, :]
u = cl.extract_strict('eastward_wind')[..., ::-1, :]
v = cl.extract_strict('northward_wind')[..., ::-1, :]
lsm = cl.extract_strict('land_binary_mask')[..., ::-1, :]

proj = 1
vert_grid = 1

nx = vo.shape[3] - 1
ny = vo.shape[2] - 1
nz = vo.shape[1]
nt = vo.shape[0]

nx1 = 10
nx2 = nx - 10
ny1 = 5
ny2 = 90

lon = vo.coord('longitude')
lat = vo.coord('latitude')
plev = vo.coord('pressure_level')
tcoord = vo.coord('time')

lons = lon.points[0]
lats = lat.points[0]
lonin = (np.diff(lon.points)).mean()
latin = (np.diff(lat.points)).mean()

levs = plev.points[::-1]

del_t = np.diff(tcoord.units.num2date(tcoord.points))[0].total_seconds()

# parameter for smoothing of vorticity
smth_type = 2
nsmth_x = 10
nsmth_y = 10
r_smth = 30.0

# parameter for detecting vortex
zeta_max0 = 2.0e-4
zeta_min0 = 1.5e-4
int_zeta_min0 = 0.02e-4
gamma = 0.25

# parameter for excluding the synoptic scale disturbances
d_cf_min = 400.0
size_synop = 40000.0
distance_ec = 300.0
del_psea_min = 0.5

# parameter for calculating steering winds
steering_type = 2
n_steering_x = 20
n_steering_y = 20
r_steering = 200.0

# parameter for linking vortex
track_type = 2
del_lon = 1.0
del_lat = 0.8
del_r = 120.0

# parameter for checking the track
period_min = 3

vort_level = 950

slp.transpose((2, 1, 0))
vo.transpose((3, 2, 1, 0))
u.transpose((3, 2, 1, 0))
v.transpose((3, 2, 1, 0))
lsm.transpose((2, 1, 0))

vo_data = vo.extract(iris.Constraint(pressure_level=950)).transpose()

slp_data = slp.data * 1e2
vo_data = vo.extract(iris.Constraint(pressure_level=950)).data
u_data = u.data[:, ::-1, ...]
v_data = v.data[:, ::-1, ...]

print(u_data.shape)
print(vo_data.shape)

tracking_main(vo_data, u_data, v_data, slp_data,
              proj, vert_grid,
              nx1, nx2, ny1, ny2,
              levs, lons, lats, lonin, latin, del_t,
              nsmth_x, nsmth_y, r_smth, smth_type,
              zeta_max0, zeta_min0, int_zeta_min0, gamma,
              n_steering_x, n_steering_y, r_steering, steering_type,
              del_lon, del_lat, del_r,
              track_type, period_min,
              d_cf_min, size_synop, del_psea_min, distance_ec,
              # nx, ny, nz, nt
              )
