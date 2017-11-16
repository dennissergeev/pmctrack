#!/usr/bin/env python
# coding: utf-8
"""
Prepare ERA5 output for tracking algorithm
"""
import argparse
import cf_units
from datetime import datetime
import iris
import numpy as np
from os import getenv
import pandas as pd
from path import Path
import sys
from tqdm import tqdm
import warnings
#
# from arke.coords import get_cube_datetimes
#
from utils import _make_dir

warnings.warn('This script is untested and soon will be deprecated!')

iris.FUTURE.netcdf_promote = True
iris.FUTURE.cell_datetime_objects = True

TOPDIR = Path(getenv('HOME')) / 'phd' / 'reanalysis' / 'ERA5'

NAMES = ('lsm', 'slp', 'u', 'v', 'vort')

CUBE_DICT = dict(
                 lsm=dict(name='land_binary_mask'),
                 slp=dict(name='air_pressure_at_sea_level',
                          scl=0.01),
                 u=dict(name='eastward_wind'),
                 v=dict(name='northward_wind'),
                 vort=dict(name='atmosphere_relative_vorticity'),
)
out_sub_dir = 'track'
out_file_mask = 'era5_slp_u_v_vo_{dt:%Y%m%d%H%M}.dat'
IN_FILES_MASK = '*.nc'

PLEVELS = (1000, 975, 950, 925, 900, 850, 800, 700)


def parse_args(args=None):
    ap = argparse.ArgumentParser(__file__,
                                 description=__doc__,
                                 formatter_class=argparse.
                                 ArgumentDefaultsHelpFormatter,
                                 epilog='Example of use:\ncoming soon...')
    ag_out_loc = ap.add_argument_group(title='Saving parameters')
    ag_out_loc.add_argument('-f', '--force', action='store_true',
                            default=False,
                            help='Overwrite output directory')
    return ap.parse_args(args)


def main(args=None):
    args = parse_args(args)

    indir = TOPDIR
    infiles = indir / IN_FILES_MASK

    cubelist = iris.load(infiles)
    print(cubelist)
    # dtrange = get_cube_datetimes(cubelist[0])
    dtrange = pd.date_range(start=datetime(2011, 1, 29),
                            end=datetime(2011, 1, 31),
                            freq='1H')

    outdir = indir / out_sub_dir  # parent
    _make_dir(outdir, args.force)

    times = []
    for dt1 in tqdm(dtrange):
        times.append(dt1)
        t_constr = iris.Constraint(time=dt1)
        cubes_slice = cubelist.extract(t_constr)

        cubes_out = iris.cube.CubeList()
        for short_name in NAMES:
            vardict = CUBE_DICT[short_name]
            if 'name' in vardict:
                cube = cubes_slice.extract(vardict['name'], strict=True)
            elif 'stash' in vardict:
                attr_constr = iris.AttributeConstraint(STASH=vardict['stash'])
                cube = cubes_slice.extract(attr_constr, strict=True)
            else:
                raise KeyError('Unable to extract cube')
            cube.attributes['short_name'] = short_name
            if 'scl' in vardict:
                cube.data = cube.data * vardict['scl']
            if 'units' in vardict:
                cube.units = cf_units.Unit(vardict['units'])
            if 'regrid_to' in vardict:
                target = cubes_slice.extract(vardict['regrid_to'], strict=True)
                cube = cube.regrid(target, iris.analysis.Linear())
            cubes_out.append(cube)

        all_data = []
        for cube in cubes_out:
            # if cube.attributes['short_name'] != 'slp':
            #     p = cube.coord('pressure')
            #     pvalues = [p[p.nearest_neighbour_index(i)].points[0]
            #                for i in PLEVELS]
            #     pconstr = iris.Constraint(pressure=lambda i: i in pvalues)
            #     cube = cube.extract(pconstr)
            ndarr = cube.data.astype(np.float32)
            if cube.attributes['short_name'] == 'slp':
                ndarr = ndarr[np.newaxis, ::-1, :]
            elif cube.attributes['short_name'] == 'lsm':
                # if dt1 == dtrange[0]:
                ndarr = ndarr[np.newaxis, ::-1, :]
                # else:
                #    ndarr = None
            else:
                ndarr = ndarr[::-1, ::-1, :]
            all_data.append(ndarr)
        all_data = np.concatenate(all_data)

        all_data.tofile(outdir / out_file_mask.format(dt=dt1))

    with (outdir / 'timecard').open('w') as fw:
        for dt1 in times:
            fw.write(f'{dt1:%Y %m %d %H %M}\n')


if __name__ == '__main__':
    sys.exit(main())
