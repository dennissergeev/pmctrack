#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot snapshots of PMC tracking results:
  * iterate over time
  * display maps of relative vorticity and SLP from the input files
  * overlay with PMC tracks
"""
# Standard packages
import argparse
# import cartopy.crs as ccrs
from datetime import datetime
import iris
from iris.experimental import equalise_cubes
# from iris.time import PartialDateTime
import daiquiri, logging  # NOQA
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import matplotlib.patheffects as PathEffects
import numpy as np
from path import Path
import pandas as pd
import sys
import warnings
# Local files
import utils

# Settings for saving figures
fmt = 'png'
svfigkw = dict(format=fmt, dpi=100, bbox_inches='tight')
imgname_mask = 'track_res_{time}.{fmt}'
# PMC settings
PMCTRACK_DIR = Path(__file__).abspath().parent.parent
conf = utils.PMCSettings(PMCTRACK_DIR / 'settings.conf')
# Directories
ORIG_DATA_DIR = (PMCTRACK_DIR / conf.datadir).expand()
TRACK_RES_DIR = (PMCTRACK_DIR / conf.outdir).expand()
# TOPDIR = Path(getenv('HOME')) / 'phd'
# ORIG_DATA_DIR = TOPDIR / 'reanalysis' / 'ERA5'
# TRACK_RES_DIR = TOPDIR / 'pmc_tracking' / 'pmctrack' / 'output'
PLOT_DIR = TRACK_RES_DIR / 'quickviews'
# File wildcards
ORIG_DATA_FILES = 'era5*2011.*.nc'
LAND_MASK_FILE = 'lsm.nc'
VORTRACK_FILES = 'vortrack*.txt'
VORMAX_FILES = 'vormax_loc_{kt:%Y%m%d%H%M}.txt'
# Other paths
LOGPATH = Path(__file__).dirname() / 'logs'
SCRIPT = Path(__file__).basename().splitext()[0]
# Subset used in tracking
lonlat = iris.Constraint()
if getattr(conf, 'lon1', None):
    lonlat &= iris.Constraint(longitude=lambda x: x > getattr(conf, 'lon1')-1)
if getattr(conf, 'lon2', None):
    lonlat &= iris.Constraint(longitude=lambda x: x < getattr(conf, 'lon2')+1)
if getattr(conf, 'lat1', None):
    lonlat &= iris.Constraint(latitude=lambda x: x > getattr(conf, 'lat1')-1)
if getattr(conf, 'lat2', None):
    lonlat &= iris.Constraint(latitude=lambda x: x < getattr(conf, 'lat2')+1)
# Column names
vor_loc_df_kw = dict(delimiter='\s+',
                     names=['lon', 'lat', 'vo', 'rad', 'vortex_type'])
vortrack_df_kw = dict(delimiter='\s+',
                      names=['lon', 'lat', 'vo', 'kt', 'rad', 'vortex_type'])
# Plotting styles
vort_scl = 1e4  # vorticity factor
vort_cmap = plt.cm.Oranges
vort_cmap.set_over('#330000')
vort_cmap.set_under('w', alpha=0)
vort_kw = dict(levels=[0.5, 1, 1.5, 2, 2.5, 3],
               cmap=vort_cmap, extend='both')
vort_kw_pc = dict(cmap=vort_cmap, vmin=0.5, vmax=3)
slp_scl = 1e-2  # pressure factor
slp_kw = dict(levels=np.arange(900, 1101, 2), colors='#FF0000', linewidths=0.5)
clab_kw = dict(inline=1, fmt='%1.0f', fontsize=10, colors='#FF0000')
vor_loc_kw = dict(marker='o', s=2**7, zorder=20)
vor_type_colors = ('C0', 'C2', 'C3', 'C4')
track_past_kw = dict(color='C0', linewidth=2)
track_future_kw = dict(color='C0', linewidth=0.5, linestyle='--')
track_start_kw = dict(marker='x', color=track_past_kw['color'])
PATH_EFFECTS_ON = True
path_effects = [PathEffects.withStroke(linewidth=2, foreground="w")]
FIGSIZE = (12, 10)
lsm_cmap = LinearSegmentedColormap.from_list('', [(1, 1, 1, 0),
                                                  (0.5, 0.5, 0.5, 0.5)])

# Logging set up
utils._make_dir(LOGPATH)
LOGFILE = LOGPATH / f'{SCRIPT}_{datetime.now():%Y%m%d%H%M}.log'
daiquiri.setup(
               level=logging.DEBUG,
               outputs=(daiquiri.output.File(LOGFILE),
                        daiquiri.output.Stream(sys.stdout))
               )
logger = daiquiri.getLogger(SCRIPT)


def parse_args(args=None):
    ap = argparse.ArgumentParser(SCRIPT,
                                 description=__doc__,
                                 formatter_class=argparse.
                                 ArgumentDefaultsHelpFormatter,
                                 epilog='Example of use: coming soon')

    ap.add_argument('-v', '--verbose', action='count',
                    default=0, help='Verbosity (-v, -vv, etc)')

    ap_time = ap.add_argument_group('Time selection')
    ap_time.add_argument('--tspan', required=True,
                         type=str,
                         help='Time span (YYYYmmddTHHMMZ:1H:YYYYmmddTHHMMZ)')

    ap_plt = ap.add_argument_group(title='Plotting settings')
    # ap_plt.add_argument('-g', '--geoax', action='store_true', default=False,
    #                     help='Plot results on the map')
    ap_plt.add_argument('-f', '--force', action='store_true',
                        default=False,
                        help='Overwrite output directory')
    ap_plt.add_argument('--gif', action='store_true', default=False,
                        help='Merge the saved figures into a .gif animation')

    return ap.parse_args(args)


def prep_canvas(anno_text='', figsize=FIGSIZE):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    ax.set_xlabel('Longitude, degrees')
    ax.set_ylabel('Latitude, degrees')
    if anno_text:
        at = AnchoredText(anno_text, prop=dict(size=10), frameon=True, loc=1)
        ax.add_artist(at)
    return fig, ax


def plot_fields(fig, ax, lons, lats, vort, slp, lsm=None):
    """ Plot vorticity and SLP in the given axes """
    # Vorticity
    # h = ax.contourf(lons, lats, vort*vort_scl, **vort_kw)
    h = ax.pcolormesh(lons, lats, vort*vort_scl, **vort_kw_pc)
    # import pdb; pdb.set_trace()
    cb = fig.colorbar(h, ax=ax)
    cb.ax.set_title(utils.unit_format(vort_scl**(-1), 's^{-1}'))

    # Sea level pressure
    h = ax.contour(lons, lats, slp*slp_scl, **slp_kw)
    ax.clabel(h, **clab_kw)

    if isinstance(lsm, np.ndarray):
        ax.pcolormesh(lons, lats, lsm, cmap=lsm_cmap, rasterized=True)


def plot_tracks(fig, ax, dt):
    """ Plot results of the tracking algorithm """
    # Plot tracks
    handles = []
    for fname in sorted(TRACK_RES_DIR.glob(VORTRACK_FILES)):
        df_track = pd.read_csv(fname, parse_dates=['kt'], **vortrack_df_kw)
        if dt in df_track['kt'].dt.to_pydatetime():  # and df_track.shape[0]>6:
            df_past = df_track[df_track['kt'].dt.to_pydatetime() <= dt]
            df_future = df_track[df_track['kt'].dt.to_pydatetime() >= dt]
            h1 = ax.plot(df_past.lon[0], df_past.lat[0], **track_start_kw)
            h2 = ax.plot(df_past.lon, df_past.lat, **track_past_kw)
            h3 = ax.plot(df_future.lon, df_future.lat, **track_future_kw)
            handles.append(h1)
            handles.append(h2)
            handles.append(h3)

    for h in handles:
        if PATH_EFFECTS_ON:
            plt.setp(h, path_effects=path_effects)

    # Plot vorticity centres at the given point
    fname = TRACK_RES_DIR / VORMAX_FILES.format(kt=dt)
    if fname.exists():
        df = pd.read_csv(fname, **vor_loc_df_kw)
        ax.scatter(df.lon, df.lat,
                   c=df.vortex_type.apply(lambda x: vor_type_colors[x]),
                   **vor_loc_kw)


def main(args=None):
    """ Main entry point of the script """
    args = parse_args(args)
    logger.setLevel((5 - min(args.verbose, 4))*10)

    # Lazy-load all the original data
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        ds = iris.load(ORIG_DATA_DIR.listdir(ORIG_DATA_FILES),
                       constraints=lonlat)
        ds += iris.load(ORIG_DATA_DIR.listdir(LAND_MASK_FILE),
                        constraints=lonlat)
    # import pdb; pdb.set_trace()
    equalise_cubes.equalise_attributes(ds)
    ds = ds.concatenate()
    vort = ds.extract_strict('atmosphere_relative_vorticity')
    slp = ds.extract_strict('air_pressure_at_sea_level')
    try:
        lsm = ds.extract_strict('land_binary_mask')
    except iris.exceptions.ConstraintMismatchError:
        lsm = ds.extract_strict('lsm')

    # Date-time span
    t0, tfreq, t1 = args.tspan.split(':')
    t0, t1 = [datetime.strptime(i, '%Y%m%dT%H%MZ') for i in [t0, t1]]

    time_span = pd.date_range(start=t0, end=t1,
                              freq=tfreq).to_pydatetime()

    # Create, if necessary, the output directory
    output_dir = PLOT_DIR
    utils._make_dir(output_dir, args.force)

    # Main time loop
    for idt in time_span:
        # pdt = PartialDateTime(**{k: getattr(idt, k)
        #                          for k in PartialDateTime.__slots__})
        logger.debug(f'Processing {idt:%Y-%m-%d %H:%M}')
        # Prepare data and coordinates
        lons, lats = iris.analysis.cartography.get_xy_grids(vort)
        lons = lons
        lats = lats
        # print(vort.coord('time'))
        # print(pdt)
        time_constr = iris.Constraint(time=idt)
        pres_constr = iris.Constraint(pressure_level=950)
        # import pdb; pdb.set_trace()
        vo_data = vort.extract(time_constr & pres_constr).data
        slp_data = slp.extract(time_constr).data
        lsm_data = lsm.data.squeeze()

        # Prepare figure
        anno_text = f'{idt: %b %d, %H%M}UTC'
        fig, ax = prep_canvas(anno_text=anno_text)
        plot_fields(fig, ax, lons, lats, vo_data, slp_data, lsm=lsm_data)
        plot_tracks(fig, ax, idt)

        imgname = imgname_mask.format(time=f'{idt:%Y%m%d%H%M}',
                                      fmt=fmt)
        fig.savefig(str(output_dir / imgname), **svfigkw)
        logger.debug(f'Saved to {output_dir / imgname}')
        plt.close()

    if args.gif:
        fnames = output_dir / f'*.{fmt}'
        gifname = imgname_mask.format(time=args.tspan.replace(':', '-'),
                                      fmt='gif')
        utils.merge_to_gif(fnames, gifname, delay=25, resize=100,
                           output_dir=output_dir)


if __name__ == '__main__':
    sys.exit(main())
