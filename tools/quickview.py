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
import daiquiri, logging  # NOQA
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import matplotlib.patheffects as PathEffects
import numpy as np
from os import getenv
from path import Path
import pandas as pd
import sys
import warnings
# Local files
import utils

# Global settings
iris.FUTURE.netcdf_promote = True
iris.FUTURE.cell_datetime_objects = True
# plt.style.use('awesome')
# Settings for saving figures
fmt = 'png'
svfigkw = dict(format=fmt, dpi=100, bbox_inches='tight')
imgname_mask = 'track_res_{time}.{fmt}'
# Directories
TOPDIR = Path(getenv('HOME')) / 'phd'
ORIG_DATA_DIR = TOPDIR / 'reanalysis' / 'ERA5'
TRACK_RES_DIR = TOPDIR / 'pmc_tracking' / 'pmctrack' / 'output'
PLOT_DIR = TRACK_RES_DIR / 'quickviews'
# File wildcards
ORIG_DATA_FILES = 'era5*.nc'
VORTRACK_FILES = 'vortrack*.txt'
VORMAX_FILES = 'vormax_loc_{kt:04d}.txt'
# Other paths
LOGPATH = Path(__file__).dirname() / 'logs'
SCRIPT = Path(__file__).basename().splitext()[0]
# Subset used in tracking
# e.g.:
#  nx1=120
#  nx2=290
#  ny1=30
#  ny2=70
# TODO: move to args?
xx = slice(1, -1)  # slice(120, 290+1)
yy = slice(1, -1)  # slice(30, 70+1)
# Column names
vor_loc_df_kw = dict(delimiter='\s+',
                     names=['lon', 'lat', 'vo', 'rad', 'vortex_type'])
vortrack_df_kw = dict(delimiter='\s+',
                      names=['lon', 'lat', 'vo', 'kt', 'rad', 'vortex_type'])
# Plotting styles
vort_scl = 1e4  # vorticity factor
vort_kw = dict(cmap=plt.cm.Oranges, vmin=0.1, vmax=5)
slp_scl = 1e-2  # pressure factor
slp_kw = dict(levels=np.arange(900, 1101, 2), colors='#FF0000', linewidths=0.5)
clab_kw = dict(inline=1, fmt='%1.0f', fontsize=10, colors='#FF0000')
vor_loc_kw = dict(linestyle='', marker='o', ms=10, mec='b', mfc='w')
track_past_kw = dict(color='b', linewidth=2)
track_future_kw = dict(color='b', linewidth=0.5, linestyle='--')
track_start_kw = dict(marker='x', color=track_past_kw['color'])
PATH_EFFECTS_ON = True
path_effects = [PathEffects.withStroke(linewidth=2, foreground="w")]
FIGSIZE = (12, 10)

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
    h = ax.pcolormesh(lons, lats, vort*vort_scl, **vort_kw)
    cb = fig.colorbar(h, ax=ax)
    cb.ax.set_title(utils.unit_format(vort_scl**(-1), 's^{-1}'))

    # Sea level pressure
    h = ax.contour(lons, lats, slp*slp_scl, **slp_kw)
    ax.clabel(h, **clab_kw)

    if isinstance(lsm, np.ndarray):
        ax.contourf(lons, lats, lsm, cmap='gray_r', alpha=0.25)


def plot_tracks(fig, ax, counter):
    """ Plot results of the tracking algorithm """
    # Plot tracks
    handles = []
    for fname in sorted(TRACK_RES_DIR.glob(VORTRACK_FILES)):
        df_track = pd.read_csv(fname, **vortrack_df_kw)
        if counter in df_track['kt'].values:
            df_past = df_track[df_track['kt'] <= counter]
            df_future = df_track[df_track['kt'] >= counter]
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
    fname = TRACK_RES_DIR / VORMAX_FILES.format(kt=counter)
    if fname.exists():
        df = pd.read_csv(fname, **vor_loc_df_kw)
        ax.plot(df.lon, df.lat, **vor_loc_kw)


def main(args=None):
    """ Main entry point of the script """
    args = parse_args(args)
    logger.setLevel((5 - min(args.verbose, 4))*10)

    # Lazy-load all the original data
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        ds = iris.load(ORIG_DATA_DIR.listdir(ORIG_DATA_FILES))
    vort = ds.extract_strict('atmosphere_relative_vorticity')
    slp = ds.extract_strict('air_pressure_at_sea_level')
    lsm = ds.extract_strict('land_binary_mask')

    # Date-time span
    t0, tfreq, t1 = args.tspan.split(':')
    t0, t1 = [datetime.strptime(i, '%Y%m%dT%H%MZ') for i in [t0, t1]]

    time_span = pd.date_range(start=t0, end=t1,
                              freq=tfreq).to_pydatetime()

    # Create, if necessary, the output directory
    output_dir = PLOT_DIR
    utils._make_dir(output_dir, args.force)

    # Main time loop
    for counter, idt in enumerate(time_span[:-1]):
        logger.debug(f'Processing {idt:%Y-%m-%d %H:%M}')
        # Prepare data and coordinates
        lons, lats = iris.analysis.cartography.get_xy_grids(vort)
        lons = lons[yy, xx]
        lats = lats[yy, xx]
        time_constr = iris.Constraint(time=idt)
        pres_constr = iris.Constraint(pressure_level=950)
        vo_data = vort.extract(time_constr & pres_constr).data[yy, xx]
        slp_data = slp.extract(time_constr).data[yy, xx]
        lsm_data = lsm.extract(time_constr).data[yy, xx]

        # Prepare figure
        anno_text = f'{idt: %b %d, %H%M}UTC'
        fig, ax = prep_canvas(anno_text=anno_text)
        plot_fields(fig, ax, lons, lats, vo_data, slp_data, lsm=lsm_data)
        plot_tracks(fig, ax, counter+1)

        imgname = imgname_mask.format(time=f'{idt:%Y%m%d%H%M}',
                                      fmt=fmt)
        fig.savefig(output_dir / imgname, **svfigkw)
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
