#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot snapshots of PMC tracking results over maps of vorticity and SLP from netCDF input files.
"""
# Standard library
import argparse
from datetime import datetime
import logging
from pathlib import Path
import sys

# Scientific stack
import daiquiri
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import matplotlib.patheffects as PathEffects
import numpy as np
from octant.core import TrackRun
from octant.parts import TrackSettings
import pandas as pd
import xarray as xr

# Local files
import utils

# Matplotlib settings
plt.style.use("pmctrack_tools.mplstyle")
# Settings for saving figures
imgname_mask = "pmctrack_result_{time}"
# PMC settings
PMCTRACK_DIR = Path(__file__).absolute().parent.parent
conf = TrackSettings(PMCTRACK_DIR / "settings.conf")
# Directories
ORIG_DATA_DIR = PMCTRACK_DIR / conf.datadir
TRACK_RES_DIR = PMCTRACK_DIR / conf.outdir
PLOT_DIR = TRACK_RES_DIR / "quickviews"
# File wildcards
ORIG_DATA_FILES = sorted(
    ORIG_DATA_DIR.glob(utils.common_wildcard(conf.fname_lvl, conf.fname_sfc))
)
ORIG_DATA_FILES.append(ORIG_DATA_DIR / "lsm.nc")
# VORMAX_FILES = "vormax_loc_{kt:%Y%m%d%H%M}.txt"
# Other paths
LOGPATH = Path(__file__).parent / "logs"
SCRIPT = Path(__file__).name
# Column names
columns = ["lon", "lat", "vo", "time", "area", "vortex_type", "slp"]
# Plotting styles
vort_scl = 1e4  # vorticity factor
vort_cmap = plt.cm.Oranges
vort_cmap.set_over("#330000")
vort_cmap.set_under("w", alpha=0)
vort_kw = dict(levels=[0.5, 1, 1.5, 2, 2.5, 3], cmap=vort_cmap, extend="both")
vort_kw_pc = dict(cmap=vort_cmap, vmin=0.5, vmax=3)
slp_scl = 1e-2  # pressure factor
slp_kw = dict(levels=np.arange(900, 1100, 2), colors="#FF0000", linewidths=0.5)
clab_kw = dict(inline=1, fmt="%1.0f", fontsize=10, colors="#FF0000")
vor_loc_kw = dict(marker="o", s=2 ** 7, zorder=20)
vor_type_colors = ("C0", "C2", "C3", "C4")
track_past_kw = dict(color="C0", linewidth=2)
track_future_kw = dict(color="C0", linewidth=0.5, linestyle="--")
track_start_kw = dict(marker="x", color=track_past_kw["color"])
PATH_EFFECTS_ON = True
path_effects = [PathEffects.withStroke(linewidth=2, foreground="w")]
lsm_cmap = LinearSegmentedColormap.from_list("", [(1, 1, 1, 0), (0.5, 0.5, 0.5, 0.5)])

# Logging set up
LOGPATH.mkdir(parents=True, exist_ok=True)
LOGFILE = LOGPATH / f"{SCRIPT}_{datetime.now():%Y%m%d%H%M}.log"
daiquiri.setup(
    level=logging.DEBUG,
    outputs=(daiquiri.output.File(LOGFILE), daiquiri.output.Stream(sys.stdout)),
)
logger = daiquiri.getLogger(SCRIPT)


def parse_args(args=None):
    ap = argparse.ArgumentParser(
        SCRIPT,
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Example of use: coming soon",
    )

    ap.add_argument(
        "-v", "--verbose", action="count", default=0, help="Verbosity (-v, -vv, etc)"
    )

    ap_time = ap.add_argument_group("Time selection")
    ap_time.add_argument(
        "--tspan",
        required=True,
        type=str,
        help="Time span (YYYY_mm_dd_HH_MM:1H:YYYY_mm_dd_HH_MM)",
    )

    ap_plt = ap.add_argument_group(title="Plotting settings")
    # ap_plt.add_argument('-g', '--geoax', action='store_true', default=False,
    #                     help='Plot results on the map')
    ap_plt.add_argument(
        "-f",
        "--force",
        action="store_true",
        default=False,
        help="Overwrite output directory",
    )
    ap_plt.add_argument(
        "--gif",
        action="store_true",
        default=False,
        help="Merge the saved figures into a .gif animation",
    )

    return ap.parse_args(args)


def prep_canvas(anno_text=""):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel("Longitude, degrees")
    ax.set_ylabel("Latitude, degrees")
    if anno_text:
        at = AnchoredText(anno_text, prop=dict(size=10), frameon=True, loc=1)
        ax.add_artist(at)
    return fig, ax


def plot_fields(fig, ax, lons, lats, vort, slp, lsm=None):
    """ Plot vorticity and SLP in the given axes """
    ax.set_xlim(lons.min() - 5, lons.max() + 5)
    ax.set_ylim(lats.min() - 1, lats.max() + 1)
    # Vorticity
    # h = ax.contourf(lons, lats, vort*vort_scl, **vort_kw)
    h = ax.pcolormesh(lons, lats, vort * vort_scl, **vort_kw_pc)
    # import pdb; pdb.set_trace()
    cb = fig.colorbar(h, ax=ax)
    cb.ax.set_title(utils.unit_format(vort_scl ** (-1), "s^{-1}"))

    # Sea level pressure
    h = ax.contour(lons, lats, slp * slp_scl, **slp_kw)
    ax.clabel(h, **clab_kw)

    if isinstance(lsm, np.ndarray):
        ax.pcolormesh(lons, lats, lsm, cmap=lsm_cmap, rasterized=True)


def plot_tracks(fig, ax, trackrun, dt):
    """ Plot results of the tracking algorithm """
    # Plot tracks
    handles = []
    for _, ot in trackrun.time_slice(end=dt).gb:
        h1 = ax.plot(ot.lon[0], ot.lat[0], **track_start_kw)
        h2 = ax.plot(ot.lon, ot.lat, **track_past_kw)
        handles.append(h1)
        handles.append(h2)
    for _, ot in trackrun.time_slice(start=dt).gb:
        h3 = ax.plot(ot.lon, ot.lat, **track_future_kw)
        handles.append(h3)

    if PATH_EFFECTS_ON:
        for h in handles:
            plt.setp(h, path_effects=path_effects)

    # # Plot vorticity centres at the given point
    # fname = TRACK_RES_DIR / VORMAX_FILES.format(kt=dt)
    # if fname.exists():
    #     df = pd.read_csv(fname, **vor_loc_df_kw)
    #     ax.scatter(df.lon, df.lat,
    #                c=df.vortex_type.apply(lambda x: vor_type_colors[x]),
    #                **vor_loc_kw)


def main(args=None):
    """ Main entry point of the script """
    args = parse_args(args)
    logger.setLevel((5 - min(args.verbose, 4)) * 10)

    # Load tracks
    tr = TrackRun(TRACK_RES_DIR, columns=columns)

    # Lazy-load all data
    # with warnings.catch_warnings():
    #     warnings.simplefilter('ignore')
    ds = xr.open_mfdataset(ORIG_DATA_FILES, combine="by_coords")
    vort = ds[conf.vort_name]
    slp = ds[conf.psea_name]
    try:
        lsm = ds["land_binary_mask"]
    except KeyError:
        lsm = ds["lsm"]
    # Prepare data and coordinates
    lons = ds[conf.x_dim]
    lats = ds[conf.y_dim]

    # Date-time span
    t0, tfreq, t1 = args.tspan.split(":")
    t0, t1 = [datetime.strptime(i, "%Y_%m_%d_%H_%M") for i in [t0, t1]]

    time_span = pd.date_range(start=t0, end=t1, freq=tfreq).to_pydatetime()

    # Create, if necessary, the output directory
    output_dir = PLOT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    # Main time loop
    for idt in time_span:
        # pdt = PartialDateTime(**{k: getattr(idt, k)
        #                          for k in PartialDateTime.__slots__})
        logger.debug(f"Processing {idt:%Y-%m-%d %H:%M}")
        vo_data = vort.sel({conf.t_dim: idt, conf.z_dim: conf.vor_lvl}).values
        slp_data = slp.sel({conf.t_dim: idt}).values
        try:
            lsm_data = lsm.sel({conf.t_dim: idt}).values
        except ValueError:
            lsm_data = lsm.values.squeeze()

        # Prepare figure
        anno_text = f"{idt: %b %d, %H%M}UTC"
        fig, ax = prep_canvas(anno_text=anno_text)
        plot_fields(fig, ax, lons, lats, vo_data, slp_data, lsm=lsm_data)
        plot_tracks(fig, ax, tr, idt)

        imgname = output_dir / imgname_mask.format(time=f"{idt:%Y%m%d%H%M}")
        fig.savefig(imgname)
        logger.debug(f"Saved to {imgname}")
        plt.close()

    if args.gif:
        fmt = plt.rcParams["savefig.format"]
        fnames = output_dir / f"*.{fmt}"
        gifname = Path(
            imgname_mask.format(time=args.tspan.replace(":", "-"))
        ).with_suffix(".gif")
        utils.merge_to_gif(fnames, gifname, delay=25, resize=100, output_dir=output_dir)


if __name__ == "__main__":
    sys.exit(main())
