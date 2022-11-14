#! /usr/bin/env python3
"""
Pysep plotting utilities
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from obspy.core.event import Catalog
from obspy.geodetics import gps2dist_azimuth

from pysep import logger
from pysep.utils.fmt import format_event_tag


def plot_geometric_spreading(distances, max_amplitudes,
                             geometric_spreading_factor, geometric_k_value,
                             station_ids=None, include=None, show=False,
                             save="geometric_spreading.png", ymax=None,
                             **kwargs):
    """
    Plot a scatter plot of peak amplitudes vs distance and a corresponding line
    representing geometric scatter. Allow ignoring stations by index.

    .. note::
        see recsec.RecordSection._calculate_geometric_spreading() for more
        information on this function

    :type distances: list or np.array
    :param distances: source-receiver distance in units of degrees
    :type max_amplitudes: list or np.array
    :param max_amplitudes: peak amplitudes, units irrelevant
    :type geometric_spreading_factor: float
    :param geometric_spreading_factor: exponent for spreading eq
    :type geometric_k_value: float
    :param geometric_k_value: constant scale factor for spreading eq
    :type station_ids: list or np.array
    :param station_ids: station ids for optional text labels
    :type include: list or np.array
    :param include: indices of stations included in equation, used for colors
    :type show: bool
    :param show: show the figure after generation
    :type save: str
    :param save: file name to save figure
    :type ymax: float
    :param ymax: set the ceiling y limit manually (incase some stations are
        scaled improperly which would throw off all the scaling). Also cuts
        off any stations above this limit and plots them differently
    """
    dpi = kwargs.get("dpi", 100)
    figsize = kwargs.get("figsize", (12, 8))

    logger.info(f"plotting geometric spreading: '{save}'")
    f, ax = plt.subplots(1, dpi=dpi, figsize=figsize)

    # Plot all stations in the list, specially mark excluded stations
    for i, (x, y) in enumerate(zip(distances, max_amplitudes)):
        if include.any() and i not in include:
            ec = "r"
        else:
            ec = "k"
        if ymax and y >= ymax:
            plt.scatter(x, ymax, marker="^", edgecolor=ec, c="w", s=50)
        else:
            plt.scatter(x, y, marker="o", edgecolor=ec, c="w", s=50)

    # Plot station ids next to markers if given
    if station_ids is not None:
        for x, y, s in zip(distances, max_amplitudes, station_ids):
            if ymax and y >= ymax:
                plt.text(x, ymax, f"{s}\n{y.max():.2E}", fontsize=8)
            else:
                plt.text(x, y, s, fontsize=8)

    # Need to generate the line that w-vector follows
    x = np.arange(1E-7, distances.max() * 1.5, 0.01)  # units: deg
    sin_d = (np.sin(np.array(x) / (180 / np.pi)))
    y = geometric_k_value / (sin_d ** geometric_spreading_factor)
    l = f"A(d)={geometric_k_value:.2E} / (sin d)**{geometric_spreading_factor}"
    w = plt.plot(x, y, "k--", lw=3, label=l)

    # Create some legend markers
    m1 = plt.scatter(-1, -1, edgecolor="k", c="w", marker="o",
                     s=50, label="Included", linestyle="None")
    m2 = plt.scatter(-1, -1, edgecolor="r", c="w", marker="o",
                     s=50, label="Excluded", linestyle="None")
    m3 = plt.scatter(-1, -1, edgecolor="r", c="w", marker="^",
                     s=50, label="Out of Bounds", linestyle="None")
    plt.legend(handles=[w[0], m1, m2, m3])

    # Final touches for 'publication-ready' look
    plt.title("Geometric Spreading Correction Factor: A(d)")
    plt.xlabel("Source-Receiver Distance (deg)")
    plt.ylabel("Peak Amplitude")
    plt.xlim([0, x.max()])
    if ymax:
        plt.ylim([0, ymax * 1.25])
    else:
        plt.ylim([0, y.max()])
    set_plot_aesthetic(ax=ax, xtick_fontsize=12, ytick_fontsize=12,
                       title_fontsize=13, label_fontsize=13,
                       xtick_minor=0.5, xtick_major=1)
    plt.grid(visible=True, which="both", axis="y", alpha=0.5, linewidth=1)

    if save:
        plt.savefig(save)
    if show:
        plt.show()


def plot_source_receiver_map(inv, event, save="./station_map.png",
                             projection=None, show=False):
    """
    Simple utility to plot station and event on a single map with Cartopy

    :type event: obspy.core.event.Event
    :param event: event to get location from
    :type inv: obspy.core.inventory.Inventory
    :param inv: inventory to get locatiosn from
    :type save: str
    :param save: filename to save station map, defaults: ./station_map.png
    :type projection: str
    :param projection: map projection, 'local, 'ortho', 'global' or None. If
        None, will be determined based on largest src-rcv distance
    :type show: bool
    :param show: show the figure in the GUI
    """
    if projection is None:
        # Calculate the maximum source receiver distance to determine what
        # projection we should be in
        dists = []
        event_latitude = event.preferred_origin().latitude
        event_longitude = event.preferred_origin().longitude
        for net in inv:
            for sta in net:
                dist_m, _, _ = gps2dist_azimuth(lat1=event_latitude,
                                                lon1=event_longitude,
                                                lat2=sta.latitude,
                                                lon2=sta.longitude
                                                )
                dists.append(dist_m * 1E-3)
        if max(dists) < 5000:  # km
            projection = "local"
        else:
            projection = "ortho"
        logger.debug(f"setting projection {projection} due to max src-rcv "
                     f"distance of {max(dists)}km")

    # Temp change the size of all text objects to get text labels small
    # and make the event size a bit smaller
    with plt.rc_context({"font.size": 6, "lines.markersize": 4}):
        inv.plot(projection=projection, resolution="i", label=True, show=False,
                 size=12, color="g", method="cartopy")
        # !!! fig=plt.gca() kwarg is a workaround for bug in
        # !!! ObsPy==1.3.0 (obspy.imaging.maps._plot_cartopy_into_axes)
        Catalog(events=[event]).plot(method="cartopy", label=None,
                                     fig=plt.gca(), show=False)

    # Hijack the ObsPy plot and make some adjustments to make it look nicer
    ax = plt.gca()
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False,
                      y_inline=False, linewidth=.25, alpha=.25, color="k")
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {"rotation": 0}
    ax.spines["geo"].set_linewidth(1.5)

    title = (f"{format_event_tag(event)}\n"
             f"NSTA: {len(inv.get_contents()['stations'])} // "
             f"DEPTH: {event.preferred_origin().depth * 1E-3}km // "
             f"{event.preferred_magnitude().magnitude_type} "
             f"{event.preferred_magnitude().mag}"
             )

    ax.set_title(title)
    if save:
        plt.savefig(save)
    if show:
        plt.show()


def set_plot_aesthetic(ax, **kwargs):
    """
    Give a nice look to the output figure by creating thick borders on the
    axis, adjusting fontsize etc. All plot aesthetics should be placed here
    so it's easiest to find.

    .. note::
        This was copy-pasted from Pyatoa.visuals.insp_plot.default_axes()
    """
    ytick_fontsize = kwargs.get("ytick_fontsize", 8)
    xtick_fontsize = kwargs.get("xtick_fontsize", 12)
    tick_linewidth = kwargs.get("tick_linewidth", 1.5)
    tick_length = kwargs.get("tick_length", 5)
    tick_direction = kwargs.get("tick_direction", "in")
    label_fontsize = kwargs.get("label_fontsize", 10)
    axis_linewidth = kwargs.get("axis_linewidth", 2.)
    title_fontsize = kwargs.get("title_fontsize", 10)
    xtick_minor = kwargs.get("xtick_minor", 25)
    xtick_major = kwargs.get("xtick_major", 100)
    spine_zorder = kwargs.get("spine_zorder", 8)

    ax.title.set_fontsize(title_fontsize)
    ax.tick_params(axis="both", which="both", width=tick_linewidth,
                        direction=tick_direction, length=tick_length)
    ax.tick_params(axis="x", labelsize=xtick_fontsize)
    ax.tick_params(axis="y", labelsize=ytick_fontsize)
    ax.xaxis.label.set_size(label_fontsize)
    ax.yaxis.label.set_size(label_fontsize)

    # Thicken up the bounding axis lines
    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(axis_linewidth)

    # Set spines above azimuth bins
    for spine in ax.spines.values():
        spine.set_zorder(spine_zorder)

    # Set xtick label major and minor which is assumed to be a time series
    ax.xaxis.set_major_locator(MultipleLocator(float(xtick_major)))
    ax.xaxis.set_minor_locator(MultipleLocator(float(xtick_minor)))

    plt.sca(ax)
    plt.grid(visible=True, which="major", axis="x", alpha=0.5, linewidth=1)
    plt.grid(visible=True, which="minor", axis="x", alpha=0.2, linewidth=.5)

    return ax
