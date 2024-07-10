#! /usr/bin/env python3
"""
Pysep plotting utilities
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from obspy.core.event import Catalog
from obspy.core.inventory import Inventory
from obspy.geodetics import gps2dist_azimuth, degrees2kilometers

from pysep import logger
from pysep.utils.fmt import format_event_tag


def plot_geometric_spreading(distances, max_amplitudes,
                             geometric_spreading_factor, geometric_k_value,
                             units="km", station_ids=None, include=None,
                             show=False, save="geometric_spreading.png",
                             ymax=None, **kwargs):
    """
    Plot a scatter plot of peak amplitudes vs distance and a corresponding line
    representing geometric scatter. Allow ignoring stations by index.

    .. note::
        see recsec.RecordSection._calculate_geometric_spreading() for more
        information on this function

    :type distances: list or np.array
    :param distances: source-receiver distance in units `units`
    :type max_amplitudes: list or np.array
    :param max_amplitudes: peak amplitudes, units irrelevant
    :type units: st
    :param units: 'km' or 'deg' for kilometers or degrees, respectively. used
        for the x-axis label only
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

    # Convert from units deg -> km if desired
    if "km" in units:
        plot_distances = degrees2kilometers(distances)
    else:
        plot_distances = distances

    # Plot all stations in the list, specially mark excluded stations
    for i, (x, y) in enumerate(zip(plot_distances, max_amplitudes)):
        # list conversion so we can check list length as a boolean
        if list(include) and i not in include:
            ec = "r"
        else:
            ec = "k"
        if ymax and y >= ymax:
            plt.scatter(x, ymax, marker="^", edgecolor=ec, c="w", s=50)
        else:
            plt.scatter(x, y, marker="o", edgecolor=ec, c="w", s=50)

    # Plot station ids next to markers if given
    if station_ids is not None:
        for x, y, s in zip(plot_distances, max_amplitudes, station_ids):
            if ymax and y >= ymax:
                plt.text(x, ymax, f"{s}\n{y.max():.2E}", fontsize=8)
            else:
                plt.text(x, y, s, fontsize=8)

    # Need to generate the line that w-vector follows with distances in unit deg
    x = np.arange(1E-7, distances.max() * 1.5, 0.01)  # units: deg
    sin_d = (np.sin(np.array(x) / (180 / np.pi)))
    y = geometric_k_value / (sin_d ** geometric_spreading_factor)
    lb = f"A(d)={geometric_k_value:.2E} / (sin d)**{geometric_spreading_factor}"
    # Convert w-vector x units deg -> km if desired
    if "km" in units:
        x = degrees2kilometers(x)
    w = plt.plot(x, y, "k--", lw=3, label=lb)

    # Create some legend markers
    m1 = plt.scatter(-1, -1, edgecolor="k", c="w", marker="o",
                     s=50, label="Included", linestyle="None")
    m2 = plt.scatter(-1, -1, edgecolor="r", c="w", marker="o",
                     s=50, label="Excluded", linestyle="None")
    m3 = plt.scatter(-1, -1, edgecolor="r", c="w", marker="^",
                     s=50, label="Out of Bounds", linestyle="None")
    plt.legend(handles=[w[0], m1, m2, m3])

    # Final touches for 'publication-ready' look
    plt.xlabel(f"Source-Receiver Distance ({units})")
    plt.ylabel("Peak Amplitude")
    plt.xlim([0, plot_distances.max() * 1.1])

    if ymax:
        plt.ylim([0, ymax])
    else:
        plt.ylim([0, max_amplitudes.max() * 1.1])

    # One last check for plot aesthetic, degrees vs km.
    if "km" in units:
        xtick_minor = 10
        xtick_major = 100
    else:  # degrees
        xtick_minor = 0.5
        xtick_major = 1

    set_plot_aesthetic(ax=ax, xtick_fontsize=12, ytick_fontsize=12,
                       title_fontsize=13, xlabel_fontsize=13,
                       ylabel_fontsize=13, xtick_minor=xtick_minor,
                       xtick_major=xtick_major)
    plt.grid(visible=True, which="both", axis="y", alpha=0.5, linewidth=1)

    if save:
        plt.savefig(save)
    if show:
        plt.show()


def plot_source_receiver_map(inv, event, subset=None, save="./station_map.png",
                             projection=None, show=False):
    """
    Simple utility to plot station and event on a single map with Cartopy

    :type event: obspy.core.event.Event
    :param event: event to get location from
    :type inv: obspy.core.inventory.Inventory
    :param inv: inventory to get locatiosn from
    :type subset: list of str
    :param subset: allow subsetting stations in the Inventory. Station IDs 
        should match the format output by function `obspy.trace.Trace.get_id()`,
        that is, NN.SSS.LL.CCC (Network, Station, Location, Channel). Location
        and channel information are not used during subsetting.
    :type save: str
    :param save: filename to save station map, defaults: ./station_map.png
    :type projection: str
    :param projection: map projection, 'local, 'ortho', 'global' or None. If
        None, will be determined based on largest src-rcv distance
    :type show: bool
    :param show: show the figure in the GUI
    """
    # Subset the Inventory to only plot certain stations (optional)
    if subset:
        netstas = []
        inv_subset = Inventory()
        for sta_id in subset:
            net, sta, *_ = sta_id.split(".")
            # Assuming NN.SSS should be unique across all networks
            netsta = f"{net}.{sta}"
            if netsta not in netstas:
                inv_subset += inv.select(network=net, station=sta)
                netstas.append(netsta)
        inv = inv_subset

    # Calculate the maximum source receiver distance to determine what 
    # projection we should be in
    if projection is None:
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
        elif max(dists) < 10000:  # ~ < circumference of Earth / 4
            projection = "ortho"
        else:
            projection = "global"
        logger.debug(f"setting projection {projection} due to max src-rcv "
                     f"distance of {max(dists)}km")

    # Use ObsPy's internal Cartopy calls to plot event and stations
    # temporarily change the size of all text objects to get text labels small
    # and make the event size a bit smaller
    with plt.rc_context({"font.size": 6, "lines.markersize": 4}):
        Catalog(events=[event]).plot(projection=projection, resolution="i", 
                                     method="cartopy", label=None, show=False,
                                     title="")
        inv.plot(label=True, show=False, size=12, color="g", method="cartopy", 
                 fig=plt.gcf())

    # Do some post-plot editing of the figure
    ax = plt.gca()

    # Get the extent all points plotted to resize the figure in local
    # since the original bounds will be confined to a small reigon around the
    # event
    if projection == "local":
        lats, lons = [], []
        for net in inv:
            for sta in net:
                lats.append(sta.latitude)
                lons.append(sta.longitude)
        lats.append(event.preferred_origin().latitude)
        lons.append(event.preferred_origin().longitude)
        lats = np.array(lats)
        lons = np.array(lons)

        # Copied verbatim from obspy.imaging.maps.plot_map(), find aspect ratio
        # based on the points plotted. Need to re-do this operation since we
        # have both event and stations
        if min(lons) < -150 and max(lons) > 150:
            max_lons = max(np.array(lons) % 360)
            min_lons = min(np.array(lons) % 360)
        else:
            max_lons = max(lons)
            min_lons = min(lons)
        lat_0 = max(lats) / 2. + min(lats) / 2.
        lon_0 = max_lons / 2. + min_lons / 2.
        if lon_0 > 180:
            lon_0 -= 360
        deg2m_lat = 2 * np.pi * 6371 * 1000 / 360
        deg2m_lon = deg2m_lat * np.cos(lat_0 / 180 * np.pi)
        if len(lats) > 1:
            height = (max(lats) - min(lats)) * deg2m_lat
            width = (max_lons - min_lons) * deg2m_lon
            margin = 0.2 * (width + height)
            height += margin
            width += margin
        else:
            height = 2.0 * deg2m_lat
            width = 5.0 * deg2m_lon
        # Do intelligent aspect calculation for local projection to adjust to
        # figure dimensions
        w, h = plt.gcf().get_size_inches()
        aspect = w / h

        if width / height < aspect:
            width = height * aspect
        else:
            height = width / aspect

        # We are ASSUMING that the event is located at the center of the figure
        # (0, 0), which it should be since we plotted it first
        x0, y0 = 0, 0  # modified from ObsPy function
        ax.set_xlim(x0 - width / 2, x0 + width / 2)
        ax.set_ylim(y0 - height / 2, y0 + height / 2)

    # Hijack the ObsPy plot and make some adjustments to make it look nicer
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False,
                      y_inline=False, linewidth=.25, alpha=.25, color="k")
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {"rotation": 0}
    ax.spines["geo"].set_linewidth(1.5)

    title = (f"{format_event_tag(event)}\n"
             f"NSTA: {len(inv.get_contents()['stations'])} // ")

    # Event may not have have depth or magnitude information, don't let that
    # break the plotting function
    if hasattr(event.preferred_origin(), "depth"):
        title += f"DEPTH: {event.preferred_origin().depth * 1E-3}km // "
    if hasattr(event.preferred_magnitude(), "mag"):
        title += (f"{event.preferred_magnitude().magnitude_type} "
                  f"{event.preferred_magnitude().mag}")

    ax.set_title(title)
    if save:
        plt.savefig(save)
    if show:
        plt.show()
    else:
        plt.close()


def set_plot_aesthetic(
        ax, ytick_fontsize=8., xtick_fontsize=12., tick_linewidth=1.5,
        tick_length=5., tick_direction="in", xlabel_fontsize=10.,
        ylabel_fontsize=10., axis_linewidth=2., spine_zorder=8, spine_top=True,
        spine_bot=True, spine_left=True, spine_right=True, title_fontsize=10.,
        xtick_minor=None, xtick_major=None, ytick_minor=None, ytick_major=None,
        xgrid_major=True, xgrid_minor=True, ygrid_major=True, ygrid_minor=True,
        **kwargs):
    """
    Give a nice look to the output figure by creating thick borders on the
    axis, adjusting fontsize etc. All plot aesthetics should be placed here
    so it's easiest to find.

    :type ax: matplotlib.axes._subplots.AxesSubplot
    :param ax: matplotlib axis figure to set plot aesthetic for
    :type ytick_fontsize: float
    :param ytick_fontsize: fontsize for labels next to Y-axis ticks
    :type xtick_fontsize: float
    :param xtick_fontsize: fontsize for labels next to X-axis ticks
    :type tick_linewidth: float
    :param tick_linewidth: thickness of tick marks for both X and Y axes
    :type tick_length: float
    :param tick_length: length of tick marks for both X and Y axes
    :type tick_direction: str
    :param tick_direction: 'in' for ticks pointing inwards, 'out' for ticks
        pointing outwards
    :type xlabel_fontsize: float
    :param xlabel_fontsize: font size for the X-axis main label (e.g., Time)
    :type ylabel_fontsize: float
    :param ylabel_fontsize: font size for the Y-axis main label (e.g., Ampli)
    :type axis_linewidth: float
    :param axis_linewidth: line thickness for the borders of the figure
    :type spine_zorder: int
    :param spine_zorder: Z order (visibility) of the axis borders (spines)
    :type spine_top: bool
    :param spine_top: toggle on/off the top axis border
    :type spine_bot: bool
    :param spine_bot: toggle on/off the bottom axis border
    :type spine_left: bool
    :param spine_left: toggle on/off the left axis border
    :type spine_right: bool
    :param spine_right: toggle on/off the right axis border
    :type title_fontsize: float
    :param title_fontsize: font size of the main title at the top of the figure
    :type xtick_minor: float
    :param xtick_minor: how often minor tick marks are drawn on X axis
    :type xtick_major: float
    :param xtick_major: how often major tick marks are drawn on X axis
    :type ytick_minor: float
    :param xtick_minor: how often minor tick marks are drawn on Y axis
    :type ytick_major: float
    :param ytick_major: how often major tick marks are drawn on Y axis
    :type xgrid_minor: bool
    :param xgrid_minor: turn on grid lines for each minor X tick
    :type xgrid_major: bool
    :param xgrid_major: turn on grid lines for each major X tick
    :type ygrid_minor: bool
    :param ygrid_minor: turn on grid lines for each minor Y tick
    :type ygrid_major: bool
    :param ygrid_major: turn on grid lines for each minor Y tick
    :rtype: matplotlib.axes._subplots.AxesSubplot
    :return: input axis with aesthetic changed
    """
    ax.title.set_fontsize(title_fontsize)
    ax.tick_params(axis="both", which="both", width=tick_linewidth,
                        direction=tick_direction, length=tick_length)
    ax.tick_params(axis="x", labelsize=xtick_fontsize)
    ax.tick_params(axis="y", labelsize=ytick_fontsize)
    ax.xaxis.label.set_size(xlabel_fontsize)
    ax.yaxis.label.set_size(ylabel_fontsize)

    # Thicken up the bounding axis lines
    for axis, flag in zip(["top", "bottom", "left", "right"],
                          [spine_top, spine_bot, spine_left, spine_right]):
        # Deal with the case where command line users are inputting strings
        if isinstance(flag, str):
            flag = bool(flag.capitalize() == "True")
        ax.spines[axis].set_visible(flag)
        ax.spines[axis].set_linewidth(axis_linewidth)

    # Set spines above azimuth bins
    for spine in ax.spines.values():
        spine.set_zorder(spine_zorder)

    # Set xtick label major and minor which is assumed to be a time series
    if xtick_major:
        ax.xaxis.set_major_locator(MultipleLocator(float(xtick_major)))
    if xtick_minor:
        ax.xaxis.set_minor_locator(MultipleLocator(float(xtick_minor)))
    if ytick_minor:
        ax.yaxis.set_major_locator(MultipleLocator(float(ytick_major)))
    if ytick_major:
        ax.yaxis.set_minor_locator(MultipleLocator(float(ytick_minor)))

    plt.sca(ax)
    if xgrid_major:
        plt.grid(visible=True, which="major", axis="x", alpha=0.5, linewidth=1)
    if xgrid_minor:
        plt.grid(visible=True, which="minor", axis="x", alpha=0.2, linewidth=.5)
    if ygrid_major:
        plt.grid(visible=True, which="major", axis="y", alpha=0.5, linewidth=1)
    if ygrid_minor:
        plt.grid(visible=True, which="minor", axis="y", alpha=0.2, linewidth=.5)

    return ax
