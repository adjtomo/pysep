#!/usr/bin/evn python3
"""
Event and station declustering and weighting algorithms for adjoint tomography.

Provides routines for downsizing an event catalog by performing a declustering
algorithm to remove nearby stations, and also implements a weighting scheme to
promote unique source receiver pairs.

:authors:
    adjTomo Dev Team (2023)
    Amanda McPherson (2021)
    Carl Tape (2018)
"""
import matplotlib.pyplot as plt
import numpy as np
import os

from obspy import Catalog
from pysep import logger


class Declust:
    """
    Declustering class in charge of declustering and source receiver weighting
    """
    def __init__(self, cat, inv=None, data_avail=None, min_lat=None,
                 max_lat=None, min_lon=None, max_lon=None):
        """
        User-input parameters to determine algorithm behavior

        :type cat: obspy.core.catalog.Catalog
        :param cat: Catalog of events to consider. Events must include origin
            information `latitude` and `longitude`
        :type inv: obspy.core.inventory.Inventory
        :param inv: Inventory of stations to consider
        :type data_avail: dict
        :param data_avail: If None, Declust assumes that all events in `cat`
            were recorded by all stations in `inv`. This is typically not the
            case however, so this dict allows the user to tell `Declust` about
            data availability. Keys of `data_avail` must match resource
            IDs, and values must be lists of station names (NN.SSSS)
        :type min_lat: float
        :param min_lat: optional, minimum latitude for bounding box defining the
            region of interest. If not given, will use the minimum latitude
            in the catalog of events
        :type max_lat: float
        :param max_lat: optional, maximum latitude for bounding box defining the
            region of interest. If not given, will use the maximum latitude
            in the catalog of events
        :type min_lon: float
        :param min_lon: optional, minimum longitude for bounding box defining
            the region of interest. If not given, will use the minimum longitude
            in the catalog of events
        :type max_lon: float
        :param max_lon: optional, maximum longitude for bounding box defining
            the region of interest. If not given, will use the maximum longitude
            in the catalog of events
        """
        self.cat = cat
        self.inv = inv
        self.data_avail = data_avail

        # Get longitude and latitude values from the Catalog and Inventory
        self.evlats = np.array(
            [event.preferred_origin().latitude for event in self.cat]
        )
        self.evlons = np.array(
            [event.preferred_origin().longitude for event in self.cat]
        )
        self.stalats, self.stalons = [], []
        if self.inv is not None:
            for net in self.inv:
                for sta in net:
                    self.stalats.append(sta.latitude)
                    self.stalons.append(sta.longitude)
        self.stalats = np.array(self.stalats)
        self.stalons = np.array(self.stalons)

        # Get additional information about events
        self.depths = np.array(
            [event.preferred_origin().depth * 1E-3 for event in self.cat]
        )
        self.mags = np.array(
            [event.preferred_magnitude().mag for event in self.cat]
        )

        # Bound domain region based on min/max lat/lons of events and stations
        self.min_lat = min_lat or np.append(self.evlats, self.stalats).min()
        self.max_lat = max_lat or np.append(self.evlats, self.stalats).max()
        self.min_lon = min_lon or np.append(self.evlons, self.stalons).min()
        self.max_lon = max_lon or np.append(self.evlons, self.stalons).max()

    def index_cat(self, idxs, cat=None):
        """
        Catalog does not allow indexing by a list of values (e.g., cat[0, 1, 3])
        so this convenience function takes care of that by forming a new
        catalog of events chosen by indices\

        :type idxs: list of int
        :param idxs: list of indices to index catalog by
        :type cat: obspy.core.catalog.Catalog
        :param cat: Catalog to index. If not given defaults to internal Cat
        """
        if cat is None:
            cat = self.cat

        # Workaround as we can't selectively index a Catalog object
        cat_out = Catalog()
        for idx in idxs:
            cat_out.append(cat[idx])
        return cat_out

    def decluster_events(self, choice="cartesian", zedges=None, min_mags=None,
                         nkeep=1, select_by="magnitude", **kwargs):
        """
        Main logic function for choosing how to decluster events

        :type choice: str
        :param choice: choice of domain partitioning, can be one of:
            - cartesian: grid the domain as a cube with `nx` by `ny` cells
            - polar: grid the domain with polar coordinates and `ntheta` bins
        :type zedges: list of float
        :param zedges: depth [km] slices to partition domain into. Each slice
            will be given equal weighting w.r.t to all other slices, independent
            of slice size. e.g., allows upweighting crustal events
        :type min_mags: list
        :param min_mags: a list of minimum magnitude thresholds for each depth
            slice. If `zedges` is None, should be a list of length==1, which
            provides minimum magnitude for entire catalog.
            Elif `zedges` is given, should be a list of len(zedges)-1, which
            defines minimum magnitude for each depth bin.
            For example if zedges=[0, 35, 400], then one example is
            min_mags=[4, 6]. Meaning between 0-34km the minimum magnitude is 4,
            and between 35-400km the minimum magnitude is 6.
        :type nkeep: int or list of int
        :param nkeep: number of events to keep per cell. If `zedges` is None,
            then this must be an integer which defines a blanket value to apply.
            If `zedges` is given, then this must be a list of length
            `zedges` - 1, defining the number of events to keep per cell, per
            depth slice. See `min_mags` definition for example.
        :type select_by: str
        :param select_by: determine how to prioritize events in the cell
            - magnitude (default): largest magnitudes prioritized
            - magnitude_r: smallest magnitudes prioritized
            - depth: shallower depths prioritized
            - depth_r: deeper depths prioritized
            - data: prioritize events which have the most data availability
        """
        acceptable_select_by = ["magnitude", "magnitude_r", "depth", "depth_r",
                                "data"]
        assert(select_by in acceptable_select_by), \
            f"`select_by` must be in {acceptable_select_by}"

        if zedges is None:
            logger.warning("`zedges` not set, all depth values will be "
                           "weighted equally")
            zedges = [min(self.depths), max(self.depths)]
            assert(isinstance(nkeep, int)), (
                f"if `zedges` not given , `nkeep` should be an integer"
            )
            nkeep = [nkeep]
        else:
            # Make `nkeep` iterable even if given as single value
            if isinstance(nkeep, int):
                nkeep = [nkeep]
            assert(len(nkeep) == len(zedges) - 1), (
                f"`nkeep` must be a list of integers with length of "
                f"`zedges` - 1 ({len(zedges)-1}"
            )

        if min_mags:
            assert(len(min_mags) == len(zedges) - 1), (
                f"`min_mags` must be a list of magnitudes with length of "
                f"`zedges` - 1 ({len(zedges)-1}"
            )

        if choice == "cartesian":
            cat = self._decluster_events_cartesian(
                zedges=zedges, min_mags=min_mags, nkeep=nkeep,
                select_by=select_by, **kwargs
            )
        elif choice == "polar":
            cat = self._decluster_events_polar(
                zedges=zedges, min_mags=min_mags, nkeep=nkeep,
                select_by=select_by, **kwargs
            )

        return cat

    def _decluster_events_cartesian(self, nx=10, ny=10, zedges=None,
                                    min_mags=None, nkeep=1,
                                    select_by="magnitude_r", plot=False,
                                    plot_dir="./", **kwargs):
        """
        Decluster event catalog by partitioning the 3D domain in the X, Y and Z
        directions, and then selecting a given number of events in each cell.

        :type nx: int
        :param nx: Number of X/longitude cells to partition domain into
        :type ny: int
        :param ny: Number of Y/latitude cells to partition domain into
        :type zedges: list of float
        :param zedges: depth [km] slices to partition domain into. Each slice
            will be given equal weighting w.r.t to all other slices, independent
            of slice size. e.g., allows upweighting crustal events
        :type min_mags: list
        :param min_mags: a list of minimum magnitude thresholds for each depth
            slice. If `zedges` is None, should be a list of length==1, which
            provides minimum magnitude for entire catalog.
            Elif `zedges` is given, should be a list of len(zedges)-1, which
            defines minimum magnitude for each depth bin.
            For example if zedges=[0, 35, 400], then one example is
            min_mags=[4, 6]. Meaning between 0-34km the minimum magnitude is 4,
            and between 35-400km the minimum magnitude is 6.
        :type nkeep: int or list of int
        :param nkeep: number of events to keep per cell. If `zedges` is None,
            then this must be an integer which defines a blanket value to apply.
            If `zedges` is given, then this must be a list of length
            `zedges` - 1, defining the number of events to keep per cell, per
            depth slice. See `min_mags` definition for example.
        :type select_by: str
        :param select_by: determine how to prioritize events in the cell
            - magnitude (default): largest magnitudes prioritized
            - magnitude_r: smallest magnitudes prioritized
            - depth: shallower depths prioritized
            - depth_r: deeper depths prioritized
            - data: prioritize events which have the most data availability
        :type plot: bool
        :param plot: create a before and after catalog scatter plot to
            compare which events were kept/removed. Plots within the cwd
        :rtype: obspy.core.catalog.Catalog
        :return: a declustered event catalog
        """
        # Figure out how to prioritize the events in each cell
        if "magnitude" in select_by:
            arr = self.mags
        elif "depth" in select_by:
            arr = self.depths
        elif select_by == "data":
            raise NotImplementedError
        else:
            raise NotImplementedError

        # Ensure depth slices are ordered
        zedges = sorted(zedges)

        # Generate equidistant bin edges for X and Y directions.
        h, xedges, yedges = np.histogram2d(
            x=self.evlons, y=self.evlats, bins=[nx, ny],
            range=((self.min_lon, self.max_lon), (self.min_lat, self.max_lat))
        )
        logger.info(f"XY partitioning (events per cell)\n{h}")

        # Use a binary data array to determine which events to return
        cat_flag = np.zeros(len(self.cat))

        # Begin brute force grid search through all possible cells
        for i, z_top in enumerate(zedges[:-1]):
            z_bot = zedges[i+1]
            for j, x_right in enumerate(xedges[:-1]):
                x_left = xedges[j+1]
                for k, y_bot in enumerate(yedges[:-1]):
                    y_top = yedges[k+1]
                    # Determine the indices that define the events in cell
                    # Use absolute values to only avoid negative signs causing
                    # indexing issues
                    idxs = np.where(
                        (np.abs(self.evlons) < np.abs(x_right)) &
                        (np.abs(self.evlons) >= np.abs(x_left)) &
                        (np.abs(self.evlats) >= np.abs(y_bot)) &
                        (np.abs(self.evlats) < np.abs(y_top)) &
                        (np.abs(self.depths) < np.abs(z_bot)) &
                        (np.abs(self.depths) >= np.abs(z_top))
                    )[0]
                    # Ignore empty cells
                    if idxs.size == 0:
                        continue

                    logger.debug(f"{len(idxs)} events found for "
                                 f"{x_left:.2f} <= lon < {x_right:.2f}; "
                                 f"{y_bot:.2f} <= lat < {y_top:.2f}; "
                                 f"{z_top:.2f} <= depth < {z_bot:.2f}"
                                 )

                    # Kick out events with magnitudes below a certain threshold
                    if min_mags:
                        _og_len = len(idxs)
                        # The minimum magnitude is tied to the given depth level
                        idxs = np.delete(
                            idxs, np.where(self.mags[idxs] < min_mags[i]),
                            axis=0
                        )
                        logger.debug(f"minimum magnitude {min_mags[i]} "
                                     f"threshold removed {_og_len - len(idxs)} "
                                     f"events")
                        if idxs.size == 0:
                            continue

                    # Sort the given events by characteristic, reverse if req.
                    sort_arr = arr[idxs].argsort()
                    if "_r" in select_by:
                        sort_arr = sort_arr[::-1]

                    # Select only the first `nkeep` events from this sorted arr
                    idxs = idxs[sort_arr][:nkeep[i]]
                    cat_flag[idxs] = 1  # flip the switch to keep these events

        if not cat_flag.any():
            logger.warning("no events found in declustered catalog")
            return None

        logger.info(f"returning {len(np.where(cat_flag == 1)[0])} events in "
                    f"declustered catalog (og={len(self.cat)})")

        cat_out = self.index_cat(idxs=np.where(cat_flag == 1)[0], cat=self.cat)

        if plot:
            # 1. Plot the original catalog
            f_og, ax_og = self.plot(
                cat=self.cat, inv=self.inv, show=False, save=None,
                title=f"Original Event Catalog N={len(self.cat)}",
                depth_min=0, depth_max=self.depths.max()
            )
            lkwargs = dict(c="k", linewidth=0.5, alpha=0.5)
            # Add gridlines to the plot
            for xe in xedges:
                ax_og.axvline(xe, **lkwargs)
            for ye in yedges:
                ax_og.axhline(ye, **lkwargs)

            plt.savefig(os.path.join(plot_dir, f"original_event_catalog.png"))
            plt.close()

            # 2. Plot the declustered catalog
            f_dc, ax_dc = self.plot(
                cat=cat_out, inv=self.inv, show=False, save=None,
                title=f"Declustered Event Catalog N={len(cat_out)}\n"
                      f"(zedges={zedges} / nkeep={nkeep} / minmags={min_mags})",
                depth_min=0, depth_max=self.depths.max()
            )
            # Add gridlines to the plot
            for xe in xedges:
                ax_dc.axvline(xe, **lkwargs)
            for ye in yedges:
                ax_dc.axhline(ye, **lkwargs)

            plt.savefig(os.path.join(plot_dir, f"decluster_event_catalog.png"))
            plt.close()

        return cat_out

    def _decluster_events_polar(self, ntheta=17, zedges=None, min_mags=None,
                                nkeep=1, select_by="magnitude_r", plot=False,
                                plot_dir="./", **kwargs):
        """
        Run the declustering agorithm but partition the domain in polar. That
        is, divide each depth slice into a pie with `ntheta` partitions and
        keep events based on events within each slice of the pie. Option to
        cut each slice of pie by radius (distance from center of domain) and
        put additional constraints (e.g., more distant events require larger
        magnitude).
        """
        # Figure out how to prioritize the events in each cell
        reverse = False  # for argsort when picking events
        if "magnitude" in select_by:
            arr = self.mags
        elif "depth" in select_by:
            arr = self.depths
        elif select_by == "data":
            raise NotImplementedError
        else:
            raise NotImplementedError

        # Partition the domain in polar with `ntheta` bins
        theta_array = np.linspace(0, 360, ntheta)  # units: deg
        dtheta = theta_array[1] - theta_array[0]

        # Ensure depth slices are ordered
        zedges = sorted(zedges)

        # Convert XY coordinates to polar by considering the middle of the
        # domain as the origin point
        mid_lat = (self.max_lat + self.min_lat) / 2
        mid_lon = (self.max_lon + self.min_lon) / 2
        # Find the distance from origin to each point
        evlats = self.evlats - mid_lat
        evlons = self.evlons - mid_lon
        # Convert from cartesian to polar coordinates
        evrad = np.sqrt(evlons**2 + evlats**2)
        evtheta = np.arctan2(evlats, evlons) * 180 / np.pi + 180  # [0, 360]

        # Use a binary data array to determine which events to return
        cat_flag = np.zeros(len(self.cat))

        # Begin brute force grid search through all possible cells
        for i, z_top in enumerate(zedges[:-1]):
            z_bot = zedges[i + 1]
            for j, theta_start in enumerate(theta_array):
                theta_end = theta_start + dtheta
                # Potentially the number of theta bins does not equally divide
                # 360 so ensure that the last bin covers all the way to 360
                if theta_start == theta_array[-2]:
                    theta_end = 360
                # Determine the indices that define the events in cell
                # Use absolute values to only avoid negative signs causing
                # indexing issues
                idxs = np.where(
                    (evtheta >= theta_start) &
                    (evtheta < theta_end) &
                    (np.abs(self.depths) < np.abs(z_bot)) &
                    (np.abs(self.depths) >= np.abs(z_top))
                )[0]
                # Ignore empty cells
                if idxs.size == 0:
                    continue

                logger.debug(f"{len(idxs)} events found for "
                             f"{theta_start:.2f} <= theta < {theta_end:.2f}")

                # Kick out events with magnitudes below a certain threshold
                if min_mags:
                    _og_len = len(idxs)
                    # The minimum magnitude is tied to the given depth level
                    idxs = np.delete(
                        idxs, np.where(self.mags[idxs] < min_mags[i]),
                        axis=0
                    )
                    logger.debug(f"minimum magnitude {min_mags[i]} "
                                 f"threshold removed {_og_len - len(idxs)} "
                                 f"events")
                    if idxs.size == 0:
                        continue

                # Sort the given events by characteristic, reverse if req.
                sort_arr = arr[idxs].argsort()
                if "_r" in select_by:
                    sort_arr = sort_arr[::-1]

                # Select only the first `nkeep` events from this sorted arr
                idxs = idxs[sort_arr][:nkeep[i]]
                cat_flag[idxs] = 1  # flip the switch to keep these events

        if not cat_flag.any():
            logger.warning("no events found in declustered catalog")
            return None

        logger.info(f"returning {len(np.where(cat_flag == 1)[0])} events in "
                    f"declustered catalog (og={len(self.cat)})")

        cat_out = self.index_cat(idxs=np.where(cat_flag == 1)[0], cat=self.cat)

        if plot:
            # 1. Plot the original catalog
            f_og, ax_og = self.plot(
                cat=self.cat, inv=self.inv, show=False, save=None,
                title=f"Original Event Catalog N={len(self.cat)}",
                depth_min=0, depth_max=self.depths.max()
            )
            ax_og.scatter(mid_lon, mid_lat, c="y", edgecolor="k", marker="o",
                          s=10, linewidth=1)

            # Plot each of the radial bins
            for theta in theta_array:
                x = 2 * evrad.max() * np.cos(theta * np.pi / 180) + mid_lon
                y = 2 * evrad.max() * np.sin(theta * np.pi / 180) + mid_lat
                ax_og.plot([mid_lon, x], [mid_lat, y], "k-", alpha=0.3)

            plt.savefig(os.path.join(plot_dir, f"original_event_catalog.png"))
            plt.close()

            # 2. Plot the declustered catalog
            f_dc, ax_dc = self.plot(
                cat=cat_out, inv=self.inv, show=False, save=None,
                title=f"Declustered Event Catalog N={len(cat_out)}\n"
                      f"(zedges={zedges} / nkeep={nkeep} / minmags={min_mags})",
                depth_min=0, depth_max=self.depths.max()
            )
            ax_dc.scatter(mid_lon, mid_lat, c="y", edgecolor="k", marker="o",
                          s=10, linewidth=1)

            # Plot each of the radial bins
            for theta in theta_array:
                x = evrad.max() * np.cos(theta * np.pi / 180) + mid_lon
                y = evrad.max() * np.sin(theta * np.pi / 180) + mid_lat
                ax_dc.plot([mid_lon, x], [mid_lat, y], "k-", alpha=0.3)

            plt.savefig(os.path.join(plot_dir, f"decluster_event_catalog.png"))
            plt.close()

        return cat_out

    def plot(self, cat=None, inv=None, show=True, save=None, depth_min=0,
             depth_max=None, title=None, cmap="inferno_r"):
        """
        Generate a simple scatter plot of events, colored by depth and sized
        by magnitude.
        """
        if cat is None:
            mags = self.mags
            depths = self.depths
            evlats = self.evlats
            evlons = self.evlons
        else:
            evlats = np.array([e.preferred_origin().latitude for e in cat])
            evlons = np.array([e.preferred_origin().longitude for e in cat])
            mags = np.array([event.preferred_magnitude().mag for event in cat])
            depths = np.array([event.preferred_origin().depth * 1E-3 for
                               event in cat])

        f, ax = plt.subplots()

        # Normalize exp of magnitudes between a and b to get good size diff.
        legend_mags = [6, 5, 4]  # these will be used for legend purposes
        mags = np.append(mags, legend_mags)
        mags = np.e ** mags
        a = 30
        b = 200
        mags = ((b - a) * (mags - mags.min()) / (mags.max() - mags.min())) + a
        sc = ax.scatter(x=evlons, y=evlats, s=mags[:-3], c=depths,
                        cmap=cmap, edgecolor="k", linewidth=1,
                        vmin=depth_min, vmax=depth_max or depths.max(),
                        zorder=10)

        plt.colorbar(sc, label="depth [km]")

        # Plot inventory if provided
        if inv is not None:
            if inv is True:
                stalats = self.stalats
                stalons = self.stalons
            else:
                stalats, stalons = [], []
                for net in inv:
                    for sta in net:
                        stalats.append(sta.latitude)
                        stalons.append(sta.longitude)
                stalats = np.array(stalats)
                stalons = np.array(stalons)
            plt.scatter(x=stalons, y=stalats, c="None", edgecolor="k",
                        linewidth=1, s=10, marker="v", zorder=5, alpha=0.5)

        # Plot attributes
        buff = 0.01
        ax.set_xlim([self.min_lon - np.abs(buff * self.min_lon),
                     self.max_lon + np.abs(buff * self.max_lon)])
        ax.set_ylim([self.min_lat - np.abs(buff * self.min_lat),
                     self.max_lat + np.abs(buff * self.max_lat)])
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        if title is not None:
            ax.set_title(title)
        # May be useful but decided things look better without
        # ax.set_aspect("equal")
        # ax.grid(True)
        plt.tight_layout()

        # Generate a scale bar for magnitude sizes. Must do this after setting
        # the x limit because we will plot based on axis display coordinates
        # Display the `legend_mags` as black markers at the top left of axis
        for i, x in enumerate([0.02, 0.05, 0.095]):
            index = -1 * (i + 1)   # -1, -2, -3,...
            ax.scatter(x, 0.95, s=mags[index], c="w", edgecolor="k",
                       linewidth=1, transform=ax.transAxes, zorder=15)
            ax.text(x, 0.9, s=f"{int(legend_mags[index]):d}", c="k",
                    fontsize=7, transform=ax.transAxes, ha="center", zorder=15)

        if save:
            plt.savefig(save)
        if show:
            plt.show()

        return f, ax
