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
from obspy.geodetics import locations2degrees
from pysep import logger
from pysep.utils.fmt import index_cat, get_data_availability
from pysep.utils.plot import set_plot_aesthetic


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

        # Allow a user-defined bounding box to define the region. Otherwise
        # will be set based on the min/max values of the Cat and Inv
        self._user_min_lat = min_lat
        self._user_max_lat = max_lat
        self._user_min_lon = min_lon
        self._user_max_lon = max_lon

        # Allow user to define data availability, otherwise it will be
        # calculated on the fly based on the `inv` object
        self._user_data_avail = data_avail

        # Instantiate Nones and then fill up internal lists with metadata getter
        self.evlats = None
        self.evlons = None
        self.evids = None
        self.stalats = None
        self.stalons = None
        self.staids = None
        self.min_lat = None
        self.max_lat = None
        self.min_lon = None
        self.max_lon = None
        self.depths = None
        self.mags = None
        self.data_avail = None
        self.navail = None
        self.update_metadata()

    def update_metadata(self, cat=None, inv=None):
        """
        Get metadata like location and event depth and magnitude from a given
        Cat and Inv and set as internal attributes.
        Needs to be as separate function as `threshold_events` will cut down
        the internal catalog representation so this will need to be re-run
        """
        if cat is None:
            cat = self.cat
        else:
            self.cat = cat

        if inv is None:
            inv = self.inv
        else:
            self.inv = inv

        # Get longitude and latitude values from the Catalog and Inventory
        self.evlats = np.array(
            [event.preferred_origin().latitude for event in cat]
        )
        self.evlons = np.array(
            [event.preferred_origin().longitude for event in cat]
        )
        # Get additional information about events
        self.depths = np.array(
            [event.preferred_origin().depth * 1E-3 for event in cat]
        )
        self.mags = np.array(
            [event.preferred_magnitude().mag for event in cat]
        )
        self.evids = np.array(
            [event.resource_id.id for event in cat]
        )

        # Get information on stations
        self.stalats, self.stalons, self.staids = [], [], []
        if inv is not None:
            for net in inv:
                for sta in net:
                    self.stalats.append(sta.latitude)
                    self.stalons.append(sta.longitude)
                    self.staids.append(f"{net.code}.{sta.code}")
        self.stalats = np.array(self.stalats)
        self.stalons = np.array(self.stalons)
        self.staids = np.array(self.staids)

        # Bound domain region based on min/max lat/lons of events and stations
        self.min_lat = \
            self._user_min_lat or np.append(self.evlats, self.stalats).min()
        self.max_lat = \
            self._user_max_lat or np.append(self.evlats, self.stalats).max()
        self.min_lon = \
            self._user_min_lon or np.append(self.evlons, self.stalons).min()
        self.max_lon = \
            self._user_max_lon or np.append(self.evlons, self.stalons).max()

        # Determine data availability based on station "on" time
        self.data_avail = \
            self._user_data_avail or get_data_availability(cat, inv)

        # Count the number of availabler stations for the entire catalog
        self.navail = np.array([len(val) for val in self.data_avail.values()])

    def threshold_catalog(self, zedges=None, min_mags=None, min_data=None):
        """
        Kick out events that fall below a given magnitude range or a given
        data availability range. Allow this to be done for various depth
        ranges or for the entire volume at once.

        .. note::

            Updates internal `cat` Catalog object and metadata in place

        :type zedges: list of float
        :param zedges: depth [km] slices to partition domain into when
            thresholding data
        :type min_mags: int or list
        :param min_mags: a list of minimum magnitude thresholds for each depth
            slice. If `zedges` is None, should be a list of length==1, which
            provides minimum magnitude for entire catalog.
            Elif `zedges` is given, should be a list of len(zedges)-1, which
            defines minimum magnitude for each depth bin.
            For example if zedges=[0, 35, 400], then one example is
            min_mags=[4, 6]. Meaning between 0-34km the minimum magnitude is 4,
            and between 35-400km the minimum magnitude is 6.
        :type min_data: int or list
        :param min_data: an integer or list of length len(zedges)-1 that defines
            the minimum number of stations on for a given event origin time,
            which allows user to prioritize stations with available data
        """
        if zedges is None:
            logger.warning("`zedges` not set, all depth values will be "
                           "weighted equally")
            zedges = [min(self.depths), max(self.depths)]
        if min_mags is not None:
            # Allow integer value of min mags as a blanket for all depths
            if isinstance(min_mags, int):
                min_mags = np.array([min_mags] * (len(zedges) - 1))
            assert(len(min_mags) == len(zedges) - 1), (
                f"`min_mags` must be a list of magnitudes with length of "
                f"`zedges` - 1 ({len(zedges)-1}"
            )
        if min_data is not None:
            # Allow integer value of min data as a blanket for all depths
            if isinstance(min_data, int):
                min_data = np.array([min_data] * (len(zedges) - 1))
            assert (len(min_data) == len(zedges) - 1), (
                f"`min_data` must be a list of available stations with "
                f"length of `zedges` - 1 ({len(zedges) - 1}"
            )

        cat_flag = np.ones(len(self.cat), dtype=int)
        for i, z_top in enumerate(zedges[:-1]):
            z_bot = zedges[i+1]

            # Deal with one depth partition at a time
            idxs = np.where(
                (np.abs(self.depths) < np.abs(z_bot)) &
                (np.abs(self.depths) >= np.abs(z_top))
            )[0]

            # Kick out events with magnitudes below a certain threshold
            if min_mags is not None:
                # Remove events that clear the threshold magnitude criteria
                idxs_remove = np.delete(
                    idxs, np.where(self.mags[idxs] >= min_mags[i]),
                    axis=0
                )

                if idxs_remove.size > 0:
                    logger.info(f"{z_top:.2f}<=Z<{z_bot:.2f} min mag "
                                f"(M{min_mags[i]}) threshold matched "
                                f"{len(idxs_remove)} events")
                    cat_flag[idxs_remove] *= 0
            # Kick out events that do not have enough stations available
            if min_data is not None:
                # Remove events that clear the threshold data availability
                idxs_remove = np.delete(
                    idxs, np.where(self.navail[idxs] >= min_data[i]),
                    axis=0
                )
                if idxs_remove.size > 0:
                    logger.info(f"{z_top:.2f}<=Z<{z_bot:.2f} min data avail "
                                f"(N={min_data[i]}) threshold matched "
                                f"{len(idxs_remove)} events")
                    cat_flag[idxs_remove] *= 0

        logger.info(f"event thresholding removed "
                    f"{len(self.cat) - cat_flag.sum()} events from cat")

        # Updates the internal Catalog and metadata in place
        self.cat = index_cat(cat=self.cat, idxs=np.where(cat_flag == 1)[0])
        self.update_metadata()

    def calculate_srcrcv_weights(self, cat=None, inv=None, write="weights.txt",
                                 plot=False, show=False, save="srcrcvwght.png"):
        """
        Calculate event and station specific weights based on event geographic
        weights, and event-dependent station geographic weights which take into
        account data availability for each event.

        Uses the internal `cat` and `inv` attributes of the class. If
        declustering was run to thin out the catalog, run `update_metadata`
        first to update internal attributes which are used for defining weights.

        :type cat: obspy.core.catalog.Catalog
        :param cat: Optional, catalog of events to consider. If none given,
            will use whatever internal catalog is available
        :type inv: obspy.core.inventory.Inventory
        :param inv: Inventory of stations to consider
        :type plot: bool
        :param plot: plot source weights and an average of station weights on
            a single figure.
        :type write: str
        :param write: filename used to write station weights to text file.
            if None, will not write
        :type show: bool
        :param show: show figure in GUI
        :type save: str
        :param save: if given, file id for the name of the output figure to save
            if not given, will not save figure
        :rtype: np.array
        :return: a 2D array where each row corresponds to a station and each
            column corresponds to an event. The value for a given row and
            column is the weight of that station w.r.t all other available
            stations. Stations that were not available (not on) are given a
            weight of 0
        """
        if cat is not None or inv is not None:
            self.update_metadata(cat, inv)

        logger.info("calculating source receiver weights")
        # Used for index lookups
        sorted_idx = np.argsort(self.staids)

        # The summation of source weights should equal the number of sources
        source_weights = self._get_weights(lons=self.evlons, lats=self.evlats,
                                           norm="avg")

        # For each source we get weights based only on available stations
        station_weights = []
        for i, eid in enumerate(self.evids):
            data_avail = np.array(self.data_avail[eid])
            logger.info(f"'{eid}' has {len(data_avail)} stations available")

            # Determine the indices of all stations available based on sta ID
            idx_avail = np.searchsorted(self.staids[sorted_idx], data_avail)

            sta_avail_weights = self._get_weights(lons=self.stalons[idx_avail],
                                                  lats=self.stalats[idx_avail],
                                                  norm="avg"
                                                  )

            # Multiply weights by the weight of the given source
            sta_avail_weights *= source_weights[i]

            # Place these weights within a full size array, fill blanks w/ 0
            # Assuming that the `get_weights` algorithm will not return 0 weight
            sta_full_weights = np.zeros(len(self.staids))  # zeros
            sta_full_weights[idx_avail] = sta_avail_weights

            station_weights.append(sta_full_weights)

        # Rows: events; Cols: stations
        station_weights = np.array(station_weights)

        # Write out a text file with weights
        if write is not None:
            header = (f"ROWS:{' '.join([_ for _ in self.evids])}\n" 
                      f"COLS:{' '.join([_ for _ in self.staids])}"
                      )
            np.savetxt(write, station_weights, header=header, fmt="%.3f")

        # Normalize by the average for each station only taking into account
        # available data (i.e. ignoring zero weight values). For plotting only
        if plot:
            avg_sta_weight = (station_weights.T.sum(axis=1) /
                              np.count_nonzero(station_weights.T, axis=1)
                              )
            title = (f"Geographic Weighting (NSTA={len(self.stalons)}, " 
                     f"NEVNT={len(self.evlons)})")
            self.plot(inv=self.inv, show=show, save=save, title=title,
                      color_by="custom", event_colors=source_weights,
                      station_colors=avg_sta_weight,
                      event_cbar_label="event weights",
                      station_cbar_label="average station weights",
                      cmap="RdBu", station_cmap="PRGn",
                      station_markersize=75, connect_data_avail=True,
                      event_size_scale=(60, 300)
                      )

        return station_weights

    def _get_weights(self, lons, lats, norm=None, plot=False,
                     save="reference_distance_scan.png"):
        """
        Given a set of coordinates (either stations or earthquakes), calculate
        a reference distance and a set of weights for each coordinate.
        Reference distance is chosen as one-third the largest possible value for
        all possible reference distances

        :type lons: np.array
        :param lons: array of longitude values
        :type lats: np.array
        :param lats: array of latitude values
        :type plot: bool
        :param plot: plot a simple scatterplot showing the weights for each
            station
        :type norm: str
        :param norm: how to normalize the weights
            - None: don't normalize, provide raw weights
            - 'max': normalize by the maximum weight
            - 'len': normalize by the length of the array
            - 'avg': normalize by the mean weight value
        :rtype: np.array
        :return: relative, normalized weights for each lon/lat pair
        """
        # Weights will be 1 if any station amount less than 2
        if len(lons) <= 2:
            weights = np.ones(len(lons))
            return weights

        # Calculate the matrix of inter-point distances (source or receiver)
        dists = []
        for i, (lat_i, lon_i) in enumerate(zip(lats, lons)):
            dists_ij = []
            for j, (lat_j, lon_j) in enumerate(zip(lats, lons)):
                if i == j:
                    dist_m = 0
                else:
                    dist_m = locations2degrees(lat1=lat_i, long1=lon_i,
                                               lat2=lat_j, long2=lon_j)
                dists_ij.append(dist_m)
            dists.append(dists_ij)
        dists = np.array(dists)

        # Calculate a range of reference distances
        min_ref_dist = 1  # km
        max_ref_dist = np.ceil(dists.max())
        ref_dists = np.arange(min_ref_dist, max_ref_dist, 0.25)

        # Loop through reference distances and calculate weights
        ratios = []
        for ref_dist in ref_dists:
            weights_try = self._covert_dist_to_weight(dists, ref_dist)
            ratios.append(max(weights_try) / min(weights_try))
        ratios = np.array(ratios)

        # Select the reference distance that provides 1/3 the max amplitude
        idx = (np.abs(ratios - ratios.max()//3)).argmin()
        optimal_ref_dist = ref_dists[idx]

        logger.debug(f"ref. dist. = {optimal_ref_dist} chosen from "
                     f"{min_ref_dist} <= d <= {max_ref_dist}")

        # Convert relative inter-source distances to a given weight
        weights = self._covert_dist_to_weight(dists, optimal_ref_dist)

        # Allow normalization for a variety of types
        if norm is None:
            normalize = 1
        elif norm == "avg":
            normalize = sum(weights) / len(weights)
        elif norm == "len":
            normalize = len(weights)
        elif norm == "max":
            normalize = max(weights)
        else:
            raise NotImplementedError(f"`norm` cannot be {norm}")
        weights /= normalize

        # Plot the reference distance scan to allow user to double-check
        if plot:
            plt.plot(dists, ratios, "ko-", linewidth=2)
            plt.axhline(dists[idx], c="r", linewidth=2, label="optimal value")
            plt.xlabel("Reference Distance (deg)")
            plt.ylabel("Condition Number")
            plt.title("Scan for Reference Distance")
            plt.legend()
            if save:
                plt.savefig(save)
            plt.show()

        return weights

    @staticmethod
    def _covert_dist_to_weight(dists, ref_dist):
        """
        Calculate distance weights from a matrix of distances and a
        reference distance

        :type dists: np.array
        :param dists: array of inter-source distances
        :type ref_dist: float
        :param ref_dist: user-defined reference distance in units km.
            larger values for reference distances increase the sensitivity to
            inter-station distances. lower values tend to reduce scatter.
        """
        weights = []
        for dists_ij in dists:
            # Calculate the weight for station i as the summation of
            # weighted distance to all other stations j (i != j)
            weight_i = 0
            for dist_ij in dists_ij:
                weight_i += np.exp(-1 * (dist_ij / ref_dist) ** 2)
            weights.append(1 / weight_i)

        return np.array(weights)

    def decluster_events(self, cat=None, inv=None, choice="cartesian",
                         zedges=None, min_mags=None, nkeep=1,
                         select_by="magnitude", **kwargs):
        """
        Main logic function for choosing how to decluster events. Allow for
        both cartesian and polar binning of the domain.

        See `_decluster_events_cartesian` and `_decluster_events_polar` for
        specific input parameters to control declustering.


        :type cat: obspy.core.catalog.Catalog
        :param cat: Optional, catalog of events to consider. If none given,
            will use whatever internal catalog is available
        :type inv: obspy.core.inventory.Inventory
        :param inv: Inventory of stations to consider
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
        if cat is not None or inv is not None:
            self.update_metadata(cat, inv)

        acceptable_select_by = ["magnitude", "magnitude_r", "depth", "depth_r",
                                "data", "data_r"]
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

    def _decluster_events_cartesian(self, nx=10, ny=10, zedges=None, nkeep=1,
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
            - data: less data availability prioritized
            - data_r: more data availability prioritized
        :type plot: bool
        :param plot: create a before and after catalog scatter plot to
            compare which events were kept/removed. Plots within the cwd
        :type plot_dir: str
        :param plot_dir: directory to save figures to. file names will be
            generated automatically
        :rtype: obspy.core.catalog.Catalog
        :return: a declustered event catalog
        """
        # Figure out how to prioritize the events in each cell
        if "magnitude" in select_by:
            arr = self.mags
        elif "depth" in select_by:
            arr = self.depths
        elif "data" in select_by:
            arr = self.navail
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

        cat_out = index_cat(cat=self.cat, idxs=np.where(cat_flag == 1)[0])

        if plot:
            # Plot the original catalog
            self._plot_cartesian(
                cat=self.cat, inv=self.inv, xedges=xedges, yedges=yedges,
                title=f"Pre-Declustered Event Catalog N={len(self.cat)}",
                save=os.path.join(plot_dir, f"pre_decluster_crtsn.png")
            )
            # Plot the declustered catalog
            self._plot_cartesian(
                cat=cat_out, inv=self.inv, xedges=xedges, yedges=yedges,
                title=f"Declustered Event Catalog N={len(cat_out)}\n"
                      f"(zedges={zedges} / nkeep={nkeep})",
                save=os.path.join(plot_dir, f"declustered_crtsn.png"),
            )
            # Plot events that were removed during declustering
            cat_rem = index_cat(cat=self.cat, idxs=np.where(cat_flag == 0)[0])

            self._plot_cartesian(
                cat=cat_rem, inv=self.inv, xedges=xedges, yedges=yedges,
                title=f"Removed Events N={len(cat_rem)}",
                save=os.path.join(plot_dir, f"removed_crtsn.png")
            )

        return cat_out

    def _decluster_events_polar(self, ntheta=16, zedges=None, nkeep=1,
                                select_by="magnitude_r", plot=False,
                                plot_dir="./", **kwargs):
        """
        Run the declustering agorithm but partition the domain in polar. That
        is, divide each depth slice into a pie with `ntheta` partitions and
        keep events based on events within each slice of the pie. Option to
        cut each slice of pie by radius (distance from center of domain) and
        put additional constraints (e.g., more distant events require larger
        magnitude).

        :type ntheta: int
        :param ntheta: Number of theta bins to break a polar search into.
            Used to break up 360 degrees, so e.g., `ntheta`==17 will return
            bins of size 22.5 degrees ([0, 22.5, 45., 67.5 .... 360.])
        :type zedges: list of float
        :param zedges: depth [km] slices to partition domain into. Each slice
            will be given equal weighting w.r.t to all other slices, independent
            of slice size. e.g., allows upweighting crustal events
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
            - data: less data availability prioritized
            - data_r: more data availability prioritized
        :type plot: bool
        :param plot: create a before and after catalog scatter plot to
            compare which events were kept/removed. Plots within the cwd
        :type plot_dir: str
        :param plot_dir: directory to save figures to. file names will be
            generated automatically
        :rtype: obspy.core.catalog.Catalog
        :return: a declustered event catalog
        """
        # Figure out how to prioritize the events in each cell
        if "magnitude" in select_by:
            arr = self.mags
        elif "depth" in select_by:
            arr = self.depths
        elif "data" in select_by:
            arr = self.navail
        else:
            raise NotImplementedError

        # Partition the domain in polar with `ntheta` bins
        theta_array = np.linspace(0, 360, ntheta + 1)  # units: deg
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

        cat_out = index_cat(cat=self.cat, idxs=np.where(cat_flag == 1)[0])

        if plot:
            self._plot_polar(
                self.cat, self.inv, evrad, theta_array, mid_lon, mid_lat,
                title=f"Pre-Declustering Event Catalog N={len(self.cat)}",
                save=os.path.join(plot_dir, f"pre_decluster_plr.png")
            )
            self._plot_polar(
                cat_out, self.inv, evrad, theta_array, mid_lon, mid_lat,
                title=f"Declustered Event Catalog N={len(cat_out)}\n"
                      f"(zedges={zedges} / nkeep={nkeep})",
                save=os.path.join(plot_dir, f"decluster_plr.png")
            )
            # Plot events that were removed during declustering
            cat_rem = index_cat(cat=self.cat, idxs=np.where(cat_flag == 0)[0])
            self._plot_polar(
                cat_rem, self.inv, evrad, theta_array, mid_lon, mid_lat,
                title=f"Removed Events N={len(cat_rem)}\n",
                save=os.path.join(plot_dir, f"removed_plr.png")
            )
        return cat_out

    def plot(self, cat=None, inv=None, color_by="depth",
             connect_data_avail=False, vmin=None, vmax=None, title=None,
             cmap="inferno_r", show=True, save=None, equal_scale=False,
             **kwargs):
        """
        Geranalized plot function used to plot an event catalog and station
        inventory.

        :type cat: obspy.core.catalog.Catalog
        :param cat: Catalog of events to consider. Events must include origin
            information `latitude` and `longitude`
        :type inv: obspy.core.inventory.Inventory
        :param inv: Inventory of stations to consider
        :type color_by: str
        :param color_by: how to color the event markers, available are
            - 'depth': color by the event depth
            - 'data': color by data availability for given event
            - 'custom': used by internal plotting routines to provide custom
                color array to the plot
        :type connect_data_avail: bool
        :param connect_data_avail: connect sources and receivers with a thin
            line based on data availability
        :type vmin: float
        :param vmin: min value for the colorbar, defaults to smallest value in
            array defined by `color_by`
        :type vmax: float:
        :param vmax: maximum value for colorbar, defaults to largest value in
            the array defined by `color_by`
        :type title: str
        :param title: custom user-input title for the figure, otherwise defaults
            to useful information about the catalog and inventory
        :type cmap: str
        :param cmap: matplotlib colormap to use for array defined by `color_by`
        :type show: bool
        :param show: show figure in GUI
        :type save: str
        :param save: if given, file id for the name of the output figure to save
            if not given, will not save figure
        :type equal_scale: bool
        :param equal_scale: set the scale of lat and lon equal, False by default
        """
        # Fine-tuning look of the plot with reasonable defaults
        event_marker = kwargs.get("event_marker", "o")
        event_legend_mags = kwargs.get("event_legend_mags", [4, 5, 6])
        event_size_scale = kwargs.get("event_size_scale", (30, 200))
        station_marker = kwargs.get("station_marker", "v")
        station_alpha = kwargs.get("station_alpha", 1)
        station_markersize = kwargs.get("station_markersize", 10)
        append_title = kwargs.get("append_title", True)

        # Kwargs for fine tuning some internal plotting routines. Not necessary
        # for general users
        event_colors = kwargs.get("event_colors", None)
        event_cbar_label = kwargs.get("event_cbar_label", None)
        station_colors = kwargs.get("station_colors", "None")
        station_cbar_label = kwargs.get("station_cbar_label", None)
        station_cmap = kwargs.get("station_cmap", "viridis")

        if cat is None:
            cat = self.cat
            mags = self.mags
            depths = self.depths
            evlats = self.evlats
            evlons = self.evlons
            data_avail = self.data_avail
        else:
            evlats = np.array([e.preferred_origin().latitude for e in cat])
            evlons = np.array([e.preferred_origin().longitude for e in cat])
            mags = np.array([event.preferred_magnitude().mag for event in cat])
            depths = np.array([event.preferred_origin().depth * 1E-3 for
                               event in cat])
            data_avail = get_data_availability(cat, inv)

        if connect_data_avail:
            assert(inv is not None), f"`connect_data_avail` requires `inv`"
            assert(data_avail is not None),  \
                f"`connect_data_avail` requires `data_avail`"

        # Calculate number of stations on for each event
        data = np.array([len(val) for val in data_avail.values()])

        f, ax = plt.subplots()

        # Normalize exp of magnitudes between `a` and `b` to get good size diff.
        legend_mags = sorted(event_legend_mags, reverse=True)
        mags = np.append(mags, legend_mags)
        mags = np.e ** mags
        a, b = event_size_scale
        mags = ((b - a) * (mags - mags.min()) / (mags.max() - mags.min())) + a

        # Choose how to color
        if color_by == "depth":
            color_arr = depths
            cbar_label = "depths [km]"
        elif color_by == "data":
            color_arr = data
            cbar_label = "data availability"
        elif color_by == "custom":
            color_arr = event_colors
            cbar_label = event_cbar_label
        else:
            raise NotImplementedError

        # Plot events as a scatter of markers, sized  by mag, colored by smthng
        sc = ax.scatter(x=evlons, y=evlats, s=mags[:-3], c=color_arr,
                        marker=event_marker, cmap=cmap, edgecolor="k",
                        linewidth=1, vmin=vmin or color_arr.min(),
                        vmax=vmax or color_arr.max(), zorder=10)
        plt.colorbar(sc, label=cbar_label)

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
            sc2 = plt.scatter(x=stalons, y=stalats, c=station_colors,
                              edgecolor="k", linewidth=1, s=station_markersize,
                              marker=station_marker, zorder=9,
                              alpha=station_alpha, cmap=station_cmap
                              )
            # Only used by internal plottign routines
            if station_cbar_label:
                plt.colorbar(sc2, label=station_cbar_label)

        # Connect sources and receivers with a straight line based on
        # data availability. Assuming that data availability order is the
        # same as the latitude and longitude
        # NOTE: Pretty brute force so this may take a while
        if connect_data_avail:
            evids = np.array([event.resource_id.id for event in cat])
            for i, (event, stalist) in enumerate(data_avail.items()):
                assert(evids[i] == event), f"incorrect event ID encountered"
                for netsta in stalist:
                    net, sta = netsta.split(".")
                    inv_ = self.inv.select(network=net, station=sta)
                    assert(len(inv_) == 1), f"too many stations found"
                    stalat = inv_[0][0].latitude
                    stalon = inv_[0][0].longitude
                    plt.plot([evlons[i], stalon], [evlats[i], stalat], c="k",
                             linewidth=0.5, alpha=0.05, zorder=8)
            if title is None:
                title = ""
            if append_title:
                title += f"\n({data.sum()} source-receiver pairs)"

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
        if equal_scale:
            ax.set_aspect("equal")

        # Generate a scale bar for magnitude sizes. Must do this after setting
        # the x limit because we will plot based on axis display coordinates
        # Display the `legend_mags` as black markers at the top left of axis
        for i, x in enumerate([0.02, 0.05, 0.095]):
            index = -1 * (i + 1)   # -1, -2, -3,...
            ax.scatter(x, 0.95, s=mags[index], c="w", edgecolor="k",
                       linewidth=1, transform=ax.transAxes, zorder=15)
            ax.text(x, 0.9, s=f"{int(legend_mags[index]):d}", c="k",
                    fontsize=7, transform=ax.transAxes, ha="center", zorder=15)

        # Adjust linewidths and fontsizes to make the figure look nicer
        set_plot_aesthetic(ax, xtick_major=5, xtick_minor=1, ytick_major=5,
                           ytick_minor=1, ytick_fontsize=12,
                           xtick_fontsize=12, xlabel_fontsize=15,
                           ylabel_fontsize=15, title_fontsize=15,
                           xgrid_major=True, xgrid_minor=False,
                           ygrid_major=True, ygrid_minor=False
                           )

        if save:
            plt.savefig(save)
        if show:
            plt.show()

        return f, ax

    def _plot_cartesian(self, cat, inv, xedges, yedges, title, save):
        """Convenience function to plot catalogs with edges"""
        f, ax = self.plot(
            cat=cat, inv=inv, show=False, save=None, title=title,
            color_by="depth", vmin=0, vmax=self.depths.max()
        )
        lkwargs = dict(c="k", linewidth=0.5, alpha=0.5)

        # Add gridlines to the plot
        for xe in xedges:
            ax.axvline(xe, **lkwargs)
        for ye in yedges:
            ax.axhline(ye, **lkwargs)

        plt.savefig(save)
        plt.close()

    def _plot_polar(self, cat, inv, evrad, theta_array, mid_lon, mid_lat,
                    title, save):
        """Convenience function to plot catalogs with radial bin lines"""
        f, ax = self.plot(
            cat=cat, inv=inv, show=False, save=None, title=title,
            color_by="depth", vmin=0, vmax=self.depths.max()
        )

        ax.scatter(mid_lon, mid_lat, c="y", edgecolor="k", marker="o", s=10,
                   linewidth=1)

        # Plot each of the radial bins
        for theta in theta_array:
            x = evrad.max() * np.cos(theta * np.pi / 180) + mid_lon
            y = evrad.max() * np.sin(theta * np.pi / 180) + mid_lat
            ax.plot([mid_lon, x], [mid_lat, y], "k-", alpha=0.3)

        plt.savefig(save)
        plt.close()
