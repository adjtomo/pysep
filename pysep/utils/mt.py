"""
Moment Tensor related functions for grabbing moment tensors, appending them
to existing event objects, and writing them out in specific formats such as the
CMTSOLUTION format required by SPECFEM
"""
import csv
import requests
import numpy as np
from obspy.core.event import source, Event
from obspy.core.event.base import Comment
from obspy.core.event.source import Tensor
from urllib.error import HTTPError
from obspy import Catalog, UTCDateTime, read_events
from obspy.clients.fdsn import Client

from pysep import logger


def seismic_moment(mt):
    """
    Return the seismic moment based on a moment tensor.
    Can take a list of tensor components, or a Tensor object from ObsPy.

    :type mt: list of floats or obspy.core.event.source.Tensor
    :param mt: the components of the moment tensor M_ij
    :rtype: float
    :return: the seismic moment, in units of N*m
    """
    if isinstance(mt, Tensor):
        # Little one liner to spit out moment tensor components into a list
        mt_temp = [getattr(mt, key) for key in mt.keys()
                   if not key.endswith("errors")]
        assert (len(mt_temp) == 6), "Moment tensor should have 6 components"
        mt = mt_temp
    return 1 / np.sqrt(2) * np.sqrt(sum([_ ** 2 for _ in mt]))


def moment_magnitude(moment, c=10.7):
    """
    Return the moment magitude, M_w, based on a seismic moment. Equation from
    Hanks & Kanamori (1979)

    :type c: float
    :param c: correction factor for conversion, 10.7 for units of N*m,
        16.1 for units of dyne*cm
    :type moment: float
    :param moment: the seismic moment, in units of N*m
    :rtype: float
    :return: moment magnitude, M_w
    """
    return 2 / 3 * np.log10(moment) - c


def half_duration_from_m0(moment):
    """
    Empirical formula for half duration used by Harvard CMT, stated in
    Dahlen and Tromp (1998, p.178).

    :type moment: float
    :param moment: seismic moment in N*m
    :rtype: float
    :return: empirically scaled half duration in unit seconds
    """
    return 2.4E-6 * moment**(1/3)


def mt_transform(mt, method):
    """
    Transform moment tensor between XYZ and RTP coordinates. Used primarily
    to transform GeoNet (John Ristau) moment tensors into the correct coordinate
    system for use in SPECFEM

    Based on Equation from 'Aki and Richards: Quantitative Seismology' book

    .. note::
        Acceptable formats for the parameter `mt`:

        1) [m11, m22, m33, m12, m13, m23]
        2) [mxx, myy, mzz, mxy, mxz, myz]
        3) [mrr, mtt, mpp, mrt, mrp, mtp]

    :type mt: dict
    :param mt: moment tensor in format above
    :type method: str
    :param method: type of conversion, "rtp2xyz" or "xyz2rtp"
    :rtype: dict
    :return: converted moment tensor dictionary
    """
    assert(method in ["xyz2rtp", "rtp2xyz"]), \
        "method must be 'xyz2rtp' or 'rtp2xyz'"

    if method == "xyz2rtp":
        if "m_xx" not in mt.keys():
            print("for xyz2rtp, dict must have keys in xyz")
        m_rr = mt["m_zz"]
        m_tt = mt["m_xx"]
        m_pp = mt["m_yy"]
        m_rt = mt["m_xz"]
        m_rp = -1 * mt["m_yz"]
        m_tp = -1 * mt["m_xy"]
        return {"m_rr": m_rr, "m_tt": m_tt, "m_pp": m_pp, "m_rt": m_rt,
                "m_rp": m_rp, "m_tp": m_tp}
    elif method == "rtp2xyz":
        if "m_tt" not in mt.keys():
            print("for rtp2xyz, dict must have keys in rtp")
        m_xx = mt["m_tt"]
        m_yy = mt["m_pp"]
        m_zz = mt["m_rr"]
        m_xy = -1 * mt["m_tp"]
        m_xz = mt["m_rt"]
        m_yz = -1 * mt["m_rp"]
        return {"m_xx": m_xx, "m_yy": m_yy, "m_zz": m_zz, "m_xy": m_xy,
                "m_xz": m_xz, "m_yz": m_yz}


def get_gcmt_moment_tensor(event=None, origintime=None, magnitude=None,
                           time_wiggle_sec=120, magnitude_wiggle=0.5):
    """
    Query online GCMT moment tensor catalog via URL access for moment tensor
    components of a given event. Searches based on origin time and magnitude
    of an event with a given amount of wiggle room for catalog mismatch of
    origin time and magnitude.

    .. note::
        input is either `event` OR `origintime` AND `magnitude`

    :type event: obspy.core.event.Event
    :param event: Event to use to query for moment tensor
    :type origintime: UTCDateTime or str
    :param origintime: event origin time
    :type magnitude: float
    :param magnitude: centroid moment magnitude for event lookup
    :type time_wiggle_sec: int
    :param time_wiggle_sec: padding on catalog filtering criteria realted to
        event origin time
    :type magnitude_wiggle: float
    :param magnitude_wiggle: padding on catalog filter for magnitude
    :rtype: obspy.core.event.Event
    :return: event object for given earthquake
    """
    if event is None:
        assert(origintime is not None and magnitude is not None), (
            "GCMT moment tensor query requires `event` or `origintime` "
            "and `magnitude"
        )
        origintime = UTCDateTime(origintime)
    else:
        origintime = event.preferred_origin().time
        magnitude = event.preferred_magnitude().mag

    # Determine filename using datetime properties
    month = origintime.strftime('%b').lower()  # e.g. 'jul'
    year_short = origintime.strftime('%y')  # e.g. '19'
    year_long = origintime.strftime('%Y')  # e.g. '2019'

    fid = f"{month}{year_short}.ndk"
    try:
        cat = read_events(
            "https://www.ldeo.columbia.edu/~gcmt/projects/CMT/"
            f"catalog/NEW_MONTHLY/{year_long}/{fid}"
        )
    except requests.HTTPError:
        cat = read_events(
            "http://www.ldeo.columbia.edu/~gcmt/projects/CMT/"
            "catalog/NEW_QUICK/qcmt.ndk"
        )
    # GCMT catalogs contain all events for a span of time
    # filter catalogs using ObsPy to find events with our specifications.
    # Magnitudes and origintimes are not always in agreement between agencies
    # So allow for some wiggle room
    cat_filt = cat.filter(f"time > {str(origintime - time_wiggle_sec)}",
                          f"time < {str(origintime + time_wiggle_sec)}",
                          f"magnitude >= {magnitude - magnitude_wiggle}",
                          f"magnitude <= {magnitude + magnitude_wiggle}",
                          )

    return cat_filt


def get_usgs_moment_tensor(event, time_wiggle_sec=120., magnitude_wiggle=.5,
                           latitude_wiggle_deg=1., longitude_wiggle_deg=1.,
                           depth_wiggle_km=5., **kwargs):
    """
    Query FDSN webservices USGS client for moment tensors using the current
    event definition, which may or may not have been collected via USGS.

    Kwargs passed to Client.get_events() for additional event constraint pars

    :type event: obspy.core.event.Event
    :param event: Event to use to query for moment tensor
    :type time_wiggle_sec: float
    :param time_wiggle_sec: padding on catalog filtering criteria realted to
        event origin time
    :type magnitude_wiggle: float
    :param magnitude_wiggle: +/- padding on magnitude search
    :type latitude_wiggle_deg: float
    :param latitude_wiggle_deg: +/- padding on latitude search
    :type longitude_wiggle_deg: float
    :param longitude_wiggle_deg: +/- padding on longitude search
    :type depth_wiggle_km: float
    :param depth_wiggle_km: +/- padding on depth search
    :rtype: obspy.core.event.Event
    :return: event object for given earthquake
    """
    c = Client("USGS")
    origintime = event.preferred_origin().time
    magnitude = event.preferred_magnitude().mag
    latitude = event.preferred_origin().latitude
    longitude = event.preferred_origin().longitude
    depth = event.preferred_origin().depth * 1E-3

    # Assuming that time, magnitude, and hypocenter are enough to uniquely
    # identify a given earthquake
    try:
        cat = c.get_events(starttime=origintime - time_wiggle_sec,
                           endtime=origintime + time_wiggle_sec,
                           minmagnitude=magnitude - magnitude_wiggle,
                           maxmagnitude=magnitude + magnitude_wiggle,
                           mindepth=depth - depth_wiggle_km,
                           maxdepth=depth + depth_wiggle_km,
                           minlatitude=latitude - latitude_wiggle_deg,
                           maxlatitude=latitude + latitude_wiggle_deg,
                           minlongitude=longitude - longitude_wiggle_deg,
                           maxlongitude=longitude + longitude_wiggle_deg,
                           includeallorigins=True, **kwargs
                           )
    # Broad failure criteria but these are usually FDSNExceptions from ObsPy
    except Exception as e:
        logger.warning(e)
        cat = None
    return cat


def get_geonet_mt(event_id, units, csv_fid=None):
    """
    Focal mechanisms created by John Ristau are written to a .csv file
    located on Github. This function will append information from the .csv file
    onto the Obspy event object so that all the information can be located in a
    single object

    :type event_id: str
    :param event_id: unique event identifier
    :type units: str
    :param units: output units of the focal mechanism, either:
        'dynecm': for dyne*cm  or
        'nm': for Newton*meter
    :type csv_fid: str
    :param csv_fid: optional local path to .csv file containing Ristau catalog
        of moment tensors. If not given, will search a default online URL
        where the catalog is assumed to be stored
    :rtype focal_mechanism: obspy.core.event.FocalMechanism
    :return focal_mechanism: generated focal mechanism
    """
    assert(units in ["dynecm", "nm"]), "units must be 'dynecm' or 'nm'"

    mtlist = query_geonet_mt_catalog(event_id, csv_fid=csv_fid)

    # Match the identifier with Goenet
    id_template = f"smi:local/geonetcsv/{mtlist['PublicID']}/{{}}"

    # Generate the Nodal Plane objects containing strike-dip-rake
    nodal_plane_1 = source.NodalPlane(
        strike=mtlist['strike1'], dip=mtlist['dip1'], rake=mtlist['rake1']
    )
    nodal_plane_2 = source.NodalPlane(
        strike=mtlist['strike2'], dip=mtlist['dip2'], rake=mtlist['rake2']
    )
    nodal_planes = source.NodalPlanes(
        nodal_plane_1, nodal_plane_2, preferred_plane=1
    )

    # Create the Principal Axes as Axis objects
    tension_axis = source.Axis(
        azimuth=mtlist['Taz'], plunge=mtlist['Tpl'], length=mtlist['Tva']
    )
    null_axis = source.Axis(
        azimuth=mtlist['Naz'], plunge=mtlist['Npl'], length=mtlist['Nva']
    )
    pressure_axis = source.Axis(
        azimuth=mtlist['Paz'], plunge=mtlist['Ppl'], length=mtlist['Pva']
    )
    principal_axes = source.PrincipalAxes(
        t_axis=tension_axis, p_axis=pressure_axis, n_axis=null_axis
    )

    # Create the Moment Tensor object with correct units and scaling
    if units == "nm":
        c = 1E-7  # conversion from dyne*cm to N*m
        logger.debug(f"GeoNet moment tensor is in units of Newton*meters")
    elif units == "dynecm":
        c = 1
        logger.debug(f"GeoNet moment tensor is in units of dyne*cm")

    # CV is the conversion from non-units to the desired output units
    cv = 1E20 * c
    seismic_moment_in_nm = mtlist['Mo'] * c

    # Convert the XYZ coordinate system of GeoNet to an RTP coordinate system
    # expected in the CMTSOLUTION file of Specfem
    rtp = mt_transform(mt={"m_xx": mtlist['Mxx']*cv, "m_yy": mtlist['Myy']*cv,
                           "m_zz": mtlist['Mzz']*cv, "m_xy": mtlist['Mxy']*cv,
                           "m_xz": mtlist['Mxz']*cv, "m_yz": mtlist['Myz']*cv
                           },
                       method="xyz2rtp"
                       )
    tensor = source.Tensor(m_rr=rtp['m_rr'], m_tt=rtp['m_tt'],
                           m_pp=rtp['m_pp'], m_rt=rtp['m_rt'],
                           m_rp=rtp['m_rp'], m_tp=rtp['m_tp']
                           )
    # Create the source time function
    source_time_function = source.SourceTimeFunction(
        duration=2 * half_duration_from_m0(seismic_moment_in_nm)
    )

    # Generate a comment for provenance
    comment = Comment(force_resource_id=True,
                      text="Automatically generated by Pyatoa via GeoNet MT CSV"
                      )

    # Fill the moment tensor object
    moment_tensor = source.MomentTensor(
        force_resource_id=True, tensor=tensor,
        source_time_function=source_time_function,
        # !!!
        # This doesn't play nice with obspy.Catalog.write(format='CMTSOLUTION')
        # so ignore the origin id
        # derived_origin_id=id_template.format('origin#ristau'),
        scalar_moment=seismic_moment_in_nm, double_couple=mtlist['DC']/100,
        variance_reduction=mtlist['VR'], comment=comment
        )

    # Finally, assemble the Focal Mechanism. Force a resource id so that
    # the event can identify its preferred focal mechanism
    focal_mechanism = source.FocalMechanism(
        force_resource_id=True, nodal_planes=nodal_planes,
        moment_tensor=moment_tensor, principal_axes=principal_axes,
        comments=[comment]
        )

    return focal_mechanism


def query_geonet_mt_catalog(event_id, csv_fid=None):
    """
    Get moment tensor information from a internal csv file,
    or from an external github repository query.
    Only relevant to the new zealand tomography problem.
    Geonet moment tensors stored with a specific column format.

    :type event_id: str
    :param event_id: unique event identifier
    :type csv_fid: str
    :param csv_fid: optional path to GeoNet CMT solution file that is stored
        locally on disk, will be accessed before querying web service
    :rtype moment_tensor: dict
    :return moment_tensor: dictionary created from rows of csv file
    """
    reader = None
    if csv_fid is not None:
        try:
            reader = csv.reader(open(csv_fid, 'r'), delimiter=',')
        except FileNotFoundError:
            pass

    if reader is None:
        # Request and open the CSV file. Assumed that GeoNet will keep their
        # moment-tensor information in their GitHub repository
        # Last accessed 23.6.19
        geonet_mt_csv = (
            "https://raw.githubusercontent.com/GeoNet/data/master/"
            "moment-tensor/GeoNet_CMT_solutions.csv"
        )
        response = requests.get(geonet_mt_csv)
        if not response.ok:
            raise FileNotFoundError(f"Response from {geonet_mt_csv} not ok")

        reader = csv.reader(response.text.splitlines(), delimiter=',')

    # Parse the CSV file
    for i, row in enumerate(reader):
        # First row contains header information
        if i == 0:
            tags = row
        # First column gives event ids
        if row[0] == event_id:
            values = []
            # Grab the relevant information from the file
            for t, v in zip(tags, row):
                if t == "Date":
                    values.append(UTCDateTime(v))
                elif t == "PublicID":
                    values.append(v)
                else:
                    values.append(float(v))

            moment_tensor = dict(zip(tags, values))
            logger.info(f"geonet moment tensor found for: {event_id}")
            return moment_tensor
    else:
        raise AttributeError(f"no geonet moment tensor found for: {event_id}")


def append_focal_mechanism_to_event(event, method="all", overwrite_focmec=False,
                                    overwrite_event=False, client=None):
    """
    Attempt to find focal mechanism information with a given ObsPy Event object.

    .. note::
        FDSN fetched events are devoid of a few bits of information that are
        useful for our applications, e.g. moment tensor, focal mechanisms.
        This function will perform the conversions and append the necessary
        information to the event located in the dataset.

    :type event: obspy.core.event.Event
    :param event: Event object to append a focal mechanism to.
    :type method: bool
    :param method: try to find correspondig focal mechanism
        using various public catalogs. Currently available:
        'all': Try all available options in order until MT is found
        'USGS': Search the USGS moment tensor catalog
        'GCMT': Search the GCMT moment tensor catalog
        False: Don't attempt to search for moment tensors
    :type client: str
    :param client: Specific `client`s come built-in with specific MT catalogs
        If matching client, will ignore other MT choices:
        'GEONET': will search John Ristau catalog for moment tensors,
    :type overwrite_focmec: bool
    :param overwrite_focmec: If the event already has a focal mechanism,
        overwrite the existing focal mechanism.
    :type overwrite_event: bool
    :param overwrite_event: A new event object is usually retrieved when
        gathering MT from USGS or GCMT. Often the locations/timing of this event
        are less accurate than the input event (which is usually sourced from
        a regional catalog). This parameter controls which event object is
        taken. If `True`, takes the USGS or GCMT catalog information, if `False`
        only takes the focal mechanism attribute.
    :rtype event: obspy.core.event.Event
    :return event: Event with a new focal mechanism if one was found
    :raises TypeError: if event is not provided as an obspy.core.event.Event
    """
    if not isinstance(event, Event):
        raise TypeError(f"`event` must be an ObsPy Event object, "
                        f"not: {type(event)}")
    # If the event already has a focal mechanism attribute, don't gather
    elif hasattr(event, "focal_mechanisms") and \
            event.focal_mechanisms and not overwrite_focmec:
        logger.debug("event already has focal mechanism, will not attempt to"
                     "append new focal mechanism")
        return event
    # Only gather moment tensors if we're already trying to do FDSN stuff
    elif client is None:
        logger.debug("client not specified, will not attempt gathering "
                     "moment tensor")
        return event

    method = method.upper()
    event_id = event.resource_id.id  # assuming datacenter tags ID with event id
    cat = Catalog()
    focal_mechanism = None
    if client.upper() == "GEONET":
        logger.info("querying GeoNet moment tensor catalog")
        focal_mechanism = get_geonet_mt(event_id=event_id, units="nm")
    else:
        # Try 1: Look at USGS catalog
        if method in ["ALL", "USGS"]:
            logger.debug("querying USGS database for moment tensor")
            cat = get_usgs_moment_tensor(event=event)
        # Try 2: Look at GCMT catalog if USGS catalog did not return
        elif (method in ["ALL", "GCMT"]) and len(cat) == 0:
            logger.debug("querying GCMT database for moment tensor")
            cat = get_gcmt_moment_tensor(event=event)
        # Try ?: Add options below for more catalog selection
        # +++++++++++++++++++++++++++++++++++++++++++++++++++
        # If multiple events found for a given set of event criteria, pick first
        if cat is not None:
            if len(cat) > 1:
                logger.warning(f"multiple ({len(cat)}) events found, "
                               f"picking zeroth index")
            # Distinguish `event_new` from `event`, sometimes you still want the
            # catalog location, not the one from USGS or GCMT. Or if nothing was
            # found, then we will return the same event
            event_new = cat[0]
            focal_mechanism = event_new.preferred_focal_mechanism()
    # Append or overwrite focal mechanism or event
    if focal_mechanism is None:
        event_out = event
    else:
        if overwrite_event:
            logger.debug("overwriting input event object with newly gathered "
                         "event containing focal mechanism")
            event_out = event_new
        else:
            logger.debug("appending gathered focal mechanism to current event")
            event_out = event.copy()
            event_out.focal_mechanisms = [focal_mechanism]
            event_out.preferred_focal_mechanism_id = focal_mechanism.resource_id

    return event_out


class Source:
    """
    A generic Source object to characterize FORCESOLUTION files and SPECFEM2D
    SOURCE files without breaking the architechture of Pyatoa which was built
    around CMTSOLUTIONs and ObsPy Event objects

    Essentially this class tries to mimic the ObsPy Event object and return
    required information that is queried throughout a Pyatoa workflow
    """
    def __init__(self, resource_id, origin_time, longitude, latitude, depth):
        """
        Only define the essential values required of a source

        :type resource_id: str
        :param resource_id: unique label for the event
        :type origin_time: str or UTCDateTime
        :param origin_time: origin time for the event
        :type longitude: float
        :param longitude: longitude or X value of the event in the domain
        :type latitude: float
        :param latitude: latitude or Y value of the event in the domain
        :type depth: float
        :param depth: depth in km, inverted Z axis, positive values means deeper
        """
        self.id = resource_id
        self.time = UTCDateTime(origin_time)
        self.longitude = float(longitude)
        self.latitude = float(latitude)
        self.depth = float(depth)

    def preferred_origin(self):
        """
        Convenience function to mimic behavior of ObsPy Event object
        """
        return self

    @property
    def resource_id(self):
        """
        Convenenience function to mimic behavior of ObsPy Event object
        """
        return self

