
"""
miscellaneous helper utilities

20160919 cralvizuri <cralvizuri@alaska.edu>
"""

from obspy.core import UTCDateTime

def otime2eid(otime):
    """
    Convert origin time to origin ID. The origin ID has format: YYYYMMDDHHMMSSsss
    See also eid2otime.

    Example
        otime = "2009-04-07T20:12:55"
        eid = otime2eid(otime)
        print(otime, eid)
    """

    yy = UTCDateTime(otime).year
    mo = UTCDateTime(otime).month
    dd = UTCDateTime(otime).day
    hh = UTCDateTime(otime).hour
    mm = UTCDateTime(otime).minute
    ss = UTCDateTime(otime).second
    ms = UTCDateTime(otime).microsecond # mili?
    eid = '%04d%02d%02d%02d%02d%02d%03d' % (yy, mo, dd, hh, mm, ss, ms)
    return eid

def eid2otime(eid):
    """
    Convert event ID to origin time. 
    See also otime2eid.

    Example
        eid = "20090407201255000"
        otime = eid2otime(eid)
        print(eid, otime)
    """

    part1 = "%s" % eid[:-3]
    part2 = "%s" % eid[-3:]
    eid = part1 + '.' + part2
    otime = UTCDateTime(eid)
    return(otime)


