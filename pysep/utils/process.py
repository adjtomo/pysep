from scipy import signal


def zerophase_chebychev_lowpass_filter(tr, freqmax):
    """
    Custom Chebychev type two zerophase lowpass filter useful for decimation
    filtering. This filter is stable up to a reduction in frequency with a
    factor of 10. If more reduction is desired, simply decimate in steps.
    Partly based on a filter in ObsPy.

    Filter parameters:
        rp: maximum ripple of passband
        rs: attenuation of stopband
        ws: stop band frequency
        wp: pass band frequency

    .. note::
        Should be replaced once ObsPy has a proper decimation filter.

    .. note::
        This function affects the trace in place!

    :type trace: obspy.core.trace.Trace
    :param trace: The trace to be filtered.
    :type freqmax: float
    :param freqmax: The desired lowpass frequency.
    """
    rp, rs, order = 1, 96, 1e99
    ws = freqmax / (tr_out.stats.sampling_rate * 0.5)
    wp = ws
    order, wn = None, None

    while True:
        if order <= 12:
            break
        wp *= 0.99
        order, wn = signal.cheb2ord(wp=wp, ws=ws, gpass=rp, gstop=rs,
                                    analog=False)

    b, a = signal.cheby2(N=order, rs=rs, Wn=wn, btype="low", analog=0,
                         output="ba")

    # Apply twice to get rid of the phase distortion.
    tr.data = signal.filtfilt(b, a, tr.data)
