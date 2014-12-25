from __future__ import print_function, division

# imports
from nmrglue.process.proc_base import zf_size, largest_power_of_2
from warnings import warn

#
# zf : zero-filling
#
def zf(source, axis=-1, **kargs):
    """
    Zero filling

    Parameters
    ----------
    source: source instance
        The source to be zero filled
    nzf: number of zero filling, optional
        if n is 0, then the data are zero filled if necessary to the power of two
        else the zero filling is double the number of times n
    axis: optional

    """
    if axis == -1:
        par = source.par
    else:
        par = source.par2

    if par.isfreq:
        warn ('This source is already transformed: zf skipped!')
        return

    # data on which we apply ft

    data = source.data.values
    if axis == 0:
        # transpose temporarily the data for indirect dimension ft
        data = data.T

    n = int(kargs.get('nzf', 0))
    si = largest_power_of_2(data.shape[-1]) * (n*2)
    data = zf_size(data, si)

    # store the new pars
    par.TD = data.shape[-1]
    par.SI = par.TD

    # untranspose the data if needed
    if axis == 0:
        data = data.T

    # change the current data with the transformed data
    source.data = data  # scale will be calculated here