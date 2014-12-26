from __future__ import print_function

from traits.api import (HasTraits, Instance)

import numpy as np
import pandas as pd

# import nmrglue as ng
# from nmrglue.process.proc_lp import lp as linear_predict
from nmrglue.process.proc_base import zf_size, largest_power_of_2

__all__ = ["em", "zf",]


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
        print('This source is already transformed: zf skipped!')
        return

    # data on which we apply ft

    data = source.data.values
    if axis == 0:
        # transpose temporarily the data for indirect dimension ft
        data = data.T

    n = int(kargs.get('nzf', 0))
    si = largest_power_of_2(data.shape[-1]) * 2**n
    data = zf_size(data, si)

    # store the new pars
    par.TD = data.shape[-1]
    par.SI = par.TD

    # un-transpose the data if needed
    if axis == 0:
        data = data.T

    # change the current data with the transformed data
    source.data = data  # scale will be calculated here

    source.history.append('Zero-filling to si:%d' % par.SI)

#
# em
#
def em(source, axis=-1, **kargs):
    """
    apodisation function multiplication

    Parameters
    ----------
    source: an instance of Source
    lb: float, optional
        Broadening factor, default is to read from the source parameters
    top: float, optional
        Position of the max of the apodisation function, default is None.
    method: string, optional, default is None
        if fullecho, top parameters is used to determine the maximum
        if k_shifted, coherence_order and satellite_order parameters are required to determine the maximums
    axis: optional, default is -1

    """
    if axis == -1:
        par = source.par
    else:
        par = source.par2

    if par.isfreq:
        print('This source is already transformed: are you sure this is what you want?')

    data = source.data
    if axis == 0:
        # transpose temporarily the data for indirect dimension ft
        data = data.T

    lb = float(kargs.get('lb', 10.))
    top = float(kargs.get('top', 0.0))
    method = kargs.get('method', None)

    func = kargs.get('func', 'em_func')

    if method == 'fullecho' and (axis == -1 or axis == 1):

        if not top:
            # id the parameter does not exist none should be return
            print("Fullecho but 'top' is 0.")
            top = 0.0

    elif method == 'k_shifted' and (axis == -1 or axis == 1) and source.is_2d:
        # in k shifted method for 2D spectra, the top of
        # the broadening function follow the top of the echoes
        # parameters p and q should be defined in parameters

        p = source.par.coherence_order # coherence
        q = source.par.satellite_order # satellite order

        if p is None or q is None:
            # id the parameter does not exist none should be return
            print("k_shifted method requires keys: 'coherence_order' and 'satellite_order' in parameters")
            return # abandon broadening

        top = np.abs(k_ratio(source, p, q))
        func = "k_shifted_%s" % func

    else:
        method = 'conventional'
        top = 0.0

    par.LB = lb


    kargs['lb'] = lb / source.get_multiplicator(axis)
    kargs['top'] = top
    kargs['states'] = True if 'STATES' in par.fnmode else False

    data = _apodization(data, func, **kargs)

    if axis == 0:
        # transpose temporarily the data for indirect dimension ft
        data = data.T

    source.data = data

#
# apodization
#
def _apodization(data, func, **kargs):
    """
    Apply a given apodization to the source data

    """
    model = ApodModels(data = data)

    _func = getattr(model, func)

    data = (data.values * _func(**kargs))

    return data

#
# ApodModels
#
class ApodModels(HasTraits):
    """
    A class to handle apodisation models

    """

    data = Instance(pd.DataFrame)

    #===============================================================================
    # em_func
    #===============================================================================
    def em_func(self, **kargs):
        """
        Calculate an exponential apodization fonction

        Parameters
        ----------
        lb: float
            apodization value in spectral units
        top: float
            k shift ratio
        axis : optional, default is -1

        Returns
        -------
        f : ndarray
            the apodisation array

        """

        lb = float(kargs.get('lb'))
        top = float(kargs.get('top'))

        x = self.data.columns.values

        return np.exp(-np.pi * np.abs(x - top) * lb)

    #===============================================================================
    # k_shifted_em_func
    #===============================================================================
    def k_shifted_em_func(self, **kargs):
        """
        Calculate a shifted exponential apodization fonction

        Parameters
        ----------
        lb: float
            apodization value in spectral units
        top: float
            k shift ratio
        axis : optional, default is -1

        Returns
        -------
        f : ndarray
            the apodisation array

        """

        data = self.data

        if not (data.shape[0] > 1):
            raise Exception('only work for two dimensional data')

        lb = float(kargs.get('lb'))
        top = float(kargs.get('top'))
        states = kargs.get('states', True)

        x = data.columns.values
        y = data.index.values

        ratio = abs(top)
        # ratio must be converted in points ratio
        r = (y[1]-y[0])/(x[1]-x[0])
        ratio = ratio * r

        f = np.fromfunction(lambda y, x, lb, ratio: np.exp(-np.pi * np.abs(x - y*ratio) * lb),
                            data.shape, dtype=float, lb=lb, ratio=ratio)
        # If we are in states mode (general case two rows correspond
        # to the same time: we need to correct f for this
        # TODO: check this from parameters

        if states:
            f[1::2] = f[::2]

        return f
