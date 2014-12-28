from __future__ import print_function

from traits.api import (HasTraits, Instance)

import numpy as np
import pandas as pd
import argparse

# import nmrglue as ng
# from nmrglue.process.proc_lp import lp as linear_predict
from nmrglue.process.proc_base import zf_size, largest_power_of_2

__all__ = ["em", "zf",]


#
# zf : zero-filling
#
def zf(source, options='', axis=-1):
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

    # options evaluation
    parser = argparse.ArgumentParser(description='ZF processing.',
                                     usage="""zf [nzf]""",)
    # positional arguments
    parser.add_argument('nzf', default=0,  nargs='?', type=int, help='number of zero-filling')
    args  = parser.parse_args(options.split())

    nzf = args.nzf
    si = largest_power_of_2(data.shape[-1]) * 2**nzf
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
def em(source, options='', axis=-1, **kargs):
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

    # options evaluation
    parser = argparse.ArgumentParser(description='EM processing.')#,
                                     ##usage="""[k_shifted_]em [lb] [--[k_]shifted [val]] [--verbose]""",)
    # positional arguments
    parser.add_argument('command', default='em', choices=['em', 'k_shifted_em'], help='broadening command')
    parser.add_argument('lb', default=10.,  nargs='?', type=float, help='broadening value in Hz')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--shifted', default=0.0, type=float, help = "position of the top in units of time")
    group.add_argument('--k_shifted', default=[3.0, 0.0], type=float, nargs=2, help = "p,q to compute the k-ratio for shifting")
    # optional arguments
    parser.add_argument('--verbose', action='store_true', help='verbose flag' )
    args  = parser.parse_args(options.split())

    lb = args.lb # float(kargs.get('lb', 10.))
    shifted = args.shifted #float(kargs.get('top', 0.0))
    k_shifted = args.k_shifted

    if args.command in ['k_shifted_em'] and (axis == -1 or axis == 1) and source.is_2d:
        # in k shifted method for 2D spectra, the top of
        # the broadening function follow the top of the echoes
        # parameters p and q should be defined in parameters

        p = k_shifted[0] # coherence
        q = k_shifted[1] # satellite order

        shifted = np.abs(k_ratio(source, p, q))
        if args.verbose:
            print('k_ratio: %.3f' % k_ratio)

    par.LB = lb

    kargs['lb'] = lb / source.get_multiplicator(axis)
    kargs['shifted'] = shifted
    kargs['states'] = True if 'STATES' in par.fnmode else False

    data = _apodization(data, args.command, **kargs)

    if axis == 0:
        # transpose temporarily the data for indirect dimension ft
        data = data.T

    source.data = data
    source.history.append('Exponential apodization lb:%.2f' % par.LB)

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
    # em
    #===============================================================================
    def em(self, **kargs):
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
        shifted = float(kargs.get('shifted'))

        x = self.data.columns.values

        return np.exp(-np.pi * np.abs(x - shifted) * lb)

    #===============================================================================
    # k_shifted_em
    #===============================================================================
    def k_shifted_em(self, **kargs):
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
        shifted = float(kargs.get('shifted'))
        states = kargs.get('states', True)

        x = data.columns.values
        y = data.index.values

        ratio = abs(shifted)
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
