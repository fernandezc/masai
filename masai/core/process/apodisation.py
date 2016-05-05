# -*- coding: utf-8 -*-
#
# ===============================================================================
# masai.core.process.apodisation
# ===============================================================================
# Copyright (C) 2015 Christian Fernandez
# Laboratoire Catalyse et Spectrochimie, Caen, France.
# christian.fernandez@ensicaen.fr
# This software is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# You can  use, modify and/ or redistribute the software under
# the terms of the CeCILL-B license as circulated by CEA, CNRS and INRIA
# at the following URL "http://www.cecill.info".
# See Licence.txt in the main masai source directory
# ===============================================================================

from __future__ import print_function

import argparse

from traits.api import (HasTraits, Instance)
import numpy as np
import pandas as pd

from masai.core.process import proc_base
from masai.core.utils.misc import position2index
from masai.core.process.mqstmas import k_ratio

__all__ = ["zf", "cs", "ls", "rs", "em", 'mc', 'pc', 'bc']

#
# bc : bc correction
#
def bc(source=None, options='', axis=-1):
    """
    bc correction

    Parameters
    ----------
    source: source instance
        The source to be zero filled
    tdeff:
    axis: optional

    """

    # options evaluation
    parser = argparse.ArgumentParser(description='BC correction.',
                                     usage="""
    bc [nav] [--tdeff] [-h]
    """, )
    # positional arguments
    parser.add_argument('nav', default=50, nargs='?', type=int, help='number of points for averaging')
    parser.add_argument('--tdeff', default=0, nargs='?', type=int, help='td effective')
    args = parser.parse_args(options.split())

    if axis == -1:
        par = source.par
    else:
        par = source.par2

    if par.isfreq:
        print('This source is already transformed: bc skipped!')
        return

    # data on which we apply zf

    datv = source.data.values
    if axis == 0:
        # transpose temporarily the data for indirect dimension ft
        datv = datv.T

    # read parameters TDeff
    tdeff = args.tdeff
    if tdeff > par.TD or tdeff < 4:
        tdeff = par.TD

    # retain only TDeff points only (if TDeff if at least 4)
    if tdeff > 4:
        tdeff = min(tdeff, datv.shape[-1])
        if datv.ndim > 1:
            datv = datv[:, :tdeff]
        else:
            datv = datv[:tdeff]

    # apply bc
    nav = args.nav
    datv = datv - datv[:,-nav:].mean()

    # store the new pars
    par.TD = datv.shape[-1]

    # un-transpose the data if needed
    if axis == 0:
        datv = datv.T

    # change the current data with the transformed data
    source.data = datv  # scale will be calculated here

    source.history.append('bc correction')


#
# pc : power calculation
#
def pc(source=None, options='', axis=-1):
    """
    Power calculation

    Parameters
    ----------
    source: source instance
        The source to be zero filled
    axis: optional

    """

    # options evaluation
    parser = argparse.ArgumentParser(description='PC processing.',
                                     usage="""
    pc npc [-h]
    """, )
    # positional arguments
    parser.add_argument('npc', default=2, nargs='?', type=int, help='power (default square)')
    args = parser.parse_args(options.split())

    if axis == -1:
        par = source.par
    else:
        par = source.par2

    if not par.isfreq:
        print('This source is not yet transformed: pc skipped!')
        return

    # data on which we apply pc

    datv = source.data.values
    if axis == 0:
        # transpose temporarily the data for indirect dimension ft
        datv = datv.T

    npc = args.npc
    datv = np.abs(datv) ** npc

    # un-transpose the data if needed
    if axis == 0:
        datv = datv.T

    # change the current data with the transformed data
    source.data = datv  # scale will be calculated here

    if npc != 1:
        source.history.append('Power spectrum calculation')
    else:
        source.history.append('Magnitude spectrum calculation')


#
# mc : magnitude calculation
#
def mc(source=None, options='', axis=-1):
    pc(source, options='1', axis=axis)


#
# zf : zero-filling
#
def zf(source=None, options='', axis=-1):
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

    # options evaluation
    parser = argparse.ArgumentParser(description='ZF processing.',
                                     usage="""
    zf [nzf] [-h]
    """, )
    # positional arguments
    parser.add_argument('nzf', default=0, nargs='?', type=int, help='number of zero-filling')
    args = parser.parse_args(options.split())

    if axis == -1:
        par = source.par
    else:
        par = source.par2

    if par.isfreq:
        print('This source is already transformed: zf skipped!')
        return

    # data on which we apply zf

    datv = source.data.values
    if axis == 0:
        # transpose temporarily the data for indirect dimension ft
        datv = datv.T

    nzf = args.nzf
    si = proc_base.largest_power_of_2(datv.shape[-1]) * 2 ** nzf
    datv = proc_base.zf_size(datv, si)

    # store the new pars
    par.TD = datv.shape[-1]
    par.SI = par.TD

    # un-transpose the data if needed
    if axis == 0:
        datv = datv.T

    # change the current data with the transformed data
    source.data = datv  # scale will be calculated here

    source.history.append('Zero-filling to si:%d' % par.SI)


#
# cs : circular shift
#
def cs(source=None, options='', axis=-1):
    """
    Circular shift

    Parameters
    ----------
    source: source instance
        The source to be zero filled
    ncs: circular shift in point, optional
        Positive for right shift
        Negative for left shift
    tcs: circular shift in unit of time
        Positive for right shift
        Negative for left shift
    axis: optional

    """

    # options evaluation
    parser = argparse.ArgumentParser(description='CS shift.',
                                     usage="""
    cs [ncs|tcs] [-h]
    """, )
    # positional arguments
    parser.add_argument('ncs', default=0, type=float, help='circular shit in points or time units')
    parser.add_argument('--time', action='store_true', help='is time')
    args = parser.parse_args(options.split())

    if axis == -1:
        par = source.par
    else:
        par = source.par2

    if par.isfreq:
        print('This source is already transformed: cs skipped!')
        return

    # data on which we apply cs

    data = source.data
    if axis == 0:
        # transpose temporarily the data for indirect dimension ft
        data = data.T

    ncs = args.ncs
    if args.time:
        sign = ncs / abs(ncs)
        ncs = sign * position2index(data, abs(ncs))
    datv = data.values
    datv = proc_base.cs(datv, int(ncs))
    #if ncs<0:
    #    datv[:, ncs:] = -np.conjugate(datv[:, ncs:])
    # store the new pars
    par.TD = datv.shape[-1]
    par.SI = par.TD

    # un-transpose the data if needed
    if axis == 0:
        datv = datv.T

    # change the current data with the transformed data
    source.data = datv  # scale will be calculated here

    source.history.append('circular shift by: %d' % ncs)


#
# ls : left shift
#
def ls(source=None, options='', axis=-1):
    """
    Left shift

    Parameters
    ----------
    source: source instance
        The source to be zero filled
    nls: int
        left shift in point, optional
    tls: float
        left shift in unit of time
    axis: optional

    """

    # options evaluation
    parser = argparse.ArgumentParser(description='LS shift.',
                                     usage="""
    ls [nls|tls] [-h]
    """, )
    # positional arguments
    parser.add_argument('nls', default=0, type=float, help='left shift in points or time units')
    parser.add_argument('--time', action='store_true', help='is time')
    args = parser.parse_args(options.split())

    if axis == -1:
        par = source.par
    else:
        par = source.par2

    if par.isfreq:
        print('This source is already transformed: ls skipped!')
        return

    # data on which we apply cs

    data = source.data
    if axis == 0:
        # transpose temporarily the data for indirect dimension ft
        data = data.T

    nls = args.nls
    if args.time:
        nls = position2index(data, nls)
    datv = data.values
    datv = proc_base.ls(datv, int(nls))

    # store the new pars
    par.TD = datv.shape[-1]
    par.SI = par.TD

    # un-transpose the data if needed
    if axis == 0:
        datv = datv.T

    # change the current data with the transformed data
    source.data = datv  # scale will be calculated here

    source.history.append('left shift by: %d' % nls)


#
# rs : circular shift
#
def rs(source=None, options='', axis=-1):
    """
    Right shift

    Parameters
    ----------
    source: source instance
        The source to be zero filled
    ncs: int
        right shift in point, optional
    tcs: float
        right shift in unit of time
    axis: optional

    """

    # options evaluation
    parser = argparse.ArgumentParser(description='RS shift.',
                                     usage="""
    rs [nrs|trs] [-h]
    """, )
    # positional arguments
    parser.add_argument('nrs', default=0, type=float, help='right shift in points or time units')
    parser.add_argument('--time', action='store_true', help='is time')
    args = parser.parse_args(options.split())

    if axis == -1:
        par = source.par
    else:
        par = source.par2

    if par.isfreq:
        print('This source is already transformed: cs skipped!')
        return

    # data on which we apply rs

    data = source.data
    if axis == 0:
        # transpose temporarily the data for indirect dimension ft
        data = data.T

    nrs = args.nrs
    if args.time:
        nrs = position2index(data, nrs)
    datv = data.values
    datv = proc_base.rs(datv, int(nrs))

    # store the new pars
    par.TD = datv.shape[-1]
    par.SI = par.TD

    # un-transpose the data if needed
    if axis == 0:
        datv = datv.T

    # change the current data with the transformed data
    source.data = datv  # scale will be calculated here

    source.history.append('right shift by: %d' % nrs)


#
# em
#
def em(source=None, options='', axis=-1, **kargs):
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
    # options evaluation
    parser = argparse.ArgumentParser(description='EM processing.',
                                     usage="""
    em [lb] [-h]
       [--shifted SHIFTED | --k_shifted P Q] [--verbose]
    """, )
    # positional arguments
    parser.add_argument('lb', default=10., nargs='?', type=float, help='broadening value in Hz')
    parser.add_argument('--shifted', default=0.0, type=float, help="position of the top in units of time")
    parser.add_argument('--k_shifted', default=None, type=int, nargs=2, help="p and q parameters")
    # optional arguments
    parser.add_argument('--verbose', action='store_true', help='verbose flag')
    args = parser.parse_args(options.split())

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

    command = 'em'

    lb = args.lb  # float(kargs.get('lb', 10.))
    shifted = args.shifted  #float(kargs.get('top', 0.0))
    k_shifted = args.k_shifted

    if k_shifted and (axis == -1 or axis == 1) and source.is_2d:
        # in k shifted method for 2D spectra, the top of
        # the broadening function follow the top of the echoes
        # parameters p and q should be defined in parameters
        command = 'k_shifted_em'

        p = k_shifted[0]  # coherence
        q = k_shifted[1]  # satellite order

        if p == 1 and q == 0:
            ratio = 1.0
        else:
            ratio = k_ratio(source, p, q)
        shifted = np.abs(ratio)
        if args.verbose:
            print('k_ratio: %.3f' % ratio)

    par.LB = lb

    kargs['lb'] = lb / source.get_multiplicator(axis)
    kargs['shifted'] = shifted
    kargs['states'] = True if 'STATES' in par.encoding else False

    datv = _apodization(data, command, **kargs)

    if axis == 0:
        # transpose temporarily the data for indirect dimension ft
        datv = datv.T

    source.data = datv
    source.history.append('Exponential apodization lb:%.2f' % par.LB)

#
# apodization
#
def _apodization(data, func, **kargs):
    """
    Apply a given apodization to the source data

    """
    model = ApodModels(data=data)

    _func = getattr(model, func)

    datv = (data.values * _func(**kargs))

    return datv


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
        states = kargs.get('states', True)

        x = self.data.columns.values

        f = np.exp(-np.pi * np.abs(x - shifted) * lb)
        if states:
            f[1::2] = f[::2]
        return f

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
        r = (y[1] - y[0]) / (x[1] - x[0])
        ratio = ratio * r

        f = np.fromfunction(lambda y, x, lb, ratio: (np.exp(-np.pi * np.abs(x - y * ratio) * lb)
                                                     + np.exp(-np.pi * np.abs(x + y * ratio) * lb)) / 2.,
                            data.shape, dtype=float, lb=lb, ratio=ratio)

        # If we are in states mode (general case two rows correspond
        # to the same time: we need to correct f for this
        # TODO: check this from parameters

        if states:
            f[1::2] = f[::2]

        return f
