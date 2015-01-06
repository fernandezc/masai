# -*- coding: utf-8 -*-
#
# ===============================================================================
# masai.core.process.sampling
# ===============================================================================
# Copyright (C) 2015 Christian Fernandez
#       Laboratoire Catalyse et Spectrochimie, Caen, France.
#       christian.fernandez@ensicaen.fr
# This software is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# You can  use, modify and/ or redistribute the software under
# the terms of the CeCILL-B license as circulated by CEA, CNRS and INRIA
# at the following URL "http://www.cecill.info".
# See Licence.txt in the main masai source directory
# ===============================================================================
import numpy as np
import argparse

__all__ = ["us", "cl" ]
#
# undersampling
#
def us(source=None, options=''):
    """
    undersampling along axis -1

    Parameters
    ----------
    n: integer, optional
        number of undersampling, default 1

    """
    # options evaluation
    parser = argparse.ArgumentParser(description='US undersampling.', usage = """
    us [-h] [nbtimes]
    """)
    # positional arguments
    parser.add_argument('nbtimes', default=1,  type=int, help='number of undersampling')
    args  = parser.parse_args(options.split())
    n = args.nbtimes

    par = source.par

    if par.isfreq:
        print('This source is already transformed: undersampling skipped!')
        return

    # data on which we apply undersampling
    data = source.data.values
    nr, nc = data.shape
    nc0 = n*int(nc/n) # we need to have nc divide by 2

    source.par.TD = nc0
    source.par.SW_h /= float(n)

    # change the current data with the transformed data
    source.data = data  # scale will be calculated here

    source.history.append('Undersampling: %d' % n)


def undersampling(source, N):
    """
    Under sampling the F2 dimension
    """
    #TODO: needs to understand why apk is not perfomring well with undersampling?
    if not source.getIsTransformed():
        source.history.append('undersampling -1')
        nr, nc = source.data.shape
        nc0 = N*int(nc/N) # we need to have nc divide by 2

        source.storePar('TD',nc0)
        source.storePar('SW_h', source.fetchPar('SW_h')/float(N))

        source.data = source.data.iloc[:, :nc0:N]

def cl(source=None, options='', axis=-1):
    """
    clipping along axis -1

    Parameters
    ----------
    limits : integer, optional
        number of undersampling, default 1

    """
    # options evaluation
    parser = argparse.ArgumentParser(description='US undersampling.', usage = """
    us [-h] [limits]
    """)

    # positional arguments
    parser.add_argument('limits', default=None, nargs=2,  type=int, help='clipping limits')
    args = parser.parse_args(options.split())

    c0, c1 = args.limits

    par = source.par

    if not par.isfreq:
        print('This source is not transformed: clipping skipped!')
        return

    # data on which we apply clipping

    if source.is_2d and axis == 0:
        data = source.data.T.copy()
        par = source.par2
    else:
        data = source.data.copy()
        axis = -1
        par = source.par

    # clip the data
    x = data.columns
    select = x[(x>=c0) & (x<=c1)]
    if select.size % 2 != 0: select = select[:-1]
    S = data.ix[:,select]

    # recompute SI, SW and new apparent SF
    SI0 = data.shape[-1]
    SI = S.shape[-1]
    par.SW_h *= float(SI) / float(SI0)
    SW_h = par.SW_h
    par.SI = SI
    BF1 = par.BF1
    SFO1 = par.SFO1
    delta = SW_h / (float(SI) * BF1)
    SF = SFO1 - (select[-1] + delta * float(SI) / 2.) * BF1 * 1.e-6
    par.SF= SF

    # change the current data with the transformed data
    if source.is_2d and axis == 0:
        source.data = S.T
    else:
        source.data = S
    source.history.append('clipping %d' % axis)

