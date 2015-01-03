# -*- coding: utf-8 -*-
#
# ===============================================================================
# masai.core.process.baseline
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

from __future__ import print_function, division
import sys

import numpy as np
import nmrglue as ng
from scipy.interpolate import pchip
import pandas as pd
import argparse

from lcslib.process.smooth import smooth

__all__ = ["ab", "basecorr"]

DEBUG = False

# =======================================================================================================================
# ab
#=======================================================================================================================
def ab(source=None, options='', axis=-1):
    """
    Automatic baseline correction - spectra are modified inline if apply is true

    Parameters
    ----------
    source: a source instance
    basetype: string, optional, default: 'linear'
        see notes - available = linear, basf, ...
    window: float/int, optional, default is 0.05
        if float <1 then the corresponding percentage ot the axis size is taken as window
    nbzone: int, optional, default is 32
        Number of zones. We will divide the size of the last axis by this number
        to determine the number of points in each zone (nw)
    mult: int
        A multiplicator. determine the number of point for the database calculation (nw*mult<n base points)
    nstd: int, optional, default is 2 times the standard error
        Another multiplicator. Multiply the standard error to determine the region in which points are from the baseline
    polynom: bool, optional, default is True
        If True a polynom is computed for the base line, else an interpolation is achieved betwwen points.
    porder: int, default is 6
        Order of the polynom to fit on the baseline points
    return_pts: bool, optional, default is False
        If True, the points abscissa used to determine the baseline are returned
    apply: bool, optional, default is True
        if apply is False, the data are not modified only the baseline is returned
    return_pts: bool, optional, default is False
        if True, the baseline reference points are returned
    axis: optional, default is -1

    Notes
    -----
    #TODO: description of these algorithms
    * linear -
    * basf -

    """
    # options evaluation
    parser = argparse.ArgumentParser(description='BC processing.', usage="""
    ab [-h] [--mode {linear,poly, svd}] [--dryrun]
                           [--window WINDOW] [--step STEP] [--nbzone NBZONE]
                           [--mult MULT] [--order ORDER] [--verbose]
    """)
    # positional arguments
    parser.add_argument('--mode', '-mo', default='linear', choices=['linear', 'poly', 'svd'], help = "mode of correction")
    parser.add_argument('--dryrun', action='store_true', help='dry flag' )

    parser.add_argument('--window', '-wi', default=0.05, type=float, help='selected window for linear and svd bc')
    parser.add_argument('--step', '-st', default=5, type=int, help='step for svd bc')
    parser.add_argument('--nbzone', '-nz', default=32, type=int, help='number of zone for poly')
    parser.add_argument('--mult', '-mt', default=4, type=int, help='multiplicator of zone for poly')
    parser.add_argument('--order', '-or', default=5, type=int, help='polynom order for poly')

    parser.add_argument('--verbose', action='store_true', help='verbose flag' )
    args  = parser.parse_args(options.split())

    source.history.append('baseline correction mode:%s' % args.mode)

    if axis == -1:
        par = source.par
    else:
        par = source.par2

    if not par.isfreq:
        pass #todo: fid correction

    # we will work on the ndarray, not DataFrame
    data = source.data

    if axis == 0:
        # transpose temporarily the data for indirect dimension ft
        data = data.T
    index=data.index
    columns=data.columns

    base, w = basecorr(data.values, args, retw=True)
    par.baseline = pd.DataFrame(base, index=index, columns=columns)
    par.basepoints = w

    if not args.dryrun:
        data = data - base

     # un-transpose the data if needed
    if axis == 0:
        data = data.T
        par.baseline = par.baseline.T
    source.data = data


# =======================================================================================================================
# basecorr
#=======================================================================================================================
def basecorr(data, args=None, retw=False):
    """
    Correct baseline

    see `autobase' for parameters definitions

    """
    if args:
        mode = args.mode
    else:
        mode = 'linear'

    if mode == 'linear':
        return _linearbase(data, args, retw)

    if mode == 'svd':
        return _svdbase(data, args, retw)

    if mode == 'poly':
        return _polybase(data, args, retw)
    else:
        raise ValueError('ab mode not known')


#####################
# private functions #
#####################
#
# _linear mode
#
def _linearbase(data, args=None, retw=False):
    """
    Apply a linear baseline correction

    """
    if not args:
        window = 0.05
    else:
        window = args.window

    if window <= 1.0:
        # percent
        window = int(data.shape[-1] * window)

    # print(window)
    if len(data.shape) == 1:
        npts = float(data.shape[-1])
        a = (data[-window:].mean() - data[:window].mean()) / (npts - 1.)
        b = data[:window].mean()
        baseline = a * np.arange(npts) + b

    else:
        npts = float(data.shape[-1])
        a = (data[:, -window:].mean(axis=-1) - data[:, :window].mean(axis=-1)) / (npts - 1.)
        b = data[:, :window].mean(axis=-1)
        baseline = (((np.ones_like(data).T * a).T * np.arange(float(npts))).T + b).T

    if retw:
        return baseline, None  #TODO: return something

    return baseline


#=======================================================================================================================
def _planeFit(points):
    """
    p, n = planeFit(points)

    Fit an n-dimensional plane to the points.
    Return a point on the plane and the normal.
    
    Parameters
    ----------
    points : 
    
    Notes
    -----
        Replace the nonlinear optimization with an SVD. 
        The following creates the moment of inertia tensor, M, and then 
        SVD's it to get the normal to the plane. 
        This should be a close approximation to the least-squares fit 
        and be much faster and more predictable. 
        It returns the point-cloud center and the normal.

    """
    from numpy.linalg import svd
    npts = points.shape[0]
    points = np.reshape(points, (npts, -1))
    assert points.shape[0] < points.shape[1]
    ctr = points.mean(axis=1)
    x = points - ctr[:, None]
    M = np.dot(x, x.T)
    return ctr, svd(M)[0][:, -1]

#=======================================================================================================================
def _svdbase(data, args=None, retw=False):
    """
    Apply a planar baseline correction to 2D data
    
    Parameters
    ----------
    data : source instance
    
    """
    if not args:
        window = 0.05
        step = 5
    else:
        window = args.window
        step = args.step

    if window <= 1.0:
        # percent
        window = int(data.shape[-1] * window)

    data = pd.DataFrame(data)  # TODO: facilitate the manipulation (but to think about further)
    a = pd.concat([data.iloc[:window], data.iloc[-window:]])
    b = pd.concat([data.iloc[window:-window, :window], data.iloc[window:-window, -window:]], axis=1)
    bs = pd.concat([a, b])
    bs = bs.stack()
    bs.sort()
    x = []
    y = []
    z = []
    for item in bs.index[::step]:
        x.append(item[0])
        y.append(item[1])
        z.append(bs[item].real)

    norm = np.max(np.abs(z))
    z = np.array(z)
    z = z / norm
    XYZ = np.array((x, y, z))
    p, n = _planeFit(XYZ)
    d = np.dot(p, n)
    if DEBUG: print(" origin baseline plane: ", p)
    if DEBUG: print(" normal vector component:", n)

    col = data.columns
    row = data.index
    X, Y = np.meshgrid(col, row)
    Z = -norm * (n[0] * X + n[1] * Y - d) / n[2]

    if retw:
        return Z, None  #TODO: return something
    return Z


#=======================================================================================================================
def _polybase(data, args=None, retw=False):
    """
    Automatic baseline correction

    """
    if isinstance(data, pd.DataFrame):
        dat = data.values.copy()
    else:
        dat = data.copy()
    if dat.ndim == 1:
        dat = np.array([dat, ])

    if not args:
        nbzone = 32
        mult = 4
        order = 5
    else:
        nbzone = args.nbzone
        mult = args.mult
        order = args.order

    npts = dat.shape[-1]
    w = np.arange(npts)

    nw = npts / nbzone
    baseline = np.ma.masked_array(dat, mask=True)

    sigma = 1.e6
    nw = int(npts / nbzone)

    # print (nw)
    # unmask extremities of the baseline
    baseline[:, :nw].mask = False
    baseline[:, -nw:].mask = False

    for j in range(nbzone):
        s = dat[:, nw * j:min(nw * (j + 1), npts + 1)]
        sigma = min(s.std(), sigma)

    nw = nw * 2  # bigger window
    nw2 = int(nw / 2)

    found = False
    nb = 0
    nstd = 2.
    while (not found) or (nb < nw * mult):
        nb = 0
        for i in range(nw2, npts - nw2 + 1, 1):
            s1 = dat[:, max(i - 1 - nw2, 0):min(i - 1 + nw2, npts + 1)]
            s2 = dat[:, max(i - nw2, 0):min(i + nw2, npts + 1)]
            s3 = dat[:, max(i + 1 - nw2, 0):min(i + 1 + nw2, npts + 1)]
            mi1, mi2, mi3 = s1.min(), s2.min(), s3.min()
            ma1, ma2, ma3 = s1.max(), s2.max(), s3.max()
            #print s1.shape
            if abs(ma1 - mi1) < float(nstd) * sigma and abs(ma2 - mi2) < float(nstd) * sigma and abs(ma3 - mi3) < float(
                    nstd) * sigma:
                found = True
                nb += 1
                baseline[:1, i].mask = False  # baseline points

        # increase nstd
        nstd = nstd * 1.1
    if DEBUG: print('basf optimized nstd: %.2F mult: %.2f' % (nstd, mult))
    
    wm = np.array(zip(*np.argwhere(~baseline[:1].mask))[1])
    bm = baseline[:, wm]  # [~baseline[0].mask]
    bm = smooth(bm.T, window_len=max(int(dat.shape[0] / 10), 3)).T
    bm = smooth(bm, window_len=max(int(dat.shape[-1] / 10), 3))

    #if not polynom:
    #    sr = pchip(wm, bm.real)
    #    si = pchip(wm, bm.imag)
    #    baseline = sr(w) + si(w) * 1.0j
    #    baseline = smooth(baseline, window_len=int(nw / 4))
    #else:
        # fit a polynom
    pf = np.polyfit(wm, bm.T, order).T
    for i, row in enumerate(pf[:]):
        poly = np.poly1d(row)
        baseline[i] = poly(w)

    if retw:
        return baseline, wm

    return baseline


if __name__ == "__main__":
    pass


