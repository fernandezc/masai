# -*- coding: utf-8 -*-
#
# ===============================================================================
# masai.core.process.mqstmas
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

import numpy as np

__all__ = ['k_ratio', 'mqmas_analysis', 'stmas_analysis']

# ===============================================================================
# k_ratio
# ===============================================================================
def k_ratio(source, p=3, q=0):
    """
    Calculation of the MQMAS or STMAS experiments shearing ratio
    see. for instance  Fernandez and Pruski, Probing Quadrupolar Nuclei
            by Solid-State NMR Spectroscopy: Recent Advances,
            Top Curr Chem (2011), DOI: 10.1007/128_2011_141

    Parameters
    ----------
    fd: an instance of NMRSource
    p: coherence order
    q: transition order

    """

    I = source.get_spin()

    n_spin = int(np.ceil(2.0 * I))
    p = float(np.abs(p))
    q = float(q)

    ratio = 0.
    if p == n_spin:
        p = -p  #
    if p > n_spin:
        print('coherence order p cannot be greater than 2*I')
        return 0
    if p > 2 or p == -n_spin:
        # MQMAS experiment assumed
        ratio = (p / 9.0) * (-17.0 * (p ** 2) + 36.0 * I * (I + 1.)
                             - 10.0) / (4.0 * I * (I + 1.0) - 3.0)
    elif p == 1 and q != 0:
        # STMAS experiment assumed
        ratio = -((17.0 / 3.0) * (q ** 2) / (4.0 * I * (I + 1.0) - 3.0) - 1.0)

    # compute ratio in terms of points
    # fres0 = 1./dic['acqu2s']['SW_h']
    # fres1 = 1./dic['acqus']['SW_h']
    # pntratio = ratio*fres0/fres1

    # Store this ratio
    source.par.k_ratio = ratio
    return ratio


# ==============================================================================
# mqmas_analysis
# ==============================================================================
def mqmas_analysis(source, rect=None):
    return quad2d_analysis(source, rect)


# ==============================================================================
# stmas_analysis
# ==============================================================================
def stmas_analysis(source, rect=None):
    return quad2d_analysis(source, rect)


# ==============================================================================
# quad2danalysis
# ==============================================================================
def quad2d_analysis(source, rect=None):
    """
    rect: list
        must be (x1, x2, y1, y2)
    """
    x = source.data.columns
    y = source.data.index
    dat = source.data.values.real

    # find gravity centers
    if rect is None:
        return "no rectangle defined"

    x1, x2, y1, y2 = rect
    x1, x2 = sorted([x1, x2])
    y1, y2 = sorted([y1, y2])

    xsel = np.argwhere((x > x1) & (x < x2))
    ysel = np.argwhere((y > y1) & (y < y2))
    xs = x[xsel]
    ys = y[ysel]
    dxp = dat.sum(axis=0)[xsel]
    dyp = dat.sum(axis=1)[ysel]
    cg2 = np.average(xs, weights=dxp)
    cgiso = np.average(ys, weights=dyp)

    # compute parameters
    cs = (17. * cgiso + 10. * cg2) / 27.
    qis = 17. * (cg2 - cgiso) / 27.
    sfo1 = source.fetchPar('SFO1')
    I = source.getSpin()
    pq = sfo1 * np.sqrt(-qis * 1.e-6 * 30. / (I * (I + 1.) - 3. / 4.)) / (3. / (2. * I * (2. * I - 1)))
    # Compute the quadripolar constant and chemical shift
    return cgiso, cg2, cs, qis, pq
