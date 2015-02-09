# -*- coding: utf-8 -*-
#
# ===============================================================================
# masai.core.analysis.picking
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
#

import numpy as np
from nmrglue.analysis.peakpick import pick
import pandas as pd

#
# picking
#
def picking(source, thres=None, index=False):
    """

    Parameters
    ----------
    source: source instance
    thres: float, optional
        percent of threshold
    index:

    Returns
    -------
    """
    if isinstance(source, pd.DataFrame):
        data = source
    else:
        data = source._data

    if thres is None:
        thres = 80.

    pthres = np.abs(data.values).max() * thres / 100.

    ar = pick(np.abs(data.values), pthres, table=False)[0]

    x = source.data.columns
    y = source.data.index
    # if x[0]>x[-1]:
    #     x = x[::-1]
    # if y[0]>y[-1]:
    #     y = y[::-1]
    ars = []
    for i, item in enumerate(ar):
        ar1 = x[item[1]]
        ar0 = y[item[0]]
        ars.append([ar0, ar1])

    if not index:
        return ars
    else:
        #print(ars)
        #print(ar)
        return ar


