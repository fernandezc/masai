# -*- coding: utf-8 -*-
#
#===============================================================================
# masai.core.utils.misc
#===============================================================================
# Copyright (C) 2012-2013 Christian Fernandez
#       Laboratoire Catalyse et Spectrochimie, Caen, France.
#       christian.fernandez@ensicaen.fr
# This software is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# You can  use, modify and/ or redistribute the software under
# the terms of the CeCILL-B license as circulated by CEA, CNRS and INRIA
# at the following URL "http://www.cecill.info".
# See Licence.txt in the main masai source directory
#===============================================================================
""".. _misc:

Various methods used in other part of the program

"""
import pandas as pd

def position2index(data, pos):
    """
    Return the index of a position along the last axis (columns) of a DataFrame

    Parameters
    ----------
        data: DataFrame
        pos: float
            a column position in the column units

    Returns
    -------
        index: int
            the index of the position

    """
    if isinstance(data, pd.DataFrame):
        a = data.columns
    else:
        a = data
    if a[0] < a[-1]:
        id = a[a<=pos].size
    else:
        id = a[a>=pos].size

    # check closest
    dpos0 = 1.e15
    idx = 0
    id0 = id
    for i in [-1,0,+1]:
        id = max(min(id0+i, a.size-1),0)
        dpos = abs(pos - a[id])
        if dpos<dpos0:
            idx=id
            dpos0 = dpos
    return idx
