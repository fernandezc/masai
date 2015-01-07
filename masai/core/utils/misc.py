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
    a = data.columns
    x0, x1 = a[0], a[-1]
    if x0 < x1:
        return a[a<=pos].size-1
    else:
        return a[a>=pos].size-1
