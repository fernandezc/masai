# -*- coding: utf-8 -*-
#
# ===============================================================================
# masai.core.source.source
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

"""
Define basic source functionalities
"""

from __future__ import print_function

# ===============================================================================
# Imports
# ===============================================================================
import os
import glob
from warnings import warn, filterwarnings

filterwarnings("ignore")

import pandas as pd
import numpy as np
from traits.api import *

from masai.core.utils.misc import position2index

class Source(HasTraits):
    # minimal plugin attribute -----------------------------------------------------
    description = Str
    type = Str

    modeldata = Instance(pd.DataFrame)

    history = List

    # *******************************************************************************
    # public functions
    # *******************************************************************************

    # -------------------------------------------------------------------------------
    def row(self, index, **kargs):
        """
        Extract a row from a two dimensional spectra
        Parameters
        ----------
        index:
        width:

        Returns
        -------
        data: DataFrame
            return a pandas Dataframe with the row

        Examples
        --------

        """
        return self._getsection(index, axis=-1, **kargs)

    # -------------------------------------------------------------------------------
    def col(self, index, **kargs):
        """
        Extract a row from a two dimensional spectra
        Parameters
        ----------
        index:
        width:

        Returns
        -------
        data: DataFrame
            return a pandas Dataframe with the row

        Examples
        --------

        """

        return self._getsection(index, axis=0, **kargs)

    # ------------------------------------------------------------------------------
    def get_limits(self, axis=-1):
        data = self.data
        lim = ([data.index[0], data.index[-1]],
               [data.columns[0], data.columns[-1]])[axis]
        return lim

    # *******************************************************************************
    # Properties
    # *******************************************************************************

    # Data origin (as read by the __init__ function): Read only --------------------
    data_orig = Property
    _data_orig = Instance(pd.DataFrame)

    def _get_data_orig(self):
        return self._data_orig.copy()

    # data output ------------------------------------------------------------------
    data = Property(depends_on='data_orig')
    _data = Instance(pd.DataFrame)

    def _get_data(self):
        if self._data is None:
            return self._data_orig
        return self._data

    def _set_data(self, value):
        if isinstance(value, pd.DataFrame):
            self._data = value
        else:
            # noinspection PyArgumentList
            self._data = self._make_dataframe(value)

    # Real data -----------------------------------------------------------------------
    real = Property

    def _get_real(self):
        r = self.data.values.real
        x = self.data.columns
        y = self.data.index
        return pd.DataFrame(r, y, x)

    # Imag data -----------------------------------------------------------------------
    imag = Property

    def _get_imag(self):
        j = self.data.values.imag
        x = self.data.columns
        y = self.data.index
        return pd.DataFrame(j, y, x)

    # ndims
    ndims = Property
    def _get_ndims(self):
        ndims = 1
        if self.is_2d:
            ndims = 2
        elif self.is_3d:
            ndims = 3
        return ndims

    # is 2d -----------------------------------------------------------------------
    is_2d = Property

    def _get_is_2d(self):
        """
        property getter: Check if the data are 2D or not
        """
        if self._data is not None:
            return len(self._data.shape) == 2 and self._data.shape[-2] > 1
        else:
            print('No data')
            return False

    # is 3d ------------------------------------------------------------------------
    is_3d = Property

    def _get_is_3d(self):
        """
        property getter: Check if the data are 2D or not
        """
        if self.data is not None:
            return len(self.data.shape) == 3 and self.data.shape[-3] > 1
        else:
            print('No data')
            return False

    def position2index(self, pos, axis=-1):
        if axis == 0:
            a = self.data.index
        else:
            a = self.data.columns

        x0, x1 = a[0], a[-1]
        if x0 < x1:
            return a[a<=pos].size-1
        else:
            return a[a>=pos].size-1

    # *******************************************************************************
    # Private functions
    # *******************************************************************************

    # -------------------------------------------------------------------------------
    def _getsection(self, index, **kargs):
        """
        Extract a section by ppm (or Hz) position or by index

        """
        if not self.is_2d :
            print("1D spectrum, so we return the input data itself!")
            return self.data

        data = self.data.copy()

        axis = kargs.get('axis', -1)
        axis = range(data.ndim)[axis]

        #print(data.ndim-1)
        if axis == data.ndim-1:
            txt1 = 'taking horizontal (axis=%d) section'%axis
        else:
            data = data.T
            txt1 = 'taking vertical (axis=%d) section'%axis

        x = data.index
        y = data.columns

        byindex = kargs.get('byindex', False)
        width = kargs.get('width', 0)

        if not byindex:
            txt2 = ' by position'
            w2 = width / 2.0
            lb = position2index(x, index-w2)
            ub = position2index(x, index+w2)
            x0 = x[lb]
            x1 = x[ub]
            if width != 0:
                txt3 = ' from %.2f to %.2f'%(x0,x1)
            else:
                txt3 = ' at %.2f'%x0
            print(txt1+txt2+txt3)

        else:
            txt2 = ' by index'
            w2 = int(width / 2)
            lb = max(index - w2, 0)
            ub = min(lb + width, x.size-1)
            ub = ub - width
            if width != 0:
                txt3 = ' from index %d to %d'%(lb,ub)
            else:
                txt3 = ' at index %d'%lb
            print(txt1+txt2+txt3)

        if lb>ub:
            lb, ub = ub, lb
        section = data.iloc[lb:ub+1].sum(axis=0)

        assert section.size == y.size

        return pd.DataFrame(section, index=y).T

if __name__ == '__main__':
    """
    test plugin
    """
    # need this to test  directly
    import pylab as pl

    import warnings
    warnings.filterwarnings("ignore")

    Bruker = Bruker()

    DATADIR = os.path.expanduser(os.path.join('~','Dropbox','storage','bruker','data'))

    user = 'christian'
    name = 'example1'
    expno = 10
    path = os.path.join(DATADIR,user,'nmr',name, str(expno))
    fd = Bruker(path, dgfilter=True)
    fd.data.T.plot(legend=False)

    user = 'aurelie'
    name = 'RbNO3'
    expno = 15
    path = os.path.join(DATADIR,user,'nmr',name, str(expno))
    fd = Bruker(path, dgfilter=True)

    row = fd.row(0, byindex=False)
    row.T.plot(legend=False)

    col = fd.col(0, width=.1, byindex=False)
    col.T.plot(legend=False)

    pl.show()