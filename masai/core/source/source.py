# -*- coding: utf-8 -*-

"""
Define basic source functionalities

"""

# Copyright (c) C. Fernandez' group @ LCS.
# Distributed under the terms of the CeCILL-B License.

from __future__ import print_function

# -----------------------------------------------------------------------------
# Globals
# -----------------------------------------------------------------------------

# stdlib
import sys

# third-party
import pandas as pd
from traits.api import HasTraits, Property, Instance, List, Str, Either

# -----------------------------------------------------------------------------
# Local utilities and constants
# -----------------------------------------------------------------------------
from masai.core.utils.misc import position2index, is_sequence
# constants
DEBUG = True


#
# Source class definition
#
class Source(HasTraits):
    """
    Abstract Source class

    """

    description = Str
    type = Str

    data_orig = Property
    _data_orig = Instance(pd.DataFrame)

    data = Property(depends_on='data_orig')
    _data = Instance(pd.DataFrame)

    real = Property

    imag = Property

    dims = Property

    ndims = Property
    is_2d = Property

    is_3d = Property

    units = Property  # units of the elements
    _units = Either(Str(), None)
    xunits = Property
    _xunits = Either(Str(), None)
    yunits = Property
    _yunits = Either(Str(), None)

    modeldata = Instance(pd.DataFrame)

    history = List

    # generic attibutes to easy access to coordinates
    X = Property
    Y = Property # for 2D data
    #Z = Property # for 3D data (max for the moment)

    loc = Property

    # deprecated function
    def row(self, *args, **kargs):
        print("'row' is deprecated, use 'get_row' instead.", file=sys.stderr)
        return self.get_row(*args, **kargs)

    def get_row(self, value, byindex=False, width=0):
        """
        Extract a row from a two dimensional spectra

        Parameters
        ----------
        value: int or float
            index or position in the axis unit of the row to extract
            depends on byindex keyword
        byindex: boolean, optional, default=False
            if False, value is a position in axis unit,
            if True it is the index of the row
        width: float, optional, default=0
            if 0, only one row is taken
            else the specified width is taken into account, the returned data
            correspond to the integration either of a number of row
            or of a range of position depending on the byindex keyword

        Returns
        -------
        data: DataFrame
            return a pandas Dataframe with the row

        """
        return self._getsection(value, axis=-1, byindex=byindex, width=width)

    # deprecated function
    def col(self, *args, **kargs):
        print("'col' is deprecated, use 'get_col' instead.", file=sys.stderr)
        return self.get_col(*args, **kargs)

    def get_col(self, value, byindex=False, width=0):
        """
        Extract a column from a two dimensional spectra

        Parameters
        ----------
        value: int or float
            index or position in the axis unit of the row to extract
            depends on byindex keyword
        byindex: boolean, optional, default=False
            if False, value is a position in axis unit,
            if True it is the index of the column
        width: float, optional, default=0
            if 0, only one row is taken
            else the specified width is taken into account, the returned data
            correspond to the integration either of a number of row
            or of a range of position depending on the byindex keyword

        Returns
        -------
        data: DataFrame
            return a pandas Dataframe with the col selected
        """

        return self._getsection(value, axis=0, byindex=byindex, width=width)

    def get_limits(self, axis=-1):
        data = self.data
        lim = ([data.index[0], data.index[-1]],
               [data.columns[0], data.columns[-1]])[axis]
        return lim

    # -------------------------------------
    #def position2index(self, pos, axis=-1):
    #    if axis == 0:
    #        a = self.data.index
    #    else:
    #        a = self.data.columns
    #
    #    x0, x1 = a[0], a[-1]
    #    if x0 < x1:
    #        return a[a <= pos].size - 1
    #    else:
    #        return a[a >= pos].size - 1

    # *******************************************************************************
    # Properties
    # *******************************************************************************

    # Data origin (as read by the __init__ function): Read only --------------------
    def _get_data_orig(self):
        return self._data_orig.copy()

    # data output ------------------------------------------------------------------
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
    def _get_real(self):
        return self.get_real()

    def get_real(self):
        r = self.data.values.real
        x = self.data.columns
        y = self.data.index
        return pd.DataFrame(r, y, x)

    # Imag data -----------------------------------------------------------------------
    def _get_imag(self):
        return self.get_imag()

    def get_imag(self):
        j = self.data.values.imag
        x = self.data.columns
        y = self.data.index
        return pd.DataFrame(j, y, x)

    # ndims
    def _get_ndims(self):
        ndims = 1
        if self.is_2d:
            ndims = 2
        elif self.is_3d:
            ndims = 3
        return ndims

    # is 2d -----------------------------------------------------------------------
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
    def _get_is_3d(self):
        """
        property getter: Check if the data are 2D or not
        """
        if self.data is not None:
            return len(self.data.shape) == 3 and self.data.shape[-3] > 1
        else:
            print('No data')
            return False

    def _get_units(self):
        if not self._units:
            self._units = 'unitless'
        return self._units

    def _set_units(self, val):
        self._units = val

    def _get_units(self):
        if not self._xunits:
            self._xunits = 'a.u.'
        return self._xunits

    def _set_xunits(self, val):
        self._xunits = val

    def _get_yunits(self):
        if not self._yunits:
            self._yunits = 'a.u.'
        return self._yunits

    def _set_yunits(self, val):
        self._yunits = val

    def _get_X(self):
        return self.coords[self.dims[-1]]

    def _get_Y(self):
        if self.ndim>2:
            return self.coords[self.dims[-2]]

    #def _get_Z(self):
    #    if self.ndim==3:
    #        return self.coords[self.dims[-3]]

    def _get_loc(self):
        """Attribute for location based indexing like pandas.
        """
        return _LocIndexer(self)

    def _get_dims(self):
        dims= ['X','Y']
        if self.ndims == 1:
            dims = ['X']
        return dims

    # *******************************************************************************
    # Private functions
    # *******************************************************************************

    # -------------------------------------------------------------------------------
    def _getsection(self, index, **kargs):
        """
        Extract a section by ppm (or Hz) position or by index

        """
        if not self.is_2d:
            print("1D data, so we return the input data itself!")
            return self.data

        data = self.data.copy()

        axis = kargs.get('axis', -1)
        axis = range(data.ndim)[axis]

        # print(data.ndim-1)
        if axis == data.ndim - 1:
            txt1 = 'taking horizontal (axis=%d) section' % axis
        else:
            data = data.T
            txt1 = 'taking vertical (axis=%d) section' % axis

        x = data.index
        y = data.columns

        byindex = kargs.get('byindex', False)
        width = kargs.get('width', 0)

        if not byindex:
            txt2 = ' by position'
            w2 = width / 2.0
            lb = position2index(x, index - w2)
            ub = position2index(x, index + w2)
            x0 = x[lb]
            x1 = x[ub]
            if width != 0:
                txt3 = ' from %.2f to %.2f' % (x0, x1)
            else:
                txt3 = ' at %.2f' % x0
            print(txt1 + txt2 + txt3)

        else:
            txt2 = ' by index'
            w2 = int(width / 2)
            lb = max(index - w2, 0)
            ub = min(lb + width, x.size - 1)
            ub -= width
            if width != 0:
                txt3 = ' from index %d to %d' % (lb, ub)
            else:
                txt3 = ' at index %d' % lb
            print(txt1 + txt2 + txt3)

        if lb > ub:
            lb, ub = ub, lb
        section = data.iloc[lb:ub + 1].sum(axis=0)

        assert section.size == y.size

        return pd.DataFrame(section, index=y).T

    def __getitem__(self, key):
        """Access coordinates or slice of data."""
        if key in self.dims:
            return self.coords[key]
        elif isinstance(key, basestring):
            print("Error: the key '%s' doesn't exist" % key,
                  file=sys.stderr)
            return
        else:
            # probably a slicing
            d = self._get_slice_by_index(key)
            if d.size < 2:
                # scalar
                return d.data
            d.name = self.name + '%d' % newid.next()
            return d

    def __setitem__(self, key, value):
        """Add an array to this dataset."""
        if key in self.dims:
            self.coords.__setitem__(key, value)

    def __delitem__(self, key):
        """Remove a variable from this dataset."""
        pass

    def _get_slice_by_index(self, inputkey):

        keys = [slice(None, None, None) for axis in self.dims]

        if not is_sequence(inputkey):  # not isinstance(key, tuple) and not isinstance(key, list):
            # by default will use a tuple
            inputkey = (inputkey, )

        for i, item in enumerate(inputkey):
            if isinstance(item, slice):
                keys[i] = item
            else:
                keys[i] = slice(item, item + 1, None)

        # slice the data
        newdata = self.data[keys]

        # act on several coordinates
        newcoords = {}
        newdims= []
        for k in sorted(self.coords.keys()):
            c = self.coords[k].copy()
            newcoords[k] = c
        for i, item in enumerate(keys):
            c = self.coords[self.dims[i]].copy()
            newcoords[self.dims[i]].data = c.data[item]  # slice it

        #preparer to the squeezing (number of dimension should be in agreement)
        for i, n in enumerate(newdata.shape):
            if n>1:
                newdims.append(self.dims[i])
            else:
                del newcoords[self.dims[i]]

        return DataArray(newdata, dims=newdims, coords=newcoords, units=self.units, attrs=self.attrs)

    def __getattr__(self, item):
        """ The item attribute was not found look for it in attrs
        """
        res = self.attrs.get(item, None)
       # if not res:
       #     # look coordinates
       #     res = [getattr(self[k], item) for k in self.dims]
       #     if len(res) == 1:
       #         res = res[0]
        return res


#
# Location Indexer
#
class _LocIndexer(object):
    def __init__(self, data):

        self.data = data


    def loc2index(self, loc, axis=-1):

        # we determine the axe from the dim key
        if axis==-1 and self.data.is_2d:
            axe = self.data.data.columns
        else:
            axe=self.data.data.index

        # if axe is a quantity but not the loc,
        # we need to work on the magnitude to assuming the index are in the units
        x0, x1 = axe[0], axe[-1]
        if x0 < x1:
            return axe[axe <= loc].size
        else:
            return axe[axe >= loc].size - 1

    def select(self, key):
        """ Here we need to transform from a slice made with location arguments, to a positional slice that
        can be used by a numpy.array

        """
        # helpers

        def _get_slice(key, i):

            if not isinstance(key, slice):
                item = self.loc2index(key, i)
                newkey = slice(item, item + 1, None)

            else:
                start, stop, step = key.start, key.stop, key.step
                start = self.loc2index(start, i) if start is not None else start
                stop = self.loc2index(stop, i) if stop is not None else stop
                step = self.loc2index(step, i) if step is not None else step
                newkey = slice(start, stop, step)

            return newkey

        def _get_slice_by_label(key, slices):
            # we have a label giving the dimension

            for axis, dim in enumerate(self.data.dims):

                if dim == key.start:

                    if not is_sequence(key.stop):
                        newkey = key.stop

                    else:
                        s, e = key.stop  # get the element of the slice if
                        newkey = slice(s, e, None)
                    ns = _get_slice(newkey, axis)
                    # we append this slice to slices
                    slices[axis] = ns

            return slices

        # we need to have a list of slice for each argument or a single slice acting on the axis=0
        # the given key can be a single argument or a single slice:

        slices = [slice(None, None, None) for axis in self.data.dims]

        if not is_sequence(key):
            keys = [key, ]
        else:
            keys = key

        for axis, key in enumerate(keys):

            if not isinstance(key, slice):
                # probably something like arr[10.7]
                slices[axis] = _get_slice(key, axis)

            else:
                # it's a slice (the labelled dimension can be in start argument
                if not isinstance(key.start, basestring): #key.start not in self.data.dims:
                    # a classical slice
                    slices[axis] = _get_slice(key, axis)
                else:
                    if key.start not in self.data.dims:
                        print("Error: the key '%s' doesn't exist" % key.start,
                              file=sys.stderr)
                        slices = [slice(0, 0, None) for axis in self.data.dims]
                        return slices
                    # a labelled  dimension is specified
                    # we need to cycle on all dimension
                    slices = _get_slice_by_label(key, slices)

        return slices

    def __getitem__(self, key):
        newkey = self.select(key)
        if len(newkey)==1:
            return self.data.data.iloc[newkey[0]]
        else:
            return self.data.data.iloc[newkey[0],newkey[1]]

    def __setitem__(self, key, value):
        newkey = self.select(key)
        self.data.data.iloc[newkey] = value
        # TODO: check if shape has changed

if __name__ == '__main__':
    pass
