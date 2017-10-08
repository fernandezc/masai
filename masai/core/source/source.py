# -*- coding: utf-8 -*-

"""
Define basic source functionalities
This is a container for the spectroscopic data or others,
and the associated parameters

"""

# Copyright (c) C. Fernandez' group @ LCS.
# Distributed under the terms of the CeCILL-B License.

from __future__ import print_function

# -----------------------------------------------------------------------------
# Globals
# -----------------------------------------------------------------------------

# stdlib
import sys
import copy

# third-party
import pandas as pd
from traits.api import HasTraits, Property, Instance, Array, \
    List, Str, Either, Bool

# -----------------------------------------------------------------------------
# Local utilities and constants
# -----------------------------------------------------------------------------
from masai.core.utils.misc import position2index, is_sequence

# constants
DEBUG = True

# =============================================================================
# Source class definition
# =============================================================================
class Source(HasTraits):
    """
    Abstract Source class

    """
    description = Property

    data_orig = Property
    _data_orig = Instance(pd.DataFrame)

    data = Property(depends_on='data_orig')
    _data = Instance(pd.DataFrame)

    real = Property

    imag = Property

    dims = Property
    _dims = List

    ndims = Property

    is_1d = Property

    is_2d = Property

    is_3d = Property

    xunits = Property
    _xunits = Either(Str(), None)
    yunits = Property
    _yunits = Either(Str(), None)
    zunits = Property
    _zunits = Either(Str(), None)
    units = Property
    _units = Either(Str(), None)

    modeldata = Instance(pd.DataFrame)

    history = List

    silent = Bool(True)

    # generic attibutes to easy access to coordinates
    X = Property
    Y = Property  # for 2D data
    Z = Property  # for 3D data (max for the moment)

    loc = Property

    sym = Property

    # *************************************************************************
    # public functions
    # *************************************************************************

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

    def get_projection(self, proj='sky', axis=1):

        data = self.data
        if proj == 'sky':
            d = data.max(axis=axis)
        if proj == 'sum':
            d = data.sum(axis=axis)
        return d

    def get_limits(self, axis=-1):
        """
        Get the data limits

        Parameters
        ----------
        axis: Int, default=-1
            The axis where to get the limits

        Returns
        -------
        limits: Tuple
            A tuple with min and max limits
        """
        data = self.data
        lim = ([data.index[0], data.index[-1]],
               [data.columns[0], data.columns[-1]])[axis]
        return lim

    def get_index(self, pos, axis=-1):
        """
        Get the index of data from a data position

        Parameters
        ----------
        pos: Float,
            The position in data units
        axis: Int, default=-1
            The axis where to select the position

        Returns
        -------
            index: Int,
                The index corresponding to the position
        """
        if axis == 0:
            a = self.data.T
        else:
            a = self.data
        return position2index(a, pos)

    def get_imag(self):
        """
        Get imaginary data

        Returns
        -------
            imag: DataFrame
                The imaginary part of the data

        Notes
        -----
        Can be advantageously replaced by the corresponding property: imag

        """
        j = self.data.values.imag
        x = self.data.columns
        y = self.data.index
        return pd.DataFrame(j, y, x)

    def get_real(self):
        """
        Get real data

        Returns
        -------
            real: DataFrame
                The imaginary part of the data

        Notes
        -----
        Can be advantageously replaced by the corresponding property: real

        """
        r = self.data.values.real
        x = self.data.columns
        y = self.data.index
        return pd.DataFrame(r, y, x)

    # *************************************************************************
    # Properties getter and setter
    # *************************************************************************

    # data
    def _get_data(self):
        """
        property getter: read data
        """
        if self._data is None:
            return self._data_orig
        return self._data

    def _set_data(self, value):
        """
        property setter: assign data
        """
        if isinstance(value, pd.DataFrame):
            self._data = value
        else:
            # noinspection PyArgumentList
            self._data = self._make_dataframe(value)

    # xunits
    def _get_xunits(self):
        """
        property getter: read units for axe X
        """
        if not self._xunits:
            self._xunits = 'a.u.'
        return self._xunits

    def _set_xunits(self, val):
        """
        property getter: set units for axe X
        """
        self._xunits = val

    # yunits
    def _get_yunits(self):
        """
        property getter: read units for axe Y
        """
        if self.is_1d:
            return None

        if not self._yunits:
            self._yunits = 'a.u.'
        return self._yunits

    # zunits
    def _get_zunits(self):
        """
        property getter: read units for axe Z
        """
        if not self.is_3d:
            return None

        if not self._zunits:
            self._zunits = 'a.u.'
        return self._zunits

    # units
    def _set_units(self, val):
        """
        property getter: set units for the data
        """
        self._units = val

    # units
    def _get_units(self):
        """
        property getter: read units of the data
        """
        if not self._units:
            self._units = 'a.u.'
        return self._units

    # dims
    def _get_dims(self):
        """
        Property getter: Return the name of dimensions
        """
        if len(self._dims)<1:
            self._dims = ['Z', 'Y', 'X']

        if self.is_1d:
            return self._dims[-1]

        elif self.is_2d:
            return self._dims[-2:]

        elif self.is_3d:
            return self._dims[-3:]

    def _set_dims(self, val):
        """
        Property setter: Set the name of dimensions
        """
        if self.is_1d:
            self._dims = val

        elif self.is_2d:
            self._dims = val[-2:]

        elif self.is_3d:
            self._dims = val[-3:]

    # Read-only properties
    # ---------------------------------------------------------------------------

    # data_orig
    def _get_data_orig(self):
        """
        property getter: read the original data
        """
        return self._data_orig.copy()

    # real
    def _get_real(self):
        """
        property getter: return the imaginary part of the original data
        """
        return self.get_real()

    # imag
    def _get_imag(self):
        """
        property getter: return the imaginary part of the original data
        """
        return self.get_imag()

    def _get_sym(self):
        w, h = self.data.shape
        if w != h:
            return "cannot symmetrize: width # height"

        r = self.data.values
        r = r + r.copy().T
        x = self.data.columns
        y = self.data.index
        newdata = pd.DataFrame(r, y, x)
        # create a new source with the new data
        newsource = self.copy()
        newsource._data_orig = newdata
        newsource.data = newdata
        return newsource

    # ndims
    def _get_ndims(self):
        """
        property getter: get the number of dimensionw
        """
        ndims = 1
        if self.is_2d:
            ndims = 2
        elif self.is_3d:
            ndims = 3
        return ndims

    # is_1d
    def _get_is_1d(self):
        """
        property getter: Check if the data are 1D or not
        """
        data = self.data
        if data is not None:
            return len(data.shape) == 2 and data.shape[-2] == 1
        else:
            print('No data')
            return False

    # is_2d
    def _get_is_2d(self):
        """
        property getter: Check if the data are 2D or not
        """
        data = self.data
        if data is not None:
            return len(data.shape) == 2 and data.shape[-2] > 1
        else:
            print('No data')
            return False

    # is_3d
    def _get_is_3d(self):
        """
        property getter: Check if the data are 2D or not
        """
        data = self.data
        if data is not None:
            return len(data.shape) == 3 and data.shape[-3] > 1
        else:
            print('No data')
            return False

    # X
    def _get_X(self):
        """
        property getter: read coords for axe X
        """
        return self.data.columns

    # Y
    def _get_Y(self):
        """
        property getter: read coords for axe Y
        """
        return self.data.index

    # Z
    def _get_Z(self):
        pass
        # not implemented

    # loc
    def _get_loc(self):
        """
        Property getter: for location based indexing like pandas.
        """
        return _LocIndexer(self)

    # *************************************************************************
    # Private functions
    # *************************************************************************

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

        if not self.silent:
            print(txt1 + txt2 + txt3)

        if lb > ub:
            lb, ub = ub, lb
        section = data.iloc[lb:ub + 1].sum(axis=0)

        assert section.size == y.size

        return pd.DataFrame(section, index=y).T

    def copy(self, deepcopy=True):

        if not deepcopy:
            return copy.copy(self)  # shallow copy
        else:
            return copy.deepcopy(self)  # deep (recursive) copy

    def __getitem__(self, key):
        """
        Access slice of data.
        """
        d = self._get_slice_by_index(key)
        return d

    def __setitem__(self, key, value):
        """Add an array to this dataset."""
        raise NotImplementedError
        # TODO: setitem

    def __delitem__(self, key):
        """Remove a variable from this dataset."""
        raise NotImplementedError
        # TODO: delitem

    def _get_slice_by_index(self, inputkey):

        keys = [slice(None, None, None) for axis in self.dims]

        if not is_sequence(
                inputkey):  # not isinstance(key, tuple) and not isinstance(
            # key, list):
            # by default will use a tuple
            inputkey = (inputkey,)

        for i, item in enumerate(inputkey):
            if isinstance(item, slice):
                keys[i] = item
            else:
                keys[i] = slice(item, item + 1, None)

        # need a correction
        # remember that 1D data are also represented by a 2D array with only
        # one row!
        if self.ndims == 1:
            keys = slice(0, 1, None), keys[0]

        # slice the data

        newdata = self.data.iloc[keys]

        # create a new source with the new data
        newsource = self.copy()

        newsource._data_orig = newdata
        newsource.data = newdata

        return newsource

    def __getattr__(self, item):
        """ The item attribute was not found
        (should be subclassed in derived classes)
        """
        raise NotImplementedError

    # events
    def _anytrait_changed(self, name, old, new):

        if 'history' in name and not self.silent:
            print(self.history[-1])


#
# Location Indexer
#
class _LocIndexer(object):
    def __init__(self, data):

        self.data = data

    def loc2index(self, loc, axis=-1):

        # we determine the axe from the dim key
        if axis == -1 and self.data.is_2d:
            axe = self.data.data.columns
        else:
            axe = self.data.data.index

        # if axe is a quantity but not the loc,
        # we need to work on the magnitude to assuming the index are in the
        # units
        x0, x1 = axe[0], axe[-1]
        if x0 < x1:
            return axe[axe <= loc].size
        else:
            return axe[axe >= loc].size - 1

    def select(self, key):
        """ Here we need to transform from a slice made with location
        arguments, to a positional slice that
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

        # we need to have a list of slice for each argument or a single slice
        #  acting on the axis=0
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
                if not isinstance(key.start, basestring):  # key.start not in
                    #  self.data.dims:
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
        if len(newkey) == 1:
            return self.data.data.iloc[newkey[0]]
        else:
            return self.data.data.iloc[newkey[0], newkey[1]]

    def __setitem__(self, key, value):
        newkey = self.select(key)
        self.data.data.iloc[newkey] = value
        # TODO: check if shape has changed

if __name__ == '__main__':
    pass
