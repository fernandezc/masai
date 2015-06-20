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
import functools
import operator

# third-party
import numpy as np
import pandas as pd
from pandas.util.testing import assert_frame_equal
from traits.api import HasTraits, Property, Instance, List, Str, Either

# -----------------------------------------------------------------------------
# Local utilities and constants
# -----------------------------------------------------------------------------
from masai.core.utils.misc import is_sequence, loc2index

# constants
DEBUG = True

#
# UFUNCS
#
_UNARY_MATH =[
    'negative', #(x[, out])	Numerical negative, element-wise.
    'absolute', #(x[, out])	Calculate the absolute value element-wise.
    'rint', #(x[, out])	Round elements of the array to the nearest integer.
    'sign', #(x[, out])	Returns an element-wise indication of the sign of a number.
    'conj', #(x[, out])	Return the complex conjugate, element-wise.
    'exp', #(x[, out])	Calculate the exponential of all elements in the input array.
    'exp2', #(x[, out])	Calculate 2**p for all p in the input array.
    'log', #(x[, out])	Natural logarithm, element-wise.
    'log2', #(x[, out])	Base-2 logarithm of x.
    'log10', #(x[, out])	Return the base 10 logarithm of the input array, element-wise.
    'expm1', #(x[, out])	Calculate exp(x) - 1 for all elements in the array.
    'log1p', #(x[, out])	Return the natural logarithm of one plus the input array, element-wise.
    'sqrt', #(x[, out])	Return the positive square-root of an array, element-wise.
    'square', #(x[, out])	Return the element-wise square of the input.
    'reciprocal', #(x[, out])	Return the reciprocal of the argument, element-wise.
    'sin', #(x[, out])	Trigonometric sine, element-wise.
    'cos', #(x[, out])	Cosine element-wise.
    'tan', #(x[, out])	Compute tangent element-wise.
    'arcsin', #(x[, out])	Inverse sine, element-wise.
    'arccos', #(x[, out])	Trigonometric inverse cosine, element-wise.
    'arctan', #(x[, out])	Trigonometric inverse tangent, element-wise.
    'sinh', #(x[, out])	Hyperbolic sine, element-wise.
    'cosh', #(x[, out])	Hyperbolic cosine, element-wise.
    'tanh', #(x[, out])	Compute hyperbolic tangent element-wise.
    'arcsinh', #(x[, out])	Inverse hyperbolic sine element-wise.
    'arccosh', #(x[, out])	Inverse hyperbolic cosine, element-wise.
    'arctanh', #(x[, out])	Inverse hyperbolic tangent element-wise.
]

_UNARY_BOOL = [  # return ndarrays
    'isreal', #(x)	Returns a bool array, where True if input element is real.
    'iscomplex', #(x)	Returns a bool array, where True if input element is complex.
]

_BINARY_MATH = [
    'add', #(x1, x2[, out])	Add arguments element-wise.
    'subtract', #(x1, x2[, out])	Subtract arguments, element-wise.
    'multiply', #(x1, x2[, out])	Multiply arguments element-wise.
    'divide', #(x1, x2[, out])	Divide arguments element-wise.
    #'logaddexp', #(x1, x2[, out])	Logarithm of the sum of exponentiations of the inputs.   # work only for real
    #'logaddexp2', #(x1, x2[, out])	Logarithm of the sum of exponentiations of the inputs in base-2. # work only for real
    'true_divide', #(x1, x2[, out])	Returns a true division of the inputs, element-wise.
    'floor_divide', #(x1, x2[, out])	Return the largest integer smaller or equal to the division of the inputs.
    'power', #(x1, x2[, out])	First array elements raised to powers from second array, element-wise.
    'fmax', #(x1, x2[, out])	Element-wise maximum of array elements.
    'fmin', #(x1, x2[, out])	Element-wise minimum of array elements.
    #'ldexp', #(x1, x2[, out])	Returns x1 * 2**x2, element-wise.
    'maximum', #(x1, x2[, out])	Element-wise maximum of array elements.
    'minimum', #(x1, x2[, out])	Element-wise minimum of array elements.
]

_BINARY_BOOL =[
    'logical_and', # (x1, x2[, out])	Compute the truth value of x1 AND x2 element-wise.
    'logical_or', #(x1, x2[, out])	Compute the truth value of x1 OR x2 element-wise.
    'logical_xor', #(x1, x2[, out])	Compute the truth value of x1 XOR x2, element-wise.
    'greater', #(x1, x2[, out])	Return the truth value of (x1 > x2) element-wise.
    'greater_equal', #(x1, x2[, out])	Return the truth value of (x1 >= x2) element-wise.
    'less', #(x1, x2[, out])	Return the truth value of (x1 < x2) element-wise.
    'less_equal', #(x1, x2[, out])	Return the truth value of (x1 =< x2) element-wise.
    'not_equal', #(x1, x2[, out])	Return (x1 != x2) element-wise.
    'equal', #(x1, x2[, out])	Return (x1 == x2) element-wise.
]

_UNARY_MATH_nd = [  # return nd array not source
    'zeros_like',
    'ones_like',
]

_UNARY_OPS = ['neg', 'pos', 'abs']
_BINARY_OPS = ['add', 'sub', 'mul', 'div', 'truediv',
                'floordiv', 'mod', 'pow', 'and', 'xor', 'or']
_METHODS = ['item', 'searchsorted']
_UNARY_METHODS = ['astype', 'clip', 'conj', 'conjugate', 'round']

# TODO: implement some of the following method directly in the class
# NUMPY_REDUCE_METHODS = ['all', 'any']
# NAN_REDUCE_METHODS = ['argmax', 'argmin', 'max', 'min', 'mean', 'sum', 'std', 'var', 'median']


# =============================================================================
# Location Indexer
# =============================================================================
class _LocIndexer(object):

    def __init__(self, source):
        self.source = source

    def loc2index(self, loc, axis=-1):
        if self.source.is_1d:
            axis = -1
        return loc2index(self.source.data, loc, axis)

    def select(self, key):
        """ Here we need to transform from a slice made with location
        arguments, to a locational slice that
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

            for axis, dim in enumerate(self.source.dims):

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

        slices = [slice(None, None, None) for axis in self.source.dims]

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
                    if key.start not in self.source.dims:
                        print("Error: the key '%s' doesn't exist" % key.start,
                              file=sys.stderr)
                        slices = [slice(0, 0, None) for axis in self.source.dims]
                        return slices
                    # a labelled  dimension is specified
                    # we need to cycle on all dimension
                    slices = _get_slice_by_label(key, slices)

        return slices

    def __getitem__(self, key):
        newkey = self.select(key)
        if len(newkey) == 1:
            return self.source[newkey[0]]
        else:
            return self.source[newkey[0], newkey[1]]

    def __setitem__(self, key, value):
        newkey = self.select(key)
        self.source[newkey] = value
        # TODO: check if shape has changed

# =============================================================================
# Source class definition
# =============================================================================
class Source(HasTraits):
    """
    Abstract Source class

    """
    name = Str
    description = Str
    author = Str

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

    shape = Property

    size = Property

    dtype = Property

    T = Property

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

    # generic attibutes to easy access to coordinates
    X = Property
    Y = Property  # for 2D data
    Z = Property  # for 3D data (max for the moment)

    loc = Property

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
            index or location in the axis unit of the row to extract
            depends on byindex keyword
        byindex: boolean, optional, default=False
            if False, value is a location in axis unit,
            if True it is the index of the row
        width: float, optional, default=0
            if 0, only one row is taken
            else the specified width is taken into account, the returned data
            correspond to the integration either of a number of row
            or of a range of location depending on the byindex keyword

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
            index or location in the axis unit of the row to extract
            depends on byindex keyword
        byindex: boolean, optional, default=False
            if False, value is a location in axis unit,
            if True it is the index of the column
        width: float, optional, default=0
            if 0, only one row is taken
            else the specified width is taken into account, the returned data
            correspond to the integration either of a number of row
            or of a range of location depending on the byindex keyword

        Returns
        -------
        data: DataFrame
            return a pandas Dataframe with the col selected
        """

        return self._getsection(value, axis=0, byindex=byindex, width=width)

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

    def get_index(self, loc, axis=-1):
        """
        Get the index of data from a data location

        Parameters
        ----------
        loc: Float,
            The location in data units
        axis: Int, default=-1
            The axis where to select the location

        Returns
        -------
            index: Int,
                The index corresponding to the location
        """
        if axis == 0:
            a = self.data.T
        else:
            a = self.data
        return loc2index(a, loc)

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

    def copy(self, deepcopy=True):
        """
        Make a copy of the source

        Parameters
        ----------
        deepcopy: boolean, default=True
            make a deep recursive copy

        Returns
        -------
        source: Source
        """
        if not deepcopy:
            return copy.copy(self)  # shallow copy
        else:
            return copy.deepcopy(self)  # deep (recursive) copy

    def transpose(self, inplace=True):
        if not self.is_1d:
            if inplace:
                self.data = self.data.T
                self._dims.reverse()
                return self
            else:
                newsource = self.copy()
                newsource._data_orig = self.data
                newsource.data = self.data.T
                newsource.dims = self.dims[:]
                newsource._dims.reverse()
                newsource.name = "copy of %s"%self.name
                return newsource
        else:
            print('Warning: can not transpose 1D data!', file=sys.stderr)
            return self

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
            raise ValueError('Must be a DataFrame')

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
            self._dims = ['Z', 'Y', 'X'][-self.ndims:]

        return self._dims

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
        newdata = self.get_real()
        # create a new source with the new data
        newsource = self.copy()
        newsource._data_orig = newdata
        newsource.data = newdata
        return newsource

    # imag
    def _get_imag(self):
        """
        property getter: return the imaginary part of the original data
        """
        newdata = self.get_imag()
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

    # shape
    def _get_shape(self):
        """
        Property getter: shape of the data array
        """
        return self.data.shape

    # size
    def _get_size(self):
        """
        Property getter: size of the data array
        """
        return self.data.size

    # dtype
    def _get_dtype(self):
        """
        Property getter: type of the data array values
        """
        return self.data.values.dtype

    # T
    def _get_T(self):
        """
        Property getter: type of the data array values
        """
        return self.transpose()

    # *************************************************************************
    # Private functions
    # *************************************************************************

    def _getsection(self, index, **kargs):
        """
        Extract a section by location or by index

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
            txt2 = ' by location'
            w2 = width / 2.0
            lb = loc2index(x, index - w2)
            ub = loc2index(x, index + w2)
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

        newdata = pd.DataFrame(section, index=y).T

        # create a new source with the new data
        newsource = self.copy()
        newsource._data_orig = newdata
        newsource.data = newdata
        return newsource

    def __nonzero__(self):
        return not (self.data.size == 0)

    def __eq__(self, other):
        try:
            assert_frame_equal(self.data, other.data)
            return True
        except:  # appeantly AssertionError doesn't catch all
            return False
        return res

    def __ne__(self, other):
        return not self.__eq__(other)


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
        if self.is_1d:
            keys = slice(0, 1, None), keys[0]

        # slice the data

        newdata = self.data.iloc[tuple(keys)]

        # create a new source with the new data
        newsource = self.copy()

        newsource._data_orig = newdata
        newsource.data = newdata

        return newsource

    def __getattr__(self, item):
        """ The item attribute was not found
        (should be subclassed in derived classes)
        """
        raise NotImplementedError('should be implemented by subclass')

    def __str__(self):

        np.set_printoptions(precision=3, threshold=5, edgeitems=3, suppress=True,
                            formatter={'all':lambda x: '%.2g' % x,
                                       'complex_kind':lambda x: '%.2g%+.2gj' % (x.real, x.imag)})

        if self:
            txt = "<< Source|name: '%s'|type: %dD|dims: %s|shape: %s >>" % (
            self.name, self.ndims, str(tuple(self.dims)), str(self.data.shape))
        else:
            txt = "<< Source >> None"
            return txt

        txt += '\n- data:\n'
        arr = '%s ' % self.data.values
        arr = arr.split('\n')
        arr = ["\t%s" % a for a in arr]
        arr = "\n".join(arr)
        txt += arr
        txt += '\n- coordinates:'
        for k in self.dims[::-1]:
            if k==self.dims[-1]:
                txt += "\n\t* %s : %s" % (k, self.data.columns.values)
            elif k==self.dims[-2]:
                txt += "\n\t* %s : %s" % (k, self.data.index.values)
            else:
                raise NotImplementedError('Z not implemented')
        return txt + '\n'

    #
    # injected  methods partly copied from xray - see github:
    #
    def _op_result(self, data):
        """Create the returned values of functions
        should be a source with the same size and coordinates
        """
        new = self.copy()
        new.data = pd.DataFrame(data.copy(), index=self.data.index, columns=self.data.columns)
        return new

    @staticmethod
    def _unary_op(f):
        @functools.wraps(f)
        def func(self, *args, **kwargs):
            data = f(self._data.values, *args, **kwargs)
            return self._op_result(data)
        return func

    @staticmethod
    def _binary_op(f, reflexive=False):
        @functools.wraps(f)
        def func(self, other):
            if isinstance(other, Source):
                other = other.data.values
            data = f(self._data.values, other) if not reflexive else f(other, self._data.values)
            return self._op_result(data)

        return func

    @staticmethod
    def _inplace_binary_op(f):
        @functools.wraps(f)
        def func(self, other):
            self.data.values = f(self._data.values, other)
            return self
        return func

    @property
    def __array_struct__(self):  # to emulate a numpy array
        struc = self.data.values.__array_struct__
        return struc

    def __array_prepare__(self, *args):  # called when ufuncs are applied to the array
        obj = self.data.values.__array_prepare__(*args)
        return obj

    def __array_wrap__(self, *args):  # called when element-wise  ufuncs are applied to the array not changing dims
        uf, objs, huh = args[1]
        nobjs = []
        for obj in objs:
            nobjs.append(obj.data)
        nargs = (args[0],(uf, tuple(nobjs), huh))
        try:
            data = self._data.values.__array_wrap__(*nargs)
        except:
            print ("Unexpected error:", sys.exc_info()[0])
            raise
        return self._op_result(data)

###### Set array operators ########
# borrowed and modified from XRAY.

def _inject_ops(cls):
    def op_str(name):
        return '__%s__' % name

    def op(name):
        return getattr(operator, op_str(name))

    def _values_method_wrapper(name):
        def func(self, *args, **kwargs):
            return getattr(self.data, name)(*args, **kwargs)
        func.__name__ = name
        func.__doc__ = getattr(np.ndarray, name).__doc__
        return func

    def _method_wrapper(name):
        def func(self, *args, **kwargs):
            f = getattr(self, name)(*args, **kwargs)
            return f
        func.__name__ = name
        func.__doc__ = getattr(np.ndarray, name).__doc__
        return func

    for name in _UNARY_OPS:
        setattr(cls, op_str(name), cls._unary_op(op(name)))

    for name in _BINARY_OPS:  #in place and reflexive variants
        setattr(cls, op_str(name), cls._binary_op(op(name)))
        setattr(cls, op_str('r' + name), cls._binary_op(op(name), reflexive=True))
        setattr(cls, op_str('i' + name), cls._inplace_binary_op(op('i' + name)))

    for name in _UNARY_METHODS:
        setattr(cls, name, cls._unary_op(_method_wrapper(name)))

    for name in _METHODS:
        setattr(cls, name, _values_method_wrapper(name))


_inject_ops(Source)


if __name__ == '__main__':

    # Test
    # ----
    import numpy as np

    source1 = Source()

    # create a simple 1D source
    x = np.arange(100.)*np.pi/50.
    y = np.sin(x)
    source1.data = pd.DataFrame([y], index=[0], columns=x )

    # original
    print(source1)

    # some features
    print('size', source1.size)
    print('shape', source1.shape)
    print('dtype', source1.dtype)

    # make a copy
    sourceb = source1.copy()

    # modify sourceb
    sourceb.data /= 2

    # compare source and source2: should be different, right?
    assert(sourceb != source1)

    # slicing on index
    print(source1[10:14])

    # slicing on location
    print(source1.loc[0.58:0.72])

    # get  a single item by index
    print(source1[10])

    # get  a single item by location (the closest)
    print(source1.loc[0.69])

    # create a simple 2D source
    source2 = Source()

    y = np.arange(20.)*np.pi/20.

    z = np.array([[np.sin(xi)*np.sin(yi) for xi in x] for yi in y])
    print(z.shape)

    source2.data = pd.DataFrame(z, index=y, columns=x )

    # on the last dimension should work as for 1d
    # original
    print(source2)

    # slicing on index
    print(source2[:,10:14])

    # slicing on location
    print(source2.loc[:,0.58:0.72])

    # get  a single item by index
    print(source2[:,10])

    # get  a single item by location (the closest)
    print(source2.loc[0.15,0.69])

    # get row by index (same as direct indexing)
    row = source2.get_row(0, byindex=True)
    print(row)
    assert(row == source2[0])

    # get a row with some width
    row = source2.get_row(0.5,width=.2)
    print(row)
    assert(row != source2[0])

    # transpose
    print("TRANSPOSE\n-------------")
    #1D
    print(source1)
    print(source1.T)

    #2D
    print(source2)
    print(source2.T)

    #2D not in place (we can not use the property T in this case
    print('NOT IN PLACE')
    print(source2)
    source3 = source2.transpose(inplace=False)
    source3.data *= 3
    print(source3)
    assert source3 != source2

    print('\nUNARY FUNCTION')
    a = source3

    for f in _UNARY_MATH:
        try:
            b = getattr(np, f)(a)
            if DEBUG : print(f, b)
        except TypeError as e:
            print("ERROR: ",e.message, file=sys.stderr)

        assert isinstance(b, Source), "%s(b) is not Source (type(b): %s)" % (f,type(b))
        assert b is not a

    print('\nUNARY MATH nd')
    for f in _UNARY_MATH_nd:
        try:
            b = getattr(np, f)(a)
            if DEBUG : print(f, b)
        except TypeError as e:
            print("ERROR: ",e.message, file=sys.stderr)

        assert isinstance(b, np.ndarray), "%s(b) is not nd.array (type(b): %s)" % (f,type(b))
        assert b is not a

    print('\nUNARY BOOL')
    for f in _UNARY_BOOL:
        try:
            b = getattr(np, f)(a)
            if DEBUG : print(f, b)
        except TypeError as e:
            print("ERROR: ",e.message, file=sys.stderr)

        assert isinstance(b, np.ndarray), "%s(b) is not nd.array (type(b): %s)" % (f,type(b))
        assert b is not a

    # Binary ufuncs
    b = source3.copy()/2.5
    for f in _BINARY_MATH:
        try:
            c = getattr(np, f)(a, b)
            if DEBUG : print(f, c)
        except TypeError as e:
            print("ERROR: ",e.message, file=sys.stderr)
        assert isinstance(b, Source), "%s(b) is not Source (type(b): %s)" % (f,type(c))

    # basic ops
    assert (a+b)[4] == a[4] + b[4]
    assert (a-b)[4] == a[4] - b[4]
    assert (a/b)[4] == a[4] / b[4]
    assert (a*b)[4] == a[4] * b[4]
    assert (a//b)[4] == a[4] // b[4]
    print("a+b")
    print (a+b)

    a = (a+b).copy()
    for f in _UNARY_OPS:
        print("try ",f)
        b = getattr(a, "__%s__" % f)()
        if DEBUG: print ('__%s__' %f, b)

    for f in _BINARY_OPS:
        try:
            c = getattr(a, "__%s__" % f)(b)
            if DEBUG: print ('__%s__' %f, c)
        except TypeError as e:
            print('__%s__' % f, e.message)

    #_UNARY_METHODS test ['astype', 'clip', 'conj', 'conjugate', 'round']

    b = a.astype(complex)
    b = b.clip(2,8)
    b = b.conj()
    b = b.conjugate()
    b = b.round(1)
    print("OK ", b)

    b = a.copy()
    b = b.T
    print(b)
    print(a)
    assert b.shape == a.shape[::-1]


    a = source1[0:10]
    b = source1[10:20]
    print("a", a)
    print("b", b)
    c = a + b / a
    print(c)
