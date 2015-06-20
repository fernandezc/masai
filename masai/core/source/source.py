# -*- coding: utf-8 -*-

"""
Define basic source functionalities

This is a container for the spectroscopic data or others,
and the associated parameters

Rich slicing is provided either by


"""

# Copyright (c) C. Fernandez' group @ LCS.
# Distributed under the terms of the CeCILL-B License.

from __future__ import print_function, division

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
from masai.core.utils.misc import is_sequence
from masai.core.source.locindexer import loc2index, LocIndexer

# constants
DEBUG = True

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
    # Constructor
    # *************************************************************************
    def __init__(self, *args, **kargs):

        ### TODO

        super(Source, self).__init__(*args, **kargs)

    # *************************************************************************
    # public functions
    # *************************************************************************

    def set(self, *args, **kargs):
        values = None
        units = None
        if len(args)>0:
            values = args[0]
        if len(args)>1:
            self.units = args[1]
        if values is None:
            values = kargs.get('values',[])
        if units is None:
            units = kargs.get('units', None)
        self.units = units
        self.data = values
        if kargs:
            super(Source, self).set(**kargs)

    def get_row(self, value, byindex=False, width=0, mean=True):
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
        mean: boolean, default=True
            if True an average of all col of the section will be performed
            else the sum is done.

        Returns
        -------
        data: DataFrame
            return a pandas Dataframe with the row

        """
        row = self.get_section(value, axis=0, byindex=byindex, width=width)
        #if mean:
        #    row = row.mean(axis=-1)
        #else:
        #    row = row.sum(axis=-1)
        return row

    def get_col(self, value, byindex=False, width=0, mean=True):
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
        mean: boolean, default=True
            if True an average of all col of the section will be performed
            else the sum is done.
        Returns
        -------
        data: DataFrame
            return a pandas Dataframe with the col selected
        """

        col = self.get_section(value, axis=-1, byindex=byindex, width=width)
        #if mean:
        #    col = col.mean(axis=-1)
        #else:
        #    col = col.sum(axis=-1)
        return col

    def get_section(self, index, **kargs):
        """
        Extract a section by location or by index

        """
        if not self.is_2d:
            print("1D data, so we return a copy of the input data itself!")
            return self.copy()

        axis = kargs.get('axis', -1)
        byindex = kargs.get('byindex', False)
        width = kargs.get('width', 0)

        if not byindex:
            w2 = width / 2.0
            lb = index - w2
            ub = index + w2
        else:
            w2 = int(width / 2)
            lb = index - w2
            ub = index + w2 + 1

        if abs(ub-lb) > 0:
            s = slice(lb,ub,None)
        else:
            s = index

        if not byindex:
            if axis==0:
                return self.loc[s]
            else:
                return self.loc[:,s]
        else:
            if axis==0:
                return self[s]
            else:
                return self[:,s]

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

        mean = kargs.get('mean', False)
        sum = kargs.get('sum', False)

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

    def copy(self):
        cls = type(self)
        new = cls()
        new.__setstate__(self.__getstate__())
        return new

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

    #def argmax(self, axis=-1):
        """
        Indices of the maximum values along an axis.

        Parameters
        ----------
        axis : int, optional
            By default, the index is into columns axis

        Returns
        -------
        index_array : ndarray of ints
            Array of indices into the array. It has the same shape as `a.shape`
            with the dimension along `axis` removed.

        Notes
        -----
        In case of multiple occurrences of the maximum values, the indices
        corresponding to the first occurrence are returned.
        """
     #   return self.data.argmax

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
            self._data = value.copy()
        else:
            self._data = pd.DataFrame(value.copy(), index=self._data.index,
                                      columns=self._data.columns)

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
    # -------------------------------------------------------------------------

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
        return LocIndexer(self)

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

    def _print_info_(self, format="str"):
        if self:
            row, col = self.shape
            tdims = ""
            if row>1:
                    tdims += "rows: '%s'|"%self.dims[0]
            if col>1:
                    tdims += "columns: '%s'"%self.dims[-1]
            tshape = '%d rows x %d columns'%self.shape
            txt = "<< Source >>\n|| name: '%s' | %s | shape: %s ||\n" % (
            self.name, tdims, tshape)
        else:
            txt = "None"
            return None
        if format=='html':
            txt = txt.replace('<<', '<tr><td><h3>')
            txt = txt.replace('>>', '</h3></td></tr>')
            txt = txt.replace('||\n', '</td></tr>')
            txt = txt.replace('||', '</td></tr><tr><td><b>')
            txt = txt.replace('|', '</td><td><b>')
            txt = txt.replace(':', '</b></td><td>')
            txt = txt.replace(':', '</b></td><td>')
        return txt

    def _print_data_(self, format="str", dim=None, val=None):
        e = np.get_printoptions()['edgeitems']
        n = 2*e+1 # number of item on a column
        n = min(self.shape[-1], n)
        if dim is None:
            txt = "|| data |"
            arr = '%s ' % self.data.values
            arr = arr.strip().split('\n')
            arr = ["\t%s" % a for a in arr]
            arr = "\n".join(arr)
            arr += '||\n'
            txt += arr.replace(' ',';')

            txt = txt.replace('\t;...,;', '\t;[%s]'%str('...,;'*n)[:-1])
        else:
            txt = "|| %s |"%dim
            arr = '%s ' % val
            arr = arr.strip().split('\n')
            arr = ["\t%s" % a for a in arr]
            arr = "\n".join(arr)
            arr += '||\n'
            txt += arr.replace(' ',';')
        txt = txt.replace(',', '')
        txt = txt.replace(';;', ';')
        if format=='html':
            colspan = (7 if self.ndims>1 else 5)
            txt = txt.replace('<<', '<tr><td><h3>')
            txt = txt.replace('>>', '</h3></td></tr>')
            txt = txt.replace('||\n', '</b></td></tr>')
            txt = txt.replace('||', '</td></tr><tr><td><b>')
            txt = txt.replace('|', '</td><td colspan=%d>'%colspan)
            txt = txt.replace(']\n\t;[', '</td></tr><tr><td>')
            txt = txt.replace(']\n\t;...,;\n\t\;[', '</td></tr><tr><td>')
            txt = txt.replace(';', '</td><td>')
            txt = txt.replace('[', '<table><tr><td>')
            txt = txt.replace(']', '</td></tr></table>')
            txt = txt.replace(';', '</td><td>')
            txt = txt.replace('<tr>', '\n<tr>')
            txt = txt.replace('<table>', '\n<table>')
            txt = txt.replace('</table>', '</table>\n')

        return txt

    def _printoptions(self, precision=3, threshold=5,
                       edgeitems=4, suppress=True,
                       formatter={
                    'all':lambda x: '%.3g' % x.real,
                    'str_kind':lambda x: '%s'%x,
                    'complex_kind':lambda x: '%.3g%+.3gj' % (x.real, x.imag)}):
        """
        Set options fotr numpy array printing
        """
        np.set_printoptions(precision=precision,
                            threshold=threshold,
                            edgeitems=edgeitems,
                            suppress=suppress,
                            formatter=formatter)

    def __str__(self):
        '''
        representation of the object as a string for prionting
        '''
        self._printoptions()

        txt = self._print_info_()
        if not txt:
            return 'None'

        txt += self._print_data_()

        row, col = self.shape
        if col>1 or self.size==1:
            txt += self._print_data_(dim=self.dims[-1],
                                 val=self.data.columns.values)
        if row>1:
            txt += self._print_data_(dim=self.dims[0],
                                     val=self.data.index.values)
        txt = txt.replace(';', ' ')
        return txt + '\n'

    def _repr_html_(self):
        """
        representation of the object as an HTML string to display in notebook
        """

        self._printoptions()

        txt = self._print_info_('html')
        if not txt:
            return 'None'
        txt += self._print_data_('html')
        row, col = self.shape
        if row>1:
            txt += self._print_data_('html', dim=self.dims[0],
                                     val=self.data.index.values)
        if col>1 or self.size==1:
            txt += self._print_data_('html', dim=self.dims[-1],
                                 val=self.data.columns.values)
        return txt


    # =========================================================================
    # injected  methods (borrowed and modified from XRAY.)
    # =========================================================================
    def _op_result(self, data):
        """Create the returned values of functions
        should be a source with the same size and coordinates
        """
        new = self.copy()
        new.data = pd.DataFrame(data, index=self.data.index,
                                columns=self.data.columns)
        return new

    @staticmethod
    def _unary_op(f):
        @functools.wraps(f)
        def func(self, *args, **kwargs):
            data = f(self._data.values.copy(), *args, **kwargs)
            return self._op_result(data)
        return func

    @staticmethod
    def _binary_op(f, reflexive=False):
        @functools.wraps(f)
        def func(self, other):
            if isinstance(other, Source):
                other = other.data.values.copy()
            data = f(self._data.values.copy(), other) \
                        if not reflexive else f(other, self._data.values.copy())
            return self._op_result(data)
        return func

    @staticmethod
    def _inplace_binary_op(f):
        @functools.wraps(f)
        def func(self, other):
            self.data = f(self._data.values, other)
            return self
        return func

    @property
    def __array_struct__(self):
        # to emulate a numpy array
        struc = self.data.values.__array_struct__
        return struc

    def __array_prepare__(self, *args):
        # called when ufuncs are applied to the array
        obj = self.data.values.__array_prepare__(*args)
        return obj

    def __array_wrap__(self, *args):
        # called when element-wise ufuncs are applied to the array
        # not changing dims
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

    @classmethod
    def _reduce_method(self, func):
        def wrapped_func(self, axis=None, **kwargs):
            return self._reduce(func, axis=axis, **kwargs)
        return wrapped_func

    @classmethod
    def _pandas_reduce_method(self, func):
        def wrapped_func(self, axis=None, **kwargs):
            return self._pandas_reduce(func, axis=axis, **kwargs)
        return wrapped_func

    def _reduce(self, func, axis=None, **kwargs):
        """Reduce this dataset by applying `func` along some dimension(s).
        Parameters
        ----------
        func : function
            Function which can be called in the form
            `f(x, axis=axis, **kwargs)` to return the result of reducing an
            np.ndarray over an integer valued axis.
        axis : int or sequence of int, optional
            Dimension(s) over which to apply `func`.  By default `func` is
            applied over all dimensions.
        **kwargs : dict
            Additional keyword arguments passed on to ``func``.
        Returns
        -------
        reduced : Dataset
            Dataset with this object's Source replaced with new Source
            of summarized data and the indicated dimension(s) removed.
        """

        #if dim is not None:
        #    axis = self.get_axis_num(dim)

        newdata = func(self.data.values, axis=axis, **kwargs)
        newdata = np.array([newdata])

        newsource = Source()
        newsource.data = pd.DataFrame(newdata)
        if (axis == -1 or axis == 1) and not self.is_1d:
            newsource.data = pd.DataFrame(newdata.T)
            newsource.data.index = self.data.index.copy()
            newsource.data.columns = [func.__name__]
        elif axis == 0:
            newsource.data = pd.DataFrame(newdata)
            newsource.data.columns = self.data.columns.copy()
            newsource.data.index = [func.__name__]
        else:
            newsource.data = pd.DataFrame(newdata)
        return newsource

    def _pandas_reduce(self, func, axis=None, **kwargs):
        """Reduce this dataset by applying `func` along some dimension(s).
        Parameters
        ----------
        func : function
            Function which can be called in the form
            `f(x, axis=axis, **kwargs)` to return the result of reducing an
            np.ndarray over an integer valued axis.
        axis : int or sequence of int, optional
            Dimension(s) over which to apply `func`.  By default `func` is
            applied over all dimensions.
        **kwargs : dict
            Additional keyword arguments passed on to ``func``.
        Returns
        -------
        reduced : Dataset
            Dataset with this object's Source replaced with new Source
            of summarized data and the indicated dimension(s) removed.
        """

        #if dim is not None:
        #    axis = self.get_axis_num(dim)
        if axis==-1:
            axis = self.ndims-1

        try:
            newdata = getattr(self.data, func)(axis=axis, **kwargs)
        except:
            axis = self.ndims-1
            newdata = getattr(self.data, func)(axis=axis, **kwargs)
        if axis is None:
            newdata = getattr(newdata, func)(axis=0, **kwargs)
            return newdata

        newsource = Source()
        newdata = pd.DataFrame(newdata)

        if axis == 1:
            newsource.data = newdata
            newsource.data.index = self.data.index.copy()
            newsource.data.columns = [func]
        elif axis == 0:
            newsource.data = newdata.T
            newsource.data.columns = self.data.columns.copy()
            newsource.data.index = [func]
        else:
            newsource.data = newdata
        return newsource

# *****************************************************************************
# Set array operators
# borrowed and modified from XRAY.
# *****************************************************************************
#
# UFUNCS
#
UNARY_MATH =[
    'negative', #(x[, out])
    # Numerical negative, element-wise.
    'absolute', #(x[, out])
    # Calculate the absolute value element-wise.
    'rint', #(x[, out])
    # Round elements of the array to the nearest integer.
    'sign', #(x[, out])
    # Returns an element-wise indication of the sign of a number.
    'conj', #(x[, out])
    # Return the complex conjugate, element-wise.
    'exp', #(x[, out])
    # Calculate the exponential of all elements in the input array.
    'exp2', #(x[, out])
    # Calculate 2**p for all p in the input array.
    'log', #(x[, out])
    # Natural logarithm, element-wise.
    'log2', #(x[, out])
    # Base-2 logarithm of x.
    'log10', #(x[, out])
    # Return the base 10 logarithm of the input array, element-wise.
    'expm1', #(x[, out])
    # Calculate exp(x) - 1 for all elements in the array.
    'log1p', #(x[, out])
    # Return the natural logarithm of one plus the input array, element-wise.
    'sqrt', #(x[, out])
    # Return the positive square-root of an array, element-wise.
    'square', #(x[, out])
    # Return the element-wise square of the input.
    'reciprocal', #(x[, out])
    # Return the reciprocal of the argument, element-wise.
    'sin', #(x[, out])
    # Trigonometric sine, element-wise.
    'cos', #(x[, out])
    # Cosine element-wise.
    'tan', #(x[, out])
    # Compute tangent element-wise.
    'arcsin', #(x[, out])
    # Inverse sine, element-wise.
    'arccos', #(x[, out])
    # Trigonometric inverse cosine, element-wise.
    'arctan', #(x[, out])
    # Trigonometric inverse tangent, element-wise.
    'sinh', #(x[, out])
    # Hyperbolic sine, element-wise.
    'cosh', #(x[, out])
    # Hyperbolic cosine, element-wise.
    'tanh', #(x[, out])
    # Compute hyperbolic tangent element-wise.
    'arcsinh', #(x[, out])
    # Inverse hyperbolic sine element-wise.
    'arccosh', #(x[, out])
    # Inverse hyperbolic cosine, element-wise.
    'arctanh', #(x[, out])
    # Inverse hyperbolic tangent element-wise.
]

UNARY_BOOL = [  # return ndarrays
    'isreal', #(x)
    # Returns a bool array, where True if input element is real.
    'iscomplex', #(x)
    # Returns a bool array, where True if input element is complex.
]

BINARY_MATH = [
    'add', #(x1, x2[, out])
    # Add arguments element-wise.
    'subtract', #(x1, x2[, out])
    # Subtract arguments, element-wise.
    'multiply', #(x1, x2[, out])
    # Multiply arguments element-wise.
    'divide', #(x1, x2[, out])
    # Divide arguments element-wise.
    'true_divide', #(x1, x2[, out])
    # Returns a true division of the inputs, element-wise.
    'floor_divide', #(x1, x2[, out])
    # Return the largest integer smaller or equal to the division of the inputs.
    'power', #(x1, x2[, out])
    # First array elements raised to powers from second array, element-wise.
    'fmax', #(x1, x2[, out])
    # Element-wise maximum of array elements.
    'fmin', #(x1, x2[, out])
    # Element-wise minimum of array elements.
    'maximum', #(x1, x2[, out])
    # Element-wise maximum of array elements.
    'minimum', #(x1, x2[, out])
    # Element-wise minimum of array elements.
]

BINARY_BOOL =[
    'logical_and', # (x1, x2[, out])
    # Compute the truth value of x1 AND x2 element-wise.
    'logical_or', #(x1, x2[, out])
    # Compute the truth value of x1 OR x2 element-wise.
    'logical_xor', #(x1, x2[, out])
    # Compute the truth value of x1 XOR x2, element-wise.
    'greater', #(x1, x2[, out])
    # Return the truth value of (x1 > x2) element-wise.
    'greater_equal', #(x1, x2[, out])
    # Return the truth value of (x1 >= x2) element-wise.
    'less', #(x1, x2[, out])
    # Return the truth value of (x1 < x2) element-wise.
    'less_equal', #(x1, x2[, out])
    # Return the truth value of (x1 =< x2) element-wise.
    'not_equal', #(x1, x2[, out])
    # Return (x1 != x2) element-wise.
    'equal', #(x1, x2[, out])
    # Return (x1 == x2) element-wise.
]

UNARY_MATH_nd = [  # return nd array not source
    'zeros_like',
    'ones_like',
]

UNARY_OPS = ['neg', 'pos', 'abs']
BINARY_OPS = ['add', 'sub', 'mul', 'div', 'truediv',
                'floordiv', 'mod', 'pow', 'and', 'xor', 'or']
METHODS = ['item', 'searchsorted']
UNARY_METHODS = ['astype', 'clip', 'conj', 'conjugate', 'round']

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

    for name in UNARY_OPS:
        setattr(cls, op_str(name), cls._unary_op(op(name)))

    for name in BINARY_OPS:  #in place and reflexive variants
        setattr(cls, op_str(name), cls._binary_op(op(name)))
        setattr(cls, op_str('r' + name), cls._binary_op(op(name),
                                                                reflexive=True))
        setattr(cls, op_str('i' + name), cls._inplace_binary_op(op('i' + name)))

    for name in UNARY_METHODS:
        setattr(cls, name, cls._unary_op(_method_wrapper(name)))

    for name in METHODS:
        setattr(cls, name, _values_method_wrapper(name))

_inject_ops(Source)

PANDAS_REDUCE_METHODS = ['all', 'any',
                         'idxmax', 'idxmin',
                         'max', 'min', 'mean',
                         'sum', 'std', 'var', 'median']

def _inject_reduce_methods(cls):
    methods = ([(name, getattr(np, name)) for name
               in PANDAS_REDUCE_METHODS])
    for name, f in methods:
        func = cls._reduce_method(f)
        func.__name__ = name
        setattr(cls, name, func)

def _inject_pandas_reduce_methods(cls):
    for name in PANDAS_REDUCE_METHODS:
        func = cls._pandas_reduce_method(name)
        setattr(cls, name, func)

_inject_pandas_reduce_methods(Source)

# TESTING ####################################################################
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
    sourceb /= 2

    # compare source and source2: should be different, right?
    print(sourceb != source1)
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
    assert(np.all(source1.loc[0.58:0.72].data.columns \
                   == source2.loc[:,0.58:0.72].data.columns))

    # get  a single item by index
    print(source2[:,10])

    # get  a single item by location (the closest)
    print(source2.loc[0.15,0.69])

    # get row by index (same as direct indexing)
    row = source2.get_row(0, byindex=True)
    print("row = source2.get_row(0, byindex=True)", row)
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
    print('NOT IN PLACE\n------------------')
    print(source2)
    source3 = source2.transpose(inplace=False)
    source3.data *= 3
    print(source3)
    assert source3 != source2

    print('\nUNARY FUNCTION\n------------------')
    a = source3

    for f in UNARY_MATH:
        try:
            b = getattr(np, f)(a)
            if DEBUG : print(f, b)
        except TypeError as e:
            print("ERROR: ",e.message, file=sys.stderr)

        assert isinstance(b, Source), \
            "%s(b) is not Source (type(b): %s)" % (f,type(b))
        assert b is not a

    print('\nUNARY MATH nd\n------------------')
    for f in UNARY_MATH_nd:
        try:
            b = getattr(np, f)(a)
            if DEBUG : print(f, b)
        except TypeError as e:
            print("ERROR: ",e.message, file=sys.stderr)

        assert isinstance(b, np.ndarray), \
            "%s(b) is not nd.array (type(b): %s)" % (f,type(b))
        assert b is not a

    print('\nUNARY BOOL\n------------------')
    for f in UNARY_BOOL:
        try:
            b = getattr(np, f)(a)
            if DEBUG : print(f, b)
        except TypeError as e:
            print("ERROR: ",e.message, file=sys.stderr)

        assert isinstance(b, np.ndarray),\
            "%s(b) is not nd.array (type(b): %s)" % (f,type(b))
        assert b is not a

    # Binary ufuncs
    print("BINARY MATH\n------------------")
    b = source3.copy()/2.5
    for f in BINARY_MATH:
        try:
            c = getattr(np, f)(a, b)
            if DEBUG : print(f, c)
        except TypeError as e:
            print("ERROR: ",e.message, file=sys.stderr)
        assert isinstance(b, Source), \
            "%s(b) is not Source (type(b): %s)" % (f,type(c))

    # Binary logical ufuncs
    b = source3.copy()/2.5
    print("BINARY BOOL\n------------------")
    for f in BINARY_BOOL:
        try:
            c = getattr(np, f)(a, b)
            if DEBUG : print(f, c)
        except TypeError as e:
            print("ERROR: ",e.message, file=sys.stderr)
        assert isinstance(b, Source), \
            "%s(b) is not Source (type(b): %s)" % (f,type(c))

    # basic ops
    print("BASIC OPS\n------------------")
    assert (a+b)[4] == a[4] + b[4]
    assert (a-b)[4] == a[4] - b[4]
    assert (a/b)[4] == a[4] / b[4]
    assert (a*b)[4] == a[4] * b[4]
    assert (a//b)[4] == a[4] // b[4]
    print("a+b")
    print (a+b)


    a = (a+b).copy()
    print("UNARY OPS\n------------------")
    for f in UNARY_OPS:
        b = getattr(a, "__%s__" % f)()
        if DEBUG: print ('__%s__' %f, b)

    print("BINARY OPS\n------------------")
    for f in BINARY_OPS:
        try:
            c = getattr(a, "__%s__" % f)(b)
            if DEBUG: print ('__%s__' %f, c)
        except TypeError as e:
            print('__%s__' % f, e.message)

    #_UNARY_METHODS test ['astype', 'clip', 'conj', 'conjugate', 'round']
    print("UNARY METHODS\n------------------")
    b = a.astype(complex)
    b = b.clip(2,8)
    b = b.conj()
    b = b.conjugate()
    b = b.round(1)

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


    print("source", a)
    print("square of source", a**2)

    a = source2.T.copy()

    print("NAN_REDUCE_METHODS+NUMPY_REDUCE_METHODS\n------------------")
    #for f in NAN_REDUCE_METHODS+NUMPY_REDUCE_METHODS:
    #    b = getattr(a, "%s" % f)(axis=-1)
    #    if DEBUG: print ('axis=-1 %s' %f, b)
    #    b = getattr(a, "%s" % f)(axis=0)
    #    if DEBUG: print ('axis=0 %s' %f, b)
    #    b = getattr(a, "%s" % f)()
    #    if DEBUG: print ('axis:None %s' %f, b)
    for f in PANDAS_REDUCE_METHODS:
        print('try %s'%f)
        b = getattr(a, "%s" % f)()
        print ('axis=None %s' %f, b)
        b = getattr(a, "%s" % f)(axis=0)
        print ('axis=0 %s' %f, b)
        b = getattr(a, "%s" % f)(axis=1)
        print ('axis=1 %s' %f, b)

    # TODO: some of these methods doesn' work for complex data type
