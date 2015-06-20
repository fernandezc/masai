# -*- coding: utf-8 -*-

"""A one-line description.

A longer description that spans multiple lines.  Explain the purpose of the
file and provide a short list of the key classes/functions it contains.  This
is the docstring shown when some does 'import foo;foo?' in IPython, so it
should be reasonably useful and informative.
"""
# Copyright (c) C. Fernandez' group @ LCS.
# Distributed under the terms of the CeCILL-B License.

from __future__ import print_function, division

# -----------------------------------------------------------------------------
# Globals
# -----------------------------------------------------------------------------

# stdlib
import sys

# third-party

# -----------------------------------------------------------------------------
# Local utilities and constants
# -----------------------------------------------------------------------------
from masai.core.utils.misc import is_sequence

# constants
DEBUG = True
__all__ = ['LocIndexer']

from pandas import DataFrame


def loc2index(data, loc, axis=-1):
    """
    Return the index of a location along the last axis (columns) of a DataFrame

    Parameters
    ----------
        data: DataFrame
        loc: float
            a column location in the column units
        axis: integer, default=-1

    Returns
    -------
        index: int
            the index of the location

    """
    if isinstance(data, DataFrame):
        if axis == -1 or axis == 1:
            a = data.columns
        else:
            a = data.index
    else:
        # should be a 1d array
        a = data

    if a[0] < a[-1]:
        id = a[a <= loc].size
    else:
        id = a[a >= loc].size-1

    # check closest
    dloc0 = 1.e15
    idx = 0
    id0 = id
    for i in [-1, 0, +1]:
        id = max(min(id0 + i, a.size - 1), 0)
        dloc = abs(loc - a[id])
        if dloc < dloc0:
            idx = id
            dloc0 = dloc
    return idx


# =============================================================================
# Location Indexer
# =============================================================================
class LocIndexer(object):

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

        def _get_slice(key, i,):

            if not isinstance(key, slice):
                item = self.loc2index(key, i)
                newkey = slice(item, item + 1, None)

            else:
                start, stop, step = key.start, key.stop, key.step
                start = self.loc2index(start, i) if start is not None else start
                stop = self.loc2index(stop, i) + 1 if stop is not None else stop
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

if __name__ == '__main__':

    # check loc2index
    # create a simple 1D source
    import numpy as np

    x = np.arange(20.)*np.pi/50.
    y = np.sin(x)
    data = DataFrame([y], columns=x)
    print(DataFrame(data.columns))

    print("\nindex", loc2index(data, 1.13))
