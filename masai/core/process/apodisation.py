from __future__ import print_function

from traits.api import (HasTraits, Instance)

import numpy as np

import nmrglue as ng
from nmrglue.process.proc_lp import lp as linear_predict


#
# em
#
def em(source, axis=-1, **kargs):
    """
    apodisation function multiplication

    Parameters
    ----------
    source: an instance of Source
    lb: float, optional
        Broadening factor, default is to read from the source parameters
    top: float, optional
        Position of the max of the apodisation function, default is None.
    method: string, optional, default is None
        if fullecho, top parameters is used to determine the maximum
        if k_shifted, coherence_order and satellite_order parameters are required to determine the maximums
    axis: optional, default is -1

    """
    source.history.append('em %s' % axis)

    lb = kargs.get('lb', source.fetchPar('LB', axis=axis))
    top = kargs.get('top', None)
    method = kargs.get('method', None)
    func = kargs.get('func', 'em_func')

    if method == 'fullecho' and (axis == -1 or axis == 1):

        if top is None:
            # id the parameter does not exist none should be return
            print("Fullecho but no 'top' arguments in em. 0.0 is used")
            top = 0.0

    elif method == 'k_shifted' and (axis == -1 or axis == 1) and source.is_2D:
        # in k shifted method for 2D spectra, the top of
        # the broadening function follow the top of the echoes
        # parameters p and q should be defined in parameters

        p = source.fetchPar('coherence_order') # coherence
        q = source.fetchPar('satellite_order') # satellite order

        if p is None or q is None:
            # id the parameter does not exist none should be return
            print("k_shifted method requires keys: 'coherence_order' and 'satellite_order' in parameters")
            return # abandon broadening

        top = np.abs(k_ratio(source, p, q))
        func = "k_shifted_%s" % func

    else:
        method = 'conventional'
        top = 0.0

    source.storePar('LB', lb, axis=axis)

    _apodization(source, func, lb, top, axis=axis)


