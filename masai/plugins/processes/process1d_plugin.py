# -*- coding: utf-8 -*-
#
# ===============================================================================
# masai.plugins.processing.process1d_plugin
# ===============================================================================
# Copyright (C) 2015 Christian Fernandez
#       Laboratoire Catalyse et Spectrochimie, Caen, France.
#       christian.fernandez@ensicaen.fr
# This software is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# You can  use, modify and/ or redistribute the software under
# the terms of the CeCILL-B license as circulated by CEA, CNRS and INRIA
# at the following URL "http://www.cecill.info".
# See Licence.txt in the main masai source directory
# ===============================================================================
"""
Define process1d plugin functionalities
"""
from __future__ import print_function, division

# ===============================================================================
#  general imports
# ===============================================================================
# noinspection PyUnresolvedReferences
from traits.api import (HasTraits, Bool, Any, provides)
import traits.has_traits

# ===============================================================================
# local imports
# ===============================================================================
from masai.plugins import IPlugin
from masai.core.api import *

traits.has_traits.CHECK_INTERFACES = 2

# ===============================================================================
# Plugin signature
# ===============================================================================
__plugin_main__ = 'Process1d'


# ===================================================================================
# Process1d plugin
# ===================================================================================
@provides(IPlugin)
class Process1d(HasTraits):

    # minimal plugin attribute -----------------------------------------------------
    description = Str("Process1d (type: processes), output a transformed 1D source.")
    type = Str("processes")

    # specific attributes ----------------------------------------------------------

    silent = Bool(True)

    _ = Any

    # *******************************************************************************
    # public functions
    # *******************************************************************************

    # -------------------------------------------------------------------------------
    def print_doc(self):
        return self.description

    # -------------------------------------------------------------------------------
    def help(self):
        return self.__call__.__doc__

    # *******************************************************************************
    # Private functions
    # *******************************************************************************

    # -------------------------------------------------------------------------------
    def __call__(self, *args, **kargs):
        """Process a 1D source

        Parameters
        ----------
        fd: source (such as NMR or IR source)

        """
        # object to apply the transform
        if args:
            # look in args
            source = args[0]
        else:
            raise ValueError("no valid source specified\nUsage: Process1d(fd)")

        # type of transform, by default em+ft+apk
        transforms = kargs.get('transform', 'em 10Hz; zf 1; ft; apk').lower().split(';')

        def __getargs(txt, **k):
            z = trans.split()
            n = None
            if len(z) > 1:
                n = z[1]
            if n:
                k[txt] = n
            return k

        for trans in transforms:

            if 'em' in trans:
                kargs = __getargs('lb', **kargs)
                em(source, **kargs)
                source.history.append('Exponential apodization lb:%.2f' % source.par.LB)

            elif 'zf' in trans:
                kargs = __getargs('nzf', **kargs)
                zf(source, **kargs)
                source.history.append('Zero-filling to si:%d' % source.par.SI)

            elif 'ft' in trans:
                kargs = __getargs('si', **kargs)
                ft(source, **kargs)
                source.history.append('Fourier transform with si:%d' % source.par.SI)

            elif 'apk' in trans:
                kargs = __getargs('auto', **kargs)
                apk(source, **kargs)
                source.history.append('Auto-phasing phc0:%.2f phc1:%.2f' % (source.par.PHC0, source.par.PHC1))

            elif 'pk' in trans:
                #kargs = __getargs('auto', **kargs)
                pk(source, **kargs)
                source.history.append('Phasing phc0:%.2f phc1:%.2f' % (source.par.PHC0, source.par.PHC1))

if __name__ == '__main__':
    pass
