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
            if 'zf' in trans:
                kargs = __getargs('nzf', **kargs)
                zf(source, **kargs)

            elif 'ft' in trans:
                kargs = __getargs('si', **kargs)
                ft(source, **kargs)
                source.history.append('Fourier transform si:%d' % source.par.SI)

if __name__ == '__main__':
    """
    test process1d plugin
    """
    # need this to test  directly
    import pylab as pl
    import os
    from masai.plugins.list import Bruker
    from masai.plugins.list import Plot1d

    Bruker = Bruker()   # necessary for direct call
    Plot1d = Plot1d()
    Process1d = Process1d()

    # ******

    DATADIR = os.path.expanduser(os.path.join('~', 'Dropbox', 'storage', 'bruker', 'data'))
    user = 'christian'
    name = 'example1'
    expno = 10
    path = os.path.join(DATADIR, user, 'nmr', name, str(expno))
    fd = Bruker(path)

    # test plot a Dataframe
    # Plot1d(fd.data)

    # transform
    Process1d(fd, transform='em; zf 2; ft; apk')

    # plot a Source
    ax = Plot1d(fd)  # ,  hold=True)
    # ax = Plot1d(fd, ax = ax, imag=True, c='orange')

    pl.show()
