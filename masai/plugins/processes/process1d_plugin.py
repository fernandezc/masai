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
import re

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

        # type of transform, by default em+zf ft+apk
        if len(args)>1:
            transforms = args[1]
        else:
            transforms = kargs.get('transform', '-h; em 10; zf; ft; pk --auto')
        transforms = transforms.lower().split(';')

        for trans in transforms:

            p = re.compile('^(\w*)\s(.*)')
            m = p.match(trans.strip())
            if m :
                options = m.group(2)
            else:
                options = ''

            if 'em' in trans:
                em(source, options=trans)

            elif 'zf' in trans:
                zf(source, options=options)

            elif 'ft' in trans:
                ft(source, options=options)

            elif 'pk' in trans:
                pk(source, options=options)

            elif 'ab' in trans:
                ab(source, options=options)

if __name__ == '__main__':
    import os
    from masai.plugins.list import Bruker
    from masai.plugins.list import Plot1d

    Bruker = Bruker()
    Plot1d = Plot1d()
    Process1d = Process1d()

    # global
    DATADIR = os.path.expanduser(os.path.join('~', 'Dropbox', 'storage', 'bruker', 'data'))
    user = 'christian'
    name = 'example1'
    expno = 10
    path = os.path.join(DATADIR, user, 'nmr', name, str(expno))
    fd = Bruker(path)
    Process1d(fd, "em 50")
    Process1d(fd, "zf 0")
    Process1d(fd, "ft")
    Process1d(fd, "pk -25 0 --auto --use_phc1 --ediff 2 --mode negmin+entropy --gamma 0.0001 --bound_phc1 200 --verbose")
    Process1d(fd, "ab --mode poly")
    Plot1d(fd)
