# -*- coding: utf-8 -*-
#
# ===============================================================================
# masai.plugins.fits.fit1d_plugin
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
Define fit1d plugin functionalities
"""
from __future__ import print_function, division

# ===============================================================================
#  general imports
# ===============================================================================
# noinspection PyUnresolvedReferences
from traits.api import (HasTraits, Bool, Any, Str, Instance, provides)
import traits.has_traits
import re

# ===============================================================================
# local imports
# ===============================================================================
from masai.plugins import IPlugin
from masai.core.fit.parameters import ParameterScript

traits.has_traits.CHECK_INTERFACES = 2

# ===============================================================================
# Plugin signature
# ===============================================================================
__plugin_main__ = 'Fit1d'


# ===================================================================================
# Process1d plugin
# ===================================================================================
@provides(IPlugin)
class Fit1d(HasTraits):

    # minimal plugin attribute -----------------------------------------------------
    description = Str("Fit1d (type: fits), fit a 1D source with some models.")
    type = Str("fits")

    # specific attributes ----------------------------------------------------------

    silent = Bool(True)

    _ = Any

    script = """

#-----------------------------------------------------------
# syntax for parameters definition:
# name: value, low_bound,  high_bound
#  * for fixed parameters
#  $ for variable parameters
#  > for reference to a parameter in the COMMON block
#    (> is forbidden in the COMMON block)
# common block parameters should not have a _ in their names
#-----------------------------------------------------------
#

COMMON:
# common parameters ex.
# $ gwidth: 1.0, 0.0, none
  $ gratio: 0.5, 0.0, 1.0

MODEL: LINE_1
shape: voigtmodel
	$ ampl:  1.0, 0.0, none
	$ pos:   0.0, -100.0, 100.0
	> ratio: gratio
	$ width: 1.0, 0, 100

"""

    paramscript = Instance(ParameterScript)

    # *******************************************************************************
    # public functions
    # *******************************************************************************

    # -------------------------------------------------------------------------------
    def print_doc(self):
        return self.description

    # -------------------------------------------------------------------------------
    def help(self):
        return """Usage:
            fit1d(fd, script)
                where fd is a valid source and script must be provided.

            a default script can be generated using :
                Fit1d.default()
            """

    def default(self):
        print(self.script)

    # *******************************************************************************
    # Private functions
    # *******************************************************************************

    # -------------------------------------------------------------------------------
    def __call__(self, *args, **kargs):
        """Process a 1D source
        """
        if args:
            # look in args
            source = args[0]
            script = args[1]
        else:
            raise ValueError("No valid source specified\n%s"%self.help())

        self.paramscript = ParameterScript(fd=source)

if __name__ == '__main__':
    import os
    from masai.plugins.list import Bruker
    from masai.plugins.list import Process1d
    from masai.plugins.list import Plot1d

    Bruker = Bruker()
    Plot1d = Plot1d()
    Process1d = Process1d()
    Fit1d = Fit1d()

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
    #Process1d(fd, "ab --mode poly")
    #Plot1d(fd)

    script = """
#-----------------------------------------------------------
# syntax for parameters definition:
# name: value, low_bound,  high_bound
#  * for fixed parameters
#  $ for variable parameters
#  > for reference to a parameter in the COMMON block
#    (> is forbidden in the COMMON block)
# common block parameters should not have a _ in their names
#-----------------------------------------------------------
#

COMMON:
# common parameters ex.
# $ gwidth: 1.0, 0.0, none
  $ gratio: 0.5, 0.0, 1.0

MODEL: LINE_1
shape: voigtmodel
    * ampl:  1.0, 0.0, none
    $ pos:   4.5, -100.0, 100.0
    > ratio: gratio
    $ width: 4.0, 0, 100

"""

    Fit1d(fd, script)
