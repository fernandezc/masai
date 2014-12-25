# -*- coding: utf-8 -*-
#
#===============================================================================
# masai.plugins.__init__
#===============================================================================
# Copyright (C) 2015 Christian Fernandez
#       Laboratoire Catalyse et Spectrochimie, Caen, France.
#       christian.fernandez@ensicaen.fr
# This software is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# You can  use, modify and/ or redistribute the software under
# the terms of the CeCILL-B license as circulated by CEA, CNRS and INRIA
# at the following URL "http://www.cecill.info".
# See Licence.txt in the main masai source directory
#===============================================================================
"""
This module provide plugin functionalities and interface
"""
from __future__ import print_function, division

#===============================================================================
# Imports
#===============================================================================
import glob      # Imported for handling command line arguments.
import sys       # Imported for file directory listing capabilities.
import os        # Imported for file path manipulation capabilities.
from addict import Dict as aDict
from traits.api import HasTraits, Interface, provides, Str

import traits.has_traits
traits.has_traits.CHECK_INTERFACES = 2

#===============================================================================
# Contants
#===============================================================================
DEBUG = False

#===================================================================================
# Bruker plugin
#===================================================================================
class IPlugin(Interface):

    #===============================================================================
    # minimal arguments
    #===============================================================================
    description = Str
    type = Str

    #-------------------------------------------------------------------------------
    def help(self):
        """display help text"""

    #===============================================================================
    # # private functions
    #===============================================================================
    def __call__(self, *args, **kargs):
        """Callable class method"""
