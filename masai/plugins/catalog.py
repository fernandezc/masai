# -*- coding: utf-8 -*-
#
# ===============================================================================
# masai.plugins.api
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
import os
import glob
from addict import Dict as aDict
import sys

#from sources import *
#from plots import *
#from processes import *
#from fits import *
#from analyses import *
#from simulations import *

x= 10
#===============================================================================
# Load all plugins
#===============================================================================
def _loadPlugins():
    """
    load existing plugins
    """

    # find the path of the masai module

    import masai
    directory = os.path.join(os.path.dirname(masai.__file__),'plugins')

    plugins = glob.glob(os.path.join(directory,"*","*_plugin.py"))

    dic = aDict()

    plugins.sort()

    for file in plugins:

        # Extract just the first part of the .py file name.
        path, basename = os.path.split(file)
        name = basename.split('_plugin')[0]

        # Dynamically set the PYTHONPATH so the user doesn't have to. It assumes
        # the plugins are contained in subdirectories where the main file lives.
        sys.path.append(path)
        # Import the plugin module temporarily long enough to instantiate an object
        # which is stored in a globally accessible dictionary.

        tempModule = __import__(name+"_plugin")
        attrs = dir(tempModule)
        if '__plugin_main__' in attrs:
            attr = getattr(tempModule,'__plugin_main__')
            cls = getattr(tempModule, attr)
            if cls:
                globals()[name.capitalize()] = cls

_loadPlugins()