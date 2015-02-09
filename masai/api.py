# -*- coding: utf-8 -*-
#
#===============================================================================
# masai.api
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
This module is the API server.
"""

# try to remove deprecation warning in preference librairies
import warnings
warnings.simplefilter("ignore",DeprecationWarning)

# import the version number and preferences
import version as v
preferences = v.preferences
version = v.get_version()

# display the version number to the end user
print("version %s"%version)

# load the catalog of plugins as modules that can be imported
from plugins.catalog import *

#===============================================================================
# Main test
#===============================================================================
if __name__ == '__main__':
    pass

# eof --------------------------------------------------------------------------
