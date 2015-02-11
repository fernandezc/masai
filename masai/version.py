# -*- coding: utf-8 -*-
#
# ===============================================================================
# masai.version
# ===============================================================================
# Copyright (C) 2015 Christian Fernandez
# Laboratoire Catalyse et Spectrochimie, Caen, France.
#       christian.fernandez@ensicaen.fr
# This software is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# You can  use, modify and/ or redistribute the software under
# the terms of the CeCILL-B license as circulated by CEA, CNRS and INRIA
# at the following URL "http://www.cecill.info".
# See Licence.txt in the main masai source directory
# ===============================================================================
"""
This module contains version information. It use git to get the last version and
revision number.
"""
from __future__ import print_function, division

import os
import commands  #TODO: may be we could use Fabric

from masai.core.utils import Preferences
#TODO: this replace the apptools version while the deprecation warning is still here
#TODO: make a pull request to apptools

#===============================================================================
# Constants
#===============================================================================
_TOOLKIT_SETUP = False
_TITLE_ = u"MASAI"
_COPYRIGHT_ = u"2015, C.Fernandez"

# set or read the preference file
_path = os.path.join(os.path.expanduser('~'), '.masai')
if not os.path.exists(_path):
    os.makedirs(_path)
_filename = os.path.join(_path, "masai.ini")
preferences = Preferences(filename=_filename)

#===============================================================================
# Standard library imports
#===============================================================================

#===============================================================================
# toolkit setup
#===============================================================================
#TODO: not sure it is useful without  GUI interface but we'll see
if _TOOLKIT_SETUP:
    from traits.etsconfig.api import ETSConfig

    ETSConfig.toolkit = "qt4"
    import matplotlib

    matplotlib.use('Qt4Agg')  # We want matplotlib to use a QT backend

#===============================================================================
# get_git_version
#===============================================================================
def get_git_version():
    """
    Get the current version  of the local git repository using
    the 'git describe' shell command.

    Returns
    =======
    gitversion: str
        The Git version string
    """
    status, output = commands.getstatusoutput('git describe')
    if status == 0:
        return output
    else:
        return None


#===============================================================================
# get_version
#===============================================================================
def get_version():
    """
    Prepare the version string and write it in the preference file

    Returns
    =======
    version: str
        The version string
    """
    verstr = get_git_version()
    if verstr:
        verstr = verstr.split('-')
        # in case of a just tagged version version str contain only one string
        if len(verstr)>=2:
            verstr = "%s.%s" % tuple(verstr[:2])
        else:
            verstr = verstr[0]

    if verstr:
        preferences.set("api.version", verstr)
        preferences.flush()
    else:
        verstr = preferences.get("api.version", "unknown")
    return verstr


#===============================================================================
# get_git_version
#===============================================================================
def get_fullversion_string():
    """
    Returns
    =======
    fullversion: str
        A complete version string
    """
    __version__ = get_version()
    txt = "*" * 40 + "\n"
    txt += "%s, version : %s\nCopyright (c)%s" % (
        _TITLE_, __version__, _COPYRIGHT_)
    txt += "\n" + "*" * 40 + "\n"
    return txt

#===============================================================================
# Main test
#===============================================================================
if __name__ == '__main__':
    print(get_fullversion_string())

# eof --------------------------------------------------------------------------
