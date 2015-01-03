# -*- coding: utf-8 -*-
#
#===============================================================================
# masai.version
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
This module contains version information. It use git to get the last version and
revision number.
"""
from __future__ import print_function, division
import string

#===============================================================================
# Constants
#===============================================================================
_TOOLKIT_SETUP = False
_TITLE_ = u"MASAI"
_COPYRIGHT_ = u"2015, C.Fernandez"

#===============================================================================
# Standard library imports
#===============================================================================
import commands

#===============================================================================
# toolkit setup
#===============================================================================
if _TOOLKIT_SETUP:
    from traits.etsconfig.api import ETSConfig
    ETSConfig.toolkit = "qt4"
    import matplotlib
    matplotlib.use('Qt4Agg')  # We want matplotlib to use a QT backend

#===============================================================================
# get_git_version
#===============================================================================
def _get_git_version():
    """Get the current version  of the local git repository using
    the 'git describe' shell command.
    """
    status, output = commands.getstatusoutput('git describe')
    if status == 0:
        return output.split('-')


#===============================================================================
# get_git_version
#===============================================================================
def get_fullversion_string():

    x = _get_git_version()
    if len(x)==3:
        _VERSION_, _REVISION_, _GITID_ = x
    else:
        _VERSION_, _REVISION_ = string.join(x,''), 0

    version = "*"*80 + "\n"
    version += "%s\nversion : %s.%s\nCopyright (c)%s" % (_TITLE_, _VERSION_,
                                 _REVISION_, _COPYRIGHT_)
    version += "\n" + "*"*80 + "\n"
    return version

#===============================================================================
# Main test
#===============================================================================
if __name__ == '__main__':
    print(get_fullversion_string())

# eof --------------------------------------------------------------------------
