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
import re
import os

#===============================================================================
# Constants
#===============================================================================
_TOOLKIT_SETUP = False
_TITLE_ = u"MASAI"
_COPYRIGHT_ = u"2015, C.Fernandez"
VERSIONFILE=os.path.join("masai","__init__.py")
if not os.path.exists(VERSIONFILE):
    VERSIONFILE=os.path.join("__init__.py")

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
def get_git_version():
    """Get the current version  of the local git repository using
    the 'git describe' shell command.
    """
    status, output = commands.getstatusoutput('git describe')
    if status == 0:
        return output
    else:
        return None

#===============================================================================
# set_version
#===============================================================================
def set_version(version):

    with open(VERSIONFILE, "rt") as f:
        verstrline = f.read()
        VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]$"
        mo = re.search(VSRE, verstrline, re.M)
        res = re.sub(VSRE, "__version__ = \"%s\""%version, verstrline, flags=re.M )
        print(res)
    with open(VERSIONFILE, "wt") as f:
        f.write(res)

#===============================================================================
# get_version
#===============================================================================
def get_version():

    verstr = get_git_version()
    if verstr:
        verstr = verstr.split('-')
        verstr = "%s.%s"%tuple(verstr[:2])

    if verstr:
        set_version(verstr) # correct the file if needed
    else:
        # if verstr was none, read from file
        with open(VERSIONFILE, "rt") as f:
            verstrline = f.read()
            VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]$"
            mo = re.search(VSRE, verstrline, flags=re.M)
            if mo:
                verstr = mo.group(1)
            else:
                raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

    return verstr


#===============================================================================
# get_git_version
#===============================================================================
def get_fullversion_string():
    __version__ = get_version()
    txt=''
    #txt = "*"*40 + "\n"
    txt += "%s, version : %s\nCopyright (c)%s" % (_TITLE_, __version__, _COPYRIGHT_)
    #txt += "\n" + "*"*40 + "\n"
    return txt

#===============================================================================
# Main test
#===============================================================================
if __name__ == '__main__':
    print(get_fullversion_string())

# eof --------------------------------------------------------------------------
