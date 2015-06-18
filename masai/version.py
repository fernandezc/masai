# -*- coding: utf-8 -*-

"""
Version.py

This module contain version functions

"""

# Copyright (c) C. Fernandez' group @ LCS.
# Distributed under the terms of the CeCILL-B License.


from __future__ import print_function, division

# -----------------------------------------------------------------------------
# Globals
# -----------------------------------------------------------------------------

# stdlib
import os
import sys
import commands

# third-party

# -----------------------------------------------------------------------------
# Local utilities and constants
# -----------------------------------------------------------------------------

from masai.core.utils import Preferences
#TODO: this replace the apptools version while the deprecation warning is still here
#TODO: make a pull request to apptools

DEBUG = False
_TITLE_ = u"MASAI"
_COPYRIGHT_ = u"2015, C.Fernandez"

# set or read the preference file
_path = os.path.join(os.path.expanduser('~'), '.masai')
if not os.path.exists(_path):
    os.makedirs(_path)
_filename = os.path.join(_path, "masai.ini")
preferences = Preferences(filename=_filename)

# get/set the current path of the application
MASAI = preferences.get("api.masai", None)
if not MASAI:
    preferences.set("api.masai", os.getcwd())
    preferences.flush()

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
    txt = "-" * 32 + "\n"
    txt += " %s, version : %s \n Copyright (c)%s " % (
        _TITLE_, __version__, _COPYRIGHT_)
    txt += "\n" + "-" * 32 + "\n"
    return txt

#===============================================================================
# Main test
#===============================================================================
if __name__ == '__main__':
    print(get_fullversion_string())

# eof --------------------------------------------------------------------------
