# -*- coding: utf-8 -*-

"""
Define TG plugin functionalities

"""

# Copyright (c) C. Fernandez' group @ LCS.
# Distributed under the terms of the CeCILL-B License.

from __future__ import print_function

# -----------------------------------------------------------------------------
# Globals
# -----------------------------------------------------------------------------

# stdlib
import os
import glob
import warnings

# third-party
import pandas as pd
import numpy as np
import addict
from traits.api import *
import traits.has_traits

# -----------------------------------------------------------------------------
# Local utilities and constants
# -----------------------------------------------------------------------------

from masai.core.source import (Source)
from masai.core.database import Isotopes
from masai.core.utils import htmldoc, is_sequence
from masai.plugins import IPlugin

# constants
DEBUG = True
traits.has_traits.CHECK_INTERFACES = 2

# settings
warnings.simplefilter("ignore", DeprecationWarning)
warnings.simplefilter("ignore", FutureWarning)


# ===============================================================================
# Plugin signature
# ===============================================================================
__plugin_main__ = 'Tg'

# ===============================================================================
# Constant
# ===============================================================================


# ===================================================================================
# Tg plugin
# ===================================================================================
@provides(IPlugin)
class Tg(Source):
    """Tg (type: source), return a source from a AGIR dataset."""

    # minimal plugin attribute -----------------------------------------------------
    description = Str("Tg (type: source), return a source from a TG dataset.")
    type = Str("source")

    # specific attributes ----------------------------------------------------------

    silent = Bool(True)

    loaded = Bool
    _ = Any

    # *******************************************************************************
    # public functions
    # *******************************************************************************

    # *******************************************************************************
    # Properties
    # *******************************************************************************

    # *******************************************************************************
    # Private functions
    # *******************************************************************************
    def __init__(self, *args, **kargs):
        """
        Tg

        Usage: source = Tg(path)

        Parameters
        ----------
        path: str,
            path of the Tg file

        """
        path = ""
        if args:
            # look in args
            path = os.path.expanduser(args[0])

        if not os.path.exists(path):
            # nothing else to try
            return
            # raise Exception('a valid path to an existing file is needed in the arguments')

        self.path = path
        self._data_orig = self.read_tg(path, *kargs)

    def read_tg(self, path, **kargs):
        """
        reading tg dataset from CSV

        We expect the data in a subdirectory TGA under the path
        """
        filename = os.path.join(path, 'TGA', 'tg.csv')
        # read each csv file as a DataFrame and store in a dict
        tgs = pd.read_csv(filename, names=['time', 'tg'],
                        header=None, index_col=0,
                        delimiter=';')

        self.loaded = True
        return tgs

    # -------------------------------------------------------------------------------
    def __repr__(self):
        """
        private function: Create a html representation of some data in the Source traits

        """
        if not self.loaded:
            return self.__init__.__doc__

        txt = "***TG Source***\n"
        txt += "- path: %s\n" % self.path

        return txt

    __str__ = __repr__

    # -------------------------------------------------------------------------------
    def _repr_html_(self):
        """
        private function: Create a html representation of some data in the Source traits

        """
        if not self.loaded:
            return htmldoc(self.__init__.__doc__)

        txt = self.__repr__().replace('***\n', '</em>')
        txt = txt.replace('***', '<em>')
        txt = txt.replace('\n', '')
        txt = txt.replace('- ', '<li>')
        return txt

if __name__ == '__main__':

    USERPATH = os.path.expanduser(os.path.join('~','Dropbox'))
    DATADIR = os.path.join(USERPATH, 'storage', 'agir', 'data')

    names = ['P350','A350', 'B350', 'P400', 'A400', 'B400']
    sources = []
    for name in names:
        path = os.path.join(DATADIR, name)
        sources.append(Tg(path))
