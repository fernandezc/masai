# -*- coding: utf-8 -*-

"""
This module is the API server.

"""

# Copyright (c) C. Fernandez' group @ LCS.
# Distributed under the terms of the CeCILL-B License.

from __future__ import print_function

# -----------------------------------------------------------------------------
# Globals
# -----------------------------------------------------------------------------

# stdlib
import os
import warnings
from functools import partial

# third-party
from IPython.core.magic import (Magics, magics_class, line_magic, cell_magic, line_cell_magic)
from IPython import get_ipython

# setup backend for ipython
ip = get_ipython()
if ip:
    # if were are running this under ipython
    ip.magic(u'matplotlib nbagg')


    #print("matplotlib and numpy namespaces are available under the name 'plt' and 'np' respectively")
    #print("matplotlib currently uses the 'nbagg' backend. \n\tYou can change it using e.g., plt.switch_backend('qt4agg') or %matplotlib inline")

from matplotlib import rcParams

# -----------------------------------------------------------------------------
# Local utilities and constants
# -----------------------------------------------------------------------------

import version as v

# constants
DEBUG = False

# settings
warnings.simplefilter("ignore", DeprecationWarning)
warnings.simplefilter("ignore", FutureWarning)
preferences = v.preferences
version = v.get_version()

# display the version number to the end user
print("MASAI api loaded: version %s" % version)

# display the current directory
print("\nThe current directory is <%s>" % os.getcwd())

# set access to some other directories
MASAI = preferences.get('api.masai')
print('MASAI root directory is <%s>' % MASAI)
myexemples = partial(os.path.join, MASAI, 'doc', 'tutorial', 'exemples')
myfigures = partial(os.path.join, MASAI, 'doc', 'tutorial', 'figures')

# plot settings
params = {
    'image.cmap': 'cubehelix',

    'figure.figsize'          : (3.4, 1.8),
    'figure.dpi'              : 150,

    'figure.facecolor'        : 'white',
    'figure.edgecolor'        : 'white',
    'figure.autolayout'       : True,
    'figure.max_open_warning' : 30,
    'figure.subplot.wspace'   :0.001,
    'figure.subplot.hspace'   :0.05,
    'figure.subplot.top'      :0.985,
    'figure.subplot.bottom'   :0.09,
    'figure.subplot.left'     :0.15,
    'figure.subplot.right'    :.96,

    'savefig.dpi'       : 150,
    'savefig.facecolor' : 'white',
    'savefig.edgecolor' : 'white',

    'font.family'       : 'serif',
    'font.serif'        : 'Times New Roman',
    'font.size'         : 10,
    'font.style'        : 'normal',
    'font.weight'       : 'normal',

    'legend.fontsize'      : 8,
    'legend.labelspacing'  : 0.2,
    'legend.handletextpad' : 0.1,
    'legend.columnspacing' : .5,

    'xtick.color'       : 'black',
    'ytick.color'       : 'black',

    'xtick.labelsize'   : 10,
    'ytick.labelsize'   : 10,

    'mathtext.fontset'  : 'stix',
    'mathtext.bf'       : 'bold',
    'mathtext.cal'      : 'cursive',
    'mathtext.default'  : 'rm',

    'text.usetex'       : False,
}
rcParams.update(params)


# load plugins
from masai.plugins.catalog import *

# ===============================================================================
# Main test
# ===============================================================================
if __name__ == '__main__':
    pass

# eof --------------------------------------------------------------------------
