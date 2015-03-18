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
import warnings

# third-party

# -----------------------------------------------------------------------------
# Local utilities and constants
# -----------------------------------------------------------------------------

import version as v

# constants
DEBUG = False

# settings
warnings.simplefilter("ignore", DeprecationWarning)
preferences = v.preferences
version = v.get_version()

# display the version number to the end user
print("MASAI, version %s" % version)

# ===============================================================================
# Main test
# ===============================================================================
if __name__ == '__main__':
    pass

# eof --------------------------------------------------------------------------
