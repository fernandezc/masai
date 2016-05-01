# -*- coding: utf-8 -*-
#
# =============================================================================
# Copyright (©) 2015-2016 LCS
# Laboratoire Catalyse et Spectrochimie, Caen, France.
#
# This software is a computer program whose purpose is to [describe
# functionalities and technical features of your software].
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software. You can use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty and the software's author, the holder of the
# economic rights, and the successive licensors have only limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading, using, modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean that it is complicated to manipulate, and that also
# therefore means that it is reserved for developers and experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and, more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
# =============================================================================
from __future__ import print_function
# -----------------------------------------------------------------------------
# Local utilities and constants
# -----------------------------------------------------------------------------
import os
from masai.api import *

DATADIR = os.path.join(MASAI, 'doc','tutorial','exemples')
assert DATADIR=="/Users/christian/Dropbox/PycharmProjects/masai/doc/tutorial/exemples"

user = 'user'
name = 'HMQC'
expno = 31
path = os.path.join(DATADIR, user,'nmr',name, str(expno))
source = Bruker(path)
print(source)

p_2 = Plot2d(source, figsize=(3, 2.7), xlim=(0, 10), ylim=(0, 5),
             nc=30, exponent=1.5, start=-1.)

d6 = source.par.D[6]*1e3 # in ms
Process(source, transform = 'em 100 --shifted %.3f'%d6)

p_2 = Plot2d(source, figsize=(3,2.7), xlim=(0, 10), ylim=(0, 5),
             nc=30, exponent=1.5, start=-1.)

pass
