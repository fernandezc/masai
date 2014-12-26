# -*- coding: utf-8 -*-
# ===============================================================================
# tests.testSource
#===============================================================================
# Copyright (C) 2012-2013 Christian Fernandez
#       Laboratoire Catalyse et Spectrochimie, Caen, France.  
#       christian.fernandez@ensicaen.fr
# This software is governed by the CeCILL-B license under French law 
# and abiding by the rules of distribution of free software.  
# You can  use, modify and/ or redistribute the software under 
# the terms of the CeCILL-B license as circulated by CEA, CNRS and INRIA
# at the following URL "http://www.cecill.info".
# See Licence.txt in the main masai source directory
#===============================================================================
import pytest

import pylab as pl
import os
import numpy as np

import warnings

warnings.filterwarnings("ignore")

from masai.plugins.list import Bruker
from masai.plugins.list import Plot1d
from masai.plugins.list import Plot2d
from masai.plugins.list import Process1d

Bruker = Bruker()
Plot1d = Plot1d()
Plot2d = Plot2d()
Process1d = Process1d()

# global
DATADIR = os.path.expanduser(os.path.join('~', 'Dropbox', 'storage', 'bruker', 'data'))

FULLTEST = False
TEST2D = False
PLOT = True


@pytest.fixture(scope='function')
def fd_fid(request):
    """
    fixture returning a fid source object
    """
    user = 'christian'
    name = 'example1'
    expno = 10
    path = os.path.join(DATADIR, user, 'nmr', name, str(expno))
    fd1 = Bruker(path)
    return fd1

@pytest.fixture(scope='function')
def fd_ser(request):
    """
    fixture returning a SER source object
    """
    user = 'aurelie'
    name = 'RbNO3'
    expno = 15
    path = os.path.join(DATADIR, user, 'nmr', name, str(expno))
    fd2= Bruker(path)
    return fd2

@pytest.mark.skipif(not FULLTEST, reason="bypassed test")
def test_open_fid(fd_fid):
    fd = fd_fid
    assert (fd.datatype == "FID")
    assert (fd.par.TD == 3120)

@pytest.mark.skipif(not FULLTEST, reason="bypassed test")
def test_open_ser(fd_ser):
    fd = fd_ser
    assert (fd.datatype == "SER")
    assert (fd.par.TD == 872)
    assert (fd.par2.TD == 128)


def test_process1d(fd_fid):
    fd = fd_fid
    Process1d(fd, transform='em 100; zf 0; ft; apk first negmin')
    Plot1d(fd)

### PLOTS #####

@pytest.mark.skipif(not FULLTEST or not PLOT, reason="bypassed test")
def test_plot1d(fd_fid):
    fd = fd_fid
    # plot a Source
    ax = Plot1d(fd, imag=True, xlim=(0, 4), ylim=(-1, 1), c='orange', hold=True)  # Hold suppress immediate output

@pytest.mark.skipif(not FULLTEST or not PLOT, reason="bypassed test")
def test_plot1d_dataframe(fd_fid):
    fd = fd_fid
    # plot a Dataframe
    Plot1d(fd.data, hold=True)

@pytest.mark.skipif(not FULLTEST or not PLOT or not TEST2D, reason="bypassed test")
def test_plot2d(fd_ser):
    fd = fd_ser
    # plot a 2D Source
    ax = Plot2d(fd, hold=True)  # Hold suppress immediate output

