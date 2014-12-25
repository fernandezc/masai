# -*- coding: utf-8 -*-
#===============================================================================
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
import os
import numpy as np
import pytest
import pylab as pl

from lcslib.api import *

# global
DATADIR = os.path.expanduser(os.path.join('~','Dropbox','storage','bruker','data'))
PLOT = True


@pytest.fixture(scope='function')
def fdfid(request):
    """
    fixture returning a fid source object
    """
    user = 'christian'
    name = 'example1'
    expno = 10
    path = os.path.join(
        DATADIR, user, 'nmr', name, str(expno))
    fd = NMRSource(path)
    return fd

def test_process_fid(fdfid):
    fd = fdfid

    reset(fd)

    em(fd)
    ft(fd)
    pk(fd)

    reset(fd)

    undersampling(fd, 2)

    lb = 10.
    store(fd, 'LB', lb)
    store(fd, 'SI', 4096)

    efp(fd)

    #clipping
    clip(fd, -90, 100)

    # automatic phasing
    apk(fd, apkmode='negmin+entropy', basemode='linear', gamma=1., use_phc1=True)

    # automatic baseline correction
    base = autobase(fd, basetype='basf', polynom=True, apply=True)

    if PLOT:
        p = Plotting(fd)
        p.plot(showz=True) #, xlim=[-40,40], xreverse=True)

    #print "Kill plot to continue!"
    #pl.show()

@pytest.fixture(scope='function')
def fdser(request):
    """
    fixture returning a SER source object
    """
    user = 'aurelie'
    name = 'RbNO3'
    expno = 15
    path = os.path.join(
        DATADIR, user, 'nmr', name, str(expno))
    fd = NMRSource(path)
    return fd

def test_k_ft_ser(fdser):

    fd = fdser

    reset(fd)
    undersampling(fd, 2)

    store(fd, 'LB', 2.)

    coh = 1
    sat = 2
    fd.storePar('coherence_order', coh) # coherence for stmas
    fd.storePar('satellite_order', sat) # satellite order for stmas
    k = np.abs(k_ratio(fd, coh, sat))

    # fourier transform + broadening and phase
    # Only the first row is phased for 2D first processing
    em(fd, method='k_shifted')
    ft(fd)
    clip(fd, -66., .0)
    store(fd, 'PHC0',0)
    store(fd, 'PHC1',0)
    apk(fd, apk='first', apkmode='negmin+entropy', use_phc1=True, bound_phc1=10.)

def test_k_efp_ser(fdser):

    fd = fdser

    reset(fd)
    undersampling(fd, 2)

    coh = 1
    sat = 2
    store(fd, 'coherence_order', coh) # coherence for stmas
    store(fd, 'satellite_order', sat) # satellite order for stmas
    k = np.abs(k_ratio(fd, coh, sat))
    store(fd, 'LB', 20.)

    # fourier transform + broadening and phase
    # Only the first row is phased for 2D first processing
    efp(fd, method='k_shifted', apk='first')

    # plot k shifted broadened
    if PLOT:
        p = Plotting(fd)
        p.alpha = 0.5
        p.colormap = 'jet'
        ax = p.contour()
        ax.text(-30,0,'efp testing')

    print "kill plots to end"
    pl.show()