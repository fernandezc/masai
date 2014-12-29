#! /usr/bin/env python
#===============================================================================
# modelling.pulsar.pypulsar
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

""" docstring """

from .forpulsar import (parameters, reset, compute)
                                    #reset_parameters, reset_operators, reset_share, reset, 
                                    #diagonalize, compute)

import logging
logger = logging.getLogger()

#def WRITE_LOG(level=None,strg=None):
#    """
#    Print a message from the FORPULSAR module
#    Callback function called by Fortran
#    
#    Parameters
#    ----------
#    level: integer 
#         Level of the debugging message
#         "1":"critical",
#         "2":"error",
#         "3":"warning",
#         "4":"info",
#         "5":"log",
#         "6":"debug"
#    strg: string 
#        message to display
#   
#    """
#    strtype={
#         "1":"critical",
#         "2":"error",
#         "3":"warning",
#         "4":"info",
#         "5":"debug",   # because log does not work
#         "6":"debug"
#    }
#
#    try:
#        LogMethod = getattr(logger, "%s" % strtype[str(level)])
#        LogMethod('%s: %s'%(strtype[str(level)],str(strg)))
#    except:
#        pass
    
if __name__ == '__main__':
    # some tests
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.WARNING)

    import numpy as np

#    for I in np.arange(1.5, 4.56, 1.):
#        reset()
#        #diagonalize.dealloc()
#        parameters.protonfrequency = 400000000.
#        parameters.idealpulse= True #False #True
#        parameters.pulse=[[1, .000001, 20000, 0, 0, 0, 0, 0]]
#        parameters.delay=[[1, 0.0, 1]]
#        parameters.chemicalshift = [[1, 0, 0, 0, 0, 0, 0]]
#        parameters.quadrupole = [[1, 4e6, 0., 0, 0, 0]]
#        parameters.t2  = [[1, 1e6]]
#        parameters.spinningangle = 54.736
#        parameters.spinningspeed = 00000.
#        parameters.qfactor = 0.000001
#        parameters.rfstep = 4
#        parameters.npts = 8192
#        parameters.nsb = 0
#        parameters.accuracy = 16
#        parameters.sw = 5000000. #parameters.spinningspeed
#        parameters.sr = 0
#        parameters.nall = 1
#        parameters.rcph = 0.
#
#        parameters.nucleus = [[1,                #index
#                           I,                #spin
#                           100000000.0,      #larmor
#                           1.0]]           #abundance
#    
#        import time
#        start = time.time()
#        done = compute(WRITE_LOG)
#        end = time.time()
#        print end-start
#        
#        import pylab as pl
#        pl.plot(parameters.spec.real)
#        
#        s=np.zeros((2*I,))
#        for i,m in enumerate(np.arange(-I, I-.1, 1.)):
#            s[i] = I*(I+1.)-m*(m+1.)
#        
#        print I, s, s.sum(), parameters.spec.real.sum()
#        print [(2*I+1)/2,(2*I+1)/2-1], parameters.ros[(2*I+1)/2,(2*I+1)/2-1]
#        
#    pl.show()
    
    reset()
    
    parameters.protonfrequency = 400000000.
    parameters.idealpulse= True #False #True
    parameters.pulse=[[1, .0, 0., 0, 0, 0, 0, 0]]
    parameters.delay=[[1, 0.0, 1]]
    parameters.t2  = [[1, 1e6]]
    parameters.spinningangle = 54.736
    sw = 50000.
    parameters.spinningspeed = sw
    parameters.qfactor = 0.000001
    parameters.rfstep = 4
    npts = 8192
    parameters.npts = npts
    parameters.nsb = 0
    parameters.accuracy = 16
    parameters.sw = sw
    parameters.sr = 0
    parameters.nall = 1
    parameters.rcph = 0.

    import pylab as pl
    import time

    CA = lambda I : (I*(I+1.)-3./4.)/(2.*I*(2.*I-1))**2
    A = lambda I,cq, w0 : CA(I)*cq**2/w0
    cq2w0fromA = lambda I,A : A/((I*(I+1.)-3./4.)/(2.*I*(2.*I-1))**2)
    
    import os
    pathdb = os.path.join(os.path.expanduser('~'), '.masai','db')

    I=1.5
    w0=100.0e6
    iter = np.arange(5000,35000,10000)
    specall = np.zeros((iter.size, npts))
    start = time.time()
    for i, Ax in enumerate(iter):
        cq = (cq2w0fromA(I,Ax)*w0)**.5
        print cq
        qis = -(I*(I+1.)-3./4.)*(3.*cq/(2.*I*(2.*I-1)))**2/w0/30
        parameters.quadrupole = [[1, cq, 0., 0, 0, 0]]
        parameters.chemicalshift = [[1, -qis, 0, 0, 0, 0, 0]]
        parameters.nucleus = [[1,      #index
                               I,      #spin
                               w0,     #larmor
                               1.0]]   #abundance
        done = compute()
        specall[i]=parameters.spec.real
    np.save('res', specall)
    end = time.time()-start
    print 'time', end
  
    start = time.time()
    done = compute()
    end = time.time()-start
    print 'time direct', end
    
    start = time.time()
    specall = np.load('res.npy')
    spec = specall[-1]
    end = time.time()-start
    print 'time file', end
    
    #pl.plot(specall.T)
    #pl.show()
    
    x = np.arange(npts)-npts/2
    pl.plot(x,specall.T)
#    pl.show()
    
    x2= (iter[1]/iter[0])*x
    s = specall[0]/2.
    pl.plot(x2,s,':')
    pl.show()
#    pl.contour(specall)
#
#    #2D interpolation
#    a,b = specall.shape 
#    print a,b   
#    from scipy import ogrid, sin, mgrid, ndimage, array
#
#    x,y = np.ogrid[0:a:1,0:b:1]
#    print x, y
#    
#    newx,newy = np.mgrid[0:a:0.1,0:b:1]
#    print "newx", newx
#    print "newy", newy
#    
#    x0 = x[0,0]
#    y0 = y[0,0]
#    dx = x[1,0] - x0
#    dy = y[0,1] - y0
#    ivals = (newx - x0)/dx
#    jvals = (newy - y0)/dy
#    coords = array([ivals, jvals])
#    
#    print coords
#    
#    newf = ndimage.map_coordinates(specall, coords)
#    pl.contour(newf)
#    pl.show()
#     
   