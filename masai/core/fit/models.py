# -*- coding: utf-8 -*-
#
#===============================================================================
# masai.core.fit.models
#===============================================================================
# Copyright (C) 2012-2015 Christian Fernandez
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
This module holds the defintions all the various models used by Masaï.

"""
#===============================================================================
# logger
#===============================================================================
import logging
logger = logging.getLogger()
import os
import numpy as np
from numba import jit

#from masai.processing.processbase import _ft
from pulsar.pypulsar import (parameters, compute, reset)


basic_1D_baselines = [  # baseline models
                       'polynomialbaseline',
                     ]
basic_2D_baselines = [  # 2Dbaseline models
                       'commonrowpolynomialbaseline',
                     ]

basic_1D_models = [  # list of sequence independant models
                    'voigtmodel',
                    'gaussianmodel',
                    'lorentzianmodel',
                   ]

basic_2D_models = [ # list of sequence independant models
                    'gaussian2dmodel',
                    'poynomialseriesofvoigtmodel',
                   ]

quadrupolar_1D_models = [
                          'quad2model',
                          'quad2ssbmodel',
                          'quadallssbmodel',
                          # 'czjzekquad2model'
                        ]


def list_of_models(nd, spin=0):
    """Return the list of models available in each cases
    depending on the spin and the number of dimensions

     Parameters
     -----------
     nd : int
         Number of dimensions
     spin : float
         Spin quantum number

     """
    # basic models
    if nd == 1:
        models = basic_1D_models[:]
        if spin > 0.5:
            models.extend(quadrupolar_1D_models)
    elif nd == 2:
        models = basic_2D_models[:]
    else:
        raise NotImplementedError('not yet implemented')
    return models


def list_of_baselines(nd):
    """Return the list of baseline models available in each cases
    depending the number of dimensions

     Parameters
     -----------
     nd : int
         Number of dimensions

     """
    # basic models
    if nd == 1:
        baselines = basic_1D_baselines[:]
    elif nd == 2:
        baselines = basic_2D_baselines[:]
    else:
        raise NotImplementedError('not yet implemented')
    return baselines


############
###      ###
###  1D  ###
###      ###
############

#===============================================================================
# PolynomialBaseline
#===============================================================================
class polynomialbaseline(object):
    """
    Arbitrary-degree polynomial (degree limited to 10, however).
    As a linear baseline is automatically calculated, this polynom is always of
    greater or equal to order 2 (parabolic function).

    .. math::
        f(x) = ampl * \\sum_{i=2}^{max} c_i*x^i

    """
    args = ['ampl']
    args.extend(['c_%d' % i for i in range(2, 11)])

    script = """MODEL: baseline%(id)d\nshape: polynomialbaseline
    # This polynom starts at the order 2
    # as a linear baseline is additionnaly fitted automatically
    # parameters must be in the form c_i where i is an integer as shown below
    $ ampl: %(scale).3g, 0.0, None
    $ c_2: 1.0, None, None
    * c_3: 0.0, None, None
    * c_4: 0.0, None, None
    # etc...
    """

    def f(self, x, ampl, *c_, **kargs):
        c = [0.0, 0.0]
        c.extend(c_)
        return ampl * np.polyval(np.array(tuple(c))[::-1], x - x[x.size / 2])


#===============================================================================
# GaussianModel
#===============================================================================
class gaussianmodel(object):
    """
    Normalized 1D gaussian function:

    .. math::
        f(x) = \\frac{ampl}{\\sqrt{2 \\pi \\sigma^2} } \\exp({\\frac{-(x-pos)^2}{2 \\sigma^2}})

    where :math:`\\sigma = \\frac{width}{2.3548}`

    """
    args = ['ampl', 'width', 'pos']
    script = """MODEL: line%(id)d\nshape: gaussianmodel
    $ ampl: %(ampl).3f, 0.0, None
    $ width: %(width).3f, 0.0, None
    $ pos: %(pos).3f, %(poslb).3f, %(poshb).3f
    """

    def f(self, x, ampl, width, pos, **kargs):
        gb = width / 2.3548
        tsq = (x - pos) * 2 ** -0.5 / gb
        w = np.exp(-tsq * tsq) * (2 * np.pi) ** -0.5 / gb
        w = w * (x[1] - x[0])
        return ampl * w


#===============================================================================
# LorentzianModel
#===============================================================================
class lorentzianmodel(object):
    """
    A standard Lorentzian function (also known as the Cauchy distribution):

    .. math::
        f(x) = \\frac{ampl * \\lambda}{\\pi [(x-pos)^2+ \\lambda^2]}

    where :math:`\\lambda = \\frac{width}{2}`

    """
    args = ['ampl', 'width', 'pos']
    script = """MODEL: line%(id)d\nshape: lorentzianmodel
    $ ampl: %(ampl).3f, 0.0, None
    $ width: %(width).3f, 0.0, None
    $ pos: %(pos).3f, %(poslb).3f, %(poshb).3f
    """

    def f(self, x, ampl, width, pos, **kargs):
        lb = width / 2.
        w = lb / np.pi / (x * x - 2 * x * pos + pos * pos + lb * lb)
        w = w * (x[1] - x[0])
        return ampl * w


#===============================================================================
# VoigtModel
#===============================================================================
class voigtmodel(object):
    """
    A Voigt model constructed as the convolution of a :class:`GaussianModel` and
    a :class:`LorentzianModel` -- commonly used for spectral line fitting.

    """
    args = ['ampl', 'width', 'ratio', 'pos']
    script = """MODEL: line%(id)d\nshape: voigtmodel
    $ ampl: %(ampl).3f, 0.0, None
    $ width: %(width).3f, 0.0, None
    $ pos: %(pos).3f, %(poslb).3f, %(poshb).3f
    $ ratio: 0.1, 0.0, 1.0
    """

    def f(self, x, ampl, width, ratio, pos, **kargs):
        from scipy.special import wofz
        gb = ratio * width / 2.3548
        lb = (1. - ratio) * width / 2.
        if gb < 1.e-16:
            w = lorentzianmodel().f(x, 1., lb * 2., pos, **kargs)
        else:
            w = wofz(((x - pos) + 1.0j * lb) * 2 ** -0.5 / gb)
            w = w.real * (2.*np.pi) ** -0.5 / gb
            w = w * (x[1] - x[0])
        return ampl * w


#===============================================================================
# Quad2model
#===============================================================================
class quad2model(object):
    """Second order quadrupolar lineshape assuming infinite MAS spinning speed.
    (work also for second-order quadrupolar lineshape of static powder.

    The lineshape is calculated using Pulsar.

    """
    args = ['ampl', 'width', 'ratio', 'pos', 'cq', 'etaq']
    script = """MODEL: line%(id)d\nshape: quad2model
    $ ampl: %(ampl).3f, 0.0, None
    $ width: %(width).3f, 0.0, None
    $ pos: %(pos).3f, %(poslb).3f, %(poshb).3f
    * ratio: 0.1, 0.0, None
    $ cq: 1000000.0, 0.0, None
    * etaq: 0.1, 0.0, 1.0
    """
    def f(self, x, ampl, width, ratio, pos, cq, etaq,
                    vr=None, czjzek=False, all=False, **kargs):
        # evaluate all parameters
        reset()
        if kargs.has_key('sequence'):
            sequence = kargs['sequence']
            if sequence.type != 'ideal_pulse':
                pass  # logger.error('sorry non-ideal pulse not yet implemented')
        fd = kargs['fd']
        nucleus = fd.NUC1[-1]

        sw_h = float(fd.SW_h[-1])
        nsb = 0
        if vr is None:  # this is different from 0 (static)
            vr = sw_h
            nsb = 0
        elif vr > 0.0:
            nsb = int(sw_h / vr / 2.)
            # here we must limit the number of ssb in order to avoid a crackh of pulsar
            # for small spinning speed
            nsb = min(30, nsb)  # TODO: let the user choose the best max number of spinning sidebands
        fratio = 5.5856912e6 / nucleus.gamma
        sfo1 = fd.SFO1[-1]
        sf = fd.SF[-1]
        sr = (sfo1 - sf) * 1.e6
        protonf = fratio * sfo1
        si = fd.ydata.shape[-1]

        # set parameters for the pulsar calculation
        logging.disable(logging.INFO)  # to avoid multiple output during calculation
        parameters.protonfrequency = protonf
        parameters.nucleus = [[1,  # index
                           nucleus.spin,  # spin
                           sfo1 * 1.e6,  # larmor
                           nucleus.abundance]]  # abundance
        parameters.idealpulse = True
        parameters.pulse = [[1, 2, 0, 0, 0, 0, 0, 0]]
        parameters.delay = [[1, 0.0, 1]]
        pos = pos * sfo1 - sr
        parameters.chemicalshift = [[1, pos, 0, 0, 0, 0, 0]]
        parameters.t2 = [[1, 1e6]]
        parameters.spinningangle = 54.74
        parameters.spinningspeed = vr
        parameters.qfactor = 0.001
        parameters.rfstep = 5
        parameters.npts = si
        parameters.nsb = nsb
        parameters.accuracy = 8
        parameters.sw = sw_h
        parameters.sr = 0
        parameters.nall = 1
        if all:
            parameters.nall = 0
        parameters.rcph = 0.  # TODO: change this in the f90 compute function

        if not czjzek:
            parameters.quadrupole = [[1, cq, etaq, 0, 0, 0]]
            compute()
            data = parameters.spec.real
        else:
            parameters.accuracy = 3
            sigmaq = cq
            # look for the max of the distribution
            vq = np.linspace(0.0, 6.*sigmaq, 128)
            mp = np.max(pczjzeck(vq, .6, sigmaq))
            # now we will use a gauss legendre quadrature (may be not the best?)

        broad = voigtmodel().f(x, 1., width, ratio, sr / sfo1)  # centered
                                                # on the middle of the window!
        data = np.convolve(data, broad, 'same')
        logging.disable(logging.NOTSET)  # to restore normal logger level

        return ampl * data  # voir pour la normalisation


#===============================================================================
# quad2modelssb
#===============================================================================
class quad2ssbmodel(object):
    """Second order quadrupolar lineshape assuming finite MAS spinning speed.
    (spinning sidebands are calculated in this model)

    The lineshape is calculated using Pulsar.

    """

    args = ['ampl', 'width', 'ratio', 'pos', 'cq', 'etaq', 'vr']
    script = """MODEL: line%(id)d\nshape: quad2ssbmodel
    $ ampl: %(ampl).3f, 0.0, None
    $ width: %(width).3f, 0.0, None
    $ pos: %(pos).3f, %(poslb).3f, %(poshb).3f
    * ratio: 0.1, 0.0, None
    $ cq: 1000000.0, 0.0, None
    * etaq: 0.1, 0.0, 1.0
    * vr: 12000., 0.0, None
    """
    def f(self, x, ampl, width, ratio, pos, cq, etaq, vr, **kargs):
        return quad2model().f(x, ampl, width, ratio, pos, cq, etaq, vr, **kargs)


#===============================================================================
# quadallmodelssb
#===============================================================================
class quadallssbmodel(object):
    """All transition lineshape calculated at the first
    and second order quadrupolar
    assuming finite MAS spinning speed.
    It therefore include the spinning sidebands.

    It is calculated using Pulsar.

    """

    args = ['ampl', 'width', 'ratio', 'pos', 'cq', 'etaq', 'vr']
    script = """MODEL: line%(id)d\nshape: quadallssbmodel
    $ ampl: %(ampl).3f, 0.0, None
    $ width: %(width).3f, 0.0, None
    $ pos: %(pos).3f, %(poslb).3f, %(poshb).3f
    * ratio: 0.1, 0.0, None
    $ cq: 1000000.0, 0.0, None
    * etaq: 0.1, 0.0, 1.0
    * vr: 12000., 0.0, None
    """

    def f(self, x, ampl, width, ratio, pos, cq, etaq, vr, **kargs):
        return quad2model().f(x, ampl, width, ratio, pos, cq, etaq,
                                                          vr, all=True, **kargs)


#===============================================================================
# #===============================================================================
# # czjzekquad2model
# #===============================================================================
# class czjzekquad2model(quad2model):
#    args = ['ampl', 'width', 'ratio', 'pos', 'sigmaq']
#    script = """MODEL: line%(id)d\nshape: czjzekquad2model
#    * ampl: %(ampl).3f, 0.0, None
#    * width: %(width).3f, 0.0, None
#    * pos: %(pos).3f, %(poslb).3f, %(poshb).3f
#    * ratio: 0.0, 0.0, None
#    * sigmaq: 0.0, 0.0, None
#    """
#    def f(self, x, ampl, width, ratio, pos, sigmaq, **kargs):
#        return quad2model.f(self, x, ampl, width, ratio, pos, sigmaq,
#                                                      0.0, czjzek=True, **kargs)
#===============================================================================
def pczjzeck(vq, eta, sigma):
    """
    Czjzek probability distribution function [#]_:

    .. math::
       P(\\nu_Q, \\eta, \\sigma) = \\frac{1}{\\sqrt{2 \\pi} \\sigma^5 } \\nu_Q^4 \\eta
       \\left(1-\\frac{\\eta^2}{9} \\right) \\exp \\left[
       -\\frac{\\nu_Q^2 \\left( 1 + \\frac{\\eta^2}{3} \\right)}{2\\sigma^2}\\right]


    where

    .. math::
        \\nu_Q = \\frac{3 C_Q}{2 I (2 I - 1)}

    .. [#] Jean-Baptiste d'Espinose de Lacaillerie, Christian Fretigny, Dominique Massiot,
        MAS NMR spectra of quadrupolar nuclei in disordered solids: The Czjzek model,
        Journal of Magnetic Resonance, Volume 192, Issue 2, June 2008,
        Pages 244-251, ISSN 1090-7807, 10.1016/j.jmr.2008.03.001.
        (`<http://www.sciencedirect.com/science/article/pii/S1090780708000888>`_)

    """
    p = 1. / (np.sqrt(2.* np.pi) * sigma ** 5)
    p = p * vq ** 4 * eta * (1. - eta ** 2 / 9.)
    p = p * np.exp(-vq ** 2 * (1. + eta ** 2 / 3.) / 2. / sigma ** 2)
    return p


#===============================================================================
# simpsonquad2model
#===============================================================================
class simpsonquad2model(object):
    """Second order quadrupolar lineshape assuming MAS spinning speed
    equal to the spectral width.

    It is calculated using Simpson. (EXPERIMENTAL... need more developpement
    for a good integration in Masai)

    It requires that Simpson is installed separately and accessible from a shell
    in the the masai main directory [#]_.

    .. [#] Simpson can be downloaded here : `Simpson <http://bionmr.chem.au.dk/bionmr/software/simpson.php>`_

    """

    args = ['ampl', 'width', 'ratio', 'pos', 'cq', 'etaq']
    script = """MODEL: line%(id)d\nshape: simpsonquad2model
    * ampl: %(ampl).3f, 0.0, None
    * width: %(width).3f, 0.0, None
    * pos: %(pos).3f, %(poslb).3f, %(poshb).3f
    * ratio: 0.0, 0.0, None
    * cq: 1.0, 0.0, None
    * etaq: 0.0, 0.0, 1.0
    """

    def f(self, x, ampl, width, ratio, pos, cq, etaq, **kargs):
        gb = width * ratio
        lb = width * (1. - ratio)
        if kargs.has_key('sequence'):
            sequence = kargs['sequence']
            if sequence.type != 'ideal_pulse':
                pass  # logger.error('sorry non-ideal pulse not yet implemented')
        fd = kargs['fd']
        nucleus = fd.NUC1[-1]
        sw_h = fd.SW_h[-1]
        vr = sw_h
        fratio = 5.5856912e6 / nucleus.gamma
        sfo1 = fd.SFO1[-1]
        sf = fd.SF[-1]
        sr = (sfo1 - sf) * 1.e6 / sfo1  # en ppm
        protonf = fratio * sfo1
        si = fd.data.shape[-1]
        lb = lb * sfo1
        gb = gb * sfo1
        r = gb / (gb + lb)

        path = os.path.join(os.path.expanduser('~'), '.masai')
        temp_in = os.path.join(path, 'temp.in')
        temp_fid = os.path.join(path, 'temp.fid')

        # in this script the vr is made exactly the spectral width
        simpson_script = \
        """
        spinsys {
              channels %(symbol)s
              nuclei   %(symbol)s
              shift 1 %(pos).3fp 0 0 0 0 0
              quadrupole 1 2 %(cq).3f %(etaq).3f 0 0 0
        }

        par {
              proton_frequency %(protonf).0f
              method          gcompute
              gamma_angles    1
              crystal_file    rep2000
              start_operator  I1x
              detect_operator I1c
              spin_rate       %(vr).3f
              np              %(td)d
              sw              spin_rate
              variable tsw    1e6/sw
              verbose         0
        }

        proc pulseq {} {
          delay 9999
        }

        proc main {} {
          global par
          set f [fsimpson]
          faddlb $f %(lb)s %(r)s
          fsave $f %(path)s -raw_bin

        }

        """ % {'symbol':nucleus, 'cq':cq, 'vr':vr, 'etaq':etaq, 'protonf': protonf,
             'td':2048, 'si':si, 'lb':lb + gb, 'r':r, 'pos':pos - sr, 'path':temp_fid}

        with open(temp_in, 'wb') as f:
            f.write(simpson_script)
        import subprocess
        subprocess.check_output(['simpson', temp_in])

        data = read_raw_bin(temp_fid)
        try:
            intensity = np.absolute(data[0])
        except IndexError:
            raise IndexError("""Cannot find the first point of the fid.
            Simpson probably failed to execute the script""")
        np.save('temp', data / intensity)
        data = _ft(data / intensity, si).real
        return ampl * data  # voir pour la normalisation


#################
#####      ######
#####  2D  ######
#####      ######
#################

#===============================================================================
# commonrowPolynomialBaseline
#===============================================================================
class commonrowpolynomialbaseline(object):
    """
    Arbitrary-degree polynomial (degree limited to 10, however).
    As a linear baseline is automatically calculated, this polynom is always of
    greater or equal to order 2 (parabolic function).
    ****same for all rows***
    .. math::
        f(x) = ampl * \\sum_{i=2}^{max} c_i*x^i

    (Note that spectra are disposed in columns, evolution are along each row)
       t_0  t_1 .... t_f
       S0_0 S0_1 .... S0_f
       S1_0 S1_1 .... S1_f
       S2_0 S2_1 .... S2_f
       .
       .
       .
       Sf_0 Sf_1 .... Sf_f)
    """
    args = ['ampl']
    args.extend(['c_%d' % i for i in range(2, 11)])

    script = """MODEL: baseline%(id)d\nshape: commonrowpolynomialbaseline
    # This polynom starts at the order 2
    # as a linear baseline is additionnaly fitted automatically
    # parameters must be in the form c_i where i is an integer as shown below
    # ****same for all rows***
    $ ampl: %(scale).3g, 0.0, None
    $ c_2: 1.0, None, None
    * c_3: 0.0, None, None
    * c_4: 0.0, None, None
    # etc...
    """

    def f(self, x, y, ampl, *c, **kargs):
        base = polynomialbaseline().f(x, ampl, *c, **kargs)
        ret = np.ones((x.size, y.size))
        ret = (ret.T * base).T
        return ret


#===============================================================================
# poynomialseriesofVoigtModel
#===============================================================================
class poynomialseriesofvoigtmodel(object):
    """
    A Voigt model constructed as the convolution of a :class:`GaussianModel` and
    a :class:`LorentzianModel` -- commonly used for spectral line fitting.
    In this 2D class we assume a series of voigt model evolving following
    following a polynomial function (order limited to 20 -- should b sufficient)
    (Note that spectra are disposed in columns, evolution are along each row)
       t_0  t_1 .... t_f
       S0_0 S0_1 .... S0_f
       S1_0 S1_1 .... S1_f
       S2_0 S2_1 .... S2_f
       .
       .
       .
       Sf_0 Sf_1 .... Sf_f)
    """
    args = ['ampl', 'width', 'ratio', 'pos']
    args.extend(['c_%d' % i for i in range(20)])
    script = """MODEL: line%(id)d\nshape: poynomialseriesofvoigtmodel
    $ ampl: %(ampl).3f, 0.0, None
    $ width: %(width).3f, 0.0, None
    $ pos: %(pos).3f, %(poslb).3f, %(poshb).3f
    $ ratio: 0.1, 0.0, 1.0
    #evolution along series
    $ c_0: %(c0).3f, None, None
    $ c_1: %(c1).3f, None, None
    $ c_2: %(c2).3f, None, None
    $ c_3: %(c3).3f, None, None
    $ c_4: %(c4).3f, None, None
    $ c_5: %(c5).3f, None, None
    * c_6: 0.0, None, None
    # etc...
    """

    def f(self, x, y, ampl, width, ratio, pos, *c, **kargs):
        vm = voigtmodel().f(x, ampl, width, ratio, pos, **kargs)
        evol = np.polyval(np.array(tuple(c))[::-1], y)
        evoln = np.argwhere(evol < 0)  # concentrations cannot be negative!
        evol[evoln] = 0.0
        ret = np.ones((x.size, y.size)) * evol
        ret = ret / np.max(ret)
        ret = (ret.T * vm).T
        return ret

#===============================================================================
# 2DGaussianModel
#===============================================================================
class gaussian2dmodel(object):
    """
    Two dimensional Gaussian model (*not* normalized - peak value is 1).
    """
    args = ['ampl','width','width1','pos','pos1', 'theta']
    def f(self, x, y, ampl, width, width1, pos, pos1, theta=0.0, **kargs):
        sigma_x = width / 2.3548
        sigma_y = width1 / 2.3548
        theta = np.deg2rad(theta)
        
        #TODO: check the constistancy with the full hm
        x_, y_ = np.meshgrid(x, y)
        xo = float(pos)
        yo = float(pos1)
        a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
        b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
        c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
        w = np.exp( - (a*((x_-xo)**2) + 2*b*(x_-xo)*(y_-yo) + c*((y_-yo)**2)))
        if 'SQ-DQ' in kargs.get('kind'):
            xo = yo - xo # find the symmetric
            w += np.exp( - (a*((x_-xo)**2) + 2*b*(x_-xo)*(y_-yo) + c*((y_-yo)**2)))    # symmetrize pattern
        #w = w * np.abs(x[1] - x[0])
        return ampl * w

#===============================================================================
# Voigt2DModel
#===============================================================================
class voigt2Dmodel(object):
    #TODO: this is not a correct lineshape.... must be modified before using
    """
    A Voigt model constructed as the convolution of a :class:`GaussianModel` and
    a :class:`LorentzianModel` -- commonly used for spectral line fitting.

    """
    args = ['ampl', 'width', 'width1', 'ratio', 'pos', 'pos1']
    script = """MODEL: line%(id)d\nshape: voigt2Qmodel
    $ ampl: %(ampl).3f, 0.01, None
    $ width: %(width).3f, 0.01, None
    $ width1: %(width1).3f, 0.01, None
    $ pos: %(pos).3f, %(poslb).3f, %(poshb).3f
    $ pos1: %(pos1).3f, %(pos1lb).3f, %(pos1hb).3f
    $ ratio: 0.1, 0.0, 1.0
    """

    def f(self, x, y, ampl, width, width1, ratio, pos, pos1,  **kargs):
        """
        ['ampl', 'width', 'width1', 'ratio', 'pos', 'pos1']
        """
        from scipy.special import wofz
        gb = ratio * width / 2.3548
        lb = (1. - ratio) * width / 2.
        gb1 = ratio * width1 / 2.3548
        incx = (x[1] - x[0])
        incy = (y[1] - y[0])
        lb1 = (1. - ratio) * width1 / 2.
        x, y = np.meshgrid(x, y)

        if gb < 1.e-16:
            w = lorentzianmodel().f(x, 1., lb * 2., pos, **kargs)
            w1 = lorentzianmodel().f(y, 1., lb1 * 2., pos1, **kargs)
        else:
            w = wofz(((x - pos) + 1.0j * lb) * 2 ** -0.5 / gb)
            w = w.real * (2.*np.pi) ** -0.5 / gb
            w = w * incx
            w1 = wofz(((y - pos1) + 1.0j * lb1) * 2 ** -0.5 / gb1)
            w1 = w1.real * (2.*np.pi) ** -0.5 / gb1
            w1 = w1 * incy
        return ampl * w * w1


#################
###           ###
###  GENERAL  ###
###           ###
#################


#===============================================================================
# getmodel
#===============================================================================
def getmodel(x, y=None, modelname=None, par=None, **kargs):
    """Get the model for a given x vector.

    Parameters
    -----------
    x : ndarray
        Array of frequency where to evaluate the model values returned by the
        f function.
    y : ndarray or None
        None for 1D, or index for the second dimension
    modelname : str
        name of the model class to use.
    par : :class:`Parameters` instance
        parameter to pass to the f function
    kargs: any
        Keywords arguments to pass the the f function

    Returns
    -------
    ndarray : float
        an array containing the calculated model.

    """
    model = par.model[modelname]
    modelcls = globals()[model]

    # take an instance of the model
    a = modelcls()

    # get the parameters for the given model
    args = []
    for p in a.args:
        try:
            args.append(par['%s_%s' % (p, modelname)])
        except KeyError as e:
            if p.startswith('c_'):
                # probably the end of the list
                # due to a limited polynomial degree
                pass
            else:
                raise KeyError(e.message)

    x = np.array(x, dtype=np.float64)
    if y is None:
        res = a.f(x, *args, **kargs)
    else:
        y = np.array(y, dtype=np.float64)
        res = a.f(x, y, *args, **kargs)
    return res

#===============================================================================
# read_raw_bin
#===============================================================================
def read_raw_bin(filename):
    """
    Read a SIMPSON file wich has been saved using::

        fsave $f filename -raw_bin

    ..note::
        adapted from Jonathan J. Helmus (jjhelmus@gmail.com)
        Version: 0.1 (2012-04-13)
        License: GPL

    Returns
    --------
        ndarray : complex
            a 1D complex array

    """
    data = np.fromfile(filename, dtype='float32')
    half = int(data.shape[0] / 2)
    cdata = np.empty((half,), dtype='complex64')
    cdata.real = data[:half]
    cdata.imag = data[half:]
    return cdata


if __name__ == '__main__':
    import pylab as pl
    x = np.linspace(0, 10., 100.)
    y = np.linspace(0, 10., 100.)
    s = gaussian2dmodel().f(x, y, 1., 1., 2., 5.,  5., 45.)
    pl.contour(s)
    pl.show()
