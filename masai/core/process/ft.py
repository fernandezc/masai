from __future__ import print_function, division

#
# imports
#
# from traits.api import (HasTraits, Instance)
import numpy as np
from nmrglue.process.proc_base import largest_power_of_2, zf_size, fft
# from nmrglue.process.proc_lp import lp as linear_predict

__all__ = ["ft", ]
#
# Fourier transform (ft)
#
def ft(source, axis=-1, **kargs):
    """
    Fourier transform along one axis (limited to 1D et 2D data)

    Parameters
    ----------
    axis: optional, default=-1
    tdeff: integer, optional
        effective TD
    si: integer, optional
        size of the transformed spectrum
    """
    if axis == -1:
        par = source.par
    else:
        par = source.par2

    if par.isfreq:
        print('This source is already transformed: ft skipped!')
        return

    # data on which we apply ft

    data = source.data.values
    if axis == 0:
        # transpose temporarily the data for indirect dimension ft
        data = data.T

    # read parameters TDeff and SI

    tdeff = kargs.get('tdeff', par.TDeff)
    if tdeff > par.TD or tdeff < 8:
        tdeff = par.TD

    # retain only TDeff points only (if TDeff if at least 2)

    if tdeff > 2:
        tdeff = min(tdeff, data.shape[-1])
        if data.ndim > 1:
            data = data[:, :tdeff]
        else:
            data = data[:tdeff]

    # zero filling (SI is taken at least taken
    # as the nearest power of two from TDeff)

    si = kargs.get("si", par.SI)
    if isinstance(si, str) and "k" in si:
        si = eval(si.replace("k", "*2**10"))

    si = largest_power_of_2(si)
    data = zf_size(data, si)

    # now perform the fft

    data = fft(data) #/float(si)

    # store the new pars

    par.TDeff = tdeff
    par.SI = data.shape[-1]
    par.REVERSE = False
    par.isfreq = True
    par.units = 'ppm'

    # un-transpose the data if needed

    if axis == 0:
        data = data.T

    # change the current data with the transformed data

    source.data = data  # scale will be calculated here

    source.history.append('Fourier transform with si:%d' % par.SI)


#
# efp
#
def efp(source, axis=-1, **kargs):
    """
    Group actions: em + ft + apk + baseline

    Parameters
    ----------
    axis: optional, default=-1

    """
    strg = 'efp'
    if kargs.get('xf2'):
        del kargs['xf2']
        strg = 'xf2'
    if kargs.get('debug', False):
        print('performing %s' % strg)

    em(source, axis=axis, **kargs)

    ft(source, axis=axis)

    apk(source, axis=axis, **kargs)


#
# xf2 - Alias of efp
#
def xf2(source, **kargs):
    """
    Performs efp along the direct axis

    see `efp' parameters

    """
    kargs['xf2'] = True
    efp(source, **kargs)


#
# xf1
#
def xf1(source, **kargs):
    """
    Fourier transform of the indirect dimension of a 2D NMR spectra

    """

    if kargs.get('debug', False): print('performing xf1')

    fnmode = source.fnmode
    tppi = True if 'TPPI' in fnmode else False
    ftppi = 0.5 if tppi else 1.

    if kargs.get('debug', False): print('... fnMODE: %s' % fnmode)

    if 'STATES' in fnmode:
        _makestates(source, tppi=tppi)
    else:
        # not yet implemented
        print("FnMODE(or MC2)= '%s' mode not yet implemented" % fnmode,
              file=sys.stderr)
        return

    sw = source.fetchPar('SW_h')
    sfo1 = source.fetchPar('SFO1')
    sfo11 = source.fetchPar('SFO1', axis=0)
    bf1 = source.fetchPar('BF1')
    bf11 = source.fetchPar('BF1', axis=0)
    sf = source.fetchPar('SF')
    sf1 = source.fetchPar('SF', axis=0)

    si = source.data.shape[-1]
    sw1 = 1. / source.fetchPar('IN')[0]
    source.storePar('SW_h', sw1, axis=0)

    p0 = 0.
    p1 = 0.
    # phase for shearing and shift in F1
    shearing = kargs.get('shearing', None)
    shift = kargs.get('shift', 0.0)

    if shearing is not None:
        p, q = shearing
        ratio = k_ratio(source, p, q)

        # shearing phase
        p1 = -2. * np.pi * ratio * sw / sw1 / float(si)

        dra = abs(ratio - p)

        # for MQMAS
        if p < 0:
            p1 = -p1

        # for TROSY
        if p == 1 and q == 0:
            dra = 1.0

        # for STMAS of inner satellite
        # if (p == 1 and q!=0): dra = abs(ratio - 1.0);

        sfo11 = sfo1 * dra
        bf11 = bf1 * dra
        sr = (sf - bf1) * 1.e6
        sr1 = sr * dra
        p0 = -2. * np.pi * shift * sfo11 / sw1
        bf11 += sr1 * 1.e-6 - shift * sfo11 * 1.e-6  # * ftppi
        sf1 = bf11

    source.storePar('SFO1', sfo11, axis=0)
    source.storePar('BF1', bf11, axis=0)
    source.storePar('SF', sf1, axis=0)
    source.storePar('O1', (sfo11 - bf11) * 1.0e6, axis=0)

    # shearing and/or shift
    if shearing is not None or shift != 0.0:
        # create a 2D array of phases and performs the phase correction
        data = source.data.T

        phases = np.fromfunction(lambda i0, i1:
                                 np.exp((p1 * (i0.astype(float) - float(si) / 2.) + p0) * i1.astype(float) * 1.0j),
                                 map(float, data.shape), dtype=float)
        data = data * phases
        source.data = data.T

    if kargs.get('lp1', True):
        # linear prediction
        lp(source, axis=0, **kargs)

    if kargs.get('lastft', True):
        # fourier transform
        ft(source, axis=0)

        # autophasing
        kargs['apk'] = 'max'

        apk(source, axis=0, **kargs)

    source.storePar('REVERSE', True, axis=0)