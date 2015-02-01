# -*- coding: utf-8 -*-
#
# ===============================================================================
# masai.core.process.ft
# ===============================================================================
# Copyright (C) 2015 Christian Fernandez
#       Laboratoire Catalyse et Spectrochimie, Caen, France.
#       christian.fernandez@ensicaen.fr
# This software is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# You can  use, modify and/ or redistribute the software under
# the terms of the CeCILL-B license as circulated by CEA, CNRS and INRIA
# at the following URL "http://www.cecill.info".
# See Licence.txt in the main masai source directory
# ===============================================================================

from __future__ import print_function, division

#
# imports
#
# from traits.api import (HasTraits, Instance)
import numpy as np
import argparse
from .proc_base import largest_power_of_2, zf_size, fft
from .mqstmas import k_ratio

__all__ = ["ft", ]
#
# Fourier transform (ft)
#
def ft(source, options='', axis=-1):
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
    # options evaluation
    parser = argparse.ArgumentParser(description='FT processing.', usage = """
    ft [-h] [--si [SI]] [--tdeff [TDEFF]]
    """)
    # positional arguments
    parser.add_argument('--si',  nargs='?', help='size of fourier transform')
    parser.add_argument('--tdeff', default=0,  nargs='?', type=int, help='size of fourier transform')
    if axis==0:
        parser.add_argument('--shear', default=None, type=int, nargs=2, help = "p and q parameters for shearing")
        parser.add_argument('--shift', default=0, type=float, help = "shift in ppm unit")

    args  = parser.parse_args(options.split())
    if source is None:
        return
    
    if axis == -1:
        par = source.par
    else:
        par = source.par2

    if par.isfreq:
        print('This source is already transformed: ft skipped!')
        return

    # data on which we apply ft
    states = ('STATES' in par.encoding)
    tppi = ('TPPI' in par.encoding)
    if states:
        _makestates(source, tppi)

    data = source.data.values

    if axis == 0:
        # transpose temporarily the data for indirect dimension ft
        data = data.T


    # read parameters TDeff and SI
    tdeff = args.tdeff #kargs.get('tdeff', par.TDeff)
    if tdeff > par.TD or tdeff < 4:
        tdeff = par.TD

    # retain only TDeff points only (if TDeff if at least 2)

    if tdeff > 4:
        tdeff = min(tdeff, data.shape[-1])
        if data.ndim > 1:
            data = data[:, :tdeff]
        else:
            data = data[:tdeff]

    # zero filling (SI is taken at least taken
    # as the nearest power of two from TDeff)

    si = args.si
    if si and isinstance(si, str) and "k" in si:
        si = eval(si.replace("k", "*2**10"))

    if si is None:
        if par.SI is None:
            si = par.TD
        else:
            si = par.SI

    si = int(si)
    si = largest_power_of_2(max(si, par.TD))
    data = zf_size(data, si)

    # shearing and shift
    if axis == 0:
        shear = args.shear
        shift = args.shift
        if (shear is not None) or (shift is not None):
            sfo12 = source.par.SFO1
            sfo11 = source.par2.SFO1
            bf12 = source.par.BF1
            bf11 = source.par2.BF1
            sf2 = source.par.SF
            sf1 = source.par2.SF
            si2 = source.par.SI
            sw2 = source.par.SW_h
            sw1 = 1. / source.par.IN[0]
            source.par2.SW_h = sw1
            O1 = (sfo12-bf12)*1.e6
            SR = (sf2-bf12)*1.e6

            p0 = 0.
            p1 = 0.

            # phase for shearing and shift in F1
            p = 0
            if shear is not None:
                p, q = shear
                p = np.abs(p)
                q = np.abs(q)
                ratio = k_ratio(source, p, q)

                I = source.get_spin()
                n_spin = int(np.ceil(2.0 * I))
                if p==n_spin:
                    p = -p

                dra = abs(ratio - float(p))

                # for TROSY
                if p == 1 and q == 0:
                    dra = 1.0

                # for STMAS of inner satellite
                # if (p == 1 and q!=0): dra = abs(ratio - 1.0);

                sfo11 = sfo12 * dra
                bf11 = bf12 * dra
                sr2 = (sf2 - bf12) * 1.e6
                sr1 = sr2 * dra

                p0 = -2. * np.pi * float(shift) * bf11 / sw1

                bf11 += sr1 * 1.e-6 - float(shift) * bf11 * 1.e-6
                sf1 = bf11

                print("dra:", dra)

            par.SFO1 = sfo11
            par.BF1 = bf11
            par.SF = sf1
            par.O1 = (sfo11 - bf11) * 1.0e6

            # for MQMAS
            sign = 1
            if p < 0:
                 sign = -1

            # shearing and/or shift
            if shear is not None or shift != 0.0:
                # create a 2D array of phases and performs the phase correction
                twopi = 2. * np.pi
                w2 = source.data.columns # it's transposed
#                _, yy = np.meshgrid(range(si2), data.shape[-1])
#                w = np.array([w2[i] for i in yy])*bf12
                def func(i2, i1):
                    w = np.array([w2[i] for i in i2])*bf12
                    return np.exp((sign * twopi * ratio * (w - O1 +SR)/ sw1 - p0) * i1.astype(float) * 1.0j)
                phases = np.fromfunction(func,
                                    data.shape, dtype=int)


                data = data * phases

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


def _makestates(source, tppi=False):
    """
    Transforms according to the States or states-tppi procedure.

    """
    data = source.data.values

    ar = data[::2]  # transposed real part
    ai = data[1::2]  # transposed imag part
    data = ar.real - ai.real * 1.0j
    #data = ar.imag + ai.imag * 1.0j

    if tppi:
         data[1::2] = -data[1::2]

    source.par2.TD = data.shape[0]

    source.data = data
    source.history.append('Data transform according to encoding %s' % source.par2.encoding)
