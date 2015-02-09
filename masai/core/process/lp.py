# -*- coding: utf-8 -*-
#
# ===============================================================================
# masai.core.process.lp
# ===============================================================================
# Copyright (C) 2015 Christian Fernandez
# Laboratoire Catalyse et Spectrochimie, Caen, France.
# christian.fernandez@ensicaen.fr
# This software is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# You can  use, modify and/ or redistribute the software under
# the terms of the CeCILL-B license as circulated by CEA, CNRS and INRIA
# at the following URL "http://www.cecill.info".
# See Licence.txt in the main masai source directory
# ===============================================================================

import argparse

from masai.core.process import proc_lp


__all__ = ["lp"]
#
# undersampling
#
def lp(source=None, options='', axis=-1):
    """
    linear prediction

    Parameters
    ----------
    n: integer, optional
        number of undersampling, default 1

    """

    # options evaluation
    parser = argparse.ArgumentParser(description='Forward linear prediction.',
                                     usage="""
    lp [ncs|tcs] [-h]
    """, )
    # positional arguments
    parser.add_argument('npred', default=0, type=int, help='points to predict')
    parser.add_argument('--tdeff', default=0, nargs='?', type=int, help='size of the valid data')
    args = parser.parse_args(options.split())

    if axis == -1:
        par = source.par
    else:
        par = source.par2

    if par.isfreq:
        print('This source is already transformed: lp skipped!')
        return

    print('Forward linear prediction started')
    # data on which we apply lp
    npred = args.npred
    tdeff = args.tdeff
    if tdeff == 0:
        tdeff = par.TD

    datv = source.data.values
    if axis == 0:
        # transpose temporarily the data for indirect dimension ft
        datv = datv.T

    datv = proc_lp.lp(datv[:, :tdeff], npred)

    if datv.shape[-1] % 2 != 0:
        datv = datv[:, :-1]  # make it even

    par.TD = datv.shape[-1]
    par.TDeff = datv.shape[-1]

    # un-transpose the data if needed
    if axis == 0:
        datv = datv.T

    # change the current data with the transformed data
    source.data = datv  # scale will be calculated here

    source.history.append('Forward linear prediction of %d points' % npred)

