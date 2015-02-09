# -*- coding: utf-8 -*-
#
# ===============================================================================
# masai.plugins.sources.multiplex_plugin
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

"""
Define multiplex plugin functionalities
"""

from __future__ import print_function

# ===============================================================================
# Imports
# ===============================================================================
import glob

import logging

logging.getLogger()

from warnings import filterwarnings

filterwarnings("ignore")

import pandas as pd
import numpy as np
import addict

from traits.api import *
from masai.core.source import read_multiplex
from masai.plugins.sources import Bruker

import traits.has_traits
from masai.plugins import IPlugin

traits.has_traits.CHECK_INTERFACES = 2
from masai.core.utils import htmldoc

# ===============================================================================
# Plugin signature
# ===============================================================================
__plugin_main__ = 'Multiplex'

# ===============================================================================
# Constant
# ===============================================================================
FnMODE = ["undefined", "QF", "QSEQ", "TPPI", "STATES",
          "STATES-TPPI", "ECHO-ANTIECHO"]


# ===================================================================================
# bruker plugin
# ===================================================================================
@provides(IPlugin)
class Multiplex(Bruker):
    # minimal plugin attribute -----------------------------------------------------
    description = Str("Multiplex (type: source), return a Multiplex source.")
    type = Str("source")

    # specific attributes ----------------------------------------------------------

    # *******************************************************************************
    # public functions
    # *******************************************************************************

    # *******************************************************************************
    # Private functions
    # *******************************************************************************

    # -------------------------------------------------------------------------------
    def __init__(self, *args, **kargs):
        """
        Multiplex

        Return a Bruker NMR multiplex source

        Parameters
        ----------
        path: str,
            path of the Bruker directory
        data_dir: str,
            main storage directory, optional
        user: str, optional
            user name of the dataset
        name: str, optional
            name of the dataset
        expno: int, optional
            experiment number
        processed: bool, optional, default is False
            should we load already bruker processed files
        procno: int
            processing number

        """

        # As this routine is subject to several call we need to restore the initial state
        if hasattr(self, 'par'):
            del self.par
            del self._data_orig
            del self._data
        if hasattr(self, 'par2'):
            del self.par2

        path = ''
        # noinspection PyUnusedLocal
        expnos = False
        data_dir = kargs.get('data_dir', '')
        user = kargs.get('user', '')
        name = kargs.get('name', '')
        # reading single expno
        expno = kargs.get('expno', None)

        if args:
            # look in args
            path = os.path.expanduser(args[0])

        if not os.path.exists(path):
            # look in kargs
            path = os.path.expanduser(kargs.get('path', ''))

        if not os.path.exists(path):
            # look alternative definitions
            path = os.path.join(str(data_dir), str(user), 'nmr', str(name), str(expno))

        if not os.path.exists(path):
            # look multiple
            expnos = kargs.get('expnos', False)
            if expnos:
                path = os.path.join(str(data_dir), str(user), 'nmr', str(name))

        if not os.path.exists(path):
            # nothing else to try
            # display help in the notebook
            return

        self._read_multiplex(path, **kargs)

    # ---------------------------------------------------------------------------------
    def _read_multiplex(self, path, **kargs):
        """
        read bruker data and parameters
        """
        # reset
        self._data = None
        self.par = None
        self.par2 = None
        self.datatype = None

        # internal funct
        def _get_par_files(_dir, _procno):
            # get all possible parameters files
            parfiles = []
            pdir = os.path.join(_dir, "pdata", _procno)

            la = glob.glob(os.path.join(_dir, "acqu*"))
            la = map(os.path.basename, la)
            for _item in la:
                parfiles.append(_item)

            lp = glob.glob(os.path.join(pdir, "proc*"))
            lp = map(os.path.basename, lp)
            for _item in lp:
                parfiles.append(os.path.join('pdata', _procno, _item))

            return parfiles

        # ------------------------

        if os.path.isdir(path) is not True:
            raise IOError("directory %s does not exist" % path)

        # Determine parameter automatically
        lowmem = kargs.get('lowmem', False)

        procno = kargs.get('procno', '1')
        datatype = None

        self.path = path

        # List of expnos?
        expnos = kargs.get('expnos', False)
        paths = None
        if expnos:
            paths = [os.path.join(path, str(expno)) for expno in expnos]
        else:
            paths = [path, ]

        # This is the final dictionary
        adic = addict.Dict()

        n = 0
        ldata = []
        lvarpars = []
        ldates = []

        for path in paths:

            # Acquisition parameters
            if os.path.isfile(os.path.join(path, "ser")):
                datatype = 'SER'
            else:
                raise IOError("No Bruker binary SER file could be found in %s" % path)

            par_files = _get_par_files(path, procno)  # we read all parameters file whatever the datatype

            dic, data = read_multiplex(path, acqus_files=par_files)
            pass
            # look the case whn the reshaping was not correct
            # for example, this happen when the number
            # of accumulated row was incomplete

            # if datatype in ['SER'] and data.ndim == 1:
            # # we must reshape using the acqu parameters
            #     td1 = dic['acqu2']['TD']
            #     try:
            #         data = data.reshape(td1,-1)
            #     except ValueError:
            #         try:
            #             td = dic['acqu']['TD']/2
            #             data = data.reshape(-1, td)
            #         except ValueError:
            #             print ("Inconsistency between TD's and data size")
            #             raise ValueError

            #reduce to td
            ntd = dic['acqus']['TD'] / 2
            data = data[..., :ntd]  # necessary for agreement with bruker data and phase

            # Clean dict for pdata keys
            for key in dic.keys():
                if key.startswith('pdata'):
                    newkey = key.split(os.path.sep)[-1]
                    dic[newkey] = dic.pop(key)

            # Eliminate the digital filter
            if kargs.get('dgfilter', True):
                data = self._remove_digital_filter(dic, data)

            if n > 0 and datatype != adic.datatype:
                raise ValueError("Experiments not compatibles. End")

            dic['datatype'] = datatype

            # we now make some rearrangement of the dic
            # to have something more user friendly
            # we assume that all experiments have similar (important) parameters so that the experiments are
            # compatibles

            # but we keep the original dictionary for further reference if necessary
            adic.fulldic['exp_%d' % n] = dic

            for key, val in dic.iteritems():
                if isinstance(val, dict):
                    dic[key] = addict.Dict(val)

            keys = sorted(dic.keys())
            for item in keys:

                newitem = item.split(os.path.sep)[-1]
                store = newitem
                if newitem[:4] in ['acqu', 'proc']:
                    store = 'par'
                    if len(newitem) > 4 and newitem[4] in ['2', '3']:
                        store += newitem[4]

                if newitem[-1] != 's':  # not a status parameters
                    if isinstance(dic[item], addict.Dict):
                        for skey in dic[item].keys():
                            adic[store][skey] = dic[item][skey]
                    else:
                        adic[store] = dic[item]
                else:
                    for skey in dic[item].keys():
                        adic[store][skey] = dic[item][skey]

            # correct some initial values
            adic.par.TD = data.shape[-1]
            adic.par.TDeff = adic.par.TD
            adic.par.units = 'ms'
            adic.par.FnMODE = 0
            adic.par.isfreq = False

            adic.par2.TD = data.shape[-2]
            adic.par2.TDeff = adic.par2.TD
            adic.par2.units = 'ms'
            adic.par2.isfreq = False
            adic.par2.REVERSE = False

            if not adic.par2.FnMODE:  # For historical reasons,
                # MC2 is interpreted when the acquisition status parameter FnMODE has the value undefined, i.e. 0
                adic.par2.FnMODE = int(adic.par2.MC2) + 1
            adic.par2.encoding = FnMODE[adic.par2.FnMODE]

        # normalised amplitudes to ns=1K and rg=1
        def _norm(dat):
            fac = float(adic.par.NS) * float(adic.par.RG)
            dat /= fac
            return dat

        data = np.conj(data * 1.0j)  # this transformation is to make data coherent with bruker processsing
        data = _norm(data)

        # make attributes from dict
        for key in adic.keys():
            setattr(self, key, adic[key])

        # we need to create a DataFrame with scales
        data = self._make_dataframe(data)

        # store temporarily
        ldata.append(data)

        # store also the varpars of the series
        varpars = kargs.get('varpars')
        if isinstance(varpars, str):
            # this should be one of the dic parameters
            lvarpars.append(adic.par[varpars])  # store the variable
            ldates.append(adic.par.DATE)  # and the date
        elif isinstance(varpars, list):
            # this should be a list of parameters
            p = []
            for var in varpars:
                p.append(adic.par[var])
            lvarpars.append(p)  # store the variable
            ldates.append(adic.par.DATE)  # and the date
        n += 1
        # ---------------------------------------

        if n > 1:
            # we need to concat the various arrays into a 2D array
            # but first check if they are all compatibles
            shapes = []
            nc = 0
            for item in ldata:
                shapes.append(item.shape)
                nc = max(item.shape[-1], nc)

            # try to make a pseudo 2D.
            data = pd.concat(ldata)
            if data.shape[-1] != nc:
                logging.warn('series of spectra not compatibles (scale different)')

            # set the index using the provided list of varpar
            # or just take the dates
            if data.shape[0] == len(lvarpars):
                data.index = lvarpars
            else:
                data.index = ldates

            data.fillna(0, inplace=True)  # replace NaN with zeros
            self.datatype = 'FIDSER'

            # correct some values
            self.par.TD = data.shape[-1]
            self.par.TDeff = self.par.TD
            self.par.units = 'ms'
            self.par.FnMODE = 0
            self.par.isfreq = False
            self.par2 = addict.Dict()
            self.par2.TD = data.shape[-2]
            self.par2.TDeff = self.par2.TD
            self.par2.label = varpars
            self.par2.isfreq = True
            self.par2.lvarpars = lvarpars

        # normalized to 1
        def _globnorm(dat):
            maxi = np.abs(dat.values).max()
            dat /= maxi
            self.norm = maxi
            return dat

        data = _globnorm(data)

        # store it as data_orig
        # we use the hidden properties because the associated properties is read_only
        self._data_orig = data

    def _repr_html_(self):
        if not self.par:
            return htmldoc(self.__init__.__doc__)
        else:
            return self.message


if __name__ == '__main__':
    """
    test plugin
    """
    # need this to test  directly
    import os
    from masai.plugins.catalog import Bruker
    from masai.plugins.catalog import Plot1d, Plot2d, Process

    Bruker = Bruker()
    Plot1d = Plot1d()
    Plot2d = Plot2d()
    Process = Process()

    Multiplex = Multiplex()

    print(__file__)
    DATADIR = os.path.expanduser(os.path.join('..', '..', '..', 'exemples'))
    FIGDIR = 'figures'  # this directory should exist in the current directory, however
    user = 'user'
    name = 'MULTIPLEX'
    expno = 4
    path = os.path.join(DATADIR, user, 'nmr', name, str(expno))
    source = Multiplex(path)
    # possibly correct the F1 encoding
    source.par2.encoding = 'STATES'

    # source.par.SF = source.par.BF1


    Process(source, transform='em 1 --k_shifted 3 0 --verbose')
    # Process(source, transform = 'us 3')
    # F2 dimension

    Process(source, transform='zf 1')  # 1 zero-filling
    Process(source, transform='ft')
    Process(source, transform='cl -40 40')  # clipping
    Process(source, transform='pk -72.558  0.214')  # 0 0 --auto --fit_phc1')

    # F1 dimension
    # Process(source, axis=0, transform='em .1')

    #Process(source, axis=0, transform='lp 128 --tdeff 106') # must be before any zf command
    Process(source, axis=0, transform='zf 1')
    Process(source, axis=0, transform='ft --shear 3 0 --shift 15')
    Process(source, axis=0, transform='pk -13.910 10.000')  #0 0 --auto --fit_phc1 --select cols')

    p_2 = Plot2d(source, figsize=(3, 2.7))  #, xlim=(-48, -20), ylim=(-40, -20) )