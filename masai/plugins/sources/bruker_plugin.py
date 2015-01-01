# -*- coding: utf-8 -*-
#
# ===============================================================================
# masai.plugins.sources.bruker_plugin
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
Define Bruker plugin functionalities
"""

from __future__ import print_function

# ===============================================================================
# Imports
# ===============================================================================
import os
import glob
from warnings import warn, filterwarnings

filterwarnings("ignore")

import pandas as pd
import numpy as np
import addict
import nmrglue as ng
from traits.api import *
from masai.core.source import Source
from masai.plugins import IPlugin
from masai.core.database import Isotopes

import traits.has_traits
traits.has_traits.CHECK_INTERFACES = 2

# ===============================================================================
# Plugin signature
# ===============================================================================
__plugin_main__ = 'Bruker'

# ===============================================================================
# Constant
# ===============================================================================
FnMODE = ["undefined", "QF", "QSEQ", "TPPI", "STATES",
          "STATES-TPPI", "ECHO-ANTIECHO"]


# ===================================================================================
# bruker plugin
# ===================================================================================
@provides(IPlugin)
class Bruker(Source):
    # minimal plugin attribute -----------------------------------------------------
    description = Str("Bruker (type: source), return a source from a Bruker path.")
    type = Str("source")

    # specific attributes ----------------------------------------------------------

    silent = Bool(True)

    par = Instance(addict.Dict)
    par2 = Instance(addict.Dict)

    # *******************************************************************************
    # public functions
    # *******************************************************************************

    # -------------------------------------------------------------------------------
    def help(self):
        return self.__call__.__doc__


    #------------------------------------
    def get_multiplicator(self, axis=-1):
        """
        Compute a unit multiplicator depending on the unit string, i.e. 'ms'--> 1.e3
        """
        if axis == -1:
            par = self.par
        else:
            par = self.par2
        unit = par.units
        if par.isfreq:
            print('this axis is already processed: multiplicator set to 1')
            return 1.
        elif unit == 'ms':
            return 1.e3
        elif unit == 'us':
            return 1.e6
        elif unit == 'sec':
            return 1.
        else:
            print('unknown unit scale')
            return 1.


    # *******************************************************************************
    # Properties
    # *******************************************************************************

    # Spin --------------------------------------------------------------------------
    def get_spin(self, axis=-1):
        if self.is_2d:
            nucleus = [self.par2.NUC2, self.par.NUC1][axis]
        else:
            nucleus = self.par.NUC1
        return Isotopes(nucleus).spin

    # is 2d -----------------------------------------------------------------------
    is_2d = Property

    def _get_is_2d(self):
        """
        property getter: Check if the data are 2D or not
        """
        if self._data is not None:
            return len(self._data.shape) == 2 and self._data.shape[-2] > 1
        else:
            # check parmode parameters
            if self.par.PARMODE == 1:
                return True
        return False

    # is 3d ------------------------------------------------------------------------
    is_3d = Property

    def _get_is_3d(self):
        """
        property getter: Check if the data are 2D or not
        """
        if self.data is not None:
            return len(self.data.shape) == 3 and self.data.shape[-3] > 1
        else:
            # check parmode parameters
            if self.par.PARMODE == 2:
                return True
        return False

    # *******************************************************************************
    # Private functions
    # *******************************************************************************

    # -------------------------------------------------------------------------------
    def __call__(self, *args, **kargs):
        """
        return the output of the bruker plugin: a valid nmr source experiment

        Parameters
        ----------
        path: str,
            path of the Bruker directory
        data_dir: str,
            main strorage directory, optional
        user: str, optional
            user name of the dataset
        name: str, optional
            name of the dataset
        expno: int, optional
            experiment number
        processed: bool, optioanl, default is False
            should we load already bruker processed files
        procno: int
            processing number

        """
        path = ''
        # noinspection PyUnusedLocal
        expnos = False
        data_dir = kargs.get('data_dir', '')
        user = kargs.get('user', '')
        name = kargs.get('name', '')
        # reading single expno
        expno = kargs.get('expno', '')

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
            raise Exception('a valid path to an existing file is needed in the arguments')

        self._read(path, **kargs)

        return self

    # ---------------------------------------------------------------------------------
    def _read(self, path, **kargs):
        """
        read bruker data and parameters
        """
        # reset
        self._data = None
        self.par = None
        self.par2 = None

        # internal funct
        def _get_par_files(_dir, _procno, _processed=False):
            # get all possible parameters files
            parfiles = []
            pdir = os.path.join(_dir, "pdata", _procno)
            la = glob.glob(os.path.join(_dir, "acqu*"))
            la = map(os.path.basename, la)
            for _item in la:
                if not _processed:
                    parfiles.append(_item)
                else:
                    parfiles.append(os.path.join('..', '..', _item))
            lp = glob.glob(os.path.join(pdir, "proc*"))
            lp = map(os.path.basename, lp)
            for _item in lp:
                if not _processed:
                    parfiles.append(os.path.join('pdata', _procno, _item))
                else:
                    parfiles.append(_item)
            if _processed:
                _dir = pdir
            return _dir, parfiles

        if os.path.isdir(path) is not True:
            raise IOError("directory %s does not exist" % path)

        # Determine parameter automatically
        processed = kargs.get('processed', False)
        lowmem = kargs.get('lowmem', False)

        procno = '1'
        p = path.split(os.path.sep)
        if 'pdata' in p:
            _pd = p.index('pdata')
            procno = p[_pd + 1]
            if not p[0]:
                p[0] = os.path.sep
            path = os.path.join(*p[:_pd])
            processed = True
        procno = kargs.get('procno', procno)
        datatype = None

        self.path = path

        # Acquisition parameters
        if os.path.isfile(os.path.join(path, "fid")):
            datatype = 'FID'
        elif os.path.isfile(os.path.join(path, "ser")):
            datatype = 'SER'
        else:
            if not processed:
                warn('No binary fid or ser found in %s.\nTry processed files...' % path)
                processed = True

        if os.path.isfile(os.path.join(path, 'pdata', procno, '1r')):
            if not datatype or processed:
                datatype = '1D'
        elif os.path.isfile(os.path.join(path, 'pdata', procno, '2rr')):
            if not datatype or processed:
                datatype = '2D'
        elif os.path.isfile(os.path.join(path, 'pdata', procno, '3rrr')):
            if not datatype or processed:
                datatype = '3D'
        else:
            if not datatype:
                raise IOError("No Bruker binary file could be found in %s" % path)
            else:
                warn("No processed Bruker binary file could be found in %s. Use fid's." % path)

        npath, par_files = _get_par_files(path, procno, processed)  # we read all parameters file whatever the datatype

        if datatype in ['FID', 'SER']:
            dic, data = ng.bruker.read(npath, acqus_files=par_files)
            # look the case whn the reshaping was not correct
            # for example, this happen when the number
            # of accumulated row was incomplete
            if datatype in ['SER'] and data.ndim == 1:
                # we must reshape using the acqu parameters
                td1 = dic['acqu2']['TD']
                data = data.reshape(td1,-1)
            #reduce to td
            ntd = dic['acqus']['TD']/2
            data = data[...,:ntd] # necessary for agreement with bruker data and phase
        else:
            dic, data = ng.bruker.read_pdata(npath, procs_files=par_files)

        dic['datatype'] = datatype

        # clean dict for pdata keys
        for key in dic.keys():
            if key.startswith('pdata'):
                newkey = key.split(os.path.sep)[-1]
                dic[newkey] = dic.pop(key)

        # eliminate the digital filter
        if datatype in ['FID', 'SER'] and kargs.get('dgfilter', True):
            data = self._remove_digital_filter(dic, data)

        # we now make some rearrangement of the dic
        # to have something more user friendly
        adic = addict.Dict()

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

        # normalised amplitudes to ns=1K and rg=1
        def _norm(dat):
            fac = float(adic.par.NS) * float(adic.par.RG)
            dat /= fac
            maxi = np.max(np.abs(dat))
            dat /= maxi
            adic.par.norm = maxi
            return dat

        data = np.conj(data * 1.0j) # this transformation is to make data coherent with bruker processsing
        data = _norm(data)

        # correct some initial values
        if datatype in ['FID', 'SER']:
            adic.par.TD = data.shape[-1]
            adic.par.TDeff = adic.par.TD
            adic.par.units = 'ms'
            adic.par.FnMODE = 0
            adic.par.isfreq = False

        if datatype in ['SER']:
            adic.par2.TD = data.shape[-2]
            adic.par2.TDeff = adic.par2.TD
            adic.par2.units = 'ms'
            adic.par2.isfreq = False
            if not adic.par2.FnMODE:  # For historical reasons,
                # MC2 is interpreted when the acquisition status parameter FnMODE has the value undefined, i.e. 0
                adic.par2.FnMODE = int(adic.par2.MC2) + 1
            adic.par2.encoding = FnMODE[adic.par2.FnMODE]

            if adic.par.PARMODE == 2:
                adic.par3.TD = data.shape[-3]
                adic.par3.TDeff = adic.par3.TD
                adic.par3.units = 'ms'
                adic.par3.isfreq = False
                if not adic.par3.FnMODE:
                    adic.par3.FnMODE = int(adic.par3.MC2) + 1
                adic.par3.encoding = FnMODE[adic.par3.FnMODE]

        if datatype in ['1D', '2D', '3D']:
            adic.par.SI = data.shape[-1]
            adic.par.units = 'ppm'
            adic.par.isfreq = True
            if datatype in ['2D', '3D']:
                adic.par2.SI = data.shape[-2]
                adic.par2.units = 'ppm'
                adic.par2.isfreq = True
                if datatype in ['3D']:
                    adic.par3.SI = data.shape[-3]
                    adic.par3.units = 'ppm'
                    adic.par2.isfreq = True

        # make attributes from dict
        for key in adic.keys():
            setattr(self, key, adic[key])

        # we need to create a DataFrame with scales
        data = self._make_dataframe(data)

        # store it as data_orig
        # we use the hidden properties because the associated properties is read_only
        self._data_orig = data

    # --------------------------------------------------------------------------------------
    def _make_dataframe(self, data):
        """
        private function: Transform an array of data to pandas dataframe
        """

        is2d = self.is_2d
        is3d = self.is_3d
        if isinstance(data, pd.DataFrame):
            data.columns = self._scale(self.par)
            if is2d:
                data.index = self._scale(self.par2)
            elif is3d:
                raise NotImplementedError('3D dataframe not yet implemented')
                # TODO: implement this with pd.Panels or multiindex DataFrame
            return data

        x = self._scale(self.par)
        if is2d:
            y = self._scale(self.par2)
            df = pd.DataFrame(data.astype(np.complex128), index=y, columns=x)
        elif is3d:
            raise NotImplementedError('3D dataframe not yet implemented')
            # TODO: implement this with pd.Panels or multiindex DataFrame
        else:
            if data.ndim == 1:
                s = pd.Series(data.astype(np.complex128), index=x)
                df = pd.DataFrame({'0': s}).T
            else:
                y = [0, ]
                df = pd.DataFrame(data.astype(np.complex128), index=y, columns=x)

        return df

    # --------------------------------------------------------------------------
    # noinspection PyUnresolvedReferences
    @staticmethod
    def _scale(dic, reverse=None):
        """
        private function: Compute scale for a given axis.
        """

        # import parameter to convert units
        sw_h = float(dic.SW_h)
        sfo1 = float(dic.SFO1)
        bf1 = float(dic.BF1)
        sf = float(dic.SF)
        si = float(dic.SI)
        td = float(dic.TD)
        sr = (sf - bf1) * 1.0e6
        o1 = (sfo1 - bf1) * 1.0e6

        # set the spectral parameters
        if dic.units.lower() in ['ppm', 'hz', 'khz', 'mhz']:
            # if the spectra is transformed we use the SI parameter
            uc = ng.fileiobase.unit_conversion(si, True, sw_h, bf1, -sr + o1)
        else:
            fact = 2.0 if dic.FnMODE in [3, 4, 5, 6] else 1.0
            uc = ng.fileiobase.unit_conversion(td, True, sw_h * fact, bf1, -sr + o1)
            # size, cplx, sw, obs, car

        scal = getattr(uc, '%s_scale' % dic.units.lower())
        if reverse is None:
            reverse = dic.REVERSE

        if reverse:
            return scal()[::-1]
        else:
            return scal()

    # ===========================================================================
    # remove_digital_filter
    # nmrglue modified Digital Filter Processing
    # ===========================================================================
    def _remove_digital_filter(self, dic, data):
        """
        Remove the digital filter from Bruker data.
        """
        if 'acqus' not in dic:
            raise ValueError("dictionary does not contain acqus parameters")

        if 'DECIM' not in dic['acqus']:
            raise ValueError("dictionary does not contain DECIM parameter")
        decim = dic['acqus']['DECIM']

        if 'DSPFVS' not in dic['acqus']:
            raise ValueError("dictionary does not contain DSPFVS parameter")
        dspfvs = dic['acqus']['DSPFVS']

        if 'GRPDLY' not in dic['acqus']:
            grpdly = 0
        else:
            grpdly = dic['acqus']['GRPDLY']

        data = self._rm_dig_filter(data, decim, dspfvs, grpdly)
        data = data[..., :dic['acqus']['TD'] / 2]

        return data

    # -------------------------------------------------------------------------------
    @staticmethod
    def _rm_dig_filter(data, decim, dspfvs, grpdly=0):
        """
        Remove the digital filter from Bruker data.
        modified from nmrglue
        """
        #
        if grpdly > 0:  # use group delay value if provided (not 0 or -1)
            phase = grpdly

        # Determine the phase correction
        else:
            if dspfvs >= 14:  # DSPFVS greater than 14 give no phase correction.
                phase = 0.
            else:  # loop up the phase in the table
                if dspfvs not in ng.bruker.bruker_dsp_table:
                    raise ValueError("dspfvs not in lookup table")
                if decim not in ng.bruker.bruker_dsp_table[dspfvs]:
                    raise ValueError("decim not in lookup table")
                phase = ng.bruker.bruker_dsp_table[dspfvs][decim]

        # fft
        si = data.shape[-1]
        pdata = np.fft.fftshift(np.fft.fft(data, si, axis=-1), -1) / float(si / 2)
        pdata = (pdata.T - pdata.T[0]).T  # this allow generally to remove Bruker smiles

        # Phasing
        si = float(pdata.shape[-1])
        ph = 2.0j * np.pi * phase * np.arange(si) / si
        pdata = pdata * np.exp(ph)

        # ifft
        ndata = np.fft.ifft(np.fft.ifftshift(pdata, -1), si, axis=-1) * float(si / 2)

        # remove last points * 2
        rp = 4*(phase//2)
        pdata = ndata[..., :-rp]

        return pdata

    # -------------------------------------------------------------------------------
    def __repr__(self):
        """
        private function: Create a html representation of some data in the Source traits

        """
        if not self.par:
            txt = """Usage:
            fd = Bruker(path, [processed=False], [lowmem=False])
              where 'path' is a valid path to a directory
            """
            return txt

        txt = "***Bruker Source***\n"
        txt += "- path: %s\n" % self.path
        txt += "- datatype: %dD NMR %s\n" % (self.par.PARMODE + 1, self.datatype)
        if self.is_2d:
            if self.datatype in ['SER']:
                txt += "- Encoding: %s\n" % self.par2.encoding
                txt += "- TD1: %d, TD2: %d\n" % (self.par2.TD, self.par.TD)
            else:
                txt += "- SI1: %d, SI2: %d\n" % (self.par2.SI, self.par.SI)
            txt += "- F1 limits: %.2f to %.3f" % tuple(self.get_limits(axis=0))
            txt += " %s\n" % self.par2.units
            txt += "- F2 limits: %.2f to %.3f" % tuple(self.get_limits())
            txt += " %s\n" % self.par.units
        else:
            if self.datatype in ['FID']:
                txt += "- TD: %d\n" % self.par.TD
            else:
                txt += "- SI: %d\n" % self.par.SI
            txt += "- Limits: %.2f to %.3f" % tuple(self.get_limits())
            txt += " %s\n" % self.par.units
        return txt

    __str__ = __repr__

    # -------------------------------------------------------------------------------
    def _repr_html_(self):
        """
        private function: Create a html representation of some data in the Source traits

        """
        if not self.par:
            txt = """Usage:<br/>
            fd = <strong>Bruker</strong>(path, <em>[processed=False]</em>, <em>[lowmem=False</em>])<br/>
              where <em>'path'</em> is a valid path to a bruker dataset directory
            """
            return txt

        txt = self.__repr__().replace('***\n', '</H6>')
        txt = txt.replace('***', '<H6>')
        txt = txt.replace('\n', '')
        txt = txt.replace('- ', '<li>')
        return txt

    # *******************************************************************************
    # Events
    # *******************************************************************************
    # -------------------------------------------------------------------------------
    def _anytrait_changed(self, name, old, new):

        if 'history' in name:
            print(self.history[-1])

        elif '_data'  in name:
            if name == '_data_orig':
                self.history = ['Loaded %s'%self.datatype]

            # set labels
            if not self.par.isfreq:
                self.par.units = 'ms'
                self.par.label = r'$\mathrm{time\,(%s)}$'%self.par.units
            else:
                self.par.units = 'ppm'
                self.par.label = r'$\mathrm{\delta\,(%s)}$'%self.par.units

            if self.is_2d or self.is_3d:
                if self.par2.units in ['us', 'ms', 's']:
                    self.par2.label = r'$\mathrm{time\,(%s)}$'%self.par.units
                else:
                    self.par2.label = r'$\mathrm{\delta\,(%s)}$'%self.par.units

        return

if __name__ == '__main__':
    """
    test plugin
    """
    # need this to test  directly
    import pylab as pl

    import warnings
    warnings.filterwarnings("ignore")

    Bruker = Bruker()

    DATADIR = os.path.expanduser(os.path.join('~','Dropbox','storage','bruker','data'))

    user = 'christian'
    name = 'example1'
    expno = 10
    path = os.path.join(DATADIR,user,'nmr',name, str(expno))
    fd = Bruker(path, dgfilter=True)
    fd.data.T.plot(legend=False)

    user = 'aurelie'
    name = 'RbNO3'
    expno = 15
    path = os.path.join(DATADIR,user,'nmr',name, str(expno))
    fd = Bruker(path, dgfilter=True)

    row = fd.row(0, byindex=False)
    row.T.plot(legend=False)

    col = fd.col(0, width=.1, byindex=False)
    col.T.plot(legend=False)

    pl.show()