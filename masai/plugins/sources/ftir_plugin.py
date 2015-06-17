# -*- coding: utf-8 -*-

"""
Define FTIR plugin functionalities

"""

# Copyright (c) C. Fernandez' group @ LCS.
# Distributed under the terms of the CeCILL-B License.

from __future__ import print_function

# -----------------------------------------------------------------------------
# Globals
# -----------------------------------------------------------------------------

# stdlib
import os
import glob
import warnings

# third-party
import pandas as pd
import numpy as np
from traits.api import *
import traits.has_traits

# -----------------------------------------------------------------------------
# Local utilities and constants
# -----------------------------------------------------------------------------

from masai.core.source import (Source)
from masai.core.database import Isotopes
from masai.core.utils import htmldoc, addict
from masai.plugins import IPlugin

# constants
DEBUG = True
traits.has_traits.CHECK_INTERFACES = 2

# settings
warnings.simplefilter("ignore", DeprecationWarning)
warnings.simplefilter("ignore", FutureWarning)


# ===============================================================================
# Plugin signature
# ===============================================================================
__plugin_main__ = 'Ftir'

# ===============================================================================
# Constant
# ===============================================================================


# ===================================================================================
# Ftir plugin
# ===================================================================================
@provides(IPlugin)
class Ftir(Source):
    """Ftir (type: source), return a source from a AGIR dataset."""

    # minimal plugin attribute -----------------------------------------------------
    description = Str("Ftir (type: source), return a source from a FTIR dataset.")
    type = Str("source")

    # specific attributes ----------------------------------------------------------

    silent = Bool(True)

    par = Instance(addict.Dict)
    par2 = Instance(addict.Dict)

    statesapplied = Bool
    _ = Any

    # *******************************************************************************
    # public functions
    # *******************************************************************************

    # *******************************************************************************
    # Properties
    # *******************************************************************************

    # *******************************************************************************
    # Private functions
    # *******************************************************************************
    def __init__(self, *args, **kargs):
        """
        Ftir

        Usage: source = Ftir(path, )

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

        # As this routine is subject to several call we need to retore the initial state
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
            return
            # raise Exception('a valid path to an existing file is needed in the arguments')

        self._read(path, **kargs)

    # ---------------------------------------------------------------------------------
    def _read(self, path, **kargs):
        """
        read bruker data and parameters
        """
        # reset
        self._data = None
        self.par = None
        self.par2 = None
        self.datatype = None

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

        # ------------------------

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
                elif processed:
                    logging.warn("No processed Bruker binary file could be found in %s. Use fid's." % path)

            npath, par_files = _get_par_files(path, procno,
                                              processed)  # we read all parameters file whatever the datatype

            if datatype in ['FID', 'SER']:
                dic, data = read(npath, acqus_files=par_files)
                # look the case whn the reshaping was not correct
                # for example, this happen when the number
                # of accumulated row was incomplete
                if datatype in ['SER'] and data.ndim == 1:
                    # we must reshape using the acqu parameters
                    td1 = dic['acqu2']['TD']
                    try:
                        data = data.reshape(td1, -1)
                    except ValueError:
                        try:
                            td = dic['acqu']['TD'] / 2
                            data = data.reshape(-1, td)
                        except ValueError:
                            print("Inconsistency between TD's and data size")
                            raise ValueError

                #reduce to td
                ntd = dic['acqus']['TD'] / 2
                data = data[..., :ntd]  # necessary for agreement with bruker data and phase
            else:
                dic, data = read_pdata(npath, procs_files=par_files)

            # Clean dict for pdata keys
            for key in dic.keys():
                if key.startswith('pdata'):
                    newkey = key.split(os.path.sep)[-1]
                    dic[newkey] = dic.pop(key)

            # Eliminate the digital filter
            if datatype in ['FID', 'SER'] and kargs.get('dgfilter', True):
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

            #store also the varpars of the series
            varpars = kargs.get('varpars')
            if isinstance(varpars, str):
                varpars = varpars.split()
                if len(varpars) == 1:
                    # this should be one of the dic parameters
                    lvarpars.append(adic.par[varpars[0]])  # store the variable
                elif len(varpars) == 2:
                    lvarpars.append(adic.par[varpars[0]][int(varpars[1])])  # store the variable
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
            if self.datatype in ['FIDSER']:
                y = self.par2.lvarpars
            else:
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
            uc = unit_conversion(si, True, sw_h, bf1, -sr + o1)
        else:
            fact = 2.0 if dic.FnMODE in [3, 4, 5, 6] else 1.0
            uc = unit_conversion(td, True, sw_h * fact, bf1, -sr + o1)
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
                if dspfvs not in bruker_dsp_table:
                    raise ValueError("dspfvs not in lookup table")
                if decim not in bruker_dsp_table[dspfvs]:
                    raise ValueError("decim not in lookup table")
                phase = bruker_dsp_table[dspfvs][decim]

        # fft
        si = data.shape[-1]
        pdata = np.fft.fftshift(np.fft.fft(data, si, axis=-1), -1) / float(si / 2)
        pdata = (pdata.T - pdata.T[0]).T  #TODO: this allow generally to remove Bruker smiles, not so sure actually

        # Phasing
        si = float(pdata.shape[-1])
        ph = 2.0j * np.pi * phase * np.arange(si) / si
        pdata = pdata * np.exp(ph)

        # ifft
        ndata = np.fft.ifft(np.fft.ifftshift(pdata, -1), si, axis=-1) * float(si / 2)

        # remove last points * 2
        rp = 4 * (phase // 2)
        pdata = ndata[..., :-rp]

        return pdata

    # -------------------------------------------------------------------------------
    def __repr__(self):
        """
        private function: Create a html representation of some data in the Source traits

        """
        if not self.par:
            return self.__init__.__doc__

        txt = "***Bruker Source***\n"
        txt += "- path: %s\n" % self.path
        txt += "- datatype: %dD NMR %s\n" % (self.par.PARMODE + 1, self.datatype)
        txt += "- pulse program: %s\n" % self.par.PULPROG
        if self.is_2d:
            if self.datatype in ['SER']:
                txt += "- Encoding: %s\n" % self.par2.encoding
                txt += "- TD1: %d, TD2: %d\n" % (self.par2.TD, self.par.TD)
            elif self.datatype in ['FIDSER']:
                if not self.par.isfreq:
                    txt += "- n: %d, TD: %d\n" % (self._data_orig.shape[0], self.par.TD)
                else:
                    txt += "- n: %d, SI: %d\n" % (self._data_orig.shape[0], self.par.SI)
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
            return htmldoc(self.__init__.__doc__)

        txt = self.__repr__().replace('***\n', '</em>')
        txt = txt.replace('***', '<em>')
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

        elif '_data' in name:
            if new is None:
                return
            if name == '_data_orig':
                self.history = ['Loaded %s' % self.datatype]

            # set labels
            if not self.par.isfreq:
                self.par.units = 'ms'
                self.par.label = r'$\mathrm{time\,(%s)}$' % self.par.units
            else:
                self.par.units = 'ppm'
                self.par.label = r'$\mathrm{\delta\,(%s)}$' % self.par.units

            if self.is_2d or self.is_3d:
                if self.par2.units in ['us', 'ms', 's']:
                    self.par2.label = r'$\mathrm{time\,(%s)}$' % self.par.units
                elif self.par2.units in ['hz', 'ppm']:
                    self.par2.label = r'$\mathrm{\delta\,(%s)}$' % self.par.units

        return
