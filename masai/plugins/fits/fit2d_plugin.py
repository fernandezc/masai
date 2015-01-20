# -*- coding: utf-8 -*-
#
# ===============================================================================
# masai.plugins.fits.fit2d_plugin
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
"""
Define fit1d plugin functionalities
"""
from __future__ import print_function, division

# ===============================================================================
#  general imports
# ===============================================================================
# noinspection PyUnresolvedReferences
from traits.api import (HasTraits, Bool, Any, ListStr,
                        Array, Str, Instance, provides, Property)
import traits.has_traits
import sys
import re
import numpy as np
import pandas as pd
from IPython import display

import logging
logger = logging.getLogger('ipython_nmr')

# ===============================================================================
# local imports
# ===============================================================================
from masai.plugins import IPlugin
from masai.core.fit.parameters import ParameterScript
from masai.core.fit.models import getmodel
from masai.core.fit.sequence import PulseSequence
from masai.core.fit.optimize import optimize

traits.has_traits.CHECK_INTERFACES = 2

# ===============================================================================
# Plugin signature
# ===============================================================================
__plugin_main__ = 'Fit2d'

all = ['Fit2d', ]

# ===================================================================================
# Process1d plugin
# ===================================================================================
@provides(IPlugin)
class Fit2d(HasTraits):

    # minimal plugin attribute -----------------------------------------------------
    description = Str("Fit2d (type: fits), fit a 2D source with some models.")
    type = Str("fits")

    # specific attributes ----------------------------------------------------------

    silent = Bool(True)

    _ = Any

    parameterscript = Instance(ParameterScript)
    sequence = Instance(PulseSequence)

    # *******************************************************************************
    # public functions
    # *******************************************************************************

    # -------------------------------------------------------------------------------
    def print_doc(self):
        return self.description

    # -------------------------------------------------------------------------------
    def help(self):
        return """Usage:
            fit2d(source, script)
                where source is a valid Bruker source and script must be provided.

            a default script can be generated using :
                Fit1d.default()
            """

    def default(self):
        s = """
#-----------------------------------------------------------
# syntax for parameters definition:
# name: value, low_bound,  high_bound
#  * for fixed parameters
#  $ for variable parameters
#  > for reference to a parameter in the COMMON block
#    (> is forbidden in the COMMON block)
# common block parameters should not have a _ in their names
#-----------------------------------------------------------
#

COMMON:
# common parameters ex.
# $ gwidth: 1.0, 0.0, none
  $ gratio: 0.5, 0.0, 1.0

MODEL: LINE_1
shape: voigtmodel
	$ ampl:  1.0, 0.0, none
	$ pos:   0.0, -100.0, 100.0
	> ratio: gratio
	$ width: 1.0, 0, 100

"""
        print(s)


    #===========================================================================
    # run
    #===========================================================================
    def run(self, maxiter=2000, maxfun=5000, every=50):
        """
        Main fitting procedure

        Parameters
        ----------
        mode: Str
            define which kind of fit is required among 'FULL2D', 'X_SKY_PROJ'
        kind: Str
            define which kind of experiment

        """
        logger.info('Entering fitting procedure')

        global niter, chi2, everyiter
        everyiter = every
        niter = 0

        # internally defined function chi2
        def chi2(params, dat):
            """
            Return sum((y - x)**2)

            """
            global chi2

            # model spectrum
            modeldata = self._get_modeldata()

            diff = dat - modeldata

            chi2 = diff**2
            chi2 = chi2.values.sum()
            return chi2

        # end chi2 function ---------------------------------------------------


        # callback function--------------------------------------------------------
        def callback(x, f=None, accepted=None):
            """
            callback print function
            """
            global niter, chi2, everyiter
            niter += 1
            if niter%everyiter != 0:
                return
            if f is not None:  # hopping
                return
            if f is None:  # simplex
                display.clear_output(wait=True)
                print("Iteration: %d (chi2: %.5f)"%(niter, chi2))
                sys.stdout.flush()

        # end callback function ---------------------------------------------------

        fp = self.fp
        x = self.source.data.index
        #  fitzone = self.fd.findregion()

        optmode = 'simplex'
        fp, fopt = optimize(chi2, fp,
                            args=(self.source.real,),
                            maxfun=maxfun,
                            maxiter=maxiter,
                            method=optmode,
                            callback=callback)

        self.parameterscript.script = str(fp)
        self.source.modeldata = self._get_modeldata()

        return

    #===========================================================================
    # _callforresults
    #===========================================================================
    def _callforresults(self, params, keys, funcall, iteration, fopt):
        """call back function called by minimize
        """
        for i in range(len(params)):
            key = keys[i]
            # restore external parameters in case of bounding
            pe = self.external(params[i], self.fp.lob[key], self.fp.upb[key])
            self.fp[key] = pe, self.fp.lob[key], self.fp.upb[key]

        # update parameter script
        self.funcalls = funcall
        self.iterations = iteration
        self.fopt = np.sqrt(fopt)

        return self.fittingthread.stopped

    # property

    fp = Property
    def _get_fp(self):
        return self.parameterscript.fp

    script = Property
    def _get_script(self):
        return self.parameterscript.script

    # *******************************************************************************
    # Private functions
    # *******************************************************************************

    # -------------------------------------------------------------------------------
    def __call__(self, *args, **kargs):
        """
        Process a 1D source

        """
        logger.debug(" Fit2d Debug mode")

        if args:
            # look in args
            self.source = source = args[0]
            script = args[1]
        else:
            raise ValueError("No valid source specified\n%s"%self.help())

        # get parameters form script
        self.parameterscript = ParameterScript(fd=source, script=script)

        #sequence = kargs.get('sequence', 'ideal_pulse')
        #self.sequence = PulseSequence(type=sequence)

        self.mode = kargs.get('mode','xproj')
        self.kind = kargs.get('kind',None)

        self.source.modeldata = self._get_modeldata()

        return self

    def _get_modeldata(self):

        # Prepare parameters
        parameters = prepare(self.fp)

        # Get the list of models
        models = self.fp.models
        nbmodels = len(models)

        # Make an array 'modeldata' with the size of the source of data
        # which will contains the data produced by the models
        # This name must always be 'modeldata'
        # which will be returned to the main program.

        #calculate automatic baseline and amplitude
        data = self.source.real

        x = data.columns.values
        y = data.index.values

        modeldata = np.zeros((y.size, x.size)).astype(float)

        if nbmodels < 1:
            return pd.DataFrame(modeldata, index=y, columns=x)

        # Calculates model data
        for model in models:
            calc = getmodel(x, y, modelname=model, par=parameters,
                            fd=self.source, sequence=self.sequence,
                            kind = self.kind)
            modeldata += calc

        # call the userfunction
        # (in case something is executed in it - nothing by default)
        # if default_userscript:
        #     modeldata = self.userfunction(modeldata, self.fp)

        return pd.DataFrame(modeldata, index=y, columns=x)

# global functions
def parsing(expr, param):
    keyword = r"\b([a-z,0-9]+)\b" #match a whole word
    pat = re.compile(keyword)
    mo = pat.findall(expr)
    for kw in mo:
        if param.has_key(kw):
            expr = expr.replace(kw, str(param[kw]))
        elif np.__dict__.has_key(kw):# check if it is a recognized math function
            expr = expr.replace(kw, "np.%s" % kw)
    return expr

def prepare(param):
    """
    This function is needed for the script related to modelfunction
    """
    new_param = param.copy()

    for key in param.iterkeys():
        if param.reference[key]:
            new_param.reference[key] = False    # important to put it here
                                                # before other instruction
            # try to interpret the given refpar expression
            refpar = param[key]
            while True:
                par = parsing(refpar, new_param)
                if par == refpar: break
                refpar = par

            try:
                new_param[key] = eval(str(refpar))
            except:
                raise ValueError('Cannot evaluate the expression: %s: %s'
                                 % (key, param[refpar]))
                return
            new_param.fixed[key] = True
            new_param.reference[key] = True     #restore it for the next call

    #for k, v in new_param.iteritems():
    #    print(k, v)
    return new_param

#===============================================================================
# automatic calculation of amplitude and baseline
#===============================================================================
def ampandbas(xi, expe, calc):
    """Automatically calculate correct amplitude A and baseline
    (baseline linear model a*i+b) by detemining the zero of the first derivative
     with respect to A, a, and b
    """
    n = xi.size
    sic = sum(xi*calc)
    sc = sum(calc)
    sie = sum(xi*expe)
    se = sum(expe)
    sce = sum(expe*calc)
    sc2 = sum(calc*calc)
    A = (
         (12*n*sic + (-6*n**2 + 18*n - 12)*sc)*sie
         + (-6*n**2 + 18*n - 12)*se*sic
         + (4*n**3 - 6*n**2 + 2*n)*sc*se
         + (-n**4 - 12*n**3 + 37*n**2 - 36*n + 12)*sce)/(12*n*sic**2
                                                + (-12*n**2 + 36*n - 24)*sc*sic
         + (-n**4 - 12*n**3 + 37*n**2 - 36*n + 12)*sc2
         + (4*n**3 - 6*n**2 + 2*n)*sc**2)
    a = -(
          (12*n*sc2 - 12*sc**2)*sie + (12*sc*se - 12*n*sce)*sic
          + (-6*n**2 + 18*n - 12)*sc2*se
          + (6*n**2 - 18*n + 12)*sc*sce)/(12*n*sic**2
                                                + (-12*n**2 + 36*n - 24)*sc*sic
          + (-n**4 - 12*n**3 + 37*n**2 - 36*n + 12)*sc2
          + (4*n**3 - 6*n**2 + 2*n)*sc**2)
    b = -(
          (12*sc*sic + (-6*n**2 + 18*n - 12)*sc2)*sie
          - 12*se*sic**2 + (6*n**2 - 18*n + 12)*sce*sic
          + (4*n**3 - 6*n**2 + 2*n)*sc2*se
          + (-4*n**3 + 6*n**2 - 2*n)*sc*sce)/(12*n*sic**2
          + (-12*n**2 + 36*n - 24)*sc*sic
          + (-n**4 - 12*n**3 + 37*n**2 - 36*n + 12)*sc2
          + (4*n**3 - 6*n**2 + 2*n)*sc**2)
    # in case the modeldata is zero, to avoid further errors
    if np.isnan(A): A=0.0
    if np.isnan(a): a=0.0
    if np.isnan(b): b=0.0
    return A, a, b

if __name__ == '__main__':
    import os
    from masai.plugins.list import Bruker
    from masai.plugins.list import Process
    from masai.plugins.list import Plot1d
    from masai.plugins.list import Plot2d

    Bruker = Bruker()
    Plot1d = Plot1d()
    Plot2d = Plot2d()
    Process = Process()
    Fit2d = Fit2d()

    from masai.core.utils.pyprind import ProgBar

    def spec2Q():

        # reading
        DATADIR = os.path.join( '..', '..', '..', 'exemples')
        user = 'user'
        name = 'VPI7'
        expno = 102
        path = os.path.join(DATADIR, user, 'nmr', name, str(expno))
        source = Bruker(path)

        # scaling preparation
        # This is a DQ experiment so:
        source.par2.SFO1 = source.par.SFO1*2.0
        source.par2.SF = source.par.SF*2.0
        # correction due to the missetted ofset
        fn = 2.0 # 2 rotation of the spectra
        source.par2.SF = source.par2.SF + fn * source.par2.SW_h*1.0e-6

        # process  1srt dimention
        Process(source, "em 1")
        Process(source, "zf 2")
        Process(source, "ft")
        Process(source, transform = 'cl -99 -80')
        Process(source, 'pk 0 0 --auto --fit_phc1')
        Process(source, transform = 'ab --mode poly')

        # process 2nd dimension
        Process(source, axis=0, transform='em 1')
        Process(source, axis=0, transform='zf 2')
        Process(source, axis=0, transform='ft')
        Process(source, axis=0, transform = 'cl -198 -160')
        Process(source, axis=0, transform='pk 180 0 --auto --fit_phc1 --select cols')
        Process(source, axis=0, transform='ab --mode linear')

        # magnitude calculation
        Process(source, transform='mc')

        # plot2D
        p_2 = Plot2d(source, figsize=(3,2.7), spacing=5, exponent=1.2, xproj='sky', hold=True)

        # plot the diagonal
        p_2.axe.plot(source.data.columns, source.data.columns*2.0, lw=.1, ls='--')

        #row = source.row(-176.6, width=.2)
        #p_1 = Plot1d(row)

        common = "COMMON:"
        common +="\n * gwidth: 0.25, 0.21, 0.3"
        common +="\n * gwidth1: 1.5, 1.0, 2.0"

        spc={
            'A': -82.,
            'B': -92.6,
            'C': -92.7,
            'D': -93.9,
            'E': -94.2,
            'F': -94.9,
            'G': -95.0,
            }

        # make a dictionary sorted by key
        from collections import OrderedDict
        species = OrderedDict(sorted(spc.items(), key=lambda t: t[0]))
        print(species)
        connectivity = np.ones((7,7))  #*source.data.values.max()

        # some correlation are known to be missing
        connectivity[0,0] = 0.0                     # A-A missing
        connectivity[0,2] = connectivity[2,0] = 0.0  # A-C
        connectivity[0,3] = connectivity[3,0] = 0.0  # A-D
        connectivity[0,4] = connectivity[4,0] = 0.0  # A-E
        connectivity[0,5] = connectivity[5,0] = 0.0  # A-F

        connectivity = connectivity*source.data.values.max()

        def make_script(common, species, connectivity):
            script = common
            for i, key in enumerate(species.keys()):
                # additional common parameters
                v = species[key]
                dv = abs(v)*.01
                script += "\n $ gpos%s: %.3f, %.3f, %.3f"%(key, v, v-dv, v+dv)
            for i, key in enumerate(species.keys()):
                for j, key1 in enumerate(species.keys()):
                    # as it is normally symmetric we can try to reduce the number of varaible
                    if j<i:
                        continue
                    script += "\n\n MODEL: LINE_"+key+key1
                    script += "\n shape: gaussian2dmodel"
                    #if i==0:
                    #    script += "\n * ampl:  1., 0.0, none"
                    #else:
                    fi = '$' if connectivity[i,j]>0.001 else '*'
                    script += "\n %s ampl:  %.3f, 0.0, none"%(fi, connectivity[i,j])
                    script += "\n > pos: gpos%s "%key
                    script += "\n > pos1: (gpos%s) + (gpos%s)"%(key, key1)
                    script += "\n > width: gwidth"
                    script += "\n > width1: gwidth1"
                    script += "\n * theta: 0.0"
            return script.lower()


        script = make_script(common, species, connectivity)
        print(script)
        f = Fit2d(source,  script, kind='SQ-DQ')

        f.run(maxiter=20000, maxfun=100000, every=100)
        x = source.data.columns
        y = source.data.index
        p_2.axe.contour(x, y, source.modeldata, 10)

        p_2.show()

        print(f.script)

    spec2Q()


