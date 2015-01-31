# -*- coding: utf-8 -*-
#
# ===============================================================================
# masai.plugins.fits.fit1d_plugin
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
from traits.api import (HasTraits, Bool, Any, ListStr, List,
                        Array, Str, Instance, provides, Property)
import traits.has_traits
import sys
import re
import numpy as np
import pandas as pd
from IPython import display

# ===============================================================================
# local imports
# ===============================================================================
from masai.plugins import IPlugin
from masai.core.fit.parameters import ParameterScript
from masai.core.fit.models import getmodel
from masai.core.fit.sequence import PulseSequence
from masai.core.fit.optimize import optimize
from masai.core.source import Source

traits.has_traits.CHECK_INTERFACES = 2

# ===============================================================================
# Plugin signature
# ===============================================================================
__plugin_main__ = 'Fit1d'

all = ['Fit1d', ]

# ===================================================================================
# Process1d plugin
# ===================================================================================
@provides(IPlugin)
class Fit1d(HasTraits):

    # minimal plugin attribute -----------------------------------------------------
    description = Str("Fit1d (type: fits), fit a 1D source with some models.")
    type = Str("fits")

    # specific attributes ----------------------------------------------------------
    sources = List(Instance(Source))
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
            fit1d(source, script)
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
experiment_variables: ampl

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
    def run(self, maxiter=5000, maxfun=None, every=100, method='simplex'):
        """Main fitting procedure
        """
        print('Entering fitting procedure')

        global niter, chi2, everyiter, ncalls
        ncalls = 0
        everyiter = every
        niter = 0

        # internally defined function chi2
        def funchi2(params, sources):
            """
            Return sum((y - x)**2)
            """
            global chi2, ncalls
            # model spectrum

            chi2 = 0
            som = 0
            ncalls +=1
            for exp_idx, source in enumerate(self.sources):
                modeldata = self._get_modeldata(source, exp_idx)
                # baseline is already summed with modeldata[-1]
                # important to work with the real part of source
                # not the complex number
                diff = source.real - modeldata.T['modelsum']
                chi2 += np.sum(diff.iloc[0]**2)
                som += np.sum(source.real.iloc[0]**2)

            chi2 = np.sqrt(chi2/som)
            return chi2

        # end chi2 function ---------------------------------------------------


        # callback function--------------------------------------------------------
        def callback(x, f=None, accepted=None):
            """
            callback print function
            """
            global niter, chi2, everyiter, ncalls
            niter += 1

            if niter%everyiter != 0:
                return
            #if f is not None:  # hopping
            #    pass
            #if f is None:  # simplex
            display.clear_output(wait=True)
            print("Iterations: %d, Calls: %d (chi2: %.5f)"%(niter, ncalls, chi2))
            sys.stdout.flush()

        # end callback function ---------------------------------------------------

        fp = self.fp # starting parameters

        fp, fopt = optimize(funchi2, fp,
                           args=(self.sources,),
                           maxfun=maxfun,
                           maxiter=maxiter,
                           method=method,
                           callback=callback)

        # replace the previous script
        self.parameterscript.script = str(fp)
        print(self.parameterscript.script)

        for exp_idx, source in enumerate(self.sources):
            source.modeldata = self._get_modeldata(source, exp_idx)

        return

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
    def __init__(self, *args, **kargs):
        """
        fit a 1D or 2D sources, or a list of sources

        """

        if args:
            # look in args
            if not isinstance(args[0], list):
                self.sources = [args[0], ]
            else:
                self.sources = args[0]
                # we create a list of source in all case
            script = args[1]

        else:
            raise ValueError("No valid source specified\n%s"%self.help())

        # get parameters form script
        self.parameterscript = ParameterScript(sources=self.sources, script=script)

        #sequence = kargs.get('sequence', 'ideal_pulse')
        #self.sequence = PulseSequence(type=sequence)

        self.mode = kargs.get('mode',None)
        self.kind = kargs.get('kind',None)

        for exp_idx, source in enumerate(self.sources):
            source.modeldata = self._get_modeldata(source, exp_idx)

    def _get_modeldata(self, source, exp_idx):

        # Prepare parameters
        parameters = prepare(self.fp, exp_idx)

        # Get the list of models
        models = self.fp.models
        nbmodels = len(models)

        # Make an array 'modeldata' with the size of the source of data
        # which will contains the data produced by the models
        # This name must always be 'modeldata'
        # which will be returned to the main program.

        data = source.real

        x =data.columns.values
        xi = np.arange(x.size)

        modeldata = np.zeros((nbmodels + 2, x.size)).astype(float)

        if nbmodels < 1:
            names = ['baseline', 'modelsum']
            return pd.DataFrame(modeldata, index=names, columns=x)

        # Calculates model data
        # The first row (i=0) of the modeldata array is the baseline,
        # so we fill the array starting at row 1
        row = 0
        names = ['baseline', ]

        for model in models:
            calc = getmodel(x, modelname=model, par=parameters,
                            fd=source, exp_idx=exp_idx)
            if not model.startswith('baseline'):
                row += 1
                modeldata[row] = calc
                names.append(model)
            else:
                modeldata[0] += calc


        # call the userfunction
        # (in case something is executed in it - nothing by default)
        # if default_userscript:
        #     modeldata = self.userfunction(modeldata, self.fp)

        # make the sum
        modeldata[row + 1] = modeldata.sum(axis=0)
        names.append('modelsum')

        #remove unused row
        modeldata = modeldata[:row + 2]

        # take into account fitzone
        #fitzone = findregion(x, self.source.zone)

        A, a, b = ampandbas(xi, data.values[-1], modeldata[-1])  #(fitzone-fitzone[0], data.take(fitzone),
                             #modeldata[-1].take(fitzone))

        modeldata = A*modeldata
        baseline = a*xi + b   #a*(xi-fitzone[0]) + b

        #update the modeldata
        modeldata[0] += baseline
        modeldata[-1] += baseline

        #self._save_modeldata = modeldata
        #self._save_names = names
        self.A = A

        # calculate intensities
        # at = 0.
        # for model in models:
        #     if model.startswith('baseline'):
        #         continue
        #     if 'ampl' in parameters.expvars:
        #         expe = "_exp%d"%exp_idx
        #     at += parameters["ampl_%s%s"%(model,expe)]
        #
        # self.intensities = 'Intensities:\n       \tppm   \tRelative (%)     \tAbsolute\n'
        # for model in models:
        #     if model.startswith('baseline'):
        #         continue
        #     if 'ampl' in parameters.expvars:
        #         expe = "_exp%d"%exp_idx
        #     self.intensities += \
        #         '%(model)7s:\t%(pos).2f\t%(ampl)g \t%(conc)g\n'%{
        #                           'model':model,
        #                           'pos':parameters["pos_%s"%model],
        #                           'ampl':parameters["ampl_%s"%(model,expe)]*100./at,
        #                  'conc':parameters["ampl_%s"%model]*self.A}

        # return dataframes
        return pd.DataFrame(modeldata, index=names, columns=x)


# global functions
def parsing(expr, param):
    keyword = r"\b([a-z]+)\b" #match a whole word
    pat = re.compile(keyword)
    mo = pat.findall(expr)
    for kw in mo:
        if param.has_key(kw):
            expr = expr.replace(kw, str(param[kw]))
        elif np.__dict__.has_key(kw):# check if it is a recognized math function
            expr = expr.replace(kw, "np.%s" % kw)
    return expr

def prepare(param, exp_idx):
    """
    This function is needed for the script related to modelfunction

    exp_idx: int
        contains the index of the experiment
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
    import os, commands
    from masai.api import Bruker, Process, Plot1d, Fit1d

    # we try here a multi source fit
    # print(commands.getoutput('pwd'))

    DATADIR = os.path.expanduser(os.path.join('..', '..', '..', 'exemples'))
    user = 'user'
    name = 'RUB-17'

    path1 = os.path.join(DATADIR, user, 'nmr', name, str(1))
    source1 = Bruker(path1)
    Process(source1, transform='em 1; zf 0; ft; cl -105 -75; pk --fit_phc1 --auto; ab --mode poly')

    path2 = os.path.join(DATADIR, user, 'nmr', name, str(2))
    source2 = Bruker(path2)
    Process(source2, transform='em 1; zf 0; ft; cl -105 -75; pk --fit_phc1 --auto; ab --mode poly')

    #p = Plot1d(source1, hold=True)
    #p = Plot1d(source2, axe=p.axe)

    # fit
    # try with variable width
    pos = [-94.9200018156, -95.0999986491, -94.3173715978, -93.9072128893,
           -92.8999998621, -92.6363220634, -82.0826836442]

    script = """
    COMMON:
      experiment_variables: ampl
      $ gratio: 0.5, 0.0, 1.0

    """
    for i, delta in enumerate(pos):
        fa = '*' if i==0 else '$'
        script +="""
            MODEL: LINE_%d
            shape: voigtmodel
                %s ampl:  0.6, 0.0, none
                $ pos: %.3f, %.3f, %.3f
                > ratio: gratio
                $ width: %.3f, %.3f, %.3f

        """%(i, fa, delta, delta-.2, delta+0.2,
                    0.57, 0.45, 0.65)

    print(script)

    # prepare fit
    f = Fit1d( source1, script)

    # initial model plot
    #p = Plot1d([source1, source2],
    #            figsize=(3.4, 2.0), dpi=180,
    #            row = 2,
    #            col = 1,
    #            lw=.5,
    #            modls = '--',
    #            plot_model=True,
    #            offset=15,
    #            xlim=(-98,-80))

    f.run(maxiter=1000, method='simplex')

    # fitted model plot
    p = Plot1d(source1,
                figsize=(3.4, 2.0), dpi=180,
                #ow = 2,
                #col = 1,
                lw=.5,
                modls = '--',
                plot_model=True,
                offset=15,
                xlim=(-98,-80))