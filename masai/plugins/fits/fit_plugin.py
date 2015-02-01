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
__plugin_main__ = 'Fit'

all = ['Fit', ]

# ===================================================================================
# Process1d plugin
# ===================================================================================
@provides(IPlugin)
class Fit(HasTraits):

    # minimal plugin attribute -----------------------------------------------------
    description = Str("Fit (type: fits), fit a source or a series of source with some models.")
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
                Fit.default()
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
    def run(self, maxiter=5000, maxfun=None, every=100,
            method='simplex', constraints=None):
        """
        Main fitting procedure
        """
        print('Entering fitting procedure')

        global niter, chi2, everyiter, ncalls
        ncalls = 0
        everyiter = every
        niter = 0

        # internally defined function chi2
        def funchi2(params, sources, *constraints):
            """
            Return sum((y - x)**2)
            """
            global chi2, ncalls
            # model spectrum

            chi2 = 0
            som = 0
            ncalls +=1

            for exp_idx, source in enumerate(sources):
                modeldata = self._get_modeldata(source, exp_idx)
                # baseline is already summed with modeldata[-1]
                # important to work with the real part of source
                # not the complex number
                data = source.real.values

                if not source.is_2d:
                    mdata = modeldata.loc['modelsum'].values
                else:
                    mdata = modeldata.values

                merror = 1.
                if source.is_2d:
                    if constraints:

                        # Case of SQ-DQ experiments
                        if self.kind == 'SQ-DQ' and \
                                    'max_connections' in constraints[0]:
                            # check connectivity numbers
                            nbconnections = {}
                            for key in params.keys():
                                if 'pos1' in key:
                                    connect = key[-2:]
                                    key = 'ampl_line_'+connect # get amplitude
                                    ki = connect[0].upper()
                                    if ki not in nbconnections.keys():
                                        nbconnections[ki]=0
                                    if int(params[key])>0:
                                        nbconnections[ki] += 1
                            for k,v in nbconnections.iteritems():
                                if v>constraints[0]['max_connections']:
                                    merror *= v*10.

                diff = data - mdata
                chi2 += np.sum(diff**2)*merror
                som += np.sum(data**2)

            chi2 = np.sqrt(chi2/som)

            return chi2

        # end chi2 function ---------------------------------------------------


        # callback function--------------------------------------------------------
        def callback(*args, **kargs):
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
                           constraints=constraints,
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

        x = data.columns.values

        if not source.is_2d: #1D
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
                calc = getmodel(x, modelname=model, par=parameters, fd=source)
                if not model.startswith('baseline'):
                    row += 1
                    modeldata[row] = calc
                    names.append(model)
                else:
                    modeldata[0] += calc

            # make the sum
            modeldata[row + 1] = modeldata.sum(axis=0)
            names.append('modelsum')

            #remove unused row
            modeldata = modeldata[:row + 2]

            xi = np.arange(float(x.size))
            A, a, b = ampbas(xi, data.values[-1], modeldata[-1])  #(fitzone-fitzone[0], data.take(fitzone),
                                 #modeldata[-1].take(fitzone))

            modeldata = A*modeldata
            baseline = a*xi + b   #a*(xi-fitzone[0]) + b

            #update the modeldata
            modeldata[0] += baseline
            modeldata[-1] += baseline

            self.A = A

            # return dataframes
            return pd.DataFrame(modeldata, index=names, columns=x)

        else: #2D
            y = data.index.values
            modeldata = np.zeros((y.size, x.size)).astype(float)

            if nbmodels < 1:
                return pd.DataFrame(modeldata, index=y, columns=x)

            # Calculates model data
            for model in models:
                calc = getmodel(x, y, modelname=model, par=parameters,
                               fd=source, kind = self.kind)
                modeldata +=  calc

            xi = np.arange(float(x.size))
            yj = np.arange(float(y.size))
            A, a, b, c = ampbas2D(xi, yj, data.values, modeldata)  #(fitzone-fitzone[0], data.take(fitzone),
                                 #modeldata[-1].take(fitzone))

            modeldata = A*modeldata
            base = np.ones_like(modeldata)
            baseline = a*xi*base + (b*yj*base.T).T +  c

            #update the modeldata
            modeldata += baseline

            self.A = A

        # call the userfunction
        # (in case something is executed in it - nothing by default)
        # if default_userscript:
        #     modeldata = self.userfunction(modeldata, self.fp)

            return pd.DataFrame(modeldata, index=y, columns=x)



# global functions
def parsing(expr, param):
    keyword = r"\b([a-z]+[0-9]*)\b" #match a whole word
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
                if par == refpar:
                    break
                refpar = par

            try:
                new_param[key] = eval(str(refpar))
            except:
                raise ValueError('Cannot evaluate the expression: %s: %s'
                                 % (key, param[refpar]))
                return

            new_param.fixed[key] = True
            new_param.reference[key] = True     #restore it for the next call

        if isinstance(new_param[key], list):
            new_param.data[key] = new_param.data[key][exp_idx]

    return new_param

#===============================================================================
# automatic calculation of amplitude and baseline
#===============================================================================
def ampbas(xi, expe, calc):
    """Automatically calculate correct amplitude A and baseline
    (baseline linear model a*i+b) by detemining the zero of the first derivative
     with respect to A, a, and b
    """

    n = xi.size
    sE = sum(expe)
    sF = sum(calc)
    sFI = sum(xi*calc)
    sFd = sum(calc*calc)
    sI = sum(xi)
    sEF = sum(expe*calc)
    sEI = sum(xi*expe)
    sId = sum(xi**2)

    a = (-sE*(sF*sFI - sFd*sI) + sEF*(n*sFI - sF*sI) - sEI*(n*sFd - sF**2)
        )/(n*sFI**2 - n*sFd*sId + sF**2*sId - 2*sF*sFI*sI + sFd*sI**2)

    A= (sE*(sF*sId - sFI*sI) - sEF*(n*sId - sI**2) + sEI*(n*sFI - sF*sI)
       )/(n*sFI**2 - n*sFd*sId + sF**2*sId - 2*sF*sFI*sI + sFd*sI**2)

    b= (sE*(sFI**2 - sFd*sId) + sEF*(sF*sId - sFI*sI) - sEI*(sF*sFI - sFd*sI)
       )/(n*sFI**2 - n*sFd*sId + sF**2*sId - 2*sF*sFI*sI + sFd*sI**2)

    # in case the modeldata is zero, to avoid further errors
    if np.isnan(A): A=0.0
    if np.isnan(a): a=0.0
    if np.isnan(b): b=0.0
    return A, a, b

def ampbas2D(xi, yj, expe, calc):

    n = float(xi.size)
    m = float(yj.size)
    sE = expe.sum()
    sF = calc.sum()
    sFI = (xi*calc).sum()
    sFJ = ((yj*calc.T)).sum()
    sFd = (calc*calc).sum()
    sI = sum(xi)
    sJ = sum(yj)
    sIJ = ((np.ones_like(calc)*xi).T*yj).sum()
    sEF = (expe*calc).sum()
    sEI = (xi*expe).sum()
    sEJ = (yj*expe.T).sum()
    sId = sum(xi**2)
    sJd = sum(yj**2)

    c= (sE*(-m*n*sFd*sId*sJd + m*sFJ**2*sId + n*sFI**2*sJd - 2*sFI*sFJ*sIJ + sFd*sIJ**2) + sEF*(m*n*sF*sId*sJd
        - m*n*sFI*sI*sJd - m*n*sFJ*sId*sJ + m*sFJ*sI*sIJ + n*sFI*sIJ*sJ - sF*sIJ**2) + sEI*(m*n*sFd*sI*sJd
        - m*sFJ**2*sI - n*sF*sFI*sJd + n*sFI*sFJ*sJ - n*sFd*sIJ*sJ + sF*sFJ*sIJ) + sEJ*(m*n*sFd*sId*sJ - m*sF*sFJ*sId
        + m*sFI*sFJ*sI - m*sFd*sI*sIJ - n*sFI**2*sJ + sF*sFI*sIJ))/(-m**2*n**2*sFd*sId*sJd + m**2*n*sFJ**2*sId
        + m**2*n*sFd*sI**2*sJd - m**2*sFJ**2*sI**2 + m*n**2*sFI**2*sJd + m*n**2*sFd*sId*sJ**2 + m*n*sF**2*sId*sJd
        - 2*m*n*sF*sFI*sI*sJd - 2*m*n*sF*sFJ*sId*sJ + 2*m*n*sFI*sFJ*sI*sJ - 2*m*n*sFI*sFJ*sIJ - 2*m*n*sFd*sI*sIJ*sJ
        + m*n*sFd*sIJ**2 + 2*m*sF*sFJ*sI*sIJ - n**2*sFI**2*sJ**2 + 2*n*sF*sFI*sIJ*sJ - sF**2*sIJ**2)

    a= (n*sEF*(m*n*sFI*sJd - m*sF*sI*sJd + m*sFJ*sI*sJ - m*sFJ*sIJ - n*sFI*sJ**2 + sF*sIJ*sJ) + n*sEI*(-m*n*sFd*sJd
        + m*sFJ**2 + n*sFd*sJ**2 + sF**2*sJd - 2*sF*sFJ*sJ) + sE*(m*n*sFd*sI*sJd - m*sFJ**2*sI - n*sF*sFI*sJd
        + n*sFI*sFJ*sJ - n*sFd*sIJ*sJ + sF*sFJ*sIJ) - sEJ*(m*n*sFI*sFJ + m*n*sFd*sI*sJ - m*n*sFd*sIJ - m*sF*sFJ*sI
        - n*sF*sFI*sJ + sF**2*sIJ))/(-m**2*n**2*sFd*sId*sJd + m**2*n*sFJ**2*sId + m**2*n*sFd*sI**2*sJd
        - m**2*sFJ**2*sI**2 + m*n**2*sFI**2*sJd + m*n**2*sFd*sId*sJ**2 + m*n*sF**2*sId*sJd - 2*m*n*sF*sFI*sI*sJd
        - 2*m*n*sF*sFJ*sId*sJ + 2*m*n*sFI*sFJ*sI*sJ - 2*m*n*sFI*sFJ*sIJ - 2*m*n*sFd*sI*sIJ*sJ + m*n*sFd*sIJ**2
        + 2*m*sF*sFJ*sI*sIJ - n**2*sFI**2*sJ**2 + 2*n*sF*sFI*sIJ*sJ - sF**2*sIJ**2)

    A= (m*n*sEF*(-m*n*sId*sJd + m*sI**2*sJd + n*sId*sJ**2 - 2*sI*sIJ*sJ + sIJ**2) + m*sEJ*(m*n*sFJ*sId - m*sFJ*sI**2
        - n*sF*sId*sJ + n*sFI*sI*sJ - n*sFI*sIJ + sF*sI*sIJ) + n*sEI*(m*n*sFI*sJd - m*sF*sI*sJd + m*sFJ*sI*sJ
        - m*sFJ*sIJ - n*sFI*sJ**2 + sF*sIJ*sJ) + sE*(m*n*sF*sId*sJd - m*n*sFI*sI*sJd - m*n*sFJ*sId*sJ + m*sFJ*sI*sIJ
        + n*sFI*sIJ*sJ - sF*sIJ**2))/(-m**2*n**2*sFd*sId*sJd + m**2*n*sFJ**2*sId + m**2*n*sFd*sI**2*sJd
        - m**2*sFJ**2*sI**2 + m*n**2*sFI**2*sJd + m*n**2*sFd*sId*sJ**2 + m*n*sF**2*sId*sJd - 2*m*n*sF*sFI*sI*sJd
        - 2*m*n*sF*sFJ*sId*sJ + 2*m*n*sFI*sFJ*sI*sJ - 2*m*n*sFI*sFJ*sIJ - 2*m*n*sFd*sI*sIJ*sJ + m*n*sFd*sIJ**2
        + 2*m*sF*sFJ*sI*sIJ - n**2*sFI**2*sJ**2 + 2*n*sF*sFI*sIJ*sJ - sF**2*sIJ**2)

    b= (m*sEF*(m*n*sFJ*sId - m*sFJ*sI**2 - n*sF*sId*sJ + n*sFI*sI*sJ - n*sFI*sIJ + sF*sI*sIJ) + m*sEJ*(-m*n*sFd*sId
        + m*sFd*sI**2 + n*sFI**2 + sF**2*sId - 2*sF*sFI*sI) + sE*(m*n*sFd*sId*sJ - m*sF*sFJ*sId + m*sFI*sFJ*sI
        - m*sFd*sI*sIJ - n*sFI**2*sJ + sF*sFI*sIJ) - sEI*(m*n*sFI*sFJ + m*n*sFd*sI*sJ - m*n*sFd*sIJ - m*sF*sFJ*sI
        - n*sF*sFI*sJ + sF**2*sIJ))/(-m**2*n**2*sFd*sId*sJd + m**2*n*sFJ**2*sId + m**2*n*sFd*sI**2*sJd
        - m**2*sFJ**2*sI**2 + m*n**2*sFI**2*sJd + m*n**2*sFd*sId*sJ**2 + m*n*sF**2*sId*sJd - 2*m*n*sF*sFI*sI*sJd
        - 2*m*n*sF*sFJ*sId*sJ + 2*m*n*sFI*sFJ*sI*sJ - 2*m*n*sFI*sFJ*sIJ - 2*m*n*sFd*sI*sIJ*sJ + m*n*sFd*sIJ**2
        + 2*m*sF*sFJ*sI*sIJ - n**2*sFI**2*sJ**2 + 2*n*sF*sFI*sIJ*sJ - sF**2*sIJ**2)

    # in case the modeldata is zero, to avoid further errors
    if np.isnan(A): A=0.0
    if np.isnan(a): a=0.0
    if np.isnan(b): b=0.0
    if np.isnan(c): c=0.0
    return A, a, b, c

if __name__ == '__main__':
    import os, commands
    from masai.api import Bruker, Process, Plot1d, Plot2d,  Fit

    # we try here a multi source fit
    # print(commands.getoutput('pwd'))

    def multifit():

        DATADIR = os.path.expanduser(os.path.join('..', '..', '..', 'exemples'))
        user = 'user'
        name = 'RUB-17'

        path1 = os.path.join(DATADIR, user, 'nmr', name, str(1))
        source1 = Bruker(path1)
        Process(source1, transform='em 1; zf 0; ft; cl -105 -75; pk --fit_phc1 --auto; ab --mode poly')

        path2 = os.path.join(DATADIR, user, 'nmr', name, str(2))
        source2 = Bruker(path2)
        Process(source2, transform='em 1; zf 0; ft; cl -105 -75; pk --fit_phc1 --auto; ab --mode poly') #')#

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
          $ gwidth: 0.6, 0.5, 0.7

        """
        for i, delta in enumerate(pos):
            fa = '*' if i==0 else '$'
            script +="""
                MODEL: LINE_%d
                shape: voigtmodel
                    %s ampl:  0.6, 0.5, 2
                    $ pos: %.3f, %.3f, %.3f
                    > ratio: gratio
                    > width: gwidth

            """%(i, fa, delta, delta-.5, delta+0.5)

        print(script)

        # prepare fit
        f = Fit( [source1, source2], script)

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

        f.run(maxiter=200, every=5, method='simplex')

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

            # fitted model plot
        p = Plot1d(source2,
                    figsize=(3.4, 2.0), dpi=180,
                    #ow = 2,
                    #col = 1,
                    lw=.5,
                    modls = '--',
                    plot_model=True,
                    offset=15,
                    xlim=(-98,-80))

    #### TEST 2D ######
    def spec2Q():

        # reading
        DATADIR = os.path.join( '..', '..', '..', 'exemples')
        user = 'user'
        name = 'RUB-17'
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
            fixdone = False
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
                    if fi=='$' and not fixdone:
                        fi = '*'
                        fixdone=True
                    script += "\n %s ampl:  %.3f, 0.0, none"%(fi, connectivity[i,j])
                    script += "\n > pos: gpos%s "%key
                    script += "\n > pos1: (gpos%s) + (gpos%s)"%(key, key1)
                    script += "\n > width: gwidth"
                    script += "\n > width1: gwidth1"
                    script += "\n * theta: 0.0"
            return script.lower()


        script = make_script(common, species, connectivity)
        print(script)
        f = Fit(source,  script, kind='SQ-DQ')

        f.run(maxiter=20000, maxfun=100000, every=100,
              constraints={'max_connections':3})
        x = source.data.columns
        y = source.data.index
        p_2.axe.contour(x, y, source.modeldata, 10)

        p_2.show()

        print(f.script)

    ###############################

    #multifit()
    spec2Q()