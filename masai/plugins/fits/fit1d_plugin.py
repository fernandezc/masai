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
from traits.api import (HasTraits, Bool, Any, ListStr,
                        Array, Str, Instance, provides, Property)
import traits.has_traits
import re
import numpy as np
import pandas as pd

import logging
LOGLEVEL = logging.DEBUG
logger = logging.getLogger()

# ===============================================================================
# local imports
# ===============================================================================
from masai.plugins import IPlugin
from masai.core.fit.parameters import ParameterScript
from masai.core.fit.models import getmodel
from masai.core.fit.sequence import PulseSequence

traits.has_traits.CHECK_INTERFACES = 2

# ===============================================================================
# Plugin signature
# ===============================================================================
__plugin_main__ = 'Fit1d'

# ===================================================================================
# Process1d plugin
# ===================================================================================
@provides(IPlugin)
class Fit1d(HasTraits):

    # minimal plugin attribute -----------------------------------------------------
    description = Str("Fit1d (type: fits), fit a 1D source with some models.")
    type = Str("fits")

    # specific attributes ----------------------------------------------------------

    silent = Bool(True)

    _ = Any

    script = """

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
        print(self.script)

    # property

    fp = Property
    def _get_fp(self):
        return self.parameterscript.fp

    # *******************************************************************************
    # Private functions
    # *******************************************************************************

    # -------------------------------------------------------------------------------
    def __call__(self, *args, **kargs):
        """Process a 1D source
        """
        logger.setLevel(LOGLEVEL)
        logging.debug(" Fit1d Debug mode")

        if args:
            # look in args
            self.source = source = args[0]
            script = args[1]
        else:
            raise ValueError("No valid source specified\n%s"%self.help())

        # get parameters form script
        self.parameterscript = ParameterScript(fd=source, script=script)

        sequence = kargs.get('sequence', 'ideal_pulse')
        self.sequence = PulseSequence(type=sequence)

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
            if not model.startswith('baseline'):
                row += 1
                modeldata[row] = getmodel(x, modelname=model, par=parameters,
                                          fd=self.source, sequence=self.sequence)
                names.append(model)
            else:
                modeldata[0] += getmodel(x, modelname=model, par=parameters,
                                         fd=self.source, sequence=self.sequence)

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
        at = 0.
        for model in models:
            if model.startswith('baseline'):
                continue
            at += parameters["ampl_%s"%model]

        self.intensities = 'Intensities:\n       \tppm   \tRelative (%)     \tAbsolute\n'
        for model in models:
            if model.startswith('baseline'):
                continue
            self.intensities += \
                '%(model)7s:\t%(pos).2f\t%(ampl)g \t%(conc)g\n'%{
                                  'model':model,
                                  'pos':parameters["pos_%s"%model],
                                  'ampl':parameters["ampl_%s"%model]*100./at,
                         'conc':parameters["ampl_%s"%model]*self.A}

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
    from masai.plugins.list import Process1d
    from masai.plugins.list import Plot1d

    Bruker = Bruker()
    Plot1d = Plot1d()
    Process1d = Process1d()
    Fit1d = Fit1d()

    # global
    DATADIR = os.path.expanduser(os.path.join('~', 'Dropbox', 'storage', 'bruker', 'data'))
    user = 'christian'
    name = 'example1'
    expno = 10
    path = os.path.join(DATADIR, user, 'nmr', name, str(expno))
    source = Bruker(path)
    Process1d(source, "em 50")
    Process1d(source, "zf 0")
    Process1d(source, "ft")
    Process1d(source, "pk -25 0 --auto --use_phc1 --ediff 2 --mode negmin+entropy --gamma 0.0001 --bound_phc1 200 --verbose")
    #Process1d(source, "ab --mode poly")
    #Plot1d(source)

    script = """
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
    * ampl:  1.0, 0.0, none
    $ pos:   4.5, -100.0, 100.0
    > ratio: gratio
    $ width: 4.0, 0, 100

MODEL: LINE_2
shape: voigtmodel
    $ ampl:  0.2, 0.0, none
    $ pos:   1.6, -100.0, 100.0
    > ratio: gratio
    $ width: 4.0, 0, 100

MODEL: LINE_3
shape: voigtmodel
    $ ampl:  0.2, 0.0, none
    $ pos:   7.9, -100.0, 100.0
    > ratio: gratio
    $ width: 4.0, 0, 100

"""

    f = Fit1d(source, script)
    ax = Plot1d(source, xlim=(-10, 30), plot_model=True)
    print("end")